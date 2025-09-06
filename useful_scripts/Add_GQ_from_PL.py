#!/usr/bin/env python3
"""
add_gq_from_pl.py  (PL-only)

Compute FORMAT/GQ from FORMAT/PL for each non-missing genotype:
  GQ = min( second_smallest(PL) - smallest(PL), 99 )
- Input: VCF (.vcf or .vcf.gz)
- Output: bgzipped VCF (.vcf.gz) if bgzip is available; otherwise plain VCF
- Multiprocessing supported, streaming (chunked) to limit memory

Notes:
- Missing or half-missing genotypes get GQ='.' (do NOT write 0).
- Only PL is considered (no GL fallback).
"""

from __future__ import annotations
import argparse
import gzip
import os
import shutil
import subprocess
import sys
import time
import multiprocessing as mp
from typing import List, Optional, Tuple, Dict, Any, Iterable


# ---------- small helpers ----------

def open_text(path: str):
    """Open text file; auto-detect gzip by magic bytes (suffix not required)."""
    try:
        with open(path, 'rb') as fh:
            if fh.read(2) == b"\x1f\x8b":
                return gzip.open(path, 'rt', encoding='utf-8', errors='replace')
    except Exception:
        pass
    return open(path, 'r', encoding='utf-8', errors='replace')


def write_text(path: str):
    return open(path, 'w', encoding='utf-8')


def parse_pl(pl: str) -> List[int]:
    if not pl or pl == '.':
        return []
    out: List[int] = []
    for tok in pl.split(','):
        tok = tok.strip()
        if not tok or tok == '.':
            continue
        try:
            out.append(int(tok))
        except Exception:
            # ignore malformed tokens
            continue
    return out


def compute_gq_from_pl(pl_vals: List[int]) -> Optional[int]:
    """GQ = min(99, second smallest PL − smallest PL)."""
    if len(pl_vals) < 2:
        return None
    s = sorted(pl_vals)
    diff = s[1] - s[0]
    if diff < 0:
        diff = 0
    return min(diff, 99)


def gt_is_missing(gt: Optional[str]) -> bool:
    """Treat '.', './.', '.|.', '0/.', './1', '0|.' etc. as missing."""
    if gt is None:
        return True
    if gt == '.':
        return True
    return '.' in gt


def add_gq_header(headers: List[str]) -> List[str]:
    """Insert GQ FORMAT header if absent (before #CHROM or after last FORMAT)."""
    if any(h.startswith('##FORMAT=<ID=GQ,') for h in headers):
        return headers
    insert_idx = 0
    for i, h in enumerate(headers):
        if h.startswith('##FORMAT='):
            insert_idx = i + 1
        if h.startswith('#CHROM'):
            break
    headers.insert(insert_idx, '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">')
    return headers


# ---------- core processing ----------

def process_chunk(lines: List[str]) -> Tuple[List[str], Dict[str, int]]:
    out_lines: List[str] = []
    stats = {'sites': 0, 'gt_total': 0, 'gt_non_missing': 0, 'non_missing_with_PL': 0}

    for ln in lines:
        parts = ln.rstrip('\n').split('\t')
        if len(parts) < 10:
            out_lines.append(ln)
            continue

        fmt_keys = parts[8].split(':') if parts[8] else []
        # Ensure GQ in FORMAT keys
        if 'GQ' not in fmt_keys:
            fmt_keys = fmt_keys + ['GQ']
            parts[8] = ':'.join(fmt_keys)

        # Indexes (may be -1 if absent)
        key_to_idx = {k: i for i, k in enumerate(fmt_keys)}
        gt_idx = key_to_idx.get('GT', -1)
        pl_idx = key_to_idx.get('PL', -1)
        gq_idx = key_to_idx['GQ']  # exists by construction

        stats['sites'] += 1
        n_samples = len(parts) - 9
        stats['gt_total'] += max(0, n_samples)

        new_samples: List[str] = []
        for s in parts[9:]:
            # Entire sample field missing
            if s == '.':
                fields = ['.'] * len(fmt_keys)
                fields[gq_idx] = '.'
                new_samples.append(':'.join(fields))
                continue

            fields = s.split(':')
            if len(fields) < len(fmt_keys):
                fields += ['.'] * (len(fmt_keys) - len(fields))

            gt_val = fields[gt_idx] if 0 <= gt_idx < len(fields) else '.'
            non_missing = not gt_is_missing(gt_val)

            if non_missing:
                stats['gt_non_missing'] += 1
                if 0 <= pl_idx < len(fields):
                    pls = parse_pl(fields[pl_idx])
                else:
                    pls = []

                if pls:
                    stats['non_missing_with_PL'] += 1
                    gq = compute_gq_from_pl(pls)
                    fields[gq_idx] = str(gq) if gq is not None else '.'
                else:
                    fields[gq_idx] = '.'
            else:
                fields[gq_idx] = '.'

            new_samples.append(':'.join(fields))

        out_parts = parts[:9] + new_samples
        out_lines.append('\t'.join(out_parts) + '\n')

    return out_lines, stats


def iter_variant_chunks(fin, first_variant_line: Optional[str], chunk_size: int) -> Iterable[List[str]]:
    """Yield lists of variant lines of size up to chunk_size."""
    chunk: List[str] = []
    if first_variant_line:
        chunk.append(first_variant_line)
    for ln in fin:
        if ln.startswith('#'):
            # Should not happen (headers already consumed), but pass through just in case
            continue
        chunk.append(ln)
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def run(input_vcf: str, output_vcf: str, processes: int = 1, chunk_size: int = 10000) -> Dict[str, Any]:
    # Read headers and capture the first variant line (stream-safe)
    with open_text(input_vcf) as fin:
        headers: List[str] = []
        first_variant: Optional[str] = None
        for ln in fin:
            if ln.startswith('#'):
                headers.append(ln.rstrip('\n'))
            else:
                first_variant = ln
                break

        if not headers:
            raise RuntimeError('No headers found in VCF')

        headers = add_gq_header(headers)

        # Prepare temp path and write headers
        tmp = output_vcf + '.tmp'
        totals = {'sites': 0, 'gt_total': 0, 'gt_non_missing': 0, 'non_missing_with_PL': 0}
        start = time.time()

        with write_text(tmp) as fo:
            for h in headers:
                fo.write(h + '\n')

            # If there are no variants at all
            if first_variant is None:
                pass
            else:
                chunks_iter = iter_variant_chunks(fin, first_variant, chunk_size)

                if processes <= 1:
                    for c in chunks_iter:
                        out_lines, stats = process_chunk(c)
                        fo.writelines(out_lines)
                        for k in totals:
                            totals[k] += stats.get(k, 0)
                else:
                    with mp.Pool(processes) as pool:
                        for out_lines, stats in pool.imap(process_chunk, chunks_iter):
                            fo.writelines(out_lines)
                            for k in totals:
                                totals[k] += stats.get(k, 0)

        # Compress/index
        bgzip_path = shutil.which('bgzip')
        tabix_path = shutil.which('tabix')

        if bgzip_path is None:
            # avoid misleading .gz extension
            if output_vcf.endswith('.gz'):
                plain_out = output_vcf[:-3]
                os.replace(tmp, plain_out)
                print(f"Warning: 'bgzip' not found. Wrote plain VCF to {plain_out}.", file=sys.stderr)
            else:
                os.replace(tmp, output_vcf)
                print("Warning: 'bgzip' not found. Wrote plain VCF (not indexed).", file=sys.stderr)
        else:
            subprocess.run([bgzip_path, '-f', tmp], check=True)
            os.replace(tmp + '.gz', output_vcf)
            if tabix_path:
                subprocess.run([tabix_path, '-f', '-p', 'vcf', output_vcf], check=False)
            else:
                print("Note: 'tabix' not found; output is bgzipped but not indexed.", file=sys.stderr)

    totals['write_time_s'] = time.time() - start
    return totals


# ---------- CLI ----------

def main() -> int:
    p = argparse.ArgumentParser(description='Compute FORMAT/GQ from PL (PL-only).')
    p.add_argument('--input_vcf', required=True, help='Input VCF (.vcf or .vcf.gz)')
    p.add_argument('--output_vcf', required=True, help='Output VCF (.vcf.gz if bgzip available)')
    p.add_argument('--processes', '-p', type=int, default=4, help='Worker processes (>=1)')
    p.add_argument('--chunk_size', '-c', type=int, default=10000, help='Variants per chunk (streaming)')
    args = p.parse_args()

    if not os.path.exists(args.input_vcf):
        print('[ERROR] Input VCF not found.', file=sys.stderr)
        return 1

    procs = max(1, args.processes)
    try:
        totals = run(args.input_vcf, args.output_vcf, processes=procs, chunk_size=max(100, args.chunk_size))
    except KeyboardInterrupt:
        print('Interrupted.', file=sys.stderr)
        return 130
    except Exception as e:
        print('[ERROR]', e, file=sys.stderr)
        return 2

    # Simple report
    nm = totals.get('gt_non_missing', 0)
    with_pl = totals.get('non_missing_with_PL', 0)
    pct = (with_pl / nm * 100.0) if nm else 0.0
    print(f"Sites: {totals.get('sites', 0)}")
    print(f"Genotypes total: {totals.get('gt_total', 0)}")
    print(f"Genotypes non-missing: {nm}")
    print(f"Non-missing with PL: {with_pl} ({pct:.2f}%)")
    print(f"Time (s): {totals.get('write_time_s', 0):.2f}")
    return 0


if __name__ == '__main__':
    sys.exit(main())

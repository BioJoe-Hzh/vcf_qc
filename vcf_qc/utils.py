"""Small utility helpers used across the vcf_qc package.

This module intentionally keeps a tiny surface area of pure-Python helpers that
are easy to unit-test and have no heavy dependencies.
"""
from typing import Dict, List, Tuple, Optional


def parse_info_field(info: str) -> Dict[str, Optional[str]]:
    """Parse a VCF INFO column (key[=value];... ) into a dict.

    Values are returned as strings; keys without value map to empty string.
    An INFO field of '.' returns an empty dict.
    """
    out: Dict[str, Optional[str]] = {}
    if not info or info == ".":
        return out
    for token in info.split(";"):
        if not token:
            continue
        if "=" in token:
            k, v = token.split("=", 1)
            out[k] = v
        else:
            out[token] = ""
    return out


def parse_format_sample(fmt: str, sample: str) -> Dict[str, Optional[str]]:
    """Parse FORMAT and a sample column into a dict mapping keys->values.

    Example: fmt='GT:AD:DP' sample='0/1:10,5:15' -> {'GT':'0/1','AD':'10,5','DP':'15'}
    Missing fields are mapped to None.
    """
    keys = fmt.split(":") if fmt else []
    vals = sample.split(":") if sample else []
    out: Dict[str, Optional[str]] = {}
    for i, k in enumerate(keys):
        out[k] = vals[i] if i < len(vals) and vals[i] != "" else None
    return out


def extract_ad_dp(sample_dict: Dict[str, Optional[str]]) -> Tuple[Optional[List[int]], Optional[int]]:
    """Return (AD_list, DP) parsed from a sample dict.

    AD is returned as a list of ints (REF then ALT depths) if present, otherwise None.
    DP is returned as int if present and parseable, otherwise None.
    """
    ad = sample_dict.get("AD") if sample_dict else None
    dp = sample_dict.get("DP") if sample_dict else None
    ad_list: Optional[List[int]] = None
    dp_int: Optional[int] = None
    if ad:
        try:
            ad_list = [int(x) for x in str(ad).split(",")]
        except Exception:
            ad_list = None
    if dp:
        try:
            dp_int = int(dp)
        except Exception:
            dp_int = None
    return ad_list, dp_int


def vaf_from_ad(ad_list: List[int], allele_index: int = 1) -> Optional[float]:
    """Compute variant allele fraction for allele_index (1-based alt index).

    ad_list is expected as [ref_count, alt1_count, alt2_count, ...]. allele_index=1 uses alt1.
    Returns None if counts are missing or denom is zero.
    """
    if not ad_list or allele_index < 1 or allele_index >= len(ad_list):
        return None
    try:
        alt = ad_list[allele_index]
        total = sum(ad_list)
        if total <= 0:
            return None
        return float(alt) / float(total)
    except Exception:
        return None


def pl_to_gq(pl_vals: List[int]) -> Optional[int]:
    """Convert a list of PL (phred-scaled likelihoods) into an approximate GQ.

    Simple heuristic: GQ = difference between smallest PL and second-smallest PL.
    According to VCF format specification, GQ should be capped at 99.
    Returns None for empty input or if fewer than 2 values.
    """
    if not pl_vals or len(pl_vals) < 2:
        return None
    try:
        sorted_pl = sorted(pl_vals)
        gq = sorted_pl[1] - sorted_pl[0]
        if gq < 0:
            gq = 0
        return int(min(gq, 99))  # Cap GQ at 99 per VCF specification
    except Exception:
        return None


def gt_is_missing(gt: Optional[str]) -> bool:
    """Check if a genotype is missing or half-missing.
    
    Treat '.', './.', '.|.', '0/.', './1', '0|.' etc. as missing.
    This follows the VCF specification for missing genotypes.
    """
    if gt is None:
        return True
    if gt == '.':
        return True
    return '.' in gt


def parse_pl_robust(pl: Optional[str]) -> List[int]:
    """Parse PL field robustly, handling missing and malformed values.
    
    Returns empty list if PL is missing or contains no valid values.
    Individual missing values (.) within the PL list are skipped.
    """
    if not pl or pl == '.':
        return []
    
    out: List[int] = []
    for tok in pl.split(','):
        tok = tok.strip()
        if not tok or tok == '.':
            continue  # Skip missing individual values
        try:
            out.append(int(tok))
        except (ValueError, TypeError):
            continue  # Skip malformed values
    return out


def normalize_chrom(chrom: Optional[str]) -> str:
    """Lightweight normalization for chromosome names.

    Examples: 'chr1' -> '1', '1' -> '1', 'MT'->'MT'
    This is intentionally conservative and only strips a leading 'chr' or 'CHR'.
    """
    if chrom is None:
        return ""
    c = str(chrom)
    if c.lower().startswith("chr"):
        return c[3:]
    return c


__all__ = [
    "parse_info_field",
    "parse_format_sample",
    "extract_ad_dp",
    "vaf_from_ad",
    "pl_to_gq",
    "gt_is_missing",
    "parse_pl_robust",
    "normalize_chrom",
]

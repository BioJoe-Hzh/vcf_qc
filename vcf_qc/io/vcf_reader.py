"""Lightweight VCF reader for QC metric extraction.

This intentionally avoids external dependencies (pysam / cyvcf2) for
an initial, streaming implementation sufficient for computing summary
statistics (missing rate, depth, heterozygosity, GQ, etc.). For large
VCFs or performance‑critical workflows you should replace this with a
backend powered by specialised libraries.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Dict, Iterator, Optional, Tuple, Iterable
import gzip

from ..utils import pl_to_gq


@dataclass
class GenotypeRecord:
	"""Container for a single locus' per-sample genotype attributes.

	Attributes
	----------
	chrom, pos, ref, alts : str
		Basic variant fields (POS is kept as string to avoid int cast cost;
		cast only if needed by caller).
	format_keys : List[str]
		FORMAT field keys in order.
	sample_fields : List[str]
		Raw colon-delimited field strings per sample.
	"""

	chrom: str
	pos: str
	ref: str
	alts: str
	qual: Optional[str]
	info: Optional[str]
	format_keys: List[str]
	sample_fields: List[str]

	def extract_sample_value(self, sample_index: int, key: str) -> Optional[str]:
		"""Extract a value (e.g., DP, GQ, GT) for a given sample.

		Returns None if key not in FORMAT or value is '.'
		"""
		try:
			fi = self.format_keys.index(key)
		except ValueError:
			return None
		raw = self.sample_fields[sample_index]
		parts = raw.split(":")
		if fi >= len(parts):
			return None
		val = parts[fi]
		return None if val == "." else val


class SimpleVCFReader:
	"""Minimal streaming VCF reader.

	Parameters
	----------
	path : str
		Path to (optionally gzipped) VCF file.
	max_records : int | None
		Optional limit for testing / faster prototyping.
	"""

	def __init__(self, path: str, max_records: Optional[int] = None):
		self.path = path
		self.max_records = max_records
		self.samples: List[str] = []

	# -- internal helpers -------------------------------------------------
	def _open(self):  # type: ignore[return-type]
		if self.path.endswith('.gz'):
			return gzip.open(self.path, 'rt')
		return open(self.path, 'rt')

	def parse(self) -> Iterator[GenotypeRecord]:
		count = 0
		with self._open() as fh:
			for line in fh:
				if not line:
					continue
				if line.startswith('##'):
					continue
				if line.startswith('#CHROM'):
					header_cols = line.rstrip().split('\t')
					# VCF fixed columns then samples from index 9
					self.samples = header_cols[9:]
					continue
				parts = line.rstrip().split('\t')
				if len(parts) < 10:  # no samples
					continue
				chrom, pos, _id, ref, alts, qual, flt, info, fmt = parts[:9]
				format_keys = fmt.split(':')
				sample_fields = parts[9:]
				yield GenotypeRecord(chrom, pos, ref, alts, qual if qual != '.' else None, info, format_keys, sample_fields)
				count += 1
				if self.max_records and count >= self.max_records:
					break

	# -- convenience metric extraction ------------------------------------
	def compute_sample_level_stats(self) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, int], Dict[str, int], Dict[str, int], Dict[str, int]]:
		"""Compute core counts needed for per-sample metrics.

		Returns tuple of dicts (total_calls, missing_calls, het_calls,
		hom_alt_calls, depth_sum, depth_count).

		Definitions
		----------
		- het_calls: heterozygous (alleles differ, both non-missing)
		- hom_alt_calls: homozygous non-reference (e.g. 1/1, 2/2)
		- depth_sum / depth_count: for computing mean depth across genotypes with numeric DP.
		"""
		if not self.samples:
			# trigger parsing until header found
			for _ in self.parse():
				break
			# reset generator for full pass
		total_calls = {s: 0 for s in self.samples}
		missing_calls = {s: 0 for s in self.samples}
		het_calls = {s: 0 for s in self.samples}
		hom_alt_calls = {s: 0 for s in self.samples}
		depth_sum = {s: 0 for s in self.samples}
		depth_count = {s: 0 for s in self.samples}

		# Re-iterate fully
		for rec in self.parse():
			for i, sample in enumerate(self.samples):
				gt = rec.extract_sample_value(i, 'GT')
				total_calls[sample] += 1
				if gt is None or gt in ('.', './.', '.|.'):
					missing_calls[sample] += 1
					continue
				# Split alleles; treat any '.' allele as missing (ignore partially called like 0/.)
				sep = '/' if '/' in gt else '|' if '|' in gt else None
				alleles = gt.split(sep) if sep else [gt]
				if any(a == '.' for a in alleles):
					missing_calls[sample] += 1
					continue
				if len(alleles) == 2:  # diploid assumptions
					a, b = alleles
					if a != b:
						het_calls[sample] += 1  # heterozygous (0/1,0/2,1/2,...)
					elif a == b and a != '0':
						hom_alt_calls[sample] += 1  # homozygous alt (1/1,2/2,...)
				# Depth
				dp = rec.extract_sample_value(i, 'DP')
				if dp and dp.isdigit():
					depth_sum[sample] += int(dp)
					depth_count[sample] += 1
		return total_calls, missing_calls, het_calls, hom_alt_calls, depth_sum, depth_count

	def iterate_gq_depth(self, compute_gq_from_pl: bool = False) -> Iterator[Tuple[str, int, int]]:
		"""Yield (sample, GQ, DP) triples for genotypes with values.

		Useful for constructing distribution plots (GQ vs depth boxplots).
		
		Args:
			compute_gq_from_pl: If True, compute GQ from PL when GQ is missing.
		"""
		if not self.samples:
			for _ in self.parse():
				break
		
		compute_count = 0
		for rec in self.parse():
			for i, sample in enumerate(self.samples):
				# Get raw values
				gt = rec.extract_sample_value(i, 'GT')
				gq = rec.extract_sample_value(i, 'GQ')
				dp = rec.extract_sample_value(i, 'DP')
				
				# Skip missing genotypes (following VCF standard)
				from ..utils import gt_is_missing
				if gt_is_missing(gt):
					continue
				
				# Try to get or compute GQ
				gq_val = None
				if gq and gq.isdigit():
					gq_val = int(gq)
				elif compute_gq_from_pl and gq in (None, '.', ''):
					# Try to compute GQ from PL
					pl_str = rec.extract_sample_value(i, 'PL')
					if pl_str and pl_str not in (None, '.', ''):
						try:
							from ..utils import parse_pl_robust, pl_to_gq
							pl_values = parse_pl_robust(pl_str)
							if pl_values:
								gq_val = pl_to_gq(pl_values)
								if gq_val is not None:
									compute_count += 1
									# Adaptive progress display: start with 1K, then increase interval
									if compute_count <= 10000:
										interval = 1000
									elif compute_count <= 100000:
										interval = 10000
									else:
										interval = 50000
									
									if compute_count % interval == 0:
										print(f"[INFO] Computing GQ from PL: {compute_count:,} genotypes processed...")
						except Exception:
							continue
				
				# Skip if we couldn't get or compute GQ
				if gq_val is None:
					continue
					
				# Parse depth
				if dp and dp.isdigit():
					dp_val = int(dp)
				else:
					dp_val = 0
					
				yield sample, gq_val, dp_val

	# -- site level -------------------------------------------------------
	def iterate_site_metrics(self) -> Iterator[Dict[str, object]]:
		"""Yield per-site summary dicts (QUAL, total depth, QD, MAC, MAF, MissingRate).

		Handles multi-allelic sites properly by leveraging INFO fields when
		present (AC, AN, MAF, DP). Falls back to per-genotype parsing when
		these are absent.
		"""
		if not self.samples:
			for _ in self.parse():
				break
		sample_count = len(self.samples)
		for rec in self.parse():
			# -- Attempt fast path via INFO ---------------------------------
			info_map: Dict[str, str] = {}
			if rec.info:
				for kv in rec.info.split(';'):
					if '=' in kv:
						k, v = kv.split('=', 1)
						info_map[k] = v
			# Missing rate requires GT parsing anyway
			missing = 0
			depth_total = 0
			depth_obs = 0  # number of samples with numeric DP (used for mean depth)
			for i, _s in enumerate(self.samples):
				gt = rec.extract_sample_value(i, 'GT')
				dp = rec.extract_sample_value(i, 'DP')
				if dp and dp.isdigit():
					depth_total += int(dp)
					depth_obs += 1
				if not gt or gt in ('.', './.', '.|.'):
					missing += 1
			# Allele counting (definition: MAF = frequency of the SECOND most frequent allele, including REF)
			# This differs from the common "minor (least frequent)" allele definition.
			# Rationale (user requirement): for multi-allelic sites ignore extremely rare alleles; take the 2nd highest.
			mac = 0  # here MAC means count of the 2nd most frequent allele
			maf = 0.0
			maf_done = False
			if 'AC' in info_map and 'AN' in info_map:
				try:
					ac_values = [int(x) for x in info_map['AC'].split(',') if x]
					an = int(info_map['AN']) if info_map['AN'].isdigit() else None
					if an and an > 0:
						ref_count = an - sum(ac_values)
						counts = []
						if ref_count > 0:
							counts.append(ref_count)
						counts.extend([c for c in ac_values if c > 0])
						if len(counts) >= 2:
							sorted_counts = sorted(counts, reverse=True)
							mac = sorted_counts[1]
							maf = mac / an
						maf_done = True
				except ValueError:
					pass
			# Fallback: derive from GTs
			if not maf_done:
				allele_counts: Dict[int, int] = {}
				for i, _s in enumerate(self.samples):
					gt = rec.extract_sample_value(i, 'GT')
					if not gt or gt in ('.', './.', '.|.'):
						continue
					sep = '/' if '/' in gt else '|' if '|' in gt else None
					alleles = gt.split(sep) if sep else [gt]
					for a in alleles:
						if a.isdigit():
							idx = int(a)
							allele_counts[idx] = allele_counts.get(idx, 0) + 1
				if allele_counts:
					counts = [c for c in allele_counts.values() if c > 0]
					an2 = sum(allele_counts.values())
					if len(counts) >= 2 and an2 > 0:
						sorted_counts = sorted(counts, reverse=True)
						mac = sorted_counts[1]
						maf = mac / an2
			qual_val = float(rec.qual) if rec.qual and rec.qual.replace('.', '', 1).isdigit() else None
			# Prefer INFO DP if present (total depth) else sum of sample DP; then convert to mean
			if 'DP' in info_map and info_map['DP'].isdigit():
				total_depth_for_mean = int(info_map['DP'])
			else:
				total_depth_for_mean = depth_total
			mean_depth = (total_depth_for_mean / depth_obs) if depth_obs > 0 else 0
			qd = (qual_val / mean_depth) if qual_val is not None and mean_depth > 0 else None
			# Count distinct alleles (including REF). ALTs may be comma separated; '.' means no ALT.
			allele_count = 1
			if rec.alts and rec.alts != '.':
				allele_count += rec.alts.count(',') + 1
			yield {
				'Chrom': rec.chrom,
				'Pos': int(rec.pos),
				'QUAL': qual_val,
				'MeanDepth': mean_depth,
				'QD': qd,
				'MAC': mac,
				'MAF': maf,
				'MissingRate': missing / sample_count if sample_count else 0.0,
				'AlleleCount': allele_count,
				'Alts': rec.alts,
			}

	# -- genotype level ---------------------------------------------------
	def iterate_genotypes(self, compute_gq_from_pl: bool = False) -> Iterator[Dict[str, object]]:
		"""Yield per-genotype dict with counts suitable for depth/GQ/VAF plots.

		If AD (allele depth) is present: AD=ref,alt; compute ALT_Count, REF_Count and VAF.
		Multi-allelic aware: if GT references alt allele index >1 (e.g. 0/2, 1/2),
		we take the sum of depths for the alt alleles *present in the genotype*.
		VAF is then defined as sum(alt depths of alleles in GT) / (ref depth (if 0 in GT) + that sum).
		This yields a "genotype allele fraction" rather than the proportion of *all* alt reads at the site.
		If you prefer denominator = sum(all AD) you can adapt easily.

		Args:
			compute_gq_from_pl: If True, compute GQ from PL when GQ is missing (slower).

		Else (no AD) VAF cannot be derived reliably; we leave it as None.
		"""
		if not self.samples:
			for _ in self.parse():
				break
		
		# Track how many GQ values were computed from PL
		gq_computed_from_pl = 0
		total_genotypes_processed = 0
		
		for rec in self.parse():
			for idx, sample in enumerate(self.samples):
				total_genotypes_processed += 1
				gt = rec.extract_sample_value(idx, 'GT')
				gq = rec.extract_sample_value(idx, 'GQ')
				dp = rec.extract_sample_value(idx, 'DP')
				ad = rec.extract_sample_value(idx, 'AD')
				pl = rec.extract_sample_value(idx, 'PL')
				
				# Skip missing genotypes (following VCF standard)
				from ..utils import gt_is_missing
				if gt_is_missing(gt):
					continue
				
				# Try to get GQ, or compute from PL if enabled and not available
				gq_val = None
				if gq and gq.isdigit():
					gq_val = int(gq)
				elif compute_gq_from_pl and pl:
					# Try to compute GQ from PL using robust parsing
					try:
						from ..utils import parse_pl_robust, pl_to_gq
						pl_vals = parse_pl_robust(pl)
						if pl_vals:
							gq_val = pl_to_gq(pl_vals)
							if gq_val is not None:
								gq_computed_from_pl += 1
								# Adaptive progress display: start with 1K, then increase interval
								if gq_computed_from_pl <= 10000:
									interval = 1000
								elif gq_computed_from_pl <= 100000:
									interval = 10000
								else:
									interval = 50000
								
								if gq_computed_from_pl % interval == 0:
									print(f"[INFO] Computing GQ from PL: {gq_computed_from_pl:,} genotypes processed...")
					except Exception:
						pass
				
				# Skip if we couldn't get or compute GQ
				if gq_val is None:
					continue
				depth_val = int(dp) if dp and dp.isdigit() else 0
				alt_count: Optional[int] = None
				ref_count: Optional[int] = None
				vaf: Optional[float] = None
				# Multi-allelic aware parsing of AD & GT
				if ad and ',' in ad:
					ad_tokens = ad.split(',')
					ad_values = []
					for tok in ad_tokens:
						if tok.isdigit():
							ad_values.append(int(tok))
						else:
							ad_values.append(0)
					if len(ad_values) >= 2:
						ref_count = ad_values[0]
						# Determine which alt allele indices are present in genotype
						allele_indices_in_gt = []
						if isinstance(gt, str) and gt not in {'.', './.', '.|.'}:
							sep = '/' if '/' in gt else '|' if '|' in gt else None
							parts = gt.split(sep) if sep else [gt]
							for a in parts:
								if a.isdigit():
									ai = int(a)
									if ai > 0:  # alt allele (1-based for first ALT)
										allele_indices_in_gt.append(ai)
						# Sum depths for alt alleles present in genotype (unique indices)
						if allele_indices_in_gt:
							unique_idx = sorted(set(allele_indices_in_gt))
							alt_sum = 0
							for ai in unique_idx:
								if ai < len(ad_values):  # ai maps: 1 -> ad_values[1]
									alt_sum += ad_values[ai]
							alt_count = alt_sum
						else:
							alt_count = 0
						# Define VAF using only alleles present in genotype (ref if 0 present)
						ref_in_gt = False
						if isinstance(gt, str) and '0' in gt.replace('|', '/').split('/'):
							ref_in_gt = True
						denom = (ref_count if ref_in_gt else 0) + (alt_count if alt_count is not None else 0)
						if denom > 0 and alt_count is not None:
							vaf = alt_count / denom
			yield {
				'Sample': sample,
				'Chrom': rec.chrom,
				'Pos': int(rec.pos),
				'GQ': gq_val,
				'Depth': depth_val,
				'ALT_Count': alt_count,
				'REF_Count': ref_count,
				'VAF': vaf,
				'GT': gt,
			}


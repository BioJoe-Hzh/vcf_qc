"""Sample-level metric computations for VCF QC.

Relies on the lightweight ``SimpleVCFReader``. Replace / extend with a
high-performance backend as needed.
"""

from __future__ import annotations

from typing import Dict, Tuple
import pandas as pd

from ..io.vcf_reader import SimpleVCFReader

__all__ = [
	"compute_sample_metrics",
	"gq_depth_long_table",
]



def compute_sample_metrics(reader: SimpleVCFReader) -> pd.DataFrame:
	"""Compute core per-sample metrics including new Heterozygosity Ratio.

	Columns returned:
		Sample, MissingRate, MeanDepth, HetRatio, CalledSites, HomAltSites

	Definition:
		HetRatio = HET / HOM_ALT  (if HOM_ALT==0 and HET>0 -> NaN; if both 0 -> 0)
	"""
	total, missing, het, hom_alt, depth_sum, depth_count = reader.compute_sample_level_stats()
	rows = []
	for s in reader.samples:
		tot = total[s] or 1  # guard division by zero
		called = tot - missing[s]
		mean_depth = depth_sum[s] / depth_count[s] if depth_count[s] else 0.0
		hom_alt_sites = hom_alt[s]
		if hom_alt_sites == 0:
			if het[s] == 0:
				het_ratio = 0.0
			else:
				het_ratio = float('nan')  # undefined ratio; could also choose large sentinel
		else:
			het_ratio = het[s] / hom_alt_sites
		rows.append({
			"Sample": s,
			"MissingRate": missing[s] / tot,
			"MeanDepth": mean_depth,
			"HetRatio": het_ratio,
			"CalledSites": called,
			"HomAltSites": hom_alt_sites,
		})
	return pd.DataFrame(rows)


def gq_depth_long_table(reader: SimpleVCFReader) -> pd.DataFrame:
	"""Return long-form table with GQ and Depth per sample genotype.

	Columns: Sample, GQ, Depth
	"""
	records = [
		{"Sample": s, "GQ": gq, "Depth": dp}
		for s, gq, dp in reader.iterate_gq_depth()
	]
	return pd.DataFrame(records)


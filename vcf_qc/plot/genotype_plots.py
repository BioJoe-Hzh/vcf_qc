"""Genotype-level QC plotting functions.

Implements:
 4.8 GQ distribution (genotype level)
 4.9 Depth distribution (genotype level)
 4.10 GQ vs depth (joint)
 4.11 VAF vs depth (heterozygous only, joint)
 4.12 ALT read count vs depth (heterozygous only, joint)
"""

from __future__ import annotations

from typing import Optional, Union, Tuple

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from .base import set_plot_style, save_figure, joint_scatter, hist_plot

__all__ = [
	"plot_gq_distribution_genotype",
	"plot_depth_distribution_genotype",
	"plot_gq_vs_depth_genotype",
	"plot_vaf_vs_depth_het",
	"plot_alt_count_vs_depth_het",
]


def plot_gq_distribution_genotype(
	gq: Union[pd.Series, np.ndarray, list],
	*,
	output_path: Optional[str] = None,
	title: str = "Genotype GQ distribution",
	bins: int = 60,
	smart_cutoff: float = 99.0,
	enable_smart_cutoff: bool = True,
) -> Optional[plt.Figure]:
	"""4.8 Histogram of genotype GQ values."""
	# Use fixed integer bins 1..99 (99 bins of width 1) regardless of user 'bins' parameter.
	# This ensures compliance with VCF standard where GQ should not exceed 99.
	fixed_bins = np.arange(1, 101)  # edges 1..100 create 1..99 bins
	return hist_plot(
		gq,
		output_path=output_path,
		title=title,
		xlabel="GQ",
		bins=fixed_bins,
		color="#1565C0",
		smart_cutoff=smart_cutoff,
		enable_smart_cutoff=enable_smart_cutoff,
	)


def plot_depth_distribution_genotype(
	depth: Union[pd.Series, np.ndarray, list],
	*,
	output_path: Optional[str] = None,
	title: str = "Genotype depth distribution",
	bins: int = 60,
	logx: bool = False,
	smart_cutoff: float = 99.5,
	enable_smart_cutoff: bool = True,
) -> Optional[plt.Figure]:
	"""4.9 Histogram of genotype depth values."""
	return hist_plot(
		depth,
		output_path=output_path,
		title=title,
		xlabel="Depth",
		bins=bins,
		logx=logx,
		color="#2E7D32",
		smart_cutoff=smart_cutoff,
		enable_smart_cutoff=enable_smart_cutoff,
	)


def plot_gq_vs_depth_genotype(
	data: pd.DataFrame,
	*,
	output_path: Optional[str] = None,
	title: str = "GQ vs depth (genotype level)",
	gq_col: str = "GQ",
	depth_col: str = "Depth",
	smart_cutoff: float = 99.0,
	enable_smart_cutoff: bool = True,
	apply_on: str = "both",
	# Distinct high-contrast colors
	color_hom_ref: str = "#1474b8",  # blue
	color_hom_alt: str = "#c70202",  # red
	color_het: str = "#f58b09",      # orange
	point_alpha: float = 0.2,
) -> Optional[plt.Figure]:
	"""4.10 Joint scatter of GQ vs depth with zygosity colouring via joint_scatter."""
	if not {gq_col, depth_col}.issubset(data.columns):
		raise ValueError("DataFrame must contain GQ and Depth columns")
	df = data[[depth_col, gq_col] + (["GT"] if "GT" in data.columns else [])].copy()
	df.rename(columns={depth_col: "Depth", gq_col: "GQ"}, inplace=True)
	if "GT" in df.columns:
		def classify(gt: object) -> str:
			if not isinstance(gt, str):
				return "Other"
			if gt in {".", "./.", ".|."}:
				return "Other"
			sep = '/' if '/' in gt else ('|' if '|' in gt else None)
			if not sep:
				return "Other"
			parts = [p for p in gt.split(sep) if p != '.']
			if len(parts) != 2:
				return "Other"
			if parts[0] != parts[1]:
				return "Heterozygous"
			# Homozygous: decide ref vs alt by allele code '0'
			return "HomRef" if parts[0] == '0' else "HomAlt"
		df['Zygosity'] = df['GT'].map(classify)
	else:
		df['Zygosity'] = "Other"
	palette = {
		"HomRef": color_hom_ref,
		"HomAlt": color_hom_alt,
		"Heterozygous": color_het,
		"Other": "#C7C5C5",
	}
	return joint_scatter(
		df,
		x="Depth",
		y="GQ",
		title=title,
		output_path=output_path,
		xlabel="Depth",
		ylabel="GQ",
		color="#00838F",
		smart_cutoff=smart_cutoff,
		enable_smart_cutoff=enable_smart_cutoff,
		apply_on=apply_on,
		hue="Zygosity",
		palette=palette,
		alpha=point_alpha,
	)


def plot_vaf_vs_depth_het(
	data: pd.DataFrame,
	*,
	output_path: Optional[str] = None,
	title: str = "VAF vs depth (heterozygous genotypes only)",
	vaf_col: str = "VAF",
	depth_col: str = "Depth",
	smart_cutoff: float = 99.5,
	enable_smart_cutoff: bool = True,
	apply_on: str = "both",
	point_alpha: float = 0.2,
) -> Optional[plt.Figure]:
	"""4.11 Joint scatter of VAF vs depth restricted to heterozygous genotypes."""
	if not {vaf_col, depth_col}.issubset(data.columns):
		raise ValueError("DataFrame must contain VAF and Depth columns")
	return joint_scatter(
		data.rename(columns={vaf_col: "VAF", depth_col: "Depth"}),
		x="Depth",
		y="VAF",
		title=title,
		output_path=output_path,
		xlabel="Depth",
		ylabel="VAF",
		color="#6A1B9A",
		smart_cutoff=smart_cutoff,
		enable_smart_cutoff=enable_smart_cutoff,
		apply_on=apply_on,
		alpha=point_alpha,
	)


def plot_alt_count_vs_depth_het(
	data: pd.DataFrame,
	*,
	output_path: Optional[str] = None,
	title: str = "ALT count vs depth (heterozygous genotypes only)",
	alt_count_col: str = "ALT_Count",
	depth_col: str = "Depth",
	smart_cutoff: float = 99.5,
	enable_smart_cutoff: bool = True,
	apply_on: str = "both",
	point_alpha: float = 0.2,
) -> Optional[plt.Figure]:
	"""4.12 Joint scatter of ALT read count vs depth (heterozygous genotypes)."""
	if not {alt_count_col, depth_col}.issubset(data.columns):
		raise ValueError("DataFrame must contain ALT_Count and Depth columns")
	return joint_scatter(
		data.rename(columns={alt_count_col: "ALT_Count", depth_col: "Depth"}),
		x="Depth",
		y="ALT_Count",
		title=title,
		output_path=output_path,
		xlabel="Depth",
		ylabel="ALT read count",
		color="#C62828",
		smart_cutoff=smart_cutoff,
		enable_smart_cutoff=enable_smart_cutoff,
		apply_on=apply_on,
		alpha=point_alpha,
	)


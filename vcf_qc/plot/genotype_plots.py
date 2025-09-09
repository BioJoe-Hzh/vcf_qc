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
	show_expected_vaf: bool = True,
) -> Optional[plt.Figure]:
	"""4.11 Joint scatter of VAF vs depth restricted to heterozygous genotypes.
	
	For heterozygous genotypes, VAF should theoretically be around 0.5.
	This plot helps identify allelic bias or other technical artifacts.
	"""
	if not {vaf_col, depth_col}.issubset(data.columns):
		raise ValueError("DataFrame must contain VAF and Depth columns")
	
	# Add sample count to title
	n_samples = len(data)
	enhanced_title = f"{title}\nn = {n_samples:,} heterozygous genotypes"
	
	# Create the base scatter plot
	fig = joint_scatter(
		data.rename(columns={vaf_col: "VAF", depth_col: "Depth"}),
		x="Depth",
		y="VAF",
		title=enhanced_title,
		output_path=None,  # We'll handle saving manually to add reference line
		xlabel="Depth",
		ylabel="VAF",
		color="#6A1B9A",
		smart_cutoff=smart_cutoff,
		enable_smart_cutoff=enable_smart_cutoff,
		apply_on=apply_on,
		alpha=point_alpha,
	)
	
	if fig is not None and show_expected_vaf:
		# Add horizontal line at VAF=0.5 (expected for heterozygotes)
		ax = fig.axes[0]  # Main scatter plot axis
		ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.8, linewidth=2, 
				   label='Expected VAF=0.5')
		ax.legend(loc='upper right')
		
		# Add some statistics text
		vaf_values = data[vaf_col].dropna()
		if len(vaf_values) > 0:
			mean_vaf = vaf_values.mean()
			median_vaf = vaf_values.median()
			stats_text = f'Mean VAF: {mean_vaf:.3f}\nMedian VAF: {median_vaf:.3f}'
			ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
					verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
	
	return save_figure(fig, output_path)


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
	show_expected_ratio: bool = True,
) -> Optional[plt.Figure]:
	"""4.12 Joint scatter of ALT read count vs depth (heterozygous genotypes).
	
	For heterozygous genotypes, ALT count should theoretically be around 50% of total depth.
	This plot helps identify allelic bias in read coverage.
	"""
	if not {alt_count_col, depth_col}.issubset(data.columns):
		raise ValueError("DataFrame must contain ALT_Count and Depth columns")
	
	# Add sample count to title
	n_samples = len(data)
	enhanced_title = f"{title}\nn = {n_samples:,} heterozygous genotypes"
	
	# Create the base scatter plot
	fig = joint_scatter(
		data.rename(columns={alt_count_col: "ALT_Count", depth_col: "Depth"}),
		x="Depth",
		y="ALT_Count",
		title=enhanced_title,
		output_path=None,  # We'll handle saving manually to add reference line
		xlabel="Depth",
		ylabel="ALT read count",
		color="#C62828",
		smart_cutoff=smart_cutoff,
		enable_smart_cutoff=enable_smart_cutoff,
		apply_on=apply_on,
		alpha=point_alpha,
	)
	
	if fig is not None and show_expected_ratio:
		# Add diagonal line showing expected 50% ratio (ALT_Count = 0.5 * Depth)
		ax = fig.axes[0]  # Main scatter plot axis
		
		# Get the current axis limits for the reference line
		xlim = ax.get_xlim()
		x_range = np.linspace(xlim[0], xlim[1], 100)
		y_expected = 0.5 * x_range
		
		ax.plot(x_range, y_expected, 'r--', alpha=0.8, linewidth=2, 
				label='Expected ALT=0.5×Depth')
		ax.legend(loc='upper left')
		
		# Add some statistics text
		alt_values = data[alt_count_col].dropna()
		depth_values = data[depth_col].dropna()
		if len(alt_values) > 0 and len(depth_values) > 0:
			# Calculate actual ALT/Depth ratio
			combined_data = data[[alt_count_col, depth_col]].dropna()
			if len(combined_data) > 0:
				actual_ratios = combined_data[alt_count_col] / combined_data[depth_col]
				mean_ratio = actual_ratios.mean()
				median_ratio = actual_ratios.median()
				stats_text = f'Mean ALT/Depth: {mean_ratio:.3f}\nMedian ALT/Depth: {median_ratio:.3f}'
				ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
						verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
	
	return save_figure(fig, output_path)


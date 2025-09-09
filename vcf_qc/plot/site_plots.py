"""Site-level QC plotting functions.

Implements:
 4.5 QUAL distribution
 4.6 Depth distribution (site-total depth)
 4.7 QUAL vs depth scatter + QD distribution
 4.13 Minor allele count vs depth (joint)
 4.14 Minor allele frequency vs depth (joint)
 4.15 Site missing rate distribution (after setGT)
"""

from __future__ import annotations

from typing import Optional, Union, Tuple, Dict

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from .base import set_plot_style, save_figure, joint_scatter, hist_plot

__all__ = [
	"plot_qual_distribution",
	"plot_site_mean_depth_distribution",
	"plot_qual_vs_mean_depth",
	"plot_qd_distribution",
	"plot_minor_allele_count_vs_depth",
	"plot_minor_allele_freq_vs_depth",
	"plot_site_missing_rate_distribution",
	"plot_allele_count_pie",
	"plot_mac_distribution",
	"plot_maf_distribution",
]


def plot_qual_distribution(
	qual: Union[pd.Series, np.ndarray, list],
	*,
	output_path: Optional[str] = None,
	title: str = "QUAL distribution (site level)",
	bins: int = 60,
	logx: bool = False,
	smart_cutoff: float = 99.5,
	enable_smart_cutoff: bool = True,
	use_smart_cutoff: bool = True,
	focus_pct: Optional[float] = None,
) -> Optional[plt.Figure]:
	"""4.5 Histogram of QUAL values for sites."""
	values = np.asarray(qual, dtype=float)
	orig_mask = ~np.isnan(values)
	focus_note = ""
	if focus_pct is not None and orig_mask.any():
		pct_val = focus_pct * 100 if 0 < focus_pct <= 1 else focus_pct
		if 0 < pct_val < 100:
			cut = np.percentile(values[orig_mask], pct_val)
			values = values[values <= cut]
			focus_note = f" (<= {pct_val:.2f}th pct, cutoff={cut:g})"
	return hist_plot(
		values,
		output_path=output_path,
		title=title + focus_note,
		xlabel="QUAL",
		bins=bins,
		logx=logx,
		color="#355C7D",
		smart_cutoff=smart_cutoff,
		enable_smart_cutoff=enable_smart_cutoff and use_smart_cutoff,
	)


def plot_site_mean_depth_distribution(
	depth: Union[pd.Series, np.ndarray, list],
	*,
	output_path: Optional[str] = None,
	title: str = "Mean depth distribution (site level)",
	bins: int = 60,
	logx: bool = False,
	smart_cutoff: float = 99.5,
	enable_smart_cutoff: bool = True,
	use_smart_cutoff: bool = True,
) -> Optional[plt.Figure]:
	"""4.6 Histogram of site mean depth values."""
	return hist_plot(
		depth,
		output_path=output_path,
		title=title,
		xlabel="Depth",
		bins=bins,
		logx=logx,
		color="#2A9D8F",
		smart_cutoff=smart_cutoff,
		enable_smart_cutoff=enable_smart_cutoff and use_smart_cutoff,
	)


def plot_qual_vs_mean_depth(
	data: pd.DataFrame,
	*,
	output_path: Optional[str] = None,
	title: str = "QUAL vs mean depth",
	qual_col: str = "QUAL",
	depth_col: str = "MeanDepth",
	use_smart_cutoff: bool = True,
) -> Optional[plt.Figure]:
	"""Scatter QUAL vs mean depth (separate from QD distribution)."""
	for c in (qual_col, depth_col):
		if c not in data.columns:
			raise ValueError(f"Missing required column '{c}'")
	
	# Apply smart cutoff if enabled (similar to joint_scatter logic)
	plot_data = data.copy()
	cutoff_note = ""
	if use_smart_cutoff and len(plot_data) > 0:
		# Apply 99.5th percentile cutoff for both dimensions
		qual_vals = plot_data[qual_col].dropna()
		depth_vals = plot_data[depth_col].dropna()
		if len(qual_vals) > 0 and len(depth_vals) > 0:
			qual_cutoff = np.percentile(qual_vals, 99.5)
			depth_cutoff = np.percentile(depth_vals, 99.5)
			plot_data = plot_data[
				(plot_data[qual_col] <= qual_cutoff) & 
				(plot_data[depth_col] <= depth_cutoff)
			]
			cutoff_note = f" (smart cutoff applied)"
	
	set_plot_style()
	fig, ax = plt.subplots(figsize=(7, 5))
	sns.scatterplot(
		data=plot_data,
		x=depth_col,
		y=qual_col,
		s=10,
		alpha=0.5,
		edgecolor="none",
		color="#264653",
		ax=ax,
	)
	ax.set_xlabel("Mean depth")
	ax.set_ylabel("QUAL")
	ax.set_title(title + cutoff_note)
	fig.tight_layout()
	return save_figure(fig, output_path)


def plot_qd_distribution(
	qd: Union[pd.Series, np.ndarray, list],
	*,
	output_path: Optional[str] = None,
	title: str = "QD distribution (site level)",
	bins: int = 60,
	smart_cutoff: float = 99.5,
	enable_smart_cutoff: bool = True,
	use_smart_cutoff: bool = True,
	focus_pct: Optional[float] = None,
) -> Optional[plt.Figure]:
	"""Histogram of QD values."""
	values = np.asarray(qd, dtype=float)
	orig_mask = ~np.isnan(values)
	focus_note = ""
	if focus_pct is not None and orig_mask.any():
		pct_val = focus_pct * 100 if 0 < focus_pct <= 1 else focus_pct
		if 0 < pct_val < 100:
			cut = np.percentile(values[orig_mask], pct_val)
			values = values[values <= cut]
			focus_note = f" (<= {pct_val:.2f}th pct, cutoff={cut:g})"
	return hist_plot(
		values,
		output_path=output_path,
		title=title + focus_note,
		xlabel="QD",
		bins=bins,
		color="#E76F51",
		smart_cutoff=smart_cutoff,
		enable_smart_cutoff=enable_smart_cutoff and use_smart_cutoff,
	)

def plot_minor_allele_count_vs_depth(
	data: pd.DataFrame,
	*,
	output_path: Optional[str] = None,
	title: str = "Minor allele count vs depth",
	mac_col: str = "MAC",
	depth_col: str = "MeanDepth",
	smart_cutoff: float = 99.5,
	enable_smart_cutoff: bool = True,
	use_smart_cutoff: bool = True,
	apply_on: str = "both",
) -> Optional[plt.Figure]:
	"""4.13 Joint scatter of minor allele count vs site depth."""
	if not {mac_col, depth_col}.issubset(data.columns):
		raise ValueError("DataFrame must contain MAC and Depth columns")
	return joint_scatter(
		data.rename(columns={mac_col: "MAC", depth_col: "Depth"}),
		x="Depth",
		y="MAC",
		title=title,
		output_path=output_path,
		xlabel="Depth",
		ylabel="Minor allele count",
		color="#1D3557",
		smart_cutoff=smart_cutoff,
		enable_smart_cutoff=enable_smart_cutoff and use_smart_cutoff,
		apply_on=apply_on,
	)

def plot_minor_allele_freq_vs_depth(
	data: pd.DataFrame,
	*,
	output_path: Optional[str] = None,
	title: str = "Minor allele frequency vs depth",
	maf_col: str = "MAF",
	depth_col: str = "MeanDepth",
	smart_cutoff: float = 99.5,
	enable_smart_cutoff: bool = True,
	use_smart_cutoff: bool = True,
	apply_on: str = "both",
) -> Optional[plt.Figure]:
	"""4.14 Joint scatter of minor allele frequency vs site depth."""
	if not {maf_col, depth_col}.issubset(data.columns):
		raise ValueError("DataFrame must contain MAF and Depth columns")
	return joint_scatter(
		data.rename(columns={maf_col: "MAF", depth_col: "Depth"}),
		x="Depth",
		y="MAF",
		title=title,
		output_path=output_path,
		xlabel="Depth",
		ylabel="Minor allele frequency",
		color="#457B9D",
		smart_cutoff=smart_cutoff,
		enable_smart_cutoff=enable_smart_cutoff and use_smart_cutoff,
		apply_on=apply_on,
	)


def plot_site_missing_rate_distribution(
	missing_rates: Union[pd.Series, np.ndarray, list, Dict[str, float]],
	*,
	output_path: Optional[str] = None,
	title: str = "Missing rate distribution (site level)",
	bins: str = "auto",
	smart_cutoff: float = 99.9,
	enable_smart_cutoff: bool = True,
	use_smart_cutoff: bool = True,
) -> Optional[plt.Figure]:
	"""4.15 Histogram of site-level missing rates after genotype setting with fixed 0-1 range."""
	if isinstance(missing_rates, dict):
		values = list(missing_rates.values())
	else:
		values = missing_rates
	
	# Convert to numpy array and remove NaN values
	arr = np.asarray(values, dtype=float)
	mask = ~np.isnan(arr)
	orig_min = float(arr[mask].min()) if mask.any() else 0.0
	orig_max = float(arr[mask].max()) if mask.any() else 1.0
	filtered = arr[mask] if mask.any() else arr
	
	# Apply smart cutoff if enabled
	cut_phrase = ""
	if enable_smart_cutoff and use_smart_cutoff and mask.any() and 0 < smart_cutoff < 100:
		pct = max(0.0, min(100.0, smart_cutoff))
		current = filtered
		for _ in range(5):  # max iterations
			cur_mask = ~np.isnan(current)
			if not cur_mask.any():
				break
			cur_values = current[cur_mask]
			cur_max = float(cur_values.max())
			p99 = float(np.percentile(cur_values, pct))
			p50 = float(np.percentile(cur_values, 50))
			if (cur_max - p99) > (p99 - p50):
				current = current[current <= p99]
			else:
				break
		filtered = current
		if len(filtered) < len(arr[mask]):
			cut_phrase = f" (smart cutoff at {pct:.2f}% | original range:{orig_min:.3f}-{orig_max:.3f})"
		else:
			cut_phrase = f" (no cutoff | original range:{orig_min:.3f}-{orig_max:.3f})"
	elif mask.any():
		cut_phrase = f" (original range:{orig_min:.3f}-{orig_max:.3f})"
	
	# Create plot with fixed range
	set_plot_style()
	fig, ax = plt.subplots(figsize=(8, 5))
	
	# Plot histogram with forced auto bins
	sns.histplot(filtered, bins="auto", kde=True, color="#BC4B51", ax=ax)
	
	# Fix X-axis range to 0-1
	ax.set_xlim(0, 1)
	ax.set_xlabel("Missing rate")
	ax.set_ylabel("Count")
	
	# Set title with cutoff information
	if cut_phrase:
		ax.set_title(f"{title}\n{cut_phrase.strip()}")
	else:
		ax.set_title(title)
	
	fig.tight_layout()
	return save_figure(fig, output_path)


def plot_allele_count_pie(
	site_df: pd.DataFrame,
	*,
	output_path: Optional[str] = None,
	title: str = "Allele count composition (aggregated sites)",
) -> Optional[plt.Figure]:
	"""Pie chart of per-site distinct allele counts after aggregating normalized multi-allelic splits.

	If a multi-allelic site was normalized (e.g. split into several biallelic records with the same Chrom+Pos
	but different ALT), we reconstruct the total allele set by grouping Chrom+Pos and unifying ALT alleles.

	Categories:
	- Biallelic (2 total alleles = REF + 1 ALT)
	- Triallelic (3 total alleles)
	- ≥4 (four or more total alleles)

	Requires columns: Chrom, Pos, Alts (string ALT field possibly comma-separated or '.')
	"""
	required = {'Chrom','Pos','Alts'}
	if not required.issubset(site_df.columns):
		raise ValueError("DataFrame must contain Chrom, Pos, Alts columns for aggregation")
	grp = site_df[['Chrom','Pos','Alts']].copy()
	def collect_alts(rows: pd.Series) -> int:
		alts_set = set()
		for a in rows:
			if not isinstance(a,str) or a=='.':
				continue
			for tok in a.split(','):
				if tok != '.':
					alts_set.add(tok)
		# AlleleCount = 1 (REF) + number of distinct ALT tokens
		return 1 + len(alts_set)
	agg_counts = grp.groupby(['Chrom','Pos'], as_index=False)['Alts'].agg(collect_alts)
	if agg_counts.empty:
		print("No site records to aggregate for allele count pie.")
		return None
	counts = agg_counts['Alts']  # now holds reconstructed allele counts
	cat_series = counts.map(lambda c: 2 if c == 2 else (3 if c == 3 else (4 if c >= 4 else c)))
	labels_map = {2: 'Biallelic (2)', 3: 'Triallelic (3)', 4: '≥4'}
	value_counts = cat_series.map(lambda v: 4 if v >= 4 else v).map(labels_map).value_counts().reindex(labels_map.values(), fill_value=0)
	total_sites = int(value_counts.sum())
	set_plot_style()
	fig, ax = plt.subplots(figsize=(5,5))
	vals = value_counts.values
	labels = value_counts.index
	if total_sites == 0:
		print("AlleleCount categories empty after aggregation.")
		return None
	percent_labels = [f"{lab}: {v} ({v/total_sites*100:.2f}%)" for lab, v in zip(labels, vals)]
	ax.pie(vals, labels=percent_labels, autopct=None, startangle=90, counterclock=False)
	ax.set_title(f"{title}\nTotal sites (aggregated): {total_sites}")
	fig.tight_layout()
	return save_figure(fig, output_path)


def plot_mac_distribution(
	mac: Union[pd.Series, np.ndarray, list],
	*,
	output_path: Optional[str] = None,
	title: str = "Minor allele count (MAC) distribution",
	bins: int = 60,
	smart_cutoff: float = 99.5,
	enable_smart_cutoff: bool = True,
	use_smart_cutoff: bool = True,
	min_value: int = 0,
) -> Optional[plt.Figure]:
	"""Histogram of MAC values for filtered sites."""
	values = np.asarray(mac, dtype=float)
	
	# Apply minimum value filtering
	if min_value > 0:
		orig_count = len(values[~np.isnan(values)])
		values = values[values >= min_value]
		filtered_count = len(values[~np.isnan(values)])
		title_suffix = f" (MAC ≥ {min_value}, {filtered_count}/{orig_count} sites)"
		title = title + title_suffix
	
	return hist_plot(
		values,
		output_path=output_path,
		title=title,
		xlabel="MAC (second-most allele count)",
		bins=bins,
		color="#455A64",
		smart_cutoff=smart_cutoff,
		enable_smart_cutoff=enable_smart_cutoff and use_smart_cutoff,
	)


def plot_maf_distribution(
	maf: Union[pd.Series, np.ndarray, list],
	*,
	output_path: Optional[str] = None,
	title: str = "Minor allele frequency (MAF) distribution",
	bins: int = 60,
	smart_cutoff: float = 99.5,
	enable_smart_cutoff: bool = True,
	use_smart_cutoff: bool = True,
	min_value: float = 0.0,
) -> Optional[plt.Figure]:
	"""Histogram of MAF values for filtered sites."""
	values = np.asarray(maf, dtype=float)
	
	# Apply minimum value filtering
	if min_value > 0.0:
		orig_count = len(values[~np.isnan(values)])
		values = values[values >= min_value]
		filtered_count = len(values[~np.isnan(values)])
		title_suffix = f" (MAF ≥ {min_value:.3f}, {filtered_count}/{orig_count} sites)"
		title = title + title_suffix
	
	return hist_plot(
		values,
		output_path=output_path,
		title=title,
		xlabel="MAF (second-most allele frequency)",
		bins=bins,
		color="#00796B",
		smart_cutoff=smart_cutoff,
		enable_smart_cutoff=enable_smart_cutoff and use_smart_cutoff,
	)


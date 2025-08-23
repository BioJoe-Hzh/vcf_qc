"""Sample-level QC plotting functions.

Contains implementations for:
 4.1 Missing rate per sample (bar)
 4.2 Mean depth per sample (bar)
 4.2b Depth vs missing rate (scatter + marginal)
 4.3 Heterozygosity fraction per sample (bar)
 4.4 GQ distribution vs depth (boxplots binned by depth or per sample)

All functions follow the convention of returning a ``matplotlib.figure.Figure``
when ``output_path`` is not provided; otherwise they save and return ``None``.
"""

from __future__ import annotations

from typing import Optional, Union, Dict, Tuple, List

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from math import ceil, ceil as _ceil

try:  # optional dependency for clustering
	from sklearn.cluster import KMeans  # type: ignore
	from sklearn.metrics import silhouette_score  # type: ignore
except Exception:  # pragma: no cover - graceful degradation
	KMeans = None  # type: ignore
	silhouette_score = None  # type: ignore

from .base import set_plot_style, save_figure, joint_scatter

__all__ = [
	"plot_missing_rate_per_sample",
	"plot_depth_per_sample",
	"plot_depth_vs_missing_rate",
	"plot_het_fraction_per_sample",
	"plot_gq_boxplot_vs_depth",
	"plot_het_ratio_vs_depth",
]


def _dict_series_frame_to_df(
	data: Union[Dict[str, float], pd.Series, pd.DataFrame],
	value_col: str,
	sample_col: str = "Sample",
) -> pd.DataFrame:
	"""Normalise different input shapes into a DataFrame.

	Accepts dict, Series or DataFrame. If a DataFrame is passed and
	already contains both columns it is returned unchanged.
	"""
	if isinstance(data, dict):
		return pd.DataFrame({sample_col: list(data.keys()), value_col: list(data.values())})
	if isinstance(data, pd.Series):
		return pd.DataFrame({sample_col: data.index.to_list(), value_col: data.values})
	if isinstance(data, pd.DataFrame):
		if {sample_col, value_col}.issubset(data.columns):
			return data[[sample_col, value_col]].copy()
		raise ValueError(f"DataFrame must contain columns: {sample_col}, {value_col}")
	raise TypeError("Input must be dict | Series | DataFrame")


def _panel_bar_plot(
	df: pd.DataFrame,
	value_col: str,
	title: str,
	ylabel: str,
	color_map: Optional[Dict[str, str]] = None,
	*,
	samples_per_panel: int = 100,
	base_color: str = "#4477AA",
	output_path: Optional[str] = None,
	rotation: int = 45,
	global_ylim: Optional[Tuple[float, float]] = None,
	custom_ticks: Optional[List[float]] = None,
	per_panel_scale: bool = False,
	cluster_mode: bool = False,
) -> Optional[plt.Figure]:
	"""Generic multi-panel bar plot for per-sample metrics.

	Each panel keeps identical height (~5 inches) to mirror original single-panel size.
	If ``global_ylim`` provided, all panels share same y-limits and ticks.
	If ``per_panel_scale`` True, each panel scales y-axis individually (overrides global).
	Pads final panel with blank samples for width consistency.
	"""
	set_plot_style()
	n = len(df)
	panels = ceil(n / samples_per_panel) if n else 1
	fig_width = max(10, min(18, samples_per_panel * 0.18))
	fig, axes = plt.subplots(
		panels, 1, figsize=(fig_width, panels * 5), sharey=not per_panel_scale and global_ylim is not None
	)
	if panels == 1:
		axes = [axes]  # type: ignore

	for pi in range(panels):
		start = pi * samples_per_panel
		end = start + samples_per_panel
		sub = df.iloc[start:end].copy()
		pad_needed = samples_per_panel - len(sub)
		if pad_needed > 0:
			sub = pd.concat([
				sub,
				pd.DataFrame({
					"Sample": ["" for _ in range(pad_needed)],
					value_col: [0 for _ in range(pad_needed)],
					"_cluster": [None for _ in range(pad_needed)],
				})
			], ignore_index=True)
		ax = axes[pi]
		if cluster_mode and color_map and "_cluster" in sub.columns:
			# add string hue with padding class
			sub["_cluster_str"] = sub["_cluster"].apply(lambda x: f"c{x}" if x is not None else "pad")
			palette_dict = {f"c{k}": v for k, v in {int(k): v for k, v in (color_map.items())}.items() if str(k).isdigit()}  # type: ignore
			palette_dict.update({k: v for k, v in color_map.items() if not str(k).isdigit()})
			if "pad" not in palette_dict:
				palette_dict["pad"] = "#DDDDDD"
			sns.barplot(
				data=sub,
				x="Sample",
				y=value_col,
				hue="_cluster_str",
				palette=palette_dict,
				dodge=False,
				ax=ax,
				legend=False,
			)
		else:
			sns.barplot(data=sub, x="Sample", y=value_col, ax=ax, color=base_color, legend=False)
		ax.set_xlabel("Sample")
		if pi == 0:
			ax.set_title(title)
		ax.set_ylabel(ylabel if pi == 0 else "")
		for label in ax.get_xticklabels():
			label.set_rotation(rotation)
			label.set_ha("right")
		# Scaling
		if per_panel_scale:
			panel_max = sub[value_col].max()
			if panel_max == 0:
				panel_max = 1
			if value_col.lower().startswith("het") or "Missing" in ylabel or value_col.lower().startswith("missing"):
				y_max = 1.0
				ticks = [0.0, 0.25, 0.5, 0.75, 1.0]
			else:
				# round up to nearest multiple of 5
				y_max = float(_ceil(panel_max / 5.0) * 5)
				ticks = list(np.linspace(0, y_max, 5))
			ax.set_ylim(0, y_max)
			ax.set_yticks(ticks)
		elif global_ylim is not None:
			ax.set_ylim(*global_ylim)
			if custom_ticks is not None:
				ax.set_yticks(custom_ticks)
		sns.despine(ax=ax)
	fig.tight_layout(h_pad=0.5)
	return save_figure(fig, output_path)


def plot_missing_rate_per_sample(
	missing_rates: Union[Dict[str, float], pd.Series, pd.DataFrame],
	*,
	output_path: Optional[str] = None,
	title: str = "Missing rate per sample",
	base_color: str = "#4477AA",
	samples_per_panel: int = 100,
) -> Optional[plt.Figure]:
	"""4.1 Multi-panel missing rate (fixed y 0..1, 5 ticks)."""
	df = _dict_series_frame_to_df(missing_rates, "MissingRate")
	df = df.sort_values("MissingRate", ascending=False)
	return _panel_bar_plot(
		df,
		"MissingRate",
		title,
		"Missing rate",
		None,
		samples_per_panel=samples_per_panel,
		base_color=base_color,
		output_path=output_path,
		global_ylim=(0, 1),
		custom_ticks=[0.0, 0.25, 0.5, 0.75, 1.0],
	)


def _cluster_depths(values: np.ndarray, max_k: int = 10) -> Tuple[int, List[int]]:
	"""Cluster depths using heuristic elbow + silhouette.

	Returns (k, labels). Falls back to single cluster if sklearn absent.
	"""
	if KMeans is None or silhouette_score is None:
		return 1, [0] * len(values)
	if len(values) < 2:
		return 1, [0] * len(values)
	max_k = min(max_k, len(values) - 1)
	if max_k < 2:
		return 1, [0] * len(values)
	X = values.reshape(-1, 1)
	wcss = []
	sil = []
	k_range = range(1, max_k + 1)
	for k in k_range:
		if k == 1:
			wcss.append(np.sum((X - np.mean(X)) ** 2))
			sil.append(0)
		else:
			km = KMeans(n_clusters=k, random_state=42, n_init=10)
			labels = km.fit_predict(X)
			wcss.append(km.inertia_)
			try:
				sil.append(silhouette_score(X, labels))
			except Exception:
				sil.append(0)
	if len(wcss) >= 3:
		deltas = np.diff(wcss)
		second = np.diff(deltas)
		elbow_idx = (np.argmax(second) + 2) if len(second) else 2
	else:
		elbow_idx = 2
	if len(sil) > 2:
		candidates = range(max(2, elbow_idx - 1), min(len(sil), elbow_idx + 3))
		best_k = max(candidates, key=lambda k: sil[k - 1])
	else:
		best_k = 2 if len(values) > 2 else 1
	if best_k == 1:
		return 1, [0] * len(values)
	km = KMeans(n_clusters=best_k, random_state=42, n_init=10)
	labels = km.fit_predict(X).tolist()
	return best_k, labels


def plot_depth_per_sample(
	depths: Union[Dict[str, float], pd.Series, pd.DataFrame],
	*,
	output_path: Optional[str] = None,
	title: str = "Mean depth per sample (clustered)",
	samples_per_panel: int = 100,
	cluster_max_k: int = 10,
	palette: Optional[List[str]] = None,
	cluster_k: Optional[int] = None,
	cluster_results_path: Optional[str] = None,
) -> Optional[plt.Figure]:
	"""4.2 Depth multi-panel with clustering & global y-scale (rounded to 5)."""
	df = _dict_series_frame_to_df(depths, "Depth")
	df = df.sort_values("Depth", ascending=False).reset_index(drop=True)
	values = df["Depth"].to_numpy(dtype=float)
	if cluster_k is not None and cluster_k > 0:
		# force a specific number of clusters
		if KMeans is not None and len(values) >= cluster_k:
			km = KMeans(n_clusters=cluster_k, random_state=42, n_init=10)
			labels = km.fit_predict(values.reshape(-1, 1)).tolist()
			k = cluster_k
		else:
			k, labels = 1, [0] * len(values)
	else:
		k, labels = _cluster_depths(values, max_k=cluster_max_k)
	df["_cluster"] = labels
	if palette is None:
		palette = sns.color_palette("tab10", n_colors=k)
	color_map = {str(c): palette[c % len(palette)] for c in sorted(set(labels))}
	if cluster_results_path:
		export_df = df[["Sample", "Depth", "_cluster"]].rename(columns={"_cluster": "Cluster"})
		export_df.to_csv(cluster_results_path, sep='\t', index=False)
	global_max = df["Depth"].max()
	if global_max == 0:
		global_max = 1
	y_top = float(_ceil(global_max / 5.0) * 5)
	ticks = list(np.linspace(0, y_top, 5))
	return _panel_bar_plot(
		df,
		"Depth",
		title + (f" (k={k})" if k else ""),
		"Depth",
		color_map,
		samples_per_panel=samples_per_panel,
		base_color="#228833",
		output_path=output_path,
		global_ylim=(0, y_top),
	custom_ticks=ticks,
	cluster_mode=True,
	)


def plot_depth_vs_missing_rate(
	data: pd.DataFrame,
	*,
	output_path: Optional[str] = None,
	title: str = "Depth vs missing rate",
	depth_col: str = "Depth",
	miss_col: str = "MissingRate",
	sample_col: str = "Sample",
) -> Optional[plt.Figure]:
	"""4.2 (extended) Scatter with marginal distributions: depth vs missing rate.

	Expects columns named via parameters. ``sample_col`` is not plotted
	directly but required for uniqueness checks if needed later.
	"""
	required = {depth_col, miss_col}
	if not required.issubset(data.columns):
		raise ValueError(f"DataFrame must contain columns: {required}")
	return joint_scatter(
		data.rename(columns={depth_col: "Depth", miss_col: "MissingRate"}),
		x="Depth",
		y="MissingRate",
		title=title,
		output_path=output_path,
		xlabel="Depth",
		ylabel="Missing rate",
		color="#5B3F95",
	)


def plot_het_ratio_vs_depth(
	data: pd.DataFrame,
	*,
	output_path: Optional[str] = None,
	title: str = "HetRatio vs depth",
	depth_col: str = "MeanDepth",
	het_ratio_col: str = "HetRatio",
	sample_col: str = "Sample",
	color: str = "#BF3985",
) -> Optional[plt.Figure]:
	"""Scatter (with marginal histograms) of HetRatio vs MeanDepth.

	Filters out NaN / infinite HetRatio. Automatically adds Pearson r to title.
	"""
	required = {depth_col, het_ratio_col}
	if not required.issubset(data.columns):
		raise ValueError(f"DataFrame must contain columns: {required}")
	df = data[[depth_col, het_ratio_col]].copy()
	df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=[het_ratio_col])
	pearson_r = None
	if len(df) > 1:
		try:
			pearson_r = df[depth_col].corr(df[het_ratio_col])
		except Exception:
			pearson_r = None
	final_title = title + (f" (r={pearson_r:.3f})" if pearson_r is not None else "")
	return joint_scatter(
		df.rename(columns={depth_col: "Depth", het_ratio_col: "HetRatio"}),
		x="Depth",
		y="HetRatio",
		title=final_title,
		output_path=output_path,
		xlabel="Depth",
		ylabel="HetRatio (Het/Hom-Alt)",
		color=color,
	)


def plot_het_fraction_per_sample(
	het_rates: Union[Dict[str, float], pd.Series, pd.DataFrame],
	*,
	output_path: Optional[str] = None,
	title: str = "Heterozygosity ratio (Het / Hom-Alt)",
	cluster_max_k: int = 8,
	palette: Optional[List[str]] = None,
	cluster_k: Optional[int] = None,
	cluster_results_path: Optional[str] = None,
) -> Optional[plt.Figure]:
	"""4.3 Heterozygosity fraction: ONE CLUSTER PER PANEL + per-panel max scaling.

	Changes per request:
	- Different clusters shown in different panels (no mixing).
	- Each panel's y-axis (0 -> panel max rounded) determined by that cluster's maximum value.
	- Sample label rotation set to +45 degrees.
	"""
	# Expect HetRatio column if DataFrame provided
	if isinstance(het_rates, pd.DataFrame):
		if "HetRatio" not in het_rates.columns:
			raise ValueError("DataFrame must contain HetRatio column")
		df = het_rates[["Sample", "HetRatio"]].rename(columns={"HetRatio": "HetMetric"}).copy()
	else:
		# dict/Series path; assume provided values are ratios
		df = _dict_series_frame_to_df(het_rates, "HetMetric")
	value_label = "HetMetric"
	values = df[value_label].to_numpy(dtype=float)
	# Handle NaNs: cluster only non-NaN values; NaNs assigned to separate cluster at end
	nan_mask = np.isnan(values)
	non_nan_values = values[~nan_mask]
	labels_full = np.full(len(values), -1, dtype=int)
	if len(non_nan_values) == 0:
		# all NaN -> single synthetic cluster 0
		labels_full[:] = 0
		k = 1
	else:
		if cluster_k is not None and cluster_k > 0:
			if KMeans is not None and len(non_nan_values) >= cluster_k:
				km = KMeans(n_clusters=cluster_k, random_state=42, n_init=10)
				labels_non_nan = km.fit_predict(non_nan_values.reshape(-1, 1)).tolist()
				k = cluster_k
			else:
				k, labels_non_nan = 1, [0] * len(non_nan_values)
		else:
			k, labels_non_nan = _cluster_depths(non_nan_values, max_k=cluster_max_k)
		labels_full[~nan_mask] = labels_non_nan
		if nan_mask.any():
			# assign NaNs to new cluster index k (separate panel)
			labels_full[nan_mask] = k
			k = k + 1
	df["_cluster"] = labels_full.tolist()
	# For plotting, replace NaN metric with 0 so bars render; keep original for export
	df["_HetOriginal"] = df[value_label]
	df[value_label] = df[value_label].fillna(0.0)
	# order clusters by mean descending
	cluster_order = df.groupby("_cluster")[value_label].mean().sort_values(ascending=False).index.tolist()
	if palette is None:
		palette = sns.color_palette("tab10", n_colors=k)
	color_map = {c: palette[i % len(palette)] for i, c in enumerate(cluster_order)}
	if cluster_results_path:
		export_cols = df[["Sample", "_HetOriginal", "_cluster"]].rename(columns={"_HetOriginal": "HetRatio", "_cluster": "Cluster"})
		export_cols.to_csv(cluster_results_path, sep='\t', index=False)
	set_plot_style()
	# Split clusters into sub-panels of at most 100 samples each
	max_per_panel = 100
	panel_specs = []  # list of (cluster_id, part_index, total_parts, dataframe)
	for c in cluster_order:
		cdf_full = df[df._cluster == c].sort_values(value_label, ascending=False).reset_index(drop=True)
		if cdf_full.empty:
			continue
		total = len(cdf_full)
		parts = int(np.ceil(total / max_per_panel))
		for part in range(parts):
			start = part * max_per_panel
			end = start + max_per_panel
			panel_specs.append((c, part + 1, parts, total, cdf_full.iloc[start:end].copy()))
	if not panel_specs:
		return None
	# Pre-compute per-cluster y scaling (use full cluster data, not per part)
	cluster_ymax: Dict[int, float] = {}
	for c in cluster_order:
		full_cluster_vals = df[df._cluster == c][value_label]
		panel_max = full_cluster_vals.max() if len(full_cluster_vals) else 0
		if panel_max <= 0:
			panel_max = 1
		cluster_ymax[c] = float(panel_max * 1.05)
	fig_height = 4 * len(panel_specs)
	max_panel_width_samples = max(len(spec[4]) for spec in panel_specs)
	fig_width = max(10, min(18, max_panel_width_samples * 0.2))
	fig, axes = plt.subplots(len(panel_specs), 1, figsize=(fig_width, fig_height), squeeze=False)
	for idx, (c, part_idx, parts_total, cluster_total, cdf) in enumerate(panel_specs):
		ax = axes[idx, 0]
		cluster_title_base = f"Cluster {c} (n={cluster_total})"
		if cdf[value_label].eq(0).all() and cdf['_HetOriginal'].isna().all():
			cluster_title_base = f"Undefined HetRatio (NaN) (n={cluster_total})"
		if parts_total > 1:
			cluster_title = f"{cluster_title_base} part {part_idx}/{parts_total}"
		else:
			cluster_title = cluster_title_base
		sns.barplot(data=cdf, x="Sample", y=value_label, ax=ax, color=color_map.get(c, '#888888'))
		ax.set_title(f"{title} - {cluster_title}" if idx == 0 else cluster_title)
		ax.set_xlabel("Sample")
		ax.set_ylabel("HetRatio (Het / Hom-Alt)" if idx == 0 else "")
		y_max = cluster_ymax.get(c, 1.0)
		ticks = list(np.linspace(0, y_max, 5))
		ax.set_ylim(0, y_max)
		ax.set_yticks(ticks)
		for lab in ax.get_xticklabels():
			lab.set_rotation(45)
			lab.set_ha("right")
		sns.despine(ax=ax)
	fig.tight_layout(h_pad=0.6)
	return save_figure(fig, output_path)


def plot_gq_boxplot_vs_depth(
	data: pd.DataFrame,
	*,
	output_path: Optional[str] = None,
	title: str = "GQ distribution across sample mean depth bins",
	depth_col: str = "Depth",
	gq_col: str = "GQ",
	figsize: Tuple[int, int] = (12, 5),
	violin_color: str = "#8888CC",
	bar_color: str = "#444444",
) -> Optional[plt.Figure]:
	"""4.4 Violin plot of genotype GQ grouped by *sample mean* depth bins.

	Fixes over-counting by assigning each sample to exactly one depth bin based on
	its mean depth (across genotypes). Secondary bar axis shows number of samples per bin.
	Grid lines removed for cleaner view.
	"""
	for col in (depth_col, gq_col):
		if col not in data.columns:
			raise ValueError(f"Column '{col}' not in DataFrame")
	set_plot_style()
	sns.set_style("white")
	# Compute mean depth per sample then map sample -> bin
	if 'Sample' in data.columns:
		sample_mean_depth = (
			data.loc[data[depth_col] > 0].groupby('Sample')[depth_col].mean().to_dict()
		)
		data = data.copy()
		data['SampleMeanDepth'] = data['Sample'].map(sample_mean_depth).fillna(0)
		bin_values = data['SampleMeanDepth']
	else:
		bin_values = data[depth_col]
	bins = [0, 1, 3, 5, 10, 15, 20, 25, 30, 35, float('inf')]
	labels = ["0-1", "1-3", "3-5", "5-10", "10-15", "15-20", "20-25", "25-30", "30-35", ">35"]
	bidx = np.digitize(bin_values, bins, right=False) - 1
	data['DepthBin'] = [labels[i] if 0 <= i < len(labels) else labels[-1] for i in bidx]
	order = labels
	fig, ax = plt.subplots(figsize=figsize)
	sns.violinplot(
		data=data,
		x='DepthBin',
		y=gq_col,
		order=order,
		ax=ax,
		inner='box',
		cut=0,
		color=violin_color,
	)
	ax.set_xlabel("Sample mean depth bin")
	ax.set_ylabel("GQ")
	ax.set_title(title)
	# Counts of samples per bin
	if 'Sample' in data.columns:
		sample_bin = data[['Sample', 'DepthBin']].drop_duplicates()
		counts = sample_bin.groupby('DepthBin')['Sample'].nunique().reindex(order).fillna(0).astype(int)
	else:
		counts = data.groupby('DepthBin')[gq_col].count().reindex(order).fillna(0).astype(int)
	ax2 = ax.twinx()
	ax2.bar(order, counts, alpha=0.25, color=bar_color)
	ax2.set_ylabel("Sample count")
	# add count labels above bars
	max_count = counts.max() if len(counts) else 0
	if max_count > 0:
		ax2.set_ylim(0, max_count * 1.15)
	for i, (lab, cnt) in enumerate(zip(order, counts)):
		ax2.text(i, cnt + max_count * 0.02 if max_count else cnt + 0.1, str(int(cnt)), ha='center', va='bottom', fontsize=8)
	for label in ax.get_xticklabels():
		label.set_rotation(-45)
		label.set_ha("right")
	# remove grids
	ax.grid(False)
	ax2.grid(False)
	fig.tight_layout()
	return save_figure(fig, output_path)


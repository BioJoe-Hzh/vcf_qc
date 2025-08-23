"""Base plotting utilities shared across QC plot modules.

This module centralises style configuration and small helper wrappers
around seaborn/matplotlib so higher‑level plot functions remain concise
and consistent. Each helper returns a matplotlib Figure when an
``output_path`` is not provided; otherwise the figure is saved and
closed (to avoid memory accumulation in batch runs) and ``None`` is
returned.
"""

from __future__ import annotations

from typing import Optional, Dict, Tuple, Union

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

__all__ = [
	"set_plot_style",
	"save_figure",
	"joint_scatter",
	"hist_plot",
]


def set_plot_style() -> None:
	"""Apply a unified visual style.

	Centralised so we can later expose style choices via configuration.
	"""
	sns.set_theme(style="whitegrid")
	plt.rcParams.update({
		"axes.titlesize": 13,
		"axes.labelsize": 11,
		"font.size": 10,
		"figure.dpi": 100,
	})


def save_figure(fig: plt.Figure, output_path: Optional[str]) -> Optional[plt.Figure]:
	"""Save figure if ``output_path`` provided else return it.

	Parameters
	----------
	fig : matplotlib.figure.Figure
		Figure to save or return.
	output_path : str | None
		Path to save. If None the figure is returned and *not* closed.
	"""
	if output_path:
		fig.savefig(output_path, bbox_inches="tight")
		plt.close(fig)
		return None
	return fig


def joint_scatter(
	data: pd.DataFrame,
	x: str,
	y: str,
	*,
	output_path: Optional[str] = None,
	title: str = "",
	xlabel: Optional[str] = None,
	ylabel: Optional[str] = None,
	color: str = "steelblue",
	height: int = 6,
	alpha: float = 0.5,
	bins: int = 50,
	kde: bool = True,
	smart_cutoff: float = 99.5,
	enable_smart_cutoff: bool = True,
	apply_on: str = "both",
	smart_cutoff_max_iter: int = 5,
	hue: Optional[str] = None,
	palette: Optional[Dict[str, str]] = None,
) -> Optional[plt.Figure]:
	"""Create a scatter plot with marginal histograms.

	Uses ``seaborn.jointplot`` so we get separate axes for marginal
	distributions. For very large points (>50k) consider down‑sampling
	before calling this helper.
	"""
	set_plot_style()
	# Retain hue column if provided so colour mapping works downstream
	cols = [x, y]
	if hue and hue in data.columns and hue not in cols:
		cols.append(hue)
	df = data[cols].copy()
	orig_xmin = float(df[x].min()) if not df.empty else 0.0
	orig_xmax = float(df[x].max()) if not df.empty else 0.0
	orig_ymin = float(df[y].min()) if not df.empty else 0.0
	orig_ymax = float(df[y].max()) if not df.empty else 0.0
	cut_note = ""
	if enable_smart_cutoff and 0 < smart_cutoff < 100 and not df.empty:
		pct = max(0.0, min(100.0, smart_cutoff))
		iters_x = 0
		iters_y = 0
		# Iteratively trim X axis first if requested
		if apply_on in ("x", "both") and df[x].notna().any():
			while iters_x < smart_cutoff_max_iter and df[x].notna().any():
				cur_max = float(df[x].max())
				p99 = float(np.percentile(df[x], pct))
				p50 = float(np.percentile(df[x], 50))
				if (cur_max - p99) > (p99 - p50):
					df = df[df[x] <= p99]
					iters_x += 1
				else:
					break
		# Then iteratively trim Y axis (after any X filtering) if requested
		if apply_on in ("y", "both") and df[y].notna().any():
			while iters_y < smart_cutoff_max_iter and df[y].notna().any():
				cur_ymax2 = float(df[y].max())
				p99y = float(np.percentile(df[y], pct))
				p50y = float(np.percentile(df[y], 50))
				if (cur_ymax2 - p99y) > (p99y - p50y):
					df = df[df[y] <= p99y]
					iters_y += 1
				else:
					break
		applied_flags = []
		if iters_x > 0:
			applied_flags.append(f"X axis ({iters_x} iter)")
		if iters_y > 0:
			applied_flags.append(f"Y axis ({iters_y} iter)")
		if applied_flags:
			axis_phrase = ", ".join(applied_flags)
			cut_note = f" (smart cutoff {pct:.2f}% on {axis_phrase} | original X:{orig_xmin:g}-{orig_xmax:g} Y:{orig_ymin:g}-{orig_ymax:g})"
		else:
			cut_note = f" (no cutoff | original X:{orig_xmin:g}-{orig_xmax:g} Y:{orig_ymin:g}-{orig_ymax:g})"
	joint_kws = {
		"data": df,
		"x": x,
		"y": y,
		"kind": "scatter",
		"height": height,
	}
	if hue and hue in df.columns:
		# Custom grid so we can control marginals without conflicting kwargs
		g = sns.JointGrid(data=df, x=x, y=y, height=height)
		sns.scatterplot(data=df, x=x, y=y, hue=hue, palette=palette, alpha=alpha, edgecolor="none", ax=g.ax_joint)
		# Marginal for x
		sns.histplot(x=df[x], bins=bins, kde=False, color="#78909C", ax=g.ax_marg_x)
		# Marginal for y: use y variable along vertical axis, counts on horizontal
		sns.histplot(y=df[y], bins=bins, kde=False, color="#78909C", ax=g.ax_marg_y)
		g.ax_marg_x.set_ylabel("")
		g.ax_marg_y.set_xlabel("")
		# Tighten y-range to actual data range (prevents count max becoming y-limit)
		if df[y].notna().any():
			g.ax_joint.set_ylim(df[y].min(), df[y].max())
			g.ax_marg_y.set_ylim(df[y].min(), df[y].max())
		# Ensure legend is present
		leg = g.ax_joint.get_legend()
		if leg is not None:
			leg.set_title(hue)
	else:
		joint_kws["color"] = color
		joint_kws["marginal_kws"] = {"bins": bins, "kde": kde, "color": color}
		# Pass alpha via joint_kws (seaborn forwards to underlying scatter)
		joint_kws["joint_kws"] = {"alpha": alpha}
		g = sns.jointplot(**joint_kws)
	if cut_note:
		g.fig.suptitle(f"{title}\n{cut_note.strip()}")
	else:
		g.fig.suptitle(title)
	g.set_axis_labels(xlabel or x, ylabel or y)
	g.fig.tight_layout()
	g.fig.subplots_adjust(top=0.92)
	return save_figure(g.fig, output_path)


def hist_plot(
	values: Union[pd.Series, np.ndarray, list],
	*,
	output_path: Optional[str] = None,
	title: str = "",
	xlabel: str = "",
	bins: int = 50,
	color: str = "steelblue",
	kde: bool = True,
	logx: bool = False,
	figsize: Tuple[int, int] = (8, 5),
	smart_cutoff: float = 99.5,
	enable_smart_cutoff: bool = True,
	smart_cutoff_max_iter: int = 5,
) -> Optional[plt.Figure]:
	"""Histogram + (optional) KDE.

	Parameters are intentionally kept minimal; callers can perform any
	value filtering / transformation before passing the array.
	"""
	set_plot_style()
	arr = np.asarray(values, dtype=float)
	mask = ~np.isnan(arr)
	orig_min = float(arr[mask].min()) if mask.any() else 0.0
	orig_max = float(arr[mask].max()) if mask.any() else 0.0
	filtered = arr
	cut_phrase = ""
	if enable_smart_cutoff and mask.any() and 0 < smart_cutoff < 100:
		pct = max(0.0, min(100.0, smart_cutoff))
		iters = 0
		current = filtered
		while iters < smart_cutoff_max_iter:
			cur_mask = ~np.isnan(current)
			if not cur_mask.any():
				break
			cur_values = current[cur_mask]
			cur_max = float(cur_values.max())
			p99 = float(np.percentile(cur_values, pct))
			p50 = float(np.percentile(cur_values, 50))
			if (cur_max - p99) > (p99 - p50):
				current = current[current <= p99]
				iters += 1
			else:
				break
		filtered = current
		if iters > 0:
			cut_phrase = f" (smart cutoff at {pct:.2f}% iter={iters} | original range:{orig_min:g}-{orig_max:g})"
		else:
			cut_phrase = f" (no cutoff | original range:{orig_min:g}-{orig_max:g})"
	elif mask.any():
		cut_phrase = f" (original range:{orig_min:g}-{orig_max:g})"
	fig, ax = plt.subplots(figsize=figsize)
	sns.histplot(filtered, bins=bins, kde=kde, color=color, ax=ax)
	if logx:
		ax.set_xscale("log")
	if cut_phrase:
		ax.set_title(f"{title}\n{cut_phrase.strip()}")
	else:
		ax.set_title(title)
	ax.set_xlabel(xlabel)
	ax.set_ylabel("Count")
	fig.tight_layout()
	return save_figure(fig, output_path)


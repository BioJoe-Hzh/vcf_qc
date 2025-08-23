"""Site-level metric assembly.

Converts the generator from ``SimpleVCFReader.iterate_site_metrics`` into a
DataFrame and offers light filtering helpers.
"""
from __future__ import annotations

from typing import Optional, Iterable, Dict, Any
import pandas as pd

from ..io import SimpleVCFReader

__all__ = ["compute_site_metrics", "filter_site_metrics"]


def compute_site_metrics(reader: SimpleVCFReader, limit: Optional[int] = None) -> pd.DataFrame:
    """Return DataFrame with columns: Chrom, Pos, QUAL, MeanDepth, QD, MAC, MAF, MissingRate."""
    rows = []
    for i, rec in enumerate(reader.iterate_site_metrics()):  # type: ignore
        rows.append(rec)
        if limit and i + 1 >= limit:
            break
    df = pd.DataFrame(rows)
    # Ensure numeric types where possible
    for col in ["QUAL", "MeanDepth", "QD", "MAC", "MAF", "MissingRate", "AlleleCount"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    # Silent by default; only warn if unexpected MAF (>0.5) are observed (should be rare with second-most allele definition)
    if "MAF" in df.columns:
        over_mask = df["MAF"] > 0.5
        count_over = int(over_mask.sum())
        if count_over > 0:
            print(f"[WARNING] Detected {count_over} sites with MAF > 0.5 (expected 0 under second-most allele definition). Showing up to 10:")
            print(df.loc[over_mask, ["Chrom", "Pos", "MAC", "MAF"]].head(10).to_string(index=False))
    return df


def filter_site_metrics(
    df: pd.DataFrame,
    *,
    min_qual: Optional[float] = None,
    max_missing: Optional[float] = None,
    min_depth: Optional[int] = None,
    min_mac: Optional[int] = None,
) -> pd.DataFrame:
    """Apply simple QC thresholds; returns filtered copy (non-destructive)."""
    out = df.copy()
    if min_qual is not None and "QUAL" in out.columns:
        out = out[out["QUAL"] >= min_qual]
    if max_missing is not None and "MissingRate" in out.columns:
        out = out[out["MissingRate"] <= max_missing]
    if min_depth is not None and "MeanDepth" in out.columns:
        out = out[out["MeanDepth"] >= min_depth]
    if min_mac is not None and "MAC" in out.columns:
        out = out[out["MAC"] >= min_mac]
    return out.reset_index(drop=True)

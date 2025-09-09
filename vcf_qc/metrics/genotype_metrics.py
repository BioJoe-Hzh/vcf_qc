"""Genotype-level metric assembly.

Converts the generator from SimpleVCFReader.iterate_genotype_* into DataFrames
and offers filtering helpers for genotype-level quality control.
"""
from __future__ import annotations

from typing import Optional, Iterable, Dict, Any
import pandas as pd

from ..io import SimpleVCFReader

__all__ = ["genotype_complete_table", "genotype_het_vaf_table", "filter_heterozygous"]


def genotype_complete_table(reader: SimpleVCFReader, compute_gq_from_pl: bool = False) -> pd.DataFrame:
    """Return DataFrame with complete genotype information.
    
    Columns: Sample, Chrom, Pos, GT, GQ, Depth, ALT_Count, REF_Count, VAF
    """
    rows = []
    for rec in reader.iterate_gq_depth(compute_gq_from_pl=compute_gq_from_pl):
        rows.append(rec)
    
    df = pd.DataFrame(rows)
    
    # Ensure numeric types where possible
    for col in ["GQ", "Depth", "ALT_Count", "REF_Count", "VAF"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    
    return df


def genotype_het_vaf_table(reader: SimpleVCFReader, compute_gq_from_pl: bool = False) -> pd.DataFrame:
    """Return DataFrame with heterozygous genotypes only for VAF analysis.
    
    Filters to only include heterozygous genotypes (0/1, 1/0, 0|1, 1|0).
    Columns: Sample, Chrom, Pos, GT, GQ, Depth, ALT_Count, REF_Count, VAF
    """
    complete_df = genotype_complete_table(reader, compute_gq_from_pl)
    return filter_heterozygous(complete_df)


def filter_heterozygous(df: pd.DataFrame) -> pd.DataFrame:
    """Filter DataFrame to include only heterozygous genotypes.
    
    Parameters
    ----------
    df : pd.DataFrame
        Genotype data with GT column
        
    Returns
    -------
    pd.DataFrame
        Filtered DataFrame containing only heterozygous genotypes
    """
    if 'GT' not in df.columns:
        raise ValueError("DataFrame must contain 'GT' column for heterozygous filtering")
    
    # Define heterozygous genotype patterns
    het_genotypes = ['0/1', '1/0', '0|1', '1|0']
    
    # Filter for heterozygous genotypes
    het_mask = df['GT'].isin(het_genotypes)
    het_df = df[het_mask].copy()
    
    # Remove rows with missing VAF data for meaningful analysis
    if 'VAF' in het_df.columns:
        het_df = het_df.dropna(subset=['VAF'])
    
    return het_df.reset_index(drop=True)
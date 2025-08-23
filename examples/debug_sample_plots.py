"""Debug script to generate sample-level QC plots from a VCF.

Usage (adjust path):
    PYTHONPATH=.. python3 examples/debug_sample_plots.py \
        --vcf /Users/zh384/Desktop/temp/Myrio190.NT.s0.100k.withGQ.vcf.gz \
        --outdir sample_qc_plots
"""
from __future__ import annotations

import argparse
from pathlib import Path

from vcf_qc.io.vcf_reader import SimpleVCFReader
from vcf_qc.metrics.sample_metrics import compute_sample_metrics, gq_depth_long_table
from vcf_qc.plot.sample_plots import (
    plot_missing_rate_per_sample,
    plot_depth_per_sample,
    plot_depth_vs_missing_rate,
    plot_het_fraction_per_sample,
    plot_gq_boxplot_vs_depth,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, help="Input VCF(.gz)")
    ap.add_argument("--outdir", required=True, help="Output directory for plots")
    ap.add_argument("--max-records", type=int, default=None, help="Limit records (debug)")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    reader = SimpleVCFReader(args.vcf, max_records=args.max_records)
    metrics_df = compute_sample_metrics(reader)

    # Basic per-sample plots
    plot_missing_rate_per_sample(
        dict(zip(metrics_df.Sample, metrics_df.MissingRate)),
        output_path=str(outdir / "sample_missing_rate.png"),
    )
    plot_depth_per_sample(
        dict(zip(metrics_df.Sample, metrics_df.MeanDepth)),
        output_path=str(outdir / "sample_mean_depth.png"),
    )
    plot_het_fraction_per_sample(
        dict(zip(metrics_df.Sample, metrics_df.HetFraction)),
        output_path=str(outdir / "sample_het_fraction.png"),
    )

    # Depth vs missing rate joint scatter
    plot_depth_vs_missing_rate(
        metrics_df.rename(columns={"MeanDepth": "Depth"}),
        output_path=str(outdir / "depth_vs_missing_rate.png"),
    )

    # GQ vs depth (boxplots) using genotype-level long data
    gq_depth_df = gq_depth_long_table(reader)
    if not gq_depth_df.empty:
        plot_gq_boxplot_vs_depth(
            gq_depth_df.rename(columns={"Sample": "Sample", "Depth": "Depth", "GQ": "GQ"}),
            output_path=str(outdir / "gq_boxplot_vs_depth.png"),
        )

    print(f"Plots written to: {outdir}")


if __name__ == "__main__":
    main()

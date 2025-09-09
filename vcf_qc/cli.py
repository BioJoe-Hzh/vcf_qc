"""Command line interface for vcf_qc.

Current subcommands:
	sample   – compute sample-level metrics and generate plots
	site     – compute site-level metrics and generate plots
	genotype – generate genotype-level plots (GQ, depth, VAF etc.)

Example:
	python -m vcf_qc.cli sample --vcf input.vcf.gz --out outdir --max-site 5000
"""

from __future__ import annotations

import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

from .io import SimpleVCFReader
from .metrics import compute_sample_metrics, gq_depth_long_table
from .metrics.genotype_metrics import genotype_complete_table, filter_heterozygous
from .metrics.site_metrics import compute_site_metrics, filter_site_metrics
from .plot import (
	plot_missing_rate_per_sample,
	plot_depth_per_sample,
	plot_depth_vs_missing_rate,
	plot_het_fraction_per_sample,
	plot_gq_boxplot_vs_depth,
)
from .plot.site_plots import (
	plot_qual_distribution,
	plot_site_mean_depth_distribution,
	plot_qual_vs_mean_depth,
	plot_qd_distribution,
	plot_minor_allele_count_vs_depth,
	plot_minor_allele_freq_vs_depth,
	plot_site_missing_rate_distribution,
	plot_allele_count_pie,
	plot_mac_distribution,
	plot_maf_distribution,
)
from .plot.genotype_plots import (
	plot_gq_distribution_genotype,
	plot_depth_distribution_genotype,
	plot_gq_vs_depth_genotype,
	plot_vaf_vs_depth_het,
	plot_alt_count_vs_depth_het,
)


def _exec_plot(func, kwargs):
	func(**kwargs)


def _run_plot_tasks(tasks, threads: int):
	if threads <= 1:
		for f, kw in tasks:
			f(**kw)
		return
	# macOS / spawn safe
	with ProcessPoolExecutor(max_workers=threads) as ex:
		futs = [ex.submit(_exec_plot, f, kw) for f, kw in tasks]
		for fut in as_completed(futs):
			_ = fut.result()


def cmd_sample(args: argparse.Namespace) -> None:
	outdir = Path(args.out)
	outdir.mkdir(parents=True, exist_ok=True)

	# First pass reader for metrics
	reader_metrics = SimpleVCFReader(args.vcf, max_records=args.max_site)
	metrics_df = compute_sample_metrics(reader_metrics)

	# Prepare plot tasks
	from .plot import plot_het_ratio_vs_depth  # local import
	tasks = [
		(plot_missing_rate_per_sample, {
			'missing_rates': dict(zip(metrics_df.Sample, metrics_df.MissingRate)),
			'output_path': str(outdir / 'sample_missing_rate.png')
		}),
		(plot_depth_per_sample, {
			'depths': dict(zip(metrics_df.Sample, metrics_df.MeanDepth)),
			'output_path': str(outdir / 'sample_mean_depth.png'),
			'cluster_k': args.cluster_k,
			'cluster_results_path': str(outdir / 'depth_clusters.tsv') if args.export_clusters else None,
		}),
		(plot_het_fraction_per_sample, {
			'het_rates': metrics_df[["Sample", "HetRatio"]].rename(columns={"HetRatio": "HetRatio"}),
			'output_path': str(outdir / 'sample_het_ratio.png'),
			'cluster_k': args.cluster_k,
			'cluster_results_path': str(outdir / 'het_clusters.tsv') if args.export_clusters else None,
		}),
		(plot_depth_vs_missing_rate, {
			'data': metrics_df.rename(columns={'MeanDepth': 'Depth'}),
			'output_path': str(outdir / 'depth_vs_missing_rate.png')
		}),
		(plot_het_ratio_vs_depth, {
			'data': metrics_df.rename(columns={'MeanDepth': 'MeanDepth'}),
			'output_path': str(outdir / 'het_ratio_vs_depth.png')
		}),
	]
	_run_plot_tasks(tasks, getattr(args, 'threads', 1))

	# Second pass reader for GQ vs depth long table (needs fresh iterator)
	reader_gq = SimpleVCFReader(args.vcf, max_records=args.max_site)
	compute_gq = not getattr(args, 'no_compute_gq_from_pl', False)
	gq_df = gq_depth_long_table(reader_gq, compute_gq_from_pl=compute_gq)
	if not gq_df.empty:
		plot_gq_boxplot_vs_depth(gq_df, output_path=str(outdir / 'gq_boxplot_vs_depth.png'))
	print(f"Sample-level plots written to {outdir}")


def cmd_site(args: argparse.Namespace) -> None:
	outdir = Path(args.out)
	outdir.mkdir(parents=True, exist_ok=True)

	reader = SimpleVCFReader(args.vcf, max_records=args.max_site)
	site_df = compute_site_metrics(reader, limit=args.max_site)

	# Optional filtering
	site_df_f = filter_site_metrics(
		site_df,
		min_qual=args.min_qual,
		max_missing=args.max_missing,
		min_depth=args.min_depth,
		max_depth=args.max_depth,
		min_mac=args.min_mac,
		min_qd=args.min_qd,
		verbose=True,
	)

	# Save raw & filtered metrics for transparency
	site_df.to_csv(outdir / 'site_metrics_raw.tsv', sep='\t', index=False)
	site_df_f.to_csv(outdir / 'site_metrics_filtered.tsv', sep='\t', index=False)

	if site_df_f.empty:
		print("No site records after filtering.")
		return

	# Determine smart cutoff settings
	use_smart_cutoff = not args.no_smart_cutoff
	no_cutoff_metrics = set(args.no_cutoff_metrics)
	
	# Helper function to determine if smart cutoff should be used for a specific metric
	def should_use_cutoff(metric_name: str) -> bool:
		if args.no_smart_cutoff:
			return False
		return metric_name not in no_cutoff_metrics

	tasks = [
		(plot_qual_distribution, {'qual': site_df_f['QUAL'], 'output_path': str(outdir / 'site_qual_distribution.png'), 'focus_pct': args.focus_qual_qd_pct, 'use_smart_cutoff': should_use_cutoff('qual')}),
		(plot_site_mean_depth_distribution, {'depth': site_df_f['MeanDepth'], 'output_path': str(outdir / 'site_mean_depth_distribution.png'), 'use_smart_cutoff': should_use_cutoff('depth')}),
		(plot_qual_vs_mean_depth, {'data': site_df_f, 'output_path': str(outdir / 'site_qual_vs_mean_depth.png'), 'use_smart_cutoff': should_use_cutoff('qual')}),
		(plot_minor_allele_count_vs_depth, {'data': site_df_f, 'output_path': str(outdir / 'site_mac_vs_mean_depth.png'), 'use_smart_cutoff': should_use_cutoff('mac')}),
		(plot_minor_allele_freq_vs_depth, {'data': site_df_f, 'output_path': str(outdir / 'site_maf_vs_mean_depth.png'), 'use_smart_cutoff': should_use_cutoff('maf')}),
		(plot_site_missing_rate_distribution, {'missing_rates': site_df_f['MissingRate'], 'output_path': str(outdir / 'site_missing_rate_distribution.png'), 'use_smart_cutoff': should_use_cutoff('missing')}),
	]
	if 'QD' in site_df_f.columns:
		tasks.append((plot_qd_distribution, {'qd': site_df_f['QD'].dropna(), 'output_path': str(outdir / 'site_qd_distribution.png'), 'focus_pct': args.focus_qual_qd_pct, 'use_smart_cutoff': should_use_cutoff('qd')}))
	if 'MAC' in site_df_f.columns:
		tasks.append((plot_mac_distribution, {'mac': site_df_f['MAC'], 'output_path': str(outdir / 'site_mac_distribution.png'), 'use_smart_cutoff': should_use_cutoff('mac'), 'min_value': args.mac_plot_min}))
	if 'MAF' in site_df_f.columns:
		tasks.append((plot_maf_distribution, {'maf': site_df_f['MAF'], 'output_path': str(outdir / 'site_maf_distribution.png'), 'use_smart_cutoff': should_use_cutoff('maf'), 'min_value': args.maf_plot_min}))
	if 'AlleleCount' in site_df_f.columns:
		tasks.append((plot_allele_count_pie, {'site_df': site_df_f, 'output_path': str(outdir / 'site_allele_count_pie.png')}))
	_run_plot_tasks(tasks, getattr(args, 'threads', 1))
	print(f"Site-level plots written to {outdir}")


def cmd_genotype(args: argparse.Namespace) -> None:
	outdir = Path(args.out)
	outdir.mkdir(parents=True, exist_ok=True)

	# Use proper genotype metrics module for complete genotype data
	reader = SimpleVCFReader(args.vcf, max_records=args.max_site)
	compute_gq = not getattr(args, 'no_compute_gq_from_pl', False)
	
	# Get complete genotype data using the new genotype_metrics module
	genotype_df = genotype_complete_table(reader, compute_gq_from_pl=compute_gq)
	
	if genotype_df.empty:
		print("No genotype records found.")
		return
	
	# Save raw genotype data to TSV
	genotype_df.to_csv(outdir / 'genotype_records_raw.tsv', sep='\t', index=False)
	print(f"Raw genotype data saved to: {outdir / 'genotype_records_raw.tsv'}")

	# Basic plotting tasks for all genotypes
	tasks = [
		(plot_gq_distribution_genotype, {
			'gq': genotype_df['GQ'], 
			'output_path': str(outdir / 'genotype_gq_distribution.png')
		}),
		(plot_depth_distribution_genotype, {
			'depth': genotype_df['Depth'], 
			'output_path': str(outdir / 'genotype_depth_distribution.png')
		}),
		(plot_gq_vs_depth_genotype, {
			'data': genotype_df[['Depth', 'GQ', 'GT']].dropna(subset=['Depth','GQ']), 
			'output_path': str(outdir / 'genotype_gq_vs_depth.png')
		}),
	]

	# Heterozygous-specific analysis for VAF and ALT count
	# Filter to only heterozygous genotypes for VAF analysis
	het_df = filter_heterozygous(genotype_df)
	
	if not het_df.empty:
		print(f"Found {len(het_df):,} heterozygous genotypes for VAF analysis")
		
		# VAF vs depth plot (heterozygous only)
		vaf_df = het_df.dropna(subset=['VAF', 'Depth'])
		if not vaf_df.empty:
			print(f"  - VAF analysis: {len(vaf_df):,} genotypes with valid VAF data")
			tasks.append((plot_vaf_vs_depth_het, {
				'data': vaf_df[['Depth', 'VAF']], 
				'output_path': str(outdir / 'genotype_vaf_vs_depth_het.png')
			}))
		
		# ALT count vs depth plot (heterozygous only)
		alt_df = het_df.dropna(subset=['ALT_Count', 'Depth'])
		if not alt_df.empty:
			print(f"  - ALT count analysis: {len(alt_df):,} genotypes with valid ALT count data")
			tasks.append((plot_alt_count_vs_depth_het, {
				'data': alt_df[['Depth', 'ALT_Count']], 
				'output_path': str(outdir / 'genotype_alt_count_vs_depth_het.png')
			}))
	else:
		print("No heterozygous genotypes found for VAF/ALT count analysis")

	_run_plot_tasks(tasks, getattr(args, 'threads', 1))
	print(f"Genotype-level plots written to {outdir}")


def build_parser() -> argparse.ArgumentParser:
	p = argparse.ArgumentParser(prog="vcf_qc", description="VCF QC toolkit")
	sub = p.add_subparsers(dest="command")
	sp = sub.add_parser("sample", help="Sample-level metrics and plots")
	sp.add_argument("--vcf", required=True, help="Input VCF or VCF.GZ file")
	sp.add_argument("--out", required=True, help="Output directory for plots")
	sp.add_argument("--max-site", type=int, default=100000, help="Limit number of variant sites parsed (debug)")
	sp.add_argument("--threads", type=int, default=1, help="Parallel plot generation processes")
	sp.add_argument('--cluster-k', type=int, default=None, help='Force specific k for depth/het clustering')
	sp.add_argument('--export-clusters', action='store_true', help='Export cluster assignment TSVs')
	sp.add_argument("--no-compute-gq-from-pl", action="store_true", help="Disable automatic GQ computation from PL when GQ is missing")
	sp.set_defaults(func=cmd_sample)

	sp2 = sub.add_parser("site", help="Site-level metrics and plots")
	sp2.add_argument("--vcf", required=True, help="Input VCF or VCF.GZ file")
	sp2.add_argument("--out", required=True, help="Output directory for plots")
	sp2.add_argument("--max-site", type=int, default=100000, help="Limit number of variant sites parsed (debug)")
	sp2.add_argument("--threads", type=int, default=1, help="Parallel plot generation processes")
	sp2.add_argument("--min-qual", type=float, default=None, help="Minimum QUAL threshold")
	sp2.add_argument("--max-missing", type=float, default=None, help="Maximum missing rate threshold")
	sp2.add_argument("--min-depth", type=int, default=None, help="Minimum site depth threshold")
	sp2.add_argument("--max-depth", type=int, default=None, help="Maximum site depth threshold")
	sp2.add_argument("--min-mac", type=int, default=None, help="Minimum minor allele count threshold")
	sp2.add_argument("--min-qd", type=float, default=None, help="Minimum QD (Quality by Depth) threshold")
	sp2.add_argument("--focus-qual-qd-pct", type=float, default=100.0, dest='focus_qual_qd_pct', help="Optional: restrict QUAL & QD to <= this percentile before smart cutoff (e.g. 95 or 0.95). Default 100 = disabled.")
	sp2.add_argument("--no-smart-cutoff", action="store_true", help="Disable smart cutoff (99.5th percentile) for all metrics")
	sp2.add_argument("--no-cutoff-metrics", type=str, nargs='*', default=[], help="Disable smart cutoff for specific metrics only. Options: qual, depth, qd, mac, maf, missing")
	sp2.add_argument("--mac-plot-min", type=int, default=0, help="Minimum MAC value to include in MAC plots (default: 0)")
	sp2.add_argument("--maf-plot-min", type=float, default=0.0, help="Minimum MAF value to include in MAF plots (default: 0.0)")
	sp2.set_defaults(func=cmd_site)

	sp3 = sub.add_parser("genotype", help="Genotype-level plots (GQ, depth, VAF)")
	sp3.add_argument("--vcf", required=True, help="Input VCF or VCF.GZ file")
	sp3.add_argument("--out", required=True, help="Output directory for plots")
	sp3.add_argument("--max-site", type=int, default=100000, help="Limit number of variant sites parsed (debug)")
	sp3.add_argument("--threads", type=int, default=1, help="Parallel plot generation processes")
	sp3.add_argument("--no-compute-gq-from-pl", action="store_true", help="Disable automatic GQ computation from PL when GQ is missing")
	sp3.set_defaults(func=cmd_genotype)
	return p


def main(argv=None):
	parser = build_parser()
	args = parser.parse_args(argv)
	if not hasattr(args, 'func'):
		parser.print_help()
		return 1
	return args.func(args)


if __name__ == "__main__":  # pragma: no cover
	main()


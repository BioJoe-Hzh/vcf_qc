## vcf_qc – VCF Quality Control & Visualization Toolkit

Light‑weight QC metrics and plots for small to medium VCF datasets. Provides sample / site / genotype level metrics, histograms, joint scatter plots, allele composition, smart tail trimming, and multi‑process rendering.

### Key Features
* Streaming line‑by‑line VCF reader.
* Sample level: missing rate, mean depth, heterozygosity ratio (Het/HomAlt), optional automatic or forced k clustering (scikit‑learn optional).
* Site level: QUAL, MeanDepth, QD, MAC, MAF (second‑most allele frequency), MissingRate, reconstructed multi‑allelic AlleleCount pie chart.
* Genotype level: GQ, per‑genotype Depth, VAF (multi‑allelic aware), ALT read count, heterozygous‑only scatter panels.
* Smart cutoff (iterative tail trimming) to emphasize the main distribution.
* Parallel plot generation via `--threads`.
* Optional install extra: `vcf_qc[cluster]` to enable clustering.

### Simplified Layout
```
vcf_qc/
    cli.py                 # CLI entry (sample/site/genotype)
    io/                    # SimpleVCFReader
    metrics/               # sample / site / genotype metrics
    plot/                  # plotting primitives + level specific functions
pyproject.toml           # packaging (PEP 621)
README.md
LICENSE
```

### Installation
Development / editable (recommended inside a virtual environment):
```bash
pip install -e .
```
With clustering extra:
```bash
pip install -e .[cluster]
```

### Command Line Usage
```
vcf_qc sample   --vcf input.vcf.gz --out sample_qc   --threads 4
vcf_qc site     --vcf input.vcf.gz --out site_qc     --threads 4 --min-qual 30 --min-mac 2
vcf_qc genotype --vcf input.vcf.gz --out genotype_qc --threads 4
```
Core common args:
* `--vcf`: input (bgzip/gzip/plain) VCF.
* `--out`: output directory (auto created).
* `--max-site`: limit parsed variant count (sampling / debug).
* `--threads`: number of parallel plotting processes.

Site subcommand filters: `--min-qual --max-missing --min-depth --min-mac --focus-qual-qd-pct`.

### Output Examples
Sample level:
* sample_missing_rate.png / sample_mean_depth.png / sample_het_ratio.png
* depth_vs_missing_rate.png / het_ratio_vs_depth.png
* gq_boxplot_vs_depth.png

Site level:
* site_qual_distribution.png / site_mean_depth_distribution.png
* site_qual_vs_mean_depth.png / site_qd_distribution.png
* site_mac_vs_mean_depth.png / site_maf_vs_mean_depth.png
* site_mac_distribution.png / site_maf_distribution.png
* site_missing_rate_distribution.png
* site_allele_count_pie.png


Genotype level:
* genotype_gq_distribution.png / genotype_depth_distribution.png
* genotype_gq_vs_depth.png
* genotype_vaf_vs_depth_het.png / genotype_alt_count_vs_depth_het.png


### License
MIT License (see LICENSE).


### Author / Contact
Zihao Huang (idea and code structure) and Copilot (coding)	
Department of Genetics, University of Cambridge	
Email: zh384@cam.ac.uk  

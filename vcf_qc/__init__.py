"""vcf_qc – A modular toolkit for VCF quality control metrics & plots.

Subpackages:
	io        – reading / lightweight parsing helpers
	metrics   – sample / site / genotype level metric computations
	plot      – visualisations organised by data granularity

The public plotting API re-exports the most common plot_* functions so
users can simply::

	from vcf_qc.plot import plot_missing_rate_per_sample

Add new functionality by extending ``metrics`` or ``plot`` with clear
single-responsibility modules.
"""

from importlib import import_module as _imp

# Re-export selected namespaces for convenience.
plot = _imp("vcf_qc.plot")  # noqa: E305
metrics = _imp("vcf_qc.metrics") if False else None  # placeholder for future

__version__ = "0.1.0"
__author__ = "Zihao Huang"
__email__ = "zh384@cam.ac.uk"
__affiliation__ = "Department of Genetics, University of Cambridge"
__all__ = ["plot", "metrics", "__version__", "__author__", "__email__", "__affiliation__"]

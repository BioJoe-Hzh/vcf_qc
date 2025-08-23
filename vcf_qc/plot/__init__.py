"""High-level plotting API for the vcf_qc package.

The submodules are separated by data granularity:
	sample_plots   – sample-level QC visualisations
	site_plots     – site (variant) level QC visualisations
	genotype_plots – genotype (sample × site) level QC visualisations

Import convenience: ``from vcf_qc.plot import plot_missing_rate_per_sample``.
"""

from .sample_plots import *  # noqa: F401,F403
from .site_plots import *  # noqa: F401,F403
from .genotype_plots import *  # noqa: F401,F403
from . import base as _base  # keep base accessible if needed

__all__ = []  # populated by star imports above

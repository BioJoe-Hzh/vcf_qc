"""I/O subpackage.

Currently exposes a lightweight streaming VCF reader. Swap / extend
with high-performance parsers as required.
"""

from .vcf_reader import SimpleVCFReader, GenotypeRecord  # noqa: F401

__all__ = ["SimpleVCFReader", "GenotypeRecord"]

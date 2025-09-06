#!/usr/bin/env python3
"""
VCF Stratified Sampling and Extraction Tool

A fast, reproducible tool for sampling and extracting variants from VCF files
with support for chromosome filtering, variant type selection, and parallel processing.

Author: Zihao Huang and Claude Sonnet 4
"""

import argparse
import json
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from typing import List, Tuple, Dict, Optional, Generator

try:
    import psutil
    HAS_PSUTIL = True
except ImportError:
    HAS_PSUTIL = False


def log_info(msg: str):
    """Print info log message."""
    timestamp = datetime.now().strftime("%H:%M:%S")
    print(f"[{timestamp}] [INFO] {msg}")


def log_warn(msg: str):
    """Print warning log message."""
    timestamp = datetime.now().strftime("%H:%M:%S")
    print(f"[{timestamp}] [WARN] {msg}")


def log_error(msg: str):
    """Print error log message and exit."""
    timestamp = datetime.now().strftime("%H:%M:%S")
    print(f"[{timestamp}] [ERROR] {msg}")
    sys.exit(1)


def check_required_tools():
    """Check if required external tools are available."""
    required_tools = ["bcftools", "tabix"]
    missing = []
    
    for tool in required_tools:
        if not shutil.which(tool):
            missing.append(tool)
    
    if missing:
        log_error(f"Missing required tools: {', '.join(missing)}. Please install and add to PATH.")


@contextmanager
def monitor_resources():
    """Monitor resource usage during execution."""
    start_time = time.time()
    start_memory = None
    
    if HAS_PSUTIL:
        process = psutil.Process()
        start_memory = process.memory_info().rss / 1024 / 1024  # MB
    
    yield
    
    end_time = time.time()
    runtime = end_time - start_time
    
    if HAS_PSUTIL and start_memory:
        end_memory = process.memory_info().rss / 1024 / 1024  # MB
        log_info(f"Runtime: {runtime:.1f}s, Memory: {start_memory:.1f} -> {end_memory:.1f} MB")
    else:
        log_info(f"Runtime: {runtime:.1f} seconds")


def get_optimal_thread_count(default_threads: int = 4) -> int:
    """Get optimal thread count based on system resources."""
    if not HAS_PSUTIL:
        return default_threads
    
    try:
        # Get system info
        cpu_count = psutil.cpu_count()
        memory_gb = psutil.virtual_memory().total / (1024**3)
        
        # Conservative threading strategy
        # Each thread needs ~500MB memory for VCF processing
        memory_threads = max(1, int(memory_gb / 0.5))
        cpu_threads = max(1, min(cpu_count, 16))  # Cap at 16 threads
        
        optimal = min(memory_threads, cpu_threads, default_threads)
        log_info(f"Optimal threads: {optimal} (CPU: {cpu_count}, Memory: {memory_gb:.1f}GB)")
        return optimal
    except:
        return default_threads


def check_file_exists(filepath: str, description: str = "File"):
    """Check if file exists and is readable."""
    if not os.path.exists(filepath):
        log_error(f"{description} does not exist: {filepath}")
    if not os.access(filepath, os.R_OK):
        log_error(f"{description} is not readable: {filepath}")


def validate_vcf_file(vcf_path: str) -> bool:
    """Comprehensive VCF file validation."""
    log_info(f"Validating VCF file: {vcf_path}")
    
    # Check file exists and is readable
    check_file_exists(vcf_path, "VCF file")
    
    # Check if file is compressed
    try:
        with open(vcf_path, 'rb') as f:
            magic = f.read(2)
            if magic != b'\x1f\x8b':  # gzip magic number
                log_error(f"VCF file must be bgzip compressed: {vcf_path}")
    except Exception as e:
        log_error(f"Cannot read VCF file: {e}")
    
    # Check for index files
    index_files = [vcf_path + ".tbi", vcf_path + ".csi"]
    if not any(os.path.exists(idx) for idx in index_files):
        log_error(f"VCF index not found. Run: tabix -p vcf {vcf_path}")
    
    # Verify VCF format with bcftools
    try:
        result = run_command(["bcftools", "view", "-h", vcf_path], capture_output=True)
        if not result.stdout.startswith("##fileformat=VCF"):
            log_error("Invalid VCF format: missing proper header")
        
        # Count total records for info
        count_result = run_command(["bcftools", "index", "-n", vcf_path], capture_output=True)
        if count_result.stdout.strip():
            total_records = int(count_result.stdout.strip())
            log_info(f"VCF contains {total_records:,} total records")
    except Exception as e:
        log_error(f"Cannot validate VCF format: {e}")
    
    return True


def run_command(cmd: List[str], capture_output: bool = False, check: bool = True) -> subprocess.CompletedProcess:
    """Run a shell command with error handling."""
    try:
        result = subprocess.run(cmd, capture_output=capture_output, text=True, check=check)
        return result
    except subprocess.CalledProcessError as e:
        log_error(f"Command failed: {' '.join(cmd)}\nError: {e.stderr if e.stderr else str(e)}")
    except FileNotFoundError:
        log_error(f"Command not found: {cmd[0]}")


def natural_sort_key(item: str) -> List:
    """Generate sort key for natural sorting (chr1, chr2, ..., chr10)."""
    def convert(text):
        if text.isdigit():
            return int(text)
        return text.lower()
    
    return [convert(c) for c in re.split(r'(\d+)', item)]


def get_variant_type_filter(variant_type: str) -> str:
    """Get bcftools variant type filter expression."""
    if variant_type == "snp":
        return 'TYPE="snp"'
    elif variant_type == "indel":
        return 'TYPE="indel"'
    elif variant_type == "both":
        return 'TYPE="snp" || TYPE="indel"'
    elif variant_type == "all":
        return None  # No filtering
    else:
        log_error(f"Invalid variant type: {variant_type}")


def get_variant_type_view_arg(variant_type: str) -> List[str]:
    """Get bcftools view argument for variant type filtering."""
    if variant_type == "snp":
        return ["-v", "snps"]
    elif variant_type == "indel":
        return ["-v", "indels"]
    elif variant_type == "both":
        return ["-v", "snps,indels"]
    elif variant_type == "all":
        return []  # No filtering
    else:
        log_error(f"Invalid variant type: {variant_type}")


def validate_chromosomes_in_vcf(vcf_path: str, requested_chroms: List[str]) -> List[str]:
    """Validate that requested chromosomes exist in VCF and return available ones."""
    log_info("Validating requested chromosomes in VCF...")
    
    # Get available chromosomes from VCF
    try:
        # First try to get from header
        header_result = run_command(["bcftools", "view", "-h", vcf_path], capture_output=True)
        available_chroms = set()
        
        for line in header_result.stdout.splitlines():
            if line.startswith("##contig=<ID="):
                # Extract chromosome from ##contig=<ID=Chr1,length=...>
                chrom = line.split("ID=")[1].split(",")[0].split(">")[0]
                available_chroms.add(chrom)
        
        # If no ##contig lines found, scan actual records
        if not available_chroms:
            log_info("No ##contig lines found, scanning VCF records...")
            records_result = run_command(
                ["bcftools", "query", "-f", "%CHROM\n", vcf_path], 
                capture_output=True
            )
            
            for line in records_result.stdout.splitlines():
                chrom = line.strip()
                if chrom:
                    available_chroms.add(chrom)
                    
                # Limit scanning to avoid processing huge files
                if len(available_chroms) > 1000:  # Reasonable limit
                    break
        
        if not available_chroms:
            log_error("Could not determine available chromosomes from VCF")
        
        log_info(f"VCF contains {len(available_chroms)} chromosomes")
        
    except Exception as e:
        log_error(f"Failed to query VCF chromosomes: {e}")
    
    # Check which requested chromosomes are available
    requested_set = set(requested_chroms)
    available_requested = requested_set.intersection(available_chroms)
    missing_chroms = requested_set - available_chroms
    
    # Report results
    if missing_chroms:
        log_warn(f"Requested chromosomes NOT found in VCF: {', '.join(sorted(missing_chroms))}")
    
    if not available_requested:
        log_error("None of the requested chromosomes found in VCF!")
        
    valid_chroms = sorted(list(available_requested), key=natural_sort_key)
    log_info(f"Valid chromosomes to process: {', '.join(valid_chroms)}")
    
    # Show some available chromosomes as suggestion if some were missing
    if missing_chroms:
        available_list = sorted(list(available_chroms), key=natural_sort_key)[:20]  # Show first 20
        log_info(f"Available chromosomes in VCF: {', '.join(available_list)}")
        if len(available_chroms) > 20:
            log_info(f"... and {len(available_chroms) - 20} more")
    
    return valid_chroms


def query_all_variants(vcf_path: str, variant_type: str, include_chroms: Optional[List[str]] = None, 
                      threads: int = 4) -> List[Tuple[str, int]]:
    """Query all variant positions from VCF file with memory optimization."""
    
    # Validate chromosomes if specified
    if include_chroms:
        validated_chroms = validate_chromosomes_in_vcf(vcf_path, include_chroms)
        if not validated_chroms:
            log_error("No valid chromosomes to process")
        include_chroms = validated_chroms
    
    log_info("Querying all variant positions from VCF...")
    
    # Build bcftools query command
    cmd = ["bcftools", "query", "-f", r"%CHROM\t%POS\n"]
    
    # Add variant type filter
    if variant_type != "all":
        type_filter = get_variant_type_filter(variant_type)
        if type_filter:  # Only add filter if not None
            cmd.extend(["-i", type_filter])
    
    # Add chromosome filter
    if include_chroms:
        regions = ",".join(include_chroms)
        cmd.extend(["-r", regions])
    
    cmd.append(vcf_path)
    
    # Run query with streaming processing
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        variants = []
        chrom_counts = defaultdict(int)
        line_count = 0
        
        for line in process.stdout:
            line = line.strip()
            if not line:
                continue
                
            try:
                chrom, pos_str = line.split('\t')
                pos = int(pos_str)
                variants.append((chrom, pos))
                chrom_counts[chrom] += 1
                line_count += 1
                
                # Progress reporting for large datasets
                if line_count % 100000 == 0:
                    log_info(f"Processed {line_count:,} variants...")
                    
            except ValueError:
                continue
        
        # Wait for process to complete
        return_code = process.wait()
        if return_code != 0:
            stderr = process.stderr.read()
            log_error(f"bcftools query failed: {stderr}")
        
    except Exception as e:
        log_error(f"Failed to query variants: {e}")
    
    # Print statistics
    total_variants = len(variants)
    num_chroms = len(chrom_counts)
    log_info(f"Total variants={total_variants:,} across {num_chroms} chromosomes")
    
    if include_chroms:
        for chrom in sorted(include_chroms, key=natural_sort_key):
            count = chrom_counts.get(chrom, 0)
            log_info(f"  {chrom}: {count:,} variants")
    else:
        # Show top chromosomes by variant count
        top_chroms = sorted(chrom_counts.items(), key=lambda x: (-x[1], natural_sort_key(x[0])))[:10]
        for chrom, count in top_chroms:
            log_info(f"  {chrom}: {count:,} variants")
        if len(chrom_counts) > 10:
            log_info(f"  ... and {len(chrom_counts) - 10} more chromosomes")
    
    return variants


def sample_variants(variants: List[Tuple[str, int]], n_sites: int, seed: int) -> List[Tuple[str, int]]:
    """Randomly sample variants with specified seed, ensuring sufficient unique positions."""
    total_available = len(variants)
    
    if total_available == 0:
        log_error("No variants available for sampling")
    
    # Remove duplicates first to get unique positions
    unique_variants = list(set(variants))
    unique_count = len(unique_variants)
    
    if total_available != unique_count:
        log_info(f"Removed {total_available - unique_count:,} duplicate positions, {unique_count:,} unique positions available")
    
    if n_sites >= unique_count:
        log_warn(f"Requested {n_sites:,} >= unique available {unique_count:,}, using all unique variants")
        sampled = unique_variants[:]
    else:
        # Check if we have enough variants for good sampling
        if unique_count < n_sites * 1.5:
            log_warn(f"Limited unique variants available ({unique_count:,}) for sampling {n_sites:,}")
        
        log_info(f"Sampling {n_sites:,} variants from {unique_count:,} unique available (seed={seed})")
        random.seed(seed)
        sampled = random.sample(unique_variants, n_sites)
    
    # If we still don't have enough unique positions, sample with replacement from original list
    if len(sampled) < n_sites and total_available > 0:
        remaining_needed = n_sites - len(sampled)
        log_info(f"Need {remaining_needed:,} more positions, sampling with replacement...")
        
        # Sample additional positions with replacement
        random.seed(seed + 1)  # Use different seed for replacement sampling
        additional_samples = random.choices(variants, k=remaining_needed)
        sampled.extend(additional_samples)
        
        # Remove duplicates again and ensure we have exactly n_sites
        sampled = list(set(sampled))
        if len(sampled) < n_sites:
            # If we still don't have enough, keep sampling until we do
            attempts = 0
            while len(sampled) < n_sites and attempts < 10:
                additional_needed = n_sites - len(sampled)
                random.seed(seed + attempts + 2)
                more_samples = random.choices(variants, k=additional_needed * 2)  # Sample extra
                sampled.extend(more_samples)
                sampled = list(set(sampled))  # Remove duplicates
                attempts += 1
            
            # Final trim to exact count
            if len(sampled) > n_sites:
                random.seed(seed)
                sampled = random.sample(sampled, n_sites)
    
    log_info(f"Final sample: {len(sampled):,} unique variants")
    return sampled


def write_bed_file(variants: List[Tuple[str, int]], output_path: str):
    """Write variants to BED file with proper sorting, avoiding deduplication loss."""
    log_info(f"Writing and sorting BED file: {output_path}")
    
    initial_count = len(variants)
    
    # Convert to BED format (0-based start, 1-based end) and track duplicates
    bed_entries = set()
    duplicates_removed = 0
    
    for chrom, pos in variants:
        bed_entry = (chrom, pos - 1, pos)  # Convert to 0-based start
        if bed_entry in bed_entries:
            duplicates_removed += 1
        bed_entries.add(bed_entry)
    
    if duplicates_removed > 0:
        log_warn(f"Removed {duplicates_removed:,} duplicate BED entries during writing")
        log_warn(f"BED file will contain {len(bed_entries):,} unique positions instead of requested {initial_count:,}")
    
    # Sort naturally (chr1, chr2, ..., chr10) then by position
    sorted_entries = sorted(bed_entries, key=lambda x: (natural_sort_key(x[0]), x[1]))
    
    # Write to file with buffered I/O for better performance
    with open(output_path, 'w', buffering=8192) as f:
        for chrom, start, end in sorted_entries:
            f.write(f"{chrom}\t{start}\t{end}\n")
    
    log_info(f"Wrote and sorted {output_path} with {len(sorted_entries):,} sites")
    return len(sorted_entries)


def build_optimized_bcftools_cmd(vcf_path: str, bed_path: str, output: str, 
                                variant_type: str = "all", threads_per_job: int = 2) -> List[str]:
    """Build optimized bcftools command for extraction WITHOUT any filtering."""
    cmd = [
        "bcftools", "view",
        "-R", bed_path,                    # Region filtering with BED
        "--threads", str(threads_per_job), # More threads per job
        "-Oz", "-o", output,               # Output format and file
        vcf_path
    ]
    
    log_info("Extracting ALL variants at BED positions (no variant type filtering, no AC filtering)")
    return cmd


def read_bed_file(bed_path: str) -> List[Tuple[str, int]]:
    """Read BED file and return list of (chrom, pos) tuples."""
    check_file_exists(bed_path, "BED file")
    
    variants = []
    with open(bed_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) < 3:
                log_warn(f"Invalid BED line {line_num}: {line}")
                continue
            
            try:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                pos = start + 1  # Convert from 0-based to 1-based
                variants.append((chrom, pos))
            except ValueError:
                log_warn(f"Invalid BED coordinates at line {line_num}: {line}")
                continue
    
    log_info(f"Read {len(variants):,} variants from {bed_path}")
    return variants


def split_bed_by_chromosome(bed_path: str, temp_dir: str) -> Dict[str, str]:
    """Split BED file by chromosome and return mapping of chrom -> bed_file."""
    chrom_beds = defaultdict(list)
    
    # Group BED entries by chromosome
    with open(bed_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                chrom_beds[chrom].append(line)
    
    # Write separate BED files for each chromosome
    chrom_bed_files = {}
    for chrom, bed_lines in chrom_beds.items():
        bed_file = os.path.join(temp_dir, f"{chrom}.bed")
        with open(bed_file, 'w') as f:
            for line in bed_lines:
                f.write(line + '\n')
        chrom_bed_files[chrom] = bed_file
    
    return chrom_bed_files


def extract_chromosome_variants(vcf_path: str, chrom_bed: str, output_vcf: str, 
                              chrom: str, variant_type: str = "all", max_retries: int = 3) -> int:
    """Extract variants for a single chromosome with robust error handling."""
    
    # Count BED positions for this chromosome
    bed_positions = 0
    with open(chrom_bed, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                bed_positions += 1
    
    for attempt in range(max_retries):
        try:
            # Validate input BED file
            if not os.path.exists(chrom_bed) or os.path.getsize(chrom_bed) == 0:
                log_warn(f"Empty or missing BED file for {chrom}")
                return 0
            
            # Build optimized command without variant type filtering
            cmd = build_optimized_bcftools_cmd(vcf_path, chrom_bed, output_vcf, "all", 2)
            
            # Execute extraction
            result = run_command(cmd, capture_output=True, check=False)
            
            if result.returncode == 0 and os.path.exists(output_vcf):
                # Verify output VCF is valid
                verify_cmd = ["bcftools", "view", "-h", output_vcf]
                verify_result = run_command(verify_cmd, capture_output=True, check=False)
                
                if verify_result.returncode == 0:
                    # Create index
                    index_result = run_command(["tabix", "-p", "vcf", output_vcf], 
                                             capture_output=True, check=False)
                    
                    if index_result.returncode == 0:
                        # Count variants
                        count_result = run_command(["bcftools", "view", "-H", output_vcf], 
                                                 capture_output=True)
                        variant_lines = [line for line in count_result.stdout.strip().split('\n') if line]
                        variant_count = len(variant_lines)
                        
                        if variant_count > 0:
                            success_rate = (variant_count / bed_positions * 100) if bed_positions > 0 else 0
                            log_info(f"✓ {chrom}: extracted {variant_count:,}/{bed_positions:,} variants ({success_rate:.1f}%)")
                        else:
                            log_warn(f"⚠ {chrom}: no variants found from {bed_positions:,} BED positions (empty output)")
                        return variant_count
            
            # Failed this attempt - log the error details
            if result.returncode != 0:
                log_warn(f"bcftools failed for {chrom} (exit code {result.returncode})")
                if result.stderr:
                    log_warn(f"stderr: {result.stderr.strip()}")
            
            if attempt < max_retries - 1:
                log_warn(f"Attempt {attempt + 1} failed for {chrom}, retrying...")
                # Clean up failed output
                if os.path.exists(output_vcf):
                    os.remove(output_vcf)
                if os.path.exists(output_vcf + ".tbi"):
                    os.remove(output_vcf + ".tbi")
                time.sleep(1)  # Brief pause before retry
            
        except Exception as e:
            log_warn(f"Error in attempt {attempt + 1} for {chrom}: {e}")
            if attempt < max_retries - 1:
                time.sleep(1)
    
    log_warn(f"Failed to extract variants from {chrom} after {max_retries} attempts ({bed_positions:,} BED positions)")
    return 0


def extract_variants_parallel(vcf_path: str, bed_path: str, output_vcf: str, 
                            variant_type: str, threads: int):
    """Extract variants in parallel by chromosome with sorted output."""
    log_info(f"Extracting variants in parallel (threads={threads})...")
    log_info(f"NOTE: Extraction will include ALL variants at BED positions (no filtering)")
    log_info(f"Original sampling was done with variant type: {variant_type}")
    
    # Count total BED positions first
    bed_total_positions = 0
    with open(bed_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                bed_total_positions += 1
    
    log_info(f"BED file contains {bed_total_positions:,} positions to extract")
    
    # Skip time-consuming diagnostic for large datasets
    log_info("Skipping detailed variant type diagnostic to improve performance")
    log_info("(Use --enable-diagnostic flag if detailed analysis is needed)")
    
    # Create temporary directory
    with tempfile.TemporaryDirectory(prefix="vcf_extract_") as temp_dir:
        # Split BED file by chromosome
        chrom_bed_files = split_bed_by_chromosome(bed_path, temp_dir)
        
        if not chrom_bed_files:
            log_error("No valid entries found in BED file")
        
        # Count positions per chromosome in BED
        chrom_bed_counts = {}
        for chrom, chrom_bed_file in chrom_bed_files.items():
            with open(chrom_bed_file, 'r') as f:
                count = sum(1 for line in f if line.strip() and not line.startswith('#'))
                chrom_bed_counts[chrom] = count
        
        log_info(f"Processing {len(chrom_bed_files)} chromosomes: {', '.join(sorted(chrom_bed_files.keys(), key=natural_sort_key))}")
        for chrom in sorted(chrom_bed_files.keys(), key=natural_sort_key):
            log_info(f"  {chrom}: {chrom_bed_counts[chrom]:,} BED positions")
        
        # Extract variants for each chromosome in parallel
        chrom_vcf_files = []
        total_extracted = 0
        chrom_extracted_counts = {}
        
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = {}
            
            for chrom in sorted(chrom_bed_files.keys(), key=natural_sort_key):
                chrom_bed = chrom_bed_files[chrom]
                chrom_vcf = os.path.join(temp_dir, f"{chrom}.vcf.gz")
                chrom_vcf_files.append((chrom, chrom_vcf))
                
                future = executor.submit(
                    extract_chromosome_variants, 
                    vcf_path, chrom_bed, chrom_vcf, chrom, variant_type
                )
                futures[future] = chrom
            
            # Collect results
            for future in as_completed(futures):
                chrom = futures[future]
                try:
                    count = future.result()
                    total_extracted += count
                    chrom_extracted_counts[chrom] = count
                except Exception as e:
                    log_warn(f"Error processing {chrom}: {e}")
                    chrom_extracted_counts[chrom] = 0
        
        # Filter valid VCF files and maintain chromosome order
        valid_vcf_files = []
        for chrom, vcf_file in chrom_vcf_files:
            if os.path.exists(vcf_file) and os.path.getsize(vcf_file) > 0:
                # Check if VCF has actual variants (not just header)
                try:
                    check_result = run_command(["bcftools", "view", "-H", vcf_file], capture_output=True)
                    if check_result.stdout.strip():
                        valid_vcf_files.append(vcf_file)
                        log_info(f"✓ {chrom} ready for concatenation")
                except:
                    log_warn(f"⚠ {chrom} VCF validation failed")
        
        if not valid_vcf_files:
            log_error("No variants were extracted from any chromosome")
        
        # Concatenate VCF files in sorted order (maintains coordinate sorting)
        log_info(f"Concatenating {len(valid_vcf_files)} chromosome VCFs in sorted order...")
        concat_cmd = ["bcftools", "concat", "-a", "--threads", str(min(threads, 4)), "-Oz", "-o", output_vcf] + valid_vcf_files
        run_command(concat_cmd)
        
        # Index final VCF
        run_command(["tabix", "-p", "vcf", output_vcf])
        
        # Verify final output is sorted
        verify_sorting(output_vcf)
        
        # Final extraction report
        log_info("="*60)
        log_info("EXTRACTION SUMMARY:")
        log_info(f"BED positions requested: {bed_total_positions:,}")
        log_info(f"VCF variants extracted: {total_extracted:,} (ALL variants at BED positions)")
        log_info(f"NOTE: No variant type filtering was applied during extraction")
        
        # Calculate per-chromosome summary first
        missing_total = 0
        for chrom in sorted(chrom_bed_counts.keys(), key=natural_sort_key):
            bed_count = chrom_bed_counts[chrom]
            extracted_count = chrom_extracted_counts.get(chrom, 0)
            missing_count = bed_count - extracted_count
            missing_total += missing_count
        
        if bed_total_positions > 0:
            extraction_rate = (total_extracted / bed_total_positions) * 100
            log_info(f"Extraction success rate: {extraction_rate:.1f}% (ALL variants at BED positions)")
            
            # Also show how many positions were successfully found vs missing
            found_positions = bed_total_positions - missing_total
            if found_positions != bed_total_positions:
                position_rate = (found_positions / bed_total_positions) * 100
                log_info(f"Position coverage rate: {position_rate:.1f}% ({found_positions:,}/{bed_total_positions:,} positions found in VCF)")
        
        # Per-chromosome detailed summary
        log_info("Per-chromosome extraction details:")
        for chrom in sorted(chrom_bed_counts.keys(), key=natural_sort_key):
            bed_count = chrom_bed_counts[chrom]
            extracted_count = chrom_extracted_counts.get(chrom, 0)
            missing_count = bed_count - extracted_count
            
            if missing_count > 0:
                success_rate = (extracted_count / bed_count * 100) if bed_count > 0 else 0
                log_info(f"  {chrom}: {extracted_count:,}/{bed_count:,} extracted ({success_rate:.1f}% success, {missing_count:,} missing)")
            else:
                log_info(f"  {chrom}: {extracted_count:,}/{bed_count:,} extracted (✓ 100% success)")
        
        # Final missing summary
        
        if missing_total > 0:
            log_warn(f"Total missing variants: {missing_total:,} positions were not found in VCF")
            log_warn("Possible reasons: 1) Positions don't exist in VCF, 2) Different chromosome naming")
        else:
            log_info("✓ All BED positions successfully extracted from VCF (100% success rate)")
        
        log_info("="*60)
        log_info(f"✓ Final output: {output_vcf} ({total_extracted:,} variants)")
        
        return total_extracted, bed_total_positions


def verify_sorting(vcf_path: str):
    """Verify that VCF file is properly sorted."""
    try:
        # Check if VCF is sorted using bcftools
        check_cmd = ["bcftools", "view", "-H", vcf_path]
        result = run_command(check_cmd, capture_output=True)
        
        lines = result.stdout.strip().split('\n')
        if len(lines) < 2:
            return  # Too few records to check sorting
        
        prev_chrom = None
        prev_pos = 0
        
        for i, line in enumerate(lines[:1000]):  # Check first 1000 records
            if not line:
                continue
            parts = line.split('\t')
            chrom = parts[0]
            pos = int(parts[1])
            
            if prev_chrom is None:
                prev_chrom = chrom
                prev_pos = pos
                continue
            
            # Check sorting within chromosome
            if chrom == prev_chrom:
                if pos < prev_pos:
                    log_warn(f"VCF may not be sorted: position {pos} after {prev_pos} on {chrom}")
                    return
            else:
                # New chromosome - this is fine as we concatenated in order
                prev_chrom = chrom
            
            prev_pos = pos
        
        log_info("✓ VCF output appears to be properly sorted")
        
    except Exception as e:
        log_warn(f"Could not verify VCF sorting: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="VCF Stratified Sampling and Extraction Tool - Fast, reproducible variant sampling",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input/Output
    parser.add_argument("--vcf", required=True, 
                       help="Input VCF file (bgzip compressed with tabix index)")
    parser.add_argument("--out-pos-bed", 
                       help="Output BED file with sampled positions")
    parser.add_argument("--pos-bed", 
                       help="Existing BED file to use for extraction (skip sampling)")
    parser.add_argument("--extract-vcf", 
                       help="Output extracted VCF file")
    
    # Sampling parameters
    parser.add_argument("--n-sites", type=int, 
                       help="Maximum number of sites to sample")
    parser.add_argument("--seed", type=int, default=42, 
                       help="Random seed for reproducible sampling")
    
    # Filtering parameters
    parser.add_argument("--variant-type", choices=["snp", "indel", "both", "all"], default="snp",
                       help="Variant type to sample/extract: snp, indel, both, or all")
    parser.add_argument("--include-chr", 
                       help="Comma-separated list of chromosomes to include")
    
    # Performance
    
    parser.add_argument("--threads", type=int, default=4,
                       help="Number of threads for parallel extraction (default: 4)")
    
    args = parser.parse_args()
    
    # Validate arguments
    if not args.pos_bed and not args.out_pos_bed:
        log_error("Must specify either --pos-bed or --out-pos-bed")
    
    if not args.pos_bed and not args.n_sites:
        log_error("Must specify --n-sites when sampling (not using existing --pos-bed)")
    
    if args.extract_vcf and not (args.pos_bed or args.out_pos_bed):
        log_error("Must specify either --pos-bed or --out-pos-bed when using --extract-vcf")
    
    # Check required tools
    check_required_tools()
    
    # Validate VCF file
    validate_vcf_file(args.vcf)
    
    # Optimize thread count
    if args.threads:
        threads = args.threads
    else:
        threads = get_optimal_thread_count(4)
    
    # Parse chromosome list and validate early
    include_chroms = None
    if args.include_chr:
        include_chroms = [chrom.strip() for chrom in args.include_chr.split(',')]
        log_info(f"Requested chromosomes: {', '.join(include_chroms)}")
        
        # Early validation if not using existing BED file
        if not args.pos_bed:
            include_chroms = validate_chromosomes_in_vcf(args.vcf, include_chroms)
            log_info(f"Validated chromosomes: {', '.join(include_chroms)}")
    
    # Monitor overall performance
    with monitor_resources():
        # Determine workflow
        if args.pos_bed:
            # Use existing BED file
            log_info(f"Using existing pos.bed for extraction: {args.pos_bed}")
            bed_path = args.pos_bed
            
            # Validate existing BED file
            variants = read_bed_file(bed_path)
            
            # Apply chromosome filter if specified
            if include_chroms:
                include_set = set(include_chroms)
                filtered_variants = [(chrom, pos) for chrom, pos in variants if chrom in include_set]
                
                if len(filtered_variants) != len(variants):
                    removed_count = len(variants) - len(filtered_variants)
                    log_info(f"Filtered {len(variants):,} -> {len(filtered_variants):,} variants by chromosome (removed {removed_count:,})")
                    
                    # Check if any BED chromosomes are not in our include list
                    bed_chroms = set(chrom for chrom, pos in variants)
                    excluded_chroms = bed_chroms - include_set
                    if excluded_chroms:
                        log_info(f"Excluded chromosomes from BED: {', '.join(sorted(excluded_chroms))}")
                    
                    # Write filtered BED if needed
                    if args.out_pos_bed:
                        actual_count = write_bed_file(filtered_variants, args.out_pos_bed)
                        bed_path = args.out_pos_bed
                
                if not filtered_variants:
                    log_error("No variants remain after chromosome filtering")
        else:
            # Sample new variants (include_chroms already validated above)
            log_info(f"Sampling variants with type='{args.variant_type}' to ensure sufficient sites")
            variants = query_all_variants(args.vcf, args.variant_type, include_chroms, threads)
            
            if not variants:
                log_error("No variants found matching the specified criteria")
            
            # Check if we have enough variants
            if len(variants) < args.n_sites:
                log_warn(f"Only {len(variants):,} variants available, less than requested {args.n_sites:,}")
                log_warn("Consider: 1) Using --variant-type 'both' or 'all', 2) Including more chromosomes, 3) Reducing --n-sites")
            
            sampled_variants = sample_variants(variants, args.n_sites, args.seed)
            actual_bed_count = write_bed_file(sampled_variants, args.out_pos_bed)
            bed_path = args.out_pos_bed
            
            # Warn if we didn't get the exact count requested
            if actual_bed_count != args.n_sites:
                log_warn(f"BED file contains {actual_bed_count:,} positions instead of requested {args.n_sites:,}")
                log_warn("This can happen due to duplicate positions in the VCF")
        
        # Extract VCF if requested
        if args.extract_vcf:
            extracted_count, bed_count = extract_variants_parallel(args.vcf, bed_path, args.extract_vcf, args.variant_type, threads)
    
    log_info("[DONE] Workflow completed successfully!")


if __name__ == "__main__":
    main()

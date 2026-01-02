"""
Common functions and utilities for mgraft_v2 Snakemake workflows.

This module provides:
- Samplesheet loading and validation
- Input detection (proteins, MGEs available?)
- Helper functions for dynamic inputs
"""
import csv
from pathlib import Path
from datetime import datetime


# ======================== Column Name Alternatives ========================
PROTEIN_COLUMNS = ("protein_fasta", "proteins", "protein", "path_to_protein_FASTA", "faa")
GENOME_COLUMNS = ("genome", "genome_fasta", "path_to_genomic_FASTA", "fna")
MGE_COLUMNS = ("mges", "mge_gff", "path_to_mge_gff3", "gff")


# ======================== Output Path Helpers ========================
# Note: OUTPUT_DIR and RETAIN_INTERMEDIATES are set in the main Snakefile
# These functions are called after those variables are defined

def out_path(sample, subdir, filename):
    """Generate output path respecting RETAIN_INTERMEDIATES.
    
    Args:
        sample: Sample ID
        subdir: Subdirectory within sample output (e.g., 'parsed', 'crispr-filtering')
        filename: Output filename
        
    Returns:
        Path string. When RETAIN_INTERMEDIATES=True, adds 'intermediates/' prefix.
    """
    if RETAIN_INTERMEDIATES:
        return f"{OUTPUT_DIR}/{sample}/intermediates/{subdir}/{filename}"
    return f"{OUTPUT_DIR}/{sample}/{subdir}/{filename}"


def global_out_path(subdir, filename):
    """Generate global (non-sample-specific) output path.
    
    Used for cross-sample outputs like clustering results.
    """
    return f"{OUTPUT_DIR}/{subdir}/{filename}"


def load_samplesheet(path):
    """
    Load samplesheet and normalize column names.
    
    Supports multiple input formats with flexible column naming.
    
    Args:
        path: Path to tab-delimited samplesheet
        
    Returns:
        dict: sample_id -> normalized row dict
    """
    with open(path, newline="") as handle:
        first_line = handle.readline().strip()
        handle.seek(0)
        
        # Detect header
        first_field = first_line.split('\t')[0] if '\t' in first_line else first_line.split()[0]
        has_header = first_field.lower() in (
            'sample', 'file_name', 'genome_id', 'id', 'name'
        ) or first_field.startswith('path') or first_field.endswith('_fasta')
        
        # Check if first field looks like a path (headerless)
        if first_field.startswith('/') or first_field.startswith('.'):
            has_header = False
            
        if has_header:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = list(reader)
        else:
            rows = []
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                fields = line.split('\t') if '\t' in line else line.split()
                if len(fields) >= 2:
                    row = {
                        'sample': fields[0],
                        'genome': fields[1],
                        'protein_fasta': fields[2] if len(fields) > 2 else None,
                        'mges': fields[3] if len(fields) > 3 else None,
                    }
                    rows.append(row)
    
    # Normalize column names
    column_mapping = {
        'file_name': 'sample',
        'path_to_genomic_FASTA': 'genome',
        'path_to_protein_FASTA': 'protein_fasta',
        'path_to_mge_gff3': 'mges',
        'fna': 'genome',
        'faa': 'protein_fasta',
    }
    
    normalized = {}
    for row in rows:
        norm_row = {}
        for old_col, new_col in column_mapping.items():
            if old_col in row and row[old_col]:
                norm_row[new_col] = row[old_col]
        
        for col, val in row.items():
            if col not in column_mapping and col not in norm_row:
                norm_row[col] = val
        
        sample_id = norm_row.get('sample') or row.get('file_name') or row.get('sample')
        if not sample_id:
            raise ValueError(f"Row missing sample identifier: {row}")
        
        normalized[sample_id] = norm_row
    
    return normalized


def has_proteins(samplesheet):
    """Check if all samples have protein FASTAs."""
    for sample, row in samplesheet.items():
        has_protein = any(row.get(col) for col in PROTEIN_COLUMNS)
        if not has_protein:
            return False
    return True


def has_mges(samplesheet):
    """Check if all samples have MGE GFF files."""
    for sample, row in samplesheet.items():
        has_mge = any(row.get(col) for col in MGE_COLUMNS)
        if not has_mge:
            return False
    return True


def genome_for(sample, samplesheet):
    """Get genome FASTA path for a sample."""
    row = samplesheet.get(sample)
    if row is None:
        raise ValueError(f"Sample {sample} missing from samplesheet")
    for column in GENOME_COLUMNS:
        path = row.get(column)
        if path:
            return path
    raise ValueError(f"No genome FASTA for {sample}")


def proteins_for(sample, samplesheet, output_dir=None):
    """
    Get protein FASTA path for a sample.
    
    If sample has protein FASTA in samplesheet, return that.
    Otherwise, return path to Pyrodigal output.
    """
    row = samplesheet.get(sample)
    if row is None:
        raise ValueError(f"Sample {sample} missing from samplesheet")
    
    for column in PROTEIN_COLUMNS:
        path = row.get(column)
        if path:
            return path
    
    # Return Pyrodigal output path
    if output_dir:
        return f"{output_dir}/{sample}/tools-raw-out/pyrodigal/{sample}.faa"
    raise ValueError(f"No protein FASTA for {sample} and output_dir not specified")


def mges_for(sample, samplesheet):
    """Get MGE GFF path for a sample, or None if not available."""
    row = samplesheet.get(sample)
    if row is None:
        raise ValueError(f"Sample {sample} missing from samplesheet")
    for column in MGE_COLUMNS:
        path = row.get(column)
        if path:
            return path
    return None


def ffn_for(sample, samplesheet, output_dir=None):
    """Get nucleotide gene FASTA (FFN) path for a sample.
    
    First checks samplesheet for user-provided FFN.
    Otherwise, return path to Pyrodigal output.
    """
    row = samplesheet.get(sample)
    if row is None:
        raise ValueError(f"Sample {sample} missing from samplesheet")
    
    # Check for user-provided FFN column
    for column in ['ffn', 'genes', 'nucleotide_genes']:
        path = row.get(column)
        if path:
            return path
    
    # Return Pyrodigal output path
    if output_dir:
        return f"{output_dir}/{sample}/tools-raw-out/pyrodigal/{sample}.ffn"
    raise ValueError(f"No FFN for {sample} and output_dir not specified")


def get_run_id():
    """Generate run ID based on date."""
    return f"v2_{datetime.now().strftime('%Y%m%d')}"


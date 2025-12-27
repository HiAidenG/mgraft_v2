"""
Common helper functions shared across mgraft_v2 workflows.
"""
import csv
from datetime import datetime

# Column name alternatives for common file types
PROTEIN_COLUMNS = ("protein_fasta", "proteins", "protein", "path_to_protein_FASTA", "faa")
GENOME_COLUMNS = ("genome", "genome_fasta", "path_to_genomic_FASTA", "fna")
GENE_COLUMNS = ("gene_fasta", "genes", "ffn", "path_to_gene_FASTA")
MGE_COLUMNS = ("mges", "mge_gff", "path_to_mge_gff3", "gff")


def load_samplesheet(path):
    """Load samplesheet and normalize column names to internal format.
    
    Supports multiple formats:
    - Header format: sample, protein_fasta, genome, mges, etc.
    - Header format: file_name, path_to_protein_FASTA, path_to_genomic_FASTA, path_to_mge_gff3
    - Headerless format (6 columns): unknown, sample, fna, faa, ffn, gff
    
    Args:
        path: Path to tab-delimited samplesheet file
        
    Returns:
        dict: Mapping of sample_id (genome_id) -> normalized row dict
        
    Raises:
        ValueError: If row missing sample identifier
    """
    with open(path, newline="") as handle:
        first_line = handle.readline().strip()
        handle.seek(0)
        
        # Detect if file has headers by checking first field
        # If first field looks like a path or known non-header value, it's headerless
        first_field = first_line.split('\t')[0] if '\t' in first_line else first_line.split()[0]
        has_header = first_field.lower() in (
            'sample', 'file_name', 'genome_id', 'id', 'name',
            'unknown', 'taxonomy', 'species'
        ) or first_field.startswith('path') or first_field.endswith('_fasta')
        
        # Check if it looks like a path (headerless data row)
        if first_field.startswith('/') or first_field.startswith('.'):
            has_header = False
        # Check if second column looks like a genome ID (headerless with unknown first col)
        fields = first_line.split('\t') if '\t' in first_line else first_line.split()
        if len(fields) >= 2 and fields[1].startswith('GUT_GENOME'):
            has_header = False
        
        if has_header:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = list(reader)
        else:
            # Headerless format: assign column names positionally
            # Format: unknown, sample, fna, faa, ffn, gff
            rows = []
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                # Split by tab, fallback to whitespace
                fields = line.split('\t') if '\t' in line else line.split()
                if len(fields) >= 4:
                    row = {
                        'taxonomy': fields[0] if len(fields) > 4 else None,
                        'sample': fields[1] if len(fields) > 4 else fields[0],
                        'fna': fields[2] if len(fields) > 4 else fields[1],
                        'faa': fields[3] if len(fields) > 4 else fields[2],
                        'ffn': fields[4] if len(fields) > 5 else (fields[3] if len(fields) > 3 else None),
                        'gff': fields[5] if len(fields) > 5 else (fields[4] if len(fields) > 4 else None),
                    }
                    rows.append(row)
    
    # Map new column names to expected internal names
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
        # Create normalized row with mapped column names
        norm_row = {}
        for old_col, new_col in column_mapping.items():
            if old_col in row and row[old_col]:
                norm_row[new_col] = row[old_col]
        
        # Copy all other columns as-is
        for col, val in row.items():
            if col not in column_mapping and col not in norm_row:
                norm_row[col] = val
        
        # Determine sample identifier
        sample_id = norm_row.get('sample') or row.get('file_name') or row.get('sample')
        if not sample_id:
            raise ValueError(f"Row missing sample identifier: {row}")
        
        normalized[sample_id] = norm_row
    
    return normalized


def proteins_for(sample, samplesheet):
    """Get protein FASTA path for a sample.
    
    Args:
        sample: Sample identifier (genome_id)
        samplesheet: Dictionary returned by load_samplesheet()
        
    Returns:
        str: Path to protein FASTA file
        
    Raises:
        ValueError: If sample not found or no protein column exists
    """
    row = samplesheet.get(sample)
    if row is None:
        raise ValueError(f"Sample {sample} is missing from the samplesheet.")
    for column in PROTEIN_COLUMNS:
        path = row.get(column)
        if path:
            return path
    raise ValueError(f"No protein FASTA column found for {sample}. Expected one of {PROTEIN_COLUMNS}.")


def genome_for(sample, samplesheet):
    """Get genome FASTA path for a sample.
    
    Args:
        sample: Sample identifier (genome_id)
        samplesheet: Dictionary returned by load_samplesheet()
        
    Returns:
        str: Path to genome FASTA file
        
    Raises:
        ValueError: If sample not found or no genome column exists
    """
    row = samplesheet.get(sample)
    if row is None:
        raise ValueError(f"Sample {sample} is missing from the samplesheet.")
    for column in GENOME_COLUMNS:
        path = row.get(column)
        if path:
            return path
    raise ValueError(f"No genome FASTA column found for {sample}. Expected one of {GENOME_COLUMNS}.")


def gene_fasta_for(sample, samplesheet):
    """Get gene nucleotide FASTA path for a sample.
    
    Args:
        sample: Sample identifier (genome_id)
        samplesheet: Dictionary returned by load_samplesheet()
        
    Returns:
        str: Path to gene FASTA file or None if not present
    """
    row = samplesheet.get(sample)
    if row is None:
        raise ValueError(f"Sample {sample} is missing from the samplesheet.")
    for column in GENE_COLUMNS:
        path = row.get(column)
        if path:
            return path
    return None


def mges_for(sample, samplesheet):
    """Get MGE GFF path for a sample.
    
    Args:
        sample: Sample identifier (genome_id)
        samplesheet: Dictionary returned by load_samplesheet()
        
    Returns:
        str: Path to MGE GFF file or None if not present
    """
    row = samplesheet.get(sample)
    if row is None:
        raise ValueError(f"Sample {sample} is missing from the samplesheet.")
    for column in MGE_COLUMNS:
        path = row.get(column)
        if path:
            return path
    return None


def get_run_id():
    """Generate a run identifier based on version and date.
    
    Returns:
        str: Run ID in format 'v2_YYYYMMDD'
    """
    date_str = datetime.now().strftime("%Y%m%d")
    return f"v2_{date_str}"

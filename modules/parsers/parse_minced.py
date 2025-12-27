#!/usr/bin/env python3
"""
Parse minced outputs to extract CRISPR arrays and spacers with deterministic ID assignment.

Outputs:
- crispr_arrays.tsv (with final IDs, ready for array-to-cas linking)
- crispr_spacers.tsv (with final IDs based on parent array)
"""
import pandas as pd
import math
from collections import Counter
from typing import Dict, List, Optional
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq


def calculate_entropy(seq):
    """Calculate 3-mer entropy for a sequence."""
    if len(seq) < 3:
        return 0.0
    kmers = [seq[i:i+3] for i in range(len(seq)-2)]
    counts = Counter(kmers)
    total = len(kmers)
    entropy = 0.0
    for count in counts.values():
        p = count / total
        entropy -= p * math.log2(p)
    return round(entropy, 6)


def load_genome_sequences(genome_fasta: str) -> Dict[str, str]:
    """Load genome sequences into a dict keyed by contig ID."""
    sequences = {}
    for record in SeqIO.parse(genome_fasta, 'fasta'):
        sequences[record.id] = str(record.seq).upper()
    return sequences


def determine_array_strand_from_spacer(
    contig_seq: str, 
    spacer_seq: str, 
    spacer_start: int, 
    spacer_end: int
) -> str:
    """Determine CRISPR array strand by comparing spacer sequence with genomic sequence.
    
    MinCED extracts spacer sequences but doesn't indicate strand. We compare the
    spacer sequence from the FASTA with the genomic sequence at that position to
    determine strand orientation.

    
    
    Args:
        contig_seq: Full contig sequence (uppercase)
        spacer_seq: Spacer sequence from MinCED FASTA (uppercase)
        spacer_start: 1-based start position of spacer
        spacer_end: 1-based end position of spacer
    
    Returns:
        '+' if forward strand, '-' if reverse strand, '.' if undetermined
    """
    if not spacer_seq or not contig_seq:
        return '.'
    
    # Convert to 0-based indexing
    start_idx = spacer_start - 1
    end_idx = spacer_end  # Already exclusive for Python slicing
    
    # Check bounds
    if start_idx < 0 or end_idx > len(contig_seq):
        return '.'
    
    # Extract genomic sequence at spacer position
    genomic_seq = contig_seq[start_idx:end_idx]
    spacer_upper = spacer_seq.upper()
    
    # Compare with forward orientation
    if genomic_seq == spacer_upper:
        return '+'
    
    # Compare with reverse complement
    spacer_revcomp = str(Seq(spacer_upper).reverse_complement())
    if genomic_seq == spacer_revcomp:
        return '-'
    
    # If neither matches exactly, return undetermined
    return '.'


def parse_arrays_gff(gff_path: str, genome_id: str):
    """Parse minced arrays GFF file and collect repeat units."""
    arrays: Dict[str, dict] = {}
    repeat_units: Dict[str, List[dict]] = {}

    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            contig, source, feature_type, start, end, score, strand, phase, attributes = fields
            start_i, end_i = int(start), int(end)

            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value

            if feature_type == 'repeat_region':
                array_id = attr_dict.get('ID', f"{contig}_{start}_{end}")
                consensus_repeat = attr_dict.get('rpt_unit_seq', '')
                arrays[array_id] = {
                    'genome_id': genome_id,
                    'feature_type': 'CRISPR_array',
                    'tool_id': array_id,
                    'ID': None,  # Will be assigned deterministically
                    'Parent': None,  # Will be set after array-to-cas linking
                    'contig': contig,
                    'start': start_i,
                    'end': end_i,
                    'strand': strand if strand in ['+', '-'] else '.',
                    'score': '.',
                    'rpt_type': attr_dict.get('rpt_type', ''),
                    'rpt_family': attr_dict.get('rpt_family', ''),
                    'repeat_sequence': consensus_repeat,
                    'repeat_length': len(consensus_repeat) if consensus_repeat else None,
                    'spacer_length_avg': None,  # Will be calculated from spacers
                    'is_orphan': True,  # Initially all orphan, updated by link_arrays_to_cas
                    'array_type': 'putative',  # Updated by link_arrays_to_cas
                    'distance_to_cas': None,  # Updated by link_arrays_to_cas
                    'mge_id': None,  # Will be assigned if overlaps with MGE
                    'repeat_units': []
                }
            elif feature_type == 'repeat_unit':
                parent_id = attr_dict.get('Parent')
                if parent_id:
                    repeat_units.setdefault(parent_id, []).append({'start': start_i, 'end': end_i})

    # Merge repeat units and derive spacer coordinates
    for array_id, array_info in arrays.items():
        units = sorted(repeat_units.get(array_id, []), key=lambda x: x['start'])
        spacers = []
        for idx in range(len(units) - 1):
            spacer_start = units[idx]['end'] + 1
            spacer_end = units[idx + 1]['start'] - 1
            spacers.append({'start': spacer_start, 'end': spacer_end})
        array_info['spacer_coords'] = spacers
        array_info['spacer_count'] = len(spacers)

    return arrays


def parse_spacers_and_determine_strands(
    fa_path: str, 
    arrays: Dict[str, dict], 
    genome_id: str, 
    array_id_map: Dict[str, str],
    genome_sequences: Dict[str, str]
):
    """Parse minced spacers FASTA file and determine array strands.
    
    Uses the first spacer of each array to determine strand by comparing
    the spacer sequence with the genomic sequence at that position.
    
    Args:
        fa_path: Path to minced spacers FASTA
        arrays: Dict of parsed arrays (will be modified to add strand)
        genome_id: Genome identifier
        array_id_map: Mapping from tool_id to deterministic ID
        genome_sequences: Dict mapping contig IDs to sequences
    
    Returns:
        pd.DataFrame with spacers including final IDs and strand
    """
    spacers = []
    
    # Track first spacer for each array to determine strand
    first_spacer_per_array: Dict[str, dict] = {}

    for record in SeqIO.parse(fa_path, 'fasta'):
        header = record.id  # e.g., MGYG000001399_1_CRISPR_1_spacer_1
        sequence = str(record.seq)

        parts = header.split('_CRISPR_')
        if len(parts) != 2:
            continue

        contig_part = parts[0]
        array_part = parts[1].split('_spacer_')
        if len(array_part) != 2:
            continue

        array_idx = array_part[0]
        spacer_idx_str = array_part[1]

        tool_id = f"CRISPR{array_idx}"
        array_info = arrays.get(tool_id, {})
        contig = array_info.get('contig', contig_part)
        spacer_coords = array_info.get('spacer_coords', [])

        try:
            spacer_num = int(spacer_idx_str)
        except ValueError:
            spacer_num = None

        coord = spacer_coords[spacer_num - 1] if spacer_num and spacer_num - 1 < len(spacer_coords) else {}
        
        # Get deterministic array ID
        det_array_id = array_id_map.get(tool_id)
        if not det_array_id:
            continue
        
        # Track first spacer for strand determination
        if tool_id not in first_spacer_per_array and spacer_num == 1:
            first_spacer_per_array[tool_id] = {
                'sequence': sequence,
                'start': coord.get('start'),
                'end': coord.get('end'),
                'contig': contig
            }
            
        # Build spacer ID
        spacer_id = f"{det_array_id}_sp{spacer_num:03d}"
        
        spacers.append({
            'genome_id': genome_id,
            'feature_type': 'CRISPR_spacer',
            'ID': spacer_id,
            'Parent': det_array_id,
            'contig': contig,
            'start': coord.get('start'),
            'end': coord.get('end'),
            'strand': '.',  # Will be updated after strand determination
            'tool_id': tool_id,
            'spacer_idx': spacer_num,
            'sequence': sequence,
            'length': len(sequence),
            'entropy': calculate_entropy(sequence)
        })
    
    # Determine strand for each array using first spacer
    for tool_id, first_spacer in first_spacer_per_array.items():
        contig = first_spacer['contig']
        spacer_seq = first_spacer['sequence']
        spacer_start = first_spacer['start']
        spacer_end = first_spacer['end']
        
        if contig and spacer_seq and spacer_start and spacer_end:
            contig_seq = genome_sequences.get(contig)
            if contig_seq:
                strand = determine_array_strand_from_spacer(
                    contig_seq, spacer_seq, spacer_start, spacer_end
                )
                # Update array info
                if tool_id in arrays:
                    arrays[tool_id]['strand'] = strand
    
    # Update spacer strands from their parent arrays
    spacers_df = pd.DataFrame(spacers)
    if not spacers_df.empty:
        for idx, row in spacers_df.iterrows():
            tool_id = row['tool_id']
            if tool_id in arrays:
                spacers_df.at[idx, 'strand'] = arrays[tool_id].get('strand', '.')

    return spacers_df


def assign_mge_overlap_arrays(arrays_df, mges_df):
    """Assign mge_id to arrays based on coordinate overlap with MGEs."""
    if arrays_df.empty or mges_df.empty:
        return arrays_df
    
    arrays_df = arrays_df.copy()
    arrays_df['mge_id'] = arrays_df['mge_id'].astype('object')
    
    for idx, arr_row in arrays_df.iterrows():
        for _, mge_row in mges_df.iterrows():
            if (arr_row['contig'] == mge_row['contig'] and
                arr_row['start'] >= mge_row['start'] and
                arr_row['end'] <= mge_row['end']):
                arrays_df.at[idx, 'mge_id'] = mge_row['ID']
                break
    
    return arrays_df


def get_empty_arrays_df():
    """Return empty DataFrame with array schema."""
    return pd.DataFrame(columns=[
        'genome_id', 'feature_type', 'ID', 'Parent', 'contig', 'start', 'end', 'strand',
        'score', 'repeat_sequence', 'repeat_length', 'spacer_length_avg',
        'is_orphan', 'array_type', 'distance_to_cas', 'mge_id', 'n_spacers', 'tool_id'
    ])


def get_empty_spacers_df():
    """Return empty DataFrame with spacer schema."""
    return pd.DataFrame(columns=[
        'genome_id', 'feature_type', 'ID', 'Parent', 'contig', 'start', 'end', 'strand',
        'tool_id', 'spacer_idx', 'sequence', 'length', 'entropy'
    ])


def main(snakemake):
    """Main function for Snakemake script."""
    arrays_gff = snakemake.input.arrays_gff
    spacers_fa = snakemake.input.spacers_fa
    genome_fasta = snakemake.input.genome_fasta
    genome_id = snakemake.params.genome_id
    arrays_tsv = snakemake.output.crispr_arrays_tsv
    spacers_tsv = snakemake.output.crispr_spacers_tsv
    
    # Get MGEs TSV if available
    mges_tsv = snakemake.input.get('mges_tsv', None)
    mges_df = pd.DataFrame()
    if mges_tsv and Path(mges_tsv).exists():
        mges_df = pd.read_csv(mges_tsv, sep='\t')
    
    # Load genome sequences for strand determination
    genome_sequences = load_genome_sequences(genome_fasta)
    
    # Parse arrays from GFF
    arrays = parse_arrays_gff(arrays_gff, genome_id)
    
    # Build arrays DataFrame and assign deterministic IDs (before parsing spacers)
    arrays_list = []
    for arr in arrays.values():
        arrays_list.append({
            'genome_id': arr['genome_id'],
            'feature_type': arr['feature_type'],
            'tool_id': arr['tool_id'],
            'ID': None,
            'Parent': arr['Parent'],
            'contig': arr['contig'],
            'start': arr['start'],
            'end': arr['end'],
            'strand': '.',  # Will be updated after spacer parsing
            'score': arr['score'],
            'repeat_sequence': arr['repeat_sequence'],
            'repeat_length': arr['repeat_length'],
            'spacer_length_avg': None,
            'is_orphan': arr['is_orphan'],
            'array_type': arr['array_type'],
            'distance_to_cas': arr['distance_to_cas'],
            'mge_id': arr['mge_id'],
            'n_spacers': arr.get('spacer_count')
        })
    
    arrays_df = pd.DataFrame(arrays_list) if arrays_list else get_empty_arrays_df()
    
    # Sort and assign deterministic IDs to arrays
    array_id_map = {}
    if not arrays_df.empty:
        arrays_df = arrays_df.sort_values(['contig', 'start']).reset_index(drop=True)
        arrays_df['ID'] = arrays_df['ID'].astype('object')
        
        for i, (idx, row) in enumerate(arrays_df.iterrows()):
            det_id = f"{genome_id}__CRA{i+1:04d}"
            array_id_map[row['tool_id']] = det_id
            arrays_df.at[idx, 'ID'] = det_id
    
    # Parse spacers and determine strand using first spacer of each array
    spacers_df = parse_spacers_and_determine_strands(
        spacers_fa, arrays, genome_id, array_id_map, genome_sequences
    )
    
    # Update arrays_df with determined strands
    if not arrays_df.empty:
        for idx, row in arrays_df.iterrows():
            tool_id = row['tool_id']
            if tool_id in arrays:
                arrays_df.at[idx, 'strand'] = arrays[tool_id].get('strand', '.')
    
    # Calculate average spacer length per array
    if not spacers_df.empty:
        spacer_avg_lengths = spacers_df.groupby('Parent')['length'].mean().to_dict()
        for idx, row in arrays_df.iterrows():
            avg_len = spacer_avg_lengths.get(row['ID'])
            if avg_len:
                arrays_df.at[idx, 'spacer_length_avg'] = round(avg_len, 2)
    
    # Assign MGE overlaps to arrays
    arrays_df = assign_mge_overlap_arrays(arrays_df, mges_df)
    
    # Write outputs
    arrays_df.to_csv(arrays_tsv, sep='\t', index=False, na_rep='NULL')
    spacers_df.to_csv(spacers_tsv, sep='\t', index=False, na_rep='NULL')


if __name__ == '__main__':
    main(snakemake)

#!/usr/bin/env python3
"""
Parse CRISPRDetect 3.0 outputs to extract CRISPR arrays and spacers with deterministic ID assignment.

Parsing strategy:
- Arrays and spacers are extracted from the TXT output (authoritative source)
- The GFF output from CRISPRDetect can contain duplicate arrays with different IDs
  for the same physical location, so we avoid using it for array parsing
- The TXT Summary line contains: ID_START_STOP_DIR, repeat sequence, spacers, score

TXT format provides:
- Array coordinates from Summary line (deduplicated by design)
- Spacer sequences (comma-separated in SPACERS field)
- Array quality score
- Questionable array flag
- Direction confidence (HIGH/MEDIUM/LOW)

Outputs:
- crispr_arrays.tsv (with final IDs, ready for array-to-cas linking)
- crispr_spacers.tsv (with final IDs based on parent array)
"""
import pandas as pd
import math
import re
from collections import Counter
from typing import Dict, List, Tuple
from pathlib import Path


def calculate_entropy(seq: str) -> float:
    """Calculate 3-mer entropy for a sequence.
    
    Low entropy indicates repetitive/low-complexity sequence.
    Spacers with entropy > 1.0 are retained for analysis.
    """
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


def parse_crisprdetect_txt(txt_path: str, genome_id: str) -> Tuple[List[dict], List[dict]]:
    """
    Parse CRISPRDetect text output to extract arrays and spacers.
    
    This is the authoritative source for array/spacer data. The GFF output from
    CRISPRDetect can contain duplicate arrays with different IDs for the same
    physical location, so we parse from txt instead.
    
    The Summary line format:
    ID_START_STOP_DIR: CONTIG-START-END-DIR; REP_DR:sequence; DR_LENGTH:N; 
    NUM_DRs:N; MUTATIONS:N; INSERTIONS:N; SCORE:N.NN; SPACERS:seq1,seq2,...
    
    Args:
        txt_path: Path to CRISPRDetect .txt output
        genome_id: Genome identifier for ID prefixing
    
    Returns:
        Tuple of (list of array dicts, list of spacer dicts)
    """
    arrays: List[dict] = []
    spacers: List[dict] = []
    
    if not Path(txt_path).exists() or Path(txt_path).stat().st_size == 0:
        return arrays, spacers
    
    with open(txt_path) as f:
        content = f.read()
    
    # Split into array blocks (separated by //)
    blocks = content.split('//')
    
    # Track unique arrays by (contig, start, end) to deduplicate
    seen_coords: set = set()
    
    for block in blocks:
        if not block.strip():
            continue
        
        # Extract contig from header line
        # Format: >MGYG000001399_1		Array_Orientation: Reverse
        header_match = re.search(r'>(\S+)\s+Array_Orientation:', block)
        if not header_match:
            continue
        contig = header_match.group(1)
        
        # Get array data from Summary line
        # Format: ID_START_STOP_DIR: MGYG000000255_1-138327-139272-F; REP_DR:...; ...
        summary_match = re.search(
            r'ID_START_STOP_DIR:\s*(\S+)-(\d+)-(\d+)-([FR]);'
            r'\s*REP_DR:([^;]+);'
            r'\s*DR_LENGTH:(\d+);'
            r'\s*NUM_DRs:(\d+);'
            r'.*?SCORE:([\d.]+);'
            r'\s*SPACERS:([^\n]*)',
            block, re.DOTALL
        )
        if not summary_match:
            continue
        
        start = int(summary_match.group(2))
        end = int(summary_match.group(3))
        direction = summary_match.group(4)  # F or R
        repeat_sequence = summary_match.group(5).strip()
        repeat_length = int(summary_match.group(6))
        # num_repeats available in group(7) but not used
        quality_score = float(summary_match.group(8))
        spacers_str = summary_match.group(9).strip()
        
        # Convert direction to strand
        strand = '+' if direction == 'F' else '-'
        
        # Deduplicate by coordinates
        coord_key = (contig, start, end)
        if coord_key in seen_coords:
            continue
        seen_coords.add(coord_key)
        
        # Extract questionable array flag
        questionable_match = re.search(r'# Questionable array\s*:\s*(YES|NO)', block)
        questionable = questionable_match.group(1) == 'YES' if questionable_match else False
        
        # Extract direction confidence
        confidence_match = re.search(r'Final direction:.*Confidence:\s*(HIGH|MEDIUM|LOW)', block)
        direction_confidence = confidence_match.group(1) if confidence_match else None
        
        # Create unique tool_id based on coordinates
        tool_id = f"{contig}_{start}_{end}"
        
        # Parse spacer sequences
        spacer_seqs = [s.strip() for s in spacers_str.split(',') if s.strip()]
        n_spacers = len(spacer_seqs)
        
        # Calculate average spacer length
        spacer_length_avg = None
        if spacer_seqs:
            spacer_length_avg = round(sum(len(s) for s in spacer_seqs) / len(spacer_seqs), 2)
        
        arrays.append({
            'genome_id': genome_id,
            'feature_type': 'CRISPR_array',
            'tool_id': tool_id,
            'ID': None,  # Will be assigned deterministically
            'Parent': None,  # Will be set after array-to-cas linking
            'contig': contig,
            'start': start,
            'end': end,
            'strand': strand,
            'score': '.',
            'repeat_sequence': repeat_sequence,
            'repeat_length': repeat_length,
            'spacer_length_avg': spacer_length_avg,
            'is_orphan': True,  # Initially all orphan, updated by link_arrays_to_cas
            'array_type': 'putative',  # Updated by link_arrays_to_cas
            'distance_to_cas': None,
            'mge_id': None,  # Will be assigned if overlaps with MGE
            'n_spacers': n_spacers,
            'quality_score': quality_score,
            'questionable_array': questionable,
            'direction_confidence': direction_confidence,
            'cas_strand': None,
            'is_revcomp': False
        })
        
        # Create spacer entries
        # For reverse strand arrays, spacers are listed 5'->3' in the txt,
        # but we need to calculate genomic coordinates
        for i, spacer_seq in enumerate(spacer_seqs, start=1):
            spacers.append({
                'genome_id': genome_id,
                'feature_type': 'CRISPR_spacer',
                'ID': None,  # Will be assigned after array ID assignment
                'Parent': None,  # Will be set to deterministic array ID
                'parent_tool_id': tool_id,
                'contig': contig,
                'start': None,  # Coordinates calculated from GFF if needed
                'end': None,
                'strand': strand,
                'tool_id': tool_id,
                'spacer_idx': i,
                'sequence': spacer_seq,
                'length': len(spacer_seq),
                'entropy': calculate_entropy(spacer_seq),
                'is_revcomp': False
            })
    
    return arrays, spacers


def parse_gff_attributes(attr_string: str) -> Dict[str, str]:
    """Parse GFF9 attribute string into dict."""
    attrs = {}
    for attr in attr_string.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attrs[key] = value
    return attrs


def parse_spacer_coordinates_from_gff(gff_path: str) -> Dict[str, List[Tuple[int, int]]]:
    """
    Parse CRISPRDetect GFF file to extract spacer coordinates only.
    
    Since the txt file doesn't provide individual spacer coordinates,
    we extract them from the GFF. We group by array coordinates to
    match with txt-parsed arrays.
    
    IMPORTANT: The GFF feature coordinates may differ slightly from the TXT Summary
    coordinates (off-by-one issues). The GFF ID attribute (e.g., CRISPR1_127935_130746)
    contains coordinates that match the TXT Summary line, so we extract coordinates
    from the ID string rather than from the GFF feature columns.
    
    Args:
        gff_path: Path to CRISPRDetect GFF output
    
    Returns:
        Dict mapping "contig_start_end" to list of (spacer_start, spacer_end) tuples
        ordered by position within the array
    """
    spacer_coords: Dict[str, List[Tuple[int, int]]] = {}
    
    if not Path(gff_path).exists() or Path(gff_path).stat().st_size == 0:
        return spacer_coords
    
    # First pass: map array GFF IDs to their coord_key (matching TXT format)
    array_id_to_coords: Dict[str, str] = {}
    
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            contig, source, feature_type, start, end, score, strand, phase, attributes = fields
            attrs = parse_gff_attributes(attributes)
            
            if feature_type == 'repeat_region':
                gff_id = attrs.get('ID', '')
                # Extract start/end from GFF ID attribute (e.g., CRISPR1_127935_130746)
                # This matches the TXT Summary line coordinates
                id_match = re.search(r'CRISPR\d+_(\d+)_(\d+)', gff_id)
                if id_match:
                    id_start = int(id_match.group(1))
                    id_end = int(id_match.group(2))
                    coord_key = f"{contig}_{id_start}_{id_end}"
                    array_id_to_coords[gff_id] = coord_key
    
    # Second pass: collect spacer coordinates per array
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            contig, source, feature_type, start, end, score, strand, phase, attributes = fields
            
            if feature_type == 'binding_site':
                attrs = parse_gff_attributes(attributes)
                parent_id = attrs.get('Parent', '')
                
                if parent_id in array_id_to_coords:
                    coord_key = array_id_to_coords[parent_id]
                    start_i, end_i = int(start), int(end)
                    
                    if coord_key not in spacer_coords:
                        spacer_coords[coord_key] = []
                    spacer_coords[coord_key].append((start_i, end_i))
    
    # Sort spacers by start position within each array
    for coord_key in spacer_coords:
        spacer_coords[coord_key].sort(key=lambda x: x[0])
    
    return spacer_coords


def assign_mge_overlap_arrays(arrays_df: pd.DataFrame, mges_df: pd.DataFrame) -> pd.DataFrame:
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


def get_empty_arrays_df() -> pd.DataFrame:
    """Return empty DataFrame with array schema."""
    return pd.DataFrame(columns=[
        'genome_id', 'feature_type', 'ID', 'Parent', 'contig', 'start', 'end', 'strand',
        'score', 'repeat_sequence', 'repeat_length', 'spacer_length_avg',
        'is_orphan', 'array_type', 'distance_to_cas', 'mge_id', 'n_spacers', 'tool_id',
        'quality_score', 'questionable_array', 'direction_confidence',
        'cas_strand', 'is_revcomp'
    ])


def get_empty_spacers_df() -> pd.DataFrame:
    """Return empty DataFrame with spacer schema."""
    return pd.DataFrame(columns=[
        'genome_id', 'feature_type', 'ID', 'Parent', 'contig', 'start', 'end', 'strand',
        'tool_id', 'spacer_idx', 'sequence', 'length', 'entropy', 'is_revcomp'
    ])


def main(snakemake):
    """Main function for Snakemake script.
    
    Parses arrays and spacers from CRISPRDetect txt output (authoritative source),
    then enriches spacers with genomic coordinates from GFF.
    """
    gff_path = snakemake.input.crisprdetect_gff
    txt_path = snakemake.input.crisprdetect_txt
    genome_id = snakemake.params.genome_id
    arrays_tsv = snakemake.output.crispr_arrays_tsv
    spacers_tsv = snakemake.output.crispr_spacers_tsv
    
    # Get MGEs TSV if available
    mges_tsv = snakemake.input.get('mges_tsv', None)
    mges_df = pd.DataFrame()
    if mges_tsv and Path(mges_tsv).exists():
        mges_df = pd.read_csv(mges_tsv, sep='\t')
    
    # Parse arrays and spacers from txt (authoritative, deduplicated)
    arrays_list, spacers_list = parse_crisprdetect_txt(txt_path, genome_id)
    
    # Get spacer coordinates from GFF to enrich txt-parsed spacers
    spacer_coords = parse_spacer_coordinates_from_gff(gff_path)
    
    # Assign coordinates to spacers from GFF data
    for spacer in spacers_list:
        tool_id = spacer['tool_id']
        spacer_idx = spacer['spacer_idx']
        
        if tool_id in spacer_coords:
            coords_list = spacer_coords[tool_id]
            # spacer_idx is 1-based, coords_list is 0-based
            if 0 < spacer_idx <= len(coords_list):
                start, end = coords_list[spacer_idx - 1]
                spacer['start'] = start
                spacer['end'] = end
    
    # Build arrays DataFrame
    arrays_df = pd.DataFrame(arrays_list) if arrays_list else get_empty_arrays_df()
    
    # Sort and assign deterministic IDs to arrays
    array_id_map = {}  # tool_id -> deterministic ID
    if not arrays_df.empty:
        arrays_df = arrays_df.sort_values(['contig', 'start']).reset_index(drop=True)
        arrays_df['ID'] = arrays_df['ID'].astype('object')
        
        for i, (idx, row) in enumerate(arrays_df.iterrows()):
            det_id = f"{genome_id}__CRA{i+1:04d}"
            array_id_map[row['tool_id']] = det_id
            arrays_df.at[idx, 'ID'] = det_id
    
    # Build spacers DataFrame with deterministic IDs
    spacers_with_ids = []
    for spacer in spacers_list:
        parent_tool_id = spacer['parent_tool_id']
        det_array_id = array_id_map.get(parent_tool_id)
        
        if det_array_id:
            spacer_id = f"{det_array_id}_sp{spacer['spacer_idx']:03d}"
            spacers_with_ids.append({
                'genome_id': spacer['genome_id'],
                'feature_type': spacer['feature_type'],
                'ID': spacer_id,
                'Parent': det_array_id,
                'contig': spacer['contig'],
                'start': spacer['start'],
                'end': spacer['end'],
                'strand': spacer['strand'],
                'tool_id': spacer['tool_id'],
                'spacer_idx': spacer['spacer_idx'],
                'sequence': spacer['sequence'],
                'length': spacer['length'],
                'entropy': spacer['entropy'],
                'is_revcomp': False  # Updated by link_arrays_to_cas if needed
            })
    
    spacers_df = pd.DataFrame(spacers_with_ids) if spacers_with_ids else get_empty_spacers_df()
    
    # Assign MGE overlaps to arrays
    arrays_df = assign_mge_overlap_arrays(arrays_df, mges_df)
    
    # Write outputs
    arrays_df.to_csv(arrays_tsv, sep='\t', index=False, na_rep='NULL')
    spacers_df.to_csv(spacers_tsv, sep='\t', index=False, na_rep='NULL')


if __name__ == '__main__':
    main(snakemake)  # noqa: F821

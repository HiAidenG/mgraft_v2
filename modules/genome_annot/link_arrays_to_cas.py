#!/usr/bin/env python3
"""
Link CRISPR arrays to Cas systems based on proximity with strand-aware processing.

This script:
1. Finds nearest cas system for each array (within max_distance)
2. Updates arrays with: Parent, is_orphan, array_type, distance_to_cas, cas_strand
3. Updates cas_systems with: array_ID
4. Optionally reverse-complements arrays/spacers when on opposite strand from cas system
5. Resolves unknown strand (.) using seqkit locate on genome

Key behaviors:
- When array strand differs from assigned cas system strand, set is_revcomp=True
- Reverse-complement repeat_sequence (arrays) and sequence (spacers)
- Update strand column to reflect corrected orientation after revcomp
- For arrays with strand='.', determine strand by locating repeat sequence in genome
"""
import pandas as pd
import subprocess
import tempfile
from pathlib import Path
from Bio.Seq import Seq


def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence using Biopython."""
    if not seq or seq == 'NULL' or pd.isna(seq):
        return seq
    return str(Seq(seq).reverse_complement())


def get_strand_from_seqkit(repeat_seq: str, genome_fasta: str, contig: str, 
                           start: int, end: int) -> str:
    """
    Determine array strand by locating repeat sequence in genome using seqkit.
    
    Args:
        repeat_seq: Consensus repeat sequence from CRISPRDetect
        genome_fasta: Path to genome FASTA file
        contig: Contig/chromosome name
        start: Array start position
        end: Array end position
    
    Returns:
        '+' or '-' based on where repeat is found, or '.' if undetermined
    """
    if not repeat_seq or repeat_seq == 'NULL' or pd.isna(repeat_seq):
        return '.'
    
    if not genome_fasta or not Path(genome_fasta).exists():
        return '.'
    
    # Create temporary file with repeat sequence as pattern
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as tmp:
            tmp.write(f">repeat\n{repeat_seq}\n")
            pattern_file = tmp.name
        
        # Run seqkit locate to find the repeat sequence
        # Using --pattern-file for exact match search
        result = subprocess.run(
            ['seqkit', 'locate', '-f', pattern_file, '-P', genome_fasta],
            capture_output=True, text=True, timeout=30
        )
        
        Path(pattern_file).unlink()  # Clean up temp file
        
        if result.returncode != 0:
            return '.'
        
        # Parse seqkit locate output (TSV format)
        # Header: seqID	patternName	pattern	strand	start	end
        lines = result.stdout.strip().split('\n')
        if len(lines) < 2:  # No matches (only header or empty)
            return '.'
        
        # Find matches within the array region
        plus_hits = 0
        minus_hits = 0
        
        for line in lines[1:]:
            fields = line.split('\t')
            if len(fields) < 6:
                continue
            
            seq_id = fields[0]
            strand = fields[3]
            match_start = int(fields[4])
            match_end = int(fields[5])
            
            # Check if match is on the correct contig and within array bounds
            if seq_id == contig and match_start >= start and match_end <= end:
                if strand == '+':
                    plus_hits += 1
                elif strand == '-':
                    minus_hits += 1
        
        # Determine strand based on majority of hits
        if plus_hits > minus_hits:
            return '+'
        elif minus_hits > plus_hits:
            return '-'
        else:
            return '.'
            
    except (subprocess.TimeoutExpired, FileNotFoundError, Exception):
        return '.'


def link_arrays_to_cas(cas_systems_df, cas_genes_df, arrays_df, spacers_df,
                       max_distance, revcomp_if_cas_opposite=True,
                       genome_fasta=None):
    """Link CRISPR arrays to Cas systems based on proximity with strand handling.
    
    IMPORTANT: Only ONE array can be assigned to each cas system. When multiple
    arrays are within range of the same cas system, we select the best one using:
    1. Closest distance to cas genes
    2. Tie-break: prefer arrays on the same strand as the cas system
    3. Second tie-break: prefer arrays with higher CRISPRDetect quality score
    
    Args:
        cas_systems_df: DataFrame with cas_system features
        cas_genes_df: DataFrame with cas_gene features  
        arrays_df: DataFrame with CRISPR_array features
        spacers_df: DataFrame with CRISPR_spacer features
        max_distance: Maximum distance (bp) to consider for assignment
        revcomp_if_cas_opposite: If True, reverse-complement arrays on opposite strand
        genome_fasta: Path to genome FASTA for strand resolution (optional)
        
    Returns:
        Tuple of (updated_arrays_df, updated_cas_systems_df, updated_spacers_df)
    """
    if arrays_df.empty:
        # Add missing columns to empty dataframe
        arrays_df['cas_strand'] = None
        arrays_df['is_revcomp'] = False
        if not spacers_df.empty:
            spacers_df = spacers_df.copy()
            spacers_df['is_revcomp'] = False
        return arrays_df, cas_systems_df, spacers_df
    
    # Ensure columns are object dtype for string IDs
    arrays_df = arrays_df.copy()
    arrays_df['Parent'] = arrays_df['Parent'].astype('object')
    arrays_df['cas_strand'] = None
    arrays_df['cas_strand'] = arrays_df['cas_strand'].astype('object')
    arrays_df['is_revcomp'] = False
    
    spacers_df = spacers_df.copy() if not spacers_df.empty else spacers_df
    if not spacers_df.empty:
        spacers_df['is_revcomp'] = False
    
    if not cas_systems_df.empty:
        cas_systems_df = cas_systems_df.copy()
        cas_systems_df['array_ID'] = None
        cas_systems_df['array_ID'] = cas_systems_df['array_ID'].astype('object')
    
    # First pass: collect all candidate (array, cas_system, distance, strand_match) tuples
    candidates = []
    
    for idx, array_row in arrays_df.iterrows():
        array_id = array_row['ID']
        array_contig = array_row['contig']
        array_start = array_row['start']
        array_end = array_row['end']
        array_midpoint = (array_start + array_end) / 2
        array_strand = array_row['strand']
        
        # Resolve unknown strand using seqkit if genome is available
        if array_strand == '.' and genome_fasta:
            repeat_seq = array_row.get('repeat_sequence')
            resolved_strand = get_strand_from_seqkit(
                repeat_seq, genome_fasta, array_contig, array_start, array_end
            )
            if resolved_strand != '.':
                arrays_df.at[idx, 'strand'] = resolved_strand
                array_strand = resolved_strand
        
        # Find nearest cas gene on same contig
        min_dist = float('inf')
        nearest_gene = None
        
        if not cas_genes_df.empty:
            same_contig_genes = cas_genes_df[cas_genes_df['contig'] == array_contig]
            
            for _, gene_row in same_contig_genes.iterrows():
                gene_start = gene_row['start']
                gene_end = gene_row['end']
                gene_midpoint = (gene_start + gene_end) / 2
                
                dist = abs(array_midpoint - gene_midpoint)
                
                if dist < min_dist:
                    min_dist = dist
                    nearest_gene = gene_row
        
        if min_dist <= max_distance and nearest_gene is not None:
            parent_sys_id = nearest_gene['Parent']
            cas_strand = nearest_gene.get('strand', '.')
            
            # Check strand match for tie-breaking
            strand_match = (array_strand == cas_strand) if array_strand != '.' and cas_strand != '.' else False
            
            # Get quality score for second tie-breaker
            quality_score = array_row.get('quality_score')
            if pd.isna(quality_score):
                quality_score = 0.0
            else:
                quality_score = float(quality_score)
            
            candidates.append({
                'array_idx': idx,
                'array_id': array_id,
                'array_strand': array_strand,
                'system_id': parent_sys_id,
                'distance': min_dist,
                'cas_strand': cas_strand,
                'strand_match': strand_match,
                'quality_score': quality_score
            })
    
    # Second pass: for each cas system, select the best array
    # Best = closest distance, tie-break by same strand, then by quality score
    system_to_best_candidate = {}
    
    for cand in candidates:
        sys_id = cand['system_id']
        
        if sys_id not in system_to_best_candidate:
            system_to_best_candidate[sys_id] = cand
        else:
            existing = system_to_best_candidate[sys_id]
            
            # Compare: prefer closer distance
            if cand['distance'] < existing['distance']:
                system_to_best_candidate[sys_id] = cand
            elif cand['distance'] == existing['distance']:
                # First tie-break: prefer same strand
                if cand['strand_match'] and not existing['strand_match']:
                    system_to_best_candidate[sys_id] = cand
                elif cand['strand_match'] == existing['strand_match']:
                    # Second tie-break: prefer higher quality score
                    if cand['quality_score'] > existing['quality_score']:
                        system_to_best_candidate[sys_id] = cand
    
    # Build set of winning arrays
    winning_arrays = {cand['array_idx']: cand for cand in system_to_best_candidate.values()}
    
    # Track which systems have assigned arrays and arrays that need revcomp
    system_to_array = {}
    arrays_to_revcomp = set()
    
    # Third pass: update arrays with assignment info
    for idx, array_row in arrays_df.iterrows():
        array_id = array_row['ID']
        
        # Initialize as orphan
        array_type = 'putative'
        is_orphan = True
        system_id = None
        distance = None
        cas_strand = None
        
        if idx in winning_arrays:
            cand = winning_arrays[idx]
            parent_sys_id = cand['system_id']
            cas_strand = cand['cas_strand']
            array_strand = cand['array_strand']
            
            if not cas_systems_df.empty:
                sys_row = cas_systems_df[cas_systems_df['ID'] == parent_sys_id]
                
                if not sys_row.empty:
                    sys_row = sys_row.iloc[0]
                    wholeness = sys_row.get('sys_wholeness')
                    
                    try:
                        wholeness_val = float(wholeness) if pd.notnull(wholeness) else None
                    except Exception:
                        wholeness_val = None
                    
                    if wholeness_val == 1.0:
                        array_type = 'canonical'
                    else:
                        array_type = 'putative'
                    
                    is_orphan = False
                    system_id = parent_sys_id
                    distance = int(cand['distance'])
                    
                    # Track assignment for updating cas_systems
                    system_to_array[parent_sys_id] = array_id
                    
                    # Check strand concordance for revcomp
                    if revcomp_if_cas_opposite and array_strand != '.':
                        if cas_strand and cas_strand != '.' and array_strand != cas_strand:
                            arrays_to_revcomp.add(array_id)
        
        # Update array
        arrays_df.at[idx, 'array_type'] = array_type
        arrays_df.at[idx, 'is_orphan'] = is_orphan
        arrays_df.at[idx, 'Parent'] = system_id
        arrays_df.at[idx, 'distance_to_cas'] = distance
        arrays_df.at[idx, 'cas_strand'] = cas_strand
    
    # Apply revcomp to arrays and spacers
    for idx, row in arrays_df.iterrows():
        array_id = row['ID']
        if array_id in arrays_to_revcomp:
            # Reverse complement the repeat sequence
            repeat_seq = row.get('repeat_sequence')
            if repeat_seq and repeat_seq != 'NULL' and pd.notna(repeat_seq):
                arrays_df.at[idx, 'repeat_sequence'] = reverse_complement(repeat_seq)
            
            # Flip the strand to match cas orientation
            current_strand = row['strand']
            if current_strand == '+':
                arrays_df.at[idx, 'strand'] = '-'
            elif current_strand == '-':
                arrays_df.at[idx, 'strand'] = '+'
            
            # Mark as revcomped
            arrays_df.at[idx, 'is_revcomp'] = True
            
            # Process spacers for this array
            if not spacers_df.empty:
                spacer_mask = spacers_df['Parent'] == array_id
                for sp_idx in spacers_df[spacer_mask].index:
                    seq = spacers_df.at[sp_idx, 'sequence']
                    if seq and seq != 'NULL' and pd.notna(seq):
                        spacers_df.at[sp_idx, 'sequence'] = reverse_complement(seq)
                    
                    # Flip spacer strand
                    sp_strand = spacers_df.at[sp_idx, 'strand']
                    if sp_strand == '+':
                        spacers_df.at[sp_idx, 'strand'] = '-'
                    elif sp_strand == '-':
                        spacers_df.at[sp_idx, 'strand'] = '+'
                    
                    spacers_df.at[sp_idx, 'is_revcomp'] = True
    
    # Update cas_systems with array_ID
    if not cas_systems_df.empty:
        for sys_id, arr_id in system_to_array.items():
            mask = cas_systems_df['ID'] == sys_id
            cas_systems_df.loc[mask, 'array_ID'] = arr_id
    
    return arrays_df, cas_systems_df, spacers_df


def main(snakemake):
    """Main function for Snakemake script."""
    # Inputs
    cas_systems_tsv = snakemake.input.cas_systems_tsv
    cas_genes_tsv = snakemake.input.cas_genes_tsv
    arrays_tsv = snakemake.input.crispr_arrays_tsv
    spacers_tsv = snakemake.input.crispr_spacers_tsv
    genome_fasta = snakemake.input.get('genome_fasta', None)
    
    # Parameters
    max_distance = snakemake.params.get('max_distance', 10000)
    revcomp_if_cas_opposite = snakemake.params.get('revcomp_if_cas_opposite', True)
    
    # Outputs
    arrays_out = snakemake.output.crispr_arrays_tsv
    cas_systems_out = snakemake.output.cas_systems_tsv
    spacers_out = snakemake.output.crispr_spacers_tsv
    
    # Load data
    cas_systems_df = pd.read_csv(cas_systems_tsv, sep='\t') if Path(cas_systems_tsv).exists() else pd.DataFrame()
    cas_genes_df = pd.read_csv(cas_genes_tsv, sep='\t') if Path(cas_genes_tsv).exists() else pd.DataFrame()
    arrays_df = pd.read_csv(arrays_tsv, sep='\t') if Path(arrays_tsv).exists() else pd.DataFrame()
    spacers_df = pd.read_csv(spacers_tsv, sep='\t') if Path(spacers_tsv).exists() else pd.DataFrame()
    
    # Link arrays to cas systems with strand handling
    arrays_df, cas_systems_df, spacers_df = link_arrays_to_cas(
        cas_systems_df, cas_genes_df, arrays_df, spacers_df,
        max_distance, revcomp_if_cas_opposite, genome_fasta
    )
    
    # Write outputs
    arrays_df.to_csv(arrays_out, sep='\t', index=False, na_rep='NULL')
    cas_systems_df.to_csv(cas_systems_out, sep='\t', index=False, na_rep='NULL')
    spacers_df.to_csv(spacers_out, sep='\t', index=False, na_rep='NULL')


if __name__ == "__main__":
    main(snakemake)  # noqa: F821

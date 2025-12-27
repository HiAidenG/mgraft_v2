#!/usr/bin/env python3
"""
Map RM (Restriction-Modification) recognition sites to target genomes.

Uses exact string matching (no mismatches) on both strands since RM sites
are typically short palindromic sequences that must be exact matches.

Outputs GFF3 directly with hits classified as mobile (within MGE) or genomic.
"""

import argparse
import csv
import re
from pathlib import Path
from typing import List, Dict, Tuple, Set


# IUPAC ambiguity codes
IUPAC_CODES = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'U': 'T',
    'R': '[AG]',
    'Y': '[CT]',
    'S': '[GC]',
    'W': '[AT]',
    'K': '[GT]',
    'M': '[AC]',
    'B': '[CGT]',
    'D': '[AGT]',
    'H': '[ACT]',
    'V': '[ACG]',
    'N': '[ACGT]',
}


def iupac_to_regex(pattern: str) -> str:
    """Convert an IUPAC pattern to a regex pattern."""
    regex_parts = []
    for char in pattern.upper():
        if char in IUPAC_CODES:
            regex_parts.append(IUPAC_CODES[char])
        else:
            regex_parts.append(re.escape(char))
    return ''.join(regex_parts)


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
                  'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
                  'D': 'H', 'H': 'D', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq.upper()))


def parse_fasta(fasta_path: str) -> Dict[str, str]:
    """
    Parse a FASTA file into a dictionary of sequences.
    
    Args:
        fasta_path: Path to the FASTA file
        
    Returns:
        dict: sequence_id -> sequence (uppercase)
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq).upper()
                current_id = line[1:].split()[0]  # First word only
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq).upper()
    
    return sequences


def load_rm_patterns(patterns_tsv: str) -> List[Dict]:
    """
    Load RM patterns from metadata TSV file.
    
    Args:
        patterns_tsv: Path to patterns_rm.tsv
        
    Returns:
        List of pattern dicts with pattern_id, sequence, rm_type, source info
    """
    patterns = []
    
    if not Path(patterns_tsv).exists():
        return patterns
    
    with open(patterns_tsv, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            patterns.append({
                'pattern_id': row['pattern_id'],
                'sequence': row['sequence'].upper(),
                'recseq': row['sequence'].upper(),  # Explicit recseq
                'rm_type': row.get('rm_type', ''),
                'source_genome_ids': row.get('source_genome_ids', ''),
                'source_system_ids': row.get('source_system_ids', '')
            })
    
    return patterns


def load_mge_bounds(mges_tsv: str, genome_id: str) -> Dict[str, List[Tuple[int, int, str]]]:
    """
    Load MGE bounds for a genome, organized by contig.
    
    Args:
        mges_tsv: Path to all_mges.tsv
        genome_id: Genome ID to filter for
        
    Returns:
        Dict: contig -> list of (start, end, mge_id) tuples
    """
    bounds = {}
    
    if not Path(mges_tsv).exists():
        return bounds
    
    with open(mges_tsv, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['genome_id'] != genome_id:
                continue
            
            contig = row['contig']
            start = int(row['start'])
            end = int(row['end'])
            mge_id = row['ID']
            
            if contig not in bounds:
                bounds[contig] = []
            bounds[contig].append((start, end, mge_id))
    
    # Sort by start position
    for contig in bounds:
        bounds[contig].sort(key=lambda x: x[0])
    
    return bounds


def find_containing_mge(
    contig: str,
    start: int,
    end: int,
    mge_bounds: Dict[str, List[Tuple[int, int, str]]]
) -> str:
    """
    Check if a hit is 100% contained within any MGE.
    
    Returns:
        MGE ID if contained, empty string otherwise
    """
    if contig not in mge_bounds:
        return ''
    
    for mge_start, mge_end, mge_id in mge_bounds[contig]:
        if start >= mge_start and end <= mge_end:
            return mge_id
    
    return ''


def find_pattern_matches(
    pattern: Dict,
    sequences: Dict[str, str],
    mge_bounds: Dict[str, List[Tuple[int, int, str]]]
) -> List[Dict]:
    """
    Find all occurrences of a pattern in sequences on both strands.
    
    Args:
        pattern: Pattern dict with sequence, recseq, rm_type, etc.
        sequences: Dict of contig_id -> sequence
        mge_bounds: MGE bounds for mobility classification
        
    Returns:
        List of hit dictionaries with all metadata
    """
    hits = []
    
    seq = pattern['sequence']
    fwd_regex = iupac_to_regex(seq)
    fwd_pattern = re.compile(fwd_regex)
    
    # Reverse complement for minus strand
    rc_seq = reverse_complement(seq)
    rc_regex = iupac_to_regex(rc_seq)
    rc_pattern = re.compile(rc_regex)
    
    for contig_id, contig_seq in sequences.items():
        # Forward strand matches
        for match in fwd_pattern.finditer(contig_seq):
            start = match.start() + 1  # 1-based
            end = match.end()
            mge_id = find_containing_mge(contig_id, start, end, mge_bounds)
            
            hits.append({
                'contig': contig_id,
                'start': start,
                'end': end,
                'strand': '+',
                'pattern_id': pattern['pattern_id'],
                'recseq': pattern['recseq'],
                'rm_type': pattern['rm_type'],
                'source_genome_ids': pattern['source_genome_ids'],
                'source_system_ids': pattern['source_system_ids'],
                'mge_id': mge_id
            })
        
        # Reverse complement matches (only if not palindromic)
        if rc_regex != fwd_regex:
            for match in rc_pattern.finditer(contig_seq):
                start = match.start() + 1
                end = match.end()
                mge_id = find_containing_mge(contig_id, start, end, mge_bounds)
                
                hits.append({
                    'contig': contig_id,
                    'start': start,
                    'end': end,
                    'strand': '-',
                    'pattern_id': pattern['pattern_id'],
                    'recseq': pattern['recseq'],
                    'rm_type': pattern['rm_type'],
                    'source_genome_ids': pattern['source_genome_ids'],
                    'source_system_ids': pattern['source_system_ids'],
                    'mge_id': mge_id
                })
    
    return hits


def format_gff_attributes(attrs: Dict) -> str:
    """Format attributes as GFF3 attribute string."""
    parts = []
    for key, value in attrs.items():
        if value:
            # URL-encode special characters
            value = str(value).replace('%', '%25').replace(';', '%3B').replace('=', '%3D').replace(',', '%2C')
            parts.append(f"{key}={value}")
    return ';'.join(parts)


def write_gff(
    hits: List[Dict],
    genome_id: str,
    output_gff: str
) -> None:
    """
    Write hits to GFF3 file.
    
    Args:
        hits: List of hit dicts
        genome_id: Target genome ID
        output_gff: Output GFF path
    """
    Path(output_gff).parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_gff, 'w') as f:
        f.write("##gff-version 3\n")
        
        for idx, hit in enumerate(hits, 1):
            attrs = {
                'ID': f"{genome_id}_rm_hit_{idx:06d}",
                'pattern_id': hit['pattern_id'],
                'recseq': hit['recseq'],
                'rm_type': hit['rm_type'],
                'source_genomes': hit['source_genome_ids'],
                'source_systems': hit['source_system_ids']
            }
            
            if hit['mge_id']:
                attrs['within_mge'] = hit['mge_id']
            
            attrs_str = format_gff_attributes(attrs)
            
            # GFF3: seqid, source, type, start, end, score, strand, phase, attributes
            line = f"{hit['contig']}\ttarget_mapping\trestriction_site\t{hit['start']}\t{hit['end']}\t.\t{hit['strand']}\t.\t{attrs_str}"
            f.write(line + '\n')


def map_rm_targets(
    genome_id: str,
    patterns_tsv: str,
    genome_fasta: str,
    mges_tsv: str,
    output_gff: str
) -> None:
    """
    Map RM recognition sites to genome and output GFF.
    
    Args:
        genome_id: Target genome ID
        patterns_tsv: Path to patterns_rm.tsv with metadata
        genome_fasta: Path to genome FASTA
        mges_tsv: Path to all_mges.tsv for mobility classification
        output_gff: Output GFF path
    """
    # Load patterns
    patterns = load_rm_patterns(patterns_tsv)
    if not patterns:
        Path(output_gff).parent.mkdir(parents=True, exist_ok=True)
        with open(output_gff, 'w') as f:
            f.write("##gff-version 3\n")
        print("No RM patterns to search")
        return
    
    # Load genome sequences
    if not Path(genome_fasta).exists():
        Path(output_gff).parent.mkdir(parents=True, exist_ok=True)
        with open(output_gff, 'w') as f:
            f.write("##gff-version 3\n")
        print(f"Genome FASTA not found: {genome_fasta}")
        return
    
    sequences = parse_fasta(genome_fasta)
    
    # Load MGE bounds for mobility classification
    mge_bounds = load_mge_bounds(mges_tsv, genome_id)
    
    # Find all matches
    all_hits = []
    for pattern in patterns:
        hits = find_pattern_matches(pattern, sequences, mge_bounds)
        all_hits.extend(hits)
    
    # Write GFF
    write_gff(all_hits, genome_id, output_gff)
    
    # Count mobile vs genomic
    mobile_count = sum(1 for h in all_hits if h['mge_id'])
    genomic_count = len(all_hits) - mobile_count
    
    print(f"Found {len(all_hits)} RM recognition site hits for {genome_id}")
    print(f"  - Mobile (within MGE): {mobile_count}")
    print(f"  - Genomic (outside MGE): {genomic_count}")


def main():
    parser = argparse.ArgumentParser(
        description='Map RM recognition sites to genome and output GFF'
    )
    parser.add_argument(
        '--genome-id', required=True,
        help='Target genome ID'
    )
    parser.add_argument(
        '--patterns-tsv', required=True,
        help='Path to patterns_rm.tsv with metadata'
    )
    parser.add_argument(
        '--genome-fasta', required=True,
        help='Path to genome FASTA'
    )
    parser.add_argument(
        '--mges-tsv', required=True,
        help='Path to all_mges.tsv for mobility classification'
    )
    parser.add_argument(
        '--output-gff', required=True,
        help='Output GFF path'
    )
    
    args = parser.parse_args()
    
    map_rm_targets(
        args.genome_id,
        args.patterns_tsv,
        args.genome_fasta,
        args.mges_tsv,
        args.output_gff
    )


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Aggregate and deduplicate defense system patterns (CRISPR spacers, RM recognition sites).

Produces deduplicated FASTA files and metadata TSVs for target mapping.
"""

import argparse
import csv
import hashlib
from collections import defaultdict
from pathlib import Path


def sequence_hash(seq: str) -> str:
    """Generate a short hash for a sequence (first 8 chars of MD5)."""
    return hashlib.md5(seq.upper().encode()).hexdigest()[:8]


def aggregate_spacer_patterns(spacers_tsv: str, output_fasta: str, output_tsv: str) -> None:
    """
    Deduplicate CRISPR spacers from all_crispr_spacers.tsv by sequence.
    
    Args:
        spacers_tsv: Path to all_crispr_spacers.tsv from main workflow
        output_fasta: Path to output patterns_spacers.fa
        output_tsv: Path to output patterns_spacers.tsv with metadata
    """
    # Group spacers by sequence (uppercase for deduplication)
    seq_to_sources = defaultdict(lambda: {
        'genome_ids': set(),
        'array_ids': set(),
        'cas_types': set(),
        'lengths': set()
    })
    
    with open(spacers_tsv, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            seq = row['sequence'].upper().strip()
            if not seq:
                continue
            
            seq_to_sources[seq]['genome_ids'].add(row['genome_id'])
            seq_to_sources[seq]['array_ids'].add(row['Parent'])  # Parent = array ID
            seq_to_sources[seq]['lengths'].add(len(seq))
    
    # Write outputs
    Path(output_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
    
    patterns = []
    for seq, sources in sorted(seq_to_sources.items(), key=lambda x: x[0]):
        pattern_id = f"spacer_{sequence_hash(seq)}"
        patterns.append({
            'pattern_id': pattern_id,
            'sequence': seq,
            'length': len(seq),
            'source_genome_ids': ','.join(sorted(sources['genome_ids'])),
            'source_array_ids': ','.join(sorted(sources['array_ids'])),
            'source_cas_types': ''  # Would need join with cas_systems
        })
    
    # Write FASTA
    with open(output_fasta, 'w') as f:
        for p in patterns:
            f.write(f">{p['pattern_id']}\n{p['sequence']}\n")
    
    # Write TSV
    with open(output_tsv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'pattern_id', 'sequence', 'length', 
            'source_genome_ids', 'source_array_ids', 'source_cas_types'
        ], delimiter='\t')
        writer.writeheader()
        writer.writerows(patterns)
    
    print(f"Aggregated {len(seq_to_sources)} unique spacer patterns from spacers TSV")


def aggregate_spacer_patterns_with_cas(
    spacers_tsv: str, 
    cas_systems_tsv: str,
    output_fasta: str, 
    output_tsv: str
) -> None:
    """
    Deduplicate CRISPR spacers with cas type information.
    
    Includes source spacer coordinates for self-hit depletion during mapping.
    
    Args:
        spacers_tsv: Path to all_crispr_spacers.tsv from main workflow
        cas_systems_tsv: Path to all_cas_systems.tsv for cas type lookup
        output_fasta: Path to output patterns_spacers.fa
        output_tsv: Path to output patterns_spacers.tsv with metadata
    """
    # Build array_id -> cas_subtype mapping
    array_to_cas = {}
    if cas_systems_tsv and Path(cas_systems_tsv).exists():
        with open(cas_systems_tsv, 'r', newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                array_id = row.get('array_ID', '')
                cas_subtype = row.get('cas_subtype', '')
                if array_id and array_id != 'NULL' and cas_subtype:
                    array_to_cas[array_id] = cas_subtype
    
    # Group spacers by sequence
    # Also track source coordinates for self-hit depletion
    seq_to_sources = defaultdict(lambda: {
        'genome_ids': set(),
        'array_ids': set(),
        'cas_types': set(),
        'lengths': set(),
        'source_coords': []  # List of (genome_id, contig, start, end) tuples
    })
    
    with open(spacers_tsv, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            seq = row['sequence'].upper().strip()
            if not seq:
                continue
            
            genome_id = row['genome_id']
            array_id = row['Parent']  # Parent = array ID
            contig = row['contig']
            start = int(row['start'])
            end = int(row['end'])
            
            seq_to_sources[seq]['genome_ids'].add(genome_id)
            seq_to_sources[seq]['array_ids'].add(array_id)
            seq_to_sources[seq]['lengths'].add(len(seq))
            
            # Store source coordinates for self-hit depletion
            seq_to_sources[seq]['source_coords'].append(
                f"{genome_id}:{contig}:{start}-{end}"
            )
            
            # Add cas type if available
            if array_id in array_to_cas:
                seq_to_sources[seq]['cas_types'].add(array_to_cas[array_id])
    
    # Write outputs
    Path(output_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
    
    patterns = []
    for seq, sources in sorted(seq_to_sources.items(), key=lambda x: x[0]):
        pattern_id = f"spacer_{sequence_hash(seq)}"
        patterns.append({
            'pattern_id': pattern_id,
            'sequence': seq,
            'length': len(seq),
            'source_genome_ids': ','.join(sorted(sources['genome_ids'])),
            'source_array_ids': ','.join(sorted(sources['array_ids'])),
            'source_cas_types': ','.join(sorted(sources['cas_types'])),
            'source_coords': ','.join(sources['source_coords'])
        })
    
    # Write FASTA
    with open(output_fasta, 'w') as f:
        for p in patterns:
            f.write(f">{p['pattern_id']}\n{p['sequence']}\n")
    
    # Write TSV
    with open(output_tsv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'pattern_id', 'sequence', 'length', 
            'source_genome_ids', 'source_array_ids', 'source_cas_types',
            'source_coords'
        ], delimiter='\t')
        writer.writeheader()
        writer.writerows(patterns)
    
    print(f"Aggregated {len(seq_to_sources)} unique spacer patterns")


def aggregate_rm_patterns(rm_systems_tsv: str, output_fasta: str, output_tsv: str) -> None:
    """
    Deduplicate RM recognition sites from all_rm_systems.tsv.
    
    Only includes systems where recseq_note == 'annotated'.
    
    Args:
        rm_systems_tsv: Path to all_rm_systems.tsv from main workflow
        output_fasta: Path to output patterns_rm.fa
        output_tsv: Path to output patterns_rm.tsv with metadata
    """
    # Group RM sites by sequence
    seq_to_sources = defaultdict(lambda: {
        'genome_ids': set(),
        'system_ids': set(),
        'rm_types': set()
    })
    
    with open(rm_systems_tsv, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # Only include annotated recognition sites
            if row.get('recseq_note', '') != 'annotated':
                continue
            
            recseq = row.get('REBASE_recseq', '').upper().strip()
            if not recseq or recseq == 'NULL':
                continue
            
            seq_to_sources[recseq]['genome_ids'].add(row['genome_id'])
            seq_to_sources[recseq]['system_ids'].add(row['ID'])
            seq_to_sources[recseq]['rm_types'].add(row.get('type', ''))
    
    # Write outputs
    Path(output_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
    
    patterns = []
    for seq, sources in sorted(seq_to_sources.items(), key=lambda x: x[0]):
        pattern_id = f"rm_{sequence_hash(seq)}"
        patterns.append({
            'pattern_id': pattern_id,
            'sequence': seq,
            'length': len(seq),
            'rm_type': ','.join(sorted(sources['rm_types'])),
            'source_genome_ids': ','.join(sorted(sources['genome_ids'])),
            'source_system_ids': ','.join(sorted(sources['system_ids']))
        })
    
    # Write FASTA
    with open(output_fasta, 'w') as f:
        for p in patterns:
            f.write(f">{p['pattern_id']}\n{p['sequence']}\n")
    
    # Write TSV
    with open(output_tsv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'pattern_id', 'sequence', 'length', 'rm_type',
            'source_genome_ids', 'source_system_ids'
        ], delimiter='\t')
        writer.writeheader()
        writer.writerows(patterns)
    
    print(f"Aggregated {len(seq_to_sources)} unique RM recognition patterns")


def main():
    parser = argparse.ArgumentParser(
        description='Aggregate and deduplicate defense system patterns'
    )
    subparsers = parser.add_subparsers(dest='command', required=True)
    
    # Spacer aggregation subcommand
    spacer_parser = subparsers.add_parser('spacers', help='Aggregate CRISPR spacers')
    spacer_parser.add_argument(
        '--spacers-tsv', required=True,
        help='Path to all_crispr_spacers.tsv'
    )
    spacer_parser.add_argument(
        '--cas-systems-tsv', required=False,
        help='Path to all_cas_systems.tsv (optional, for cas type annotation)'
    )
    spacer_parser.add_argument(
        '--output-fasta', required=True,
        help='Output path for patterns_spacers.fa'
    )
    spacer_parser.add_argument(
        '--output-tsv', required=True,
        help='Output path for patterns_spacers.tsv'
    )
    
    # RM aggregation subcommand
    rm_parser = subparsers.add_parser('rm', help='Aggregate RM recognition sites')
    rm_parser.add_argument(
        '--rm-systems-tsv', required=True,
        help='Path to all_rm_systems.tsv'
    )
    rm_parser.add_argument(
        '--output-fasta', required=True,
        help='Output path for patterns_rm.fa'
    )
    rm_parser.add_argument(
        '--output-tsv', required=True,
        help='Output path for patterns_rm.tsv'
    )
    
    args = parser.parse_args()
    
    if args.command == 'spacers':
        if args.cas_systems_tsv:
            aggregate_spacer_patterns_with_cas(
                args.spacers_tsv,
                args.cas_systems_tsv,
                args.output_fasta,
                args.output_tsv
            )
        else:
            aggregate_spacer_patterns(
                args.spacers_tsv,
                args.output_fasta,
                args.output_tsv
            )
    elif args.command == 'rm':
        aggregate_rm_patterns(
            args.rm_systems_tsv,
            args.output_fasta,
            args.output_tsv
        )


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Extract MGE (Mobile Genetic Element) sequences from genome FASTAs.

For each genome, reads the MGE annotations and extracts subsequences,
creating a FASTA file of MGE sequences for target mapping.
"""

import argparse
import csv
from pathlib import Path


def parse_fasta(fasta_path: str) -> dict:
    """
    Parse a FASTA file into a dictionary of sequences.
    
    Args:
        fasta_path: Path to the FASTA file
        
    Returns:
        dict: contig_id -> sequence (uppercase)
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
                # Get ID (first word after >)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq).upper()
    
    return sequences


def extract_mge_sequences(
    mges_tsv: str,
    genome_fasta: str,
    genome_id: str,
    output_fasta: str
) -> None:
    """
    Extract MGE sequences from a genome FASTA based on MGE annotations.
    
    Args:
        mges_tsv: Path to all_mges.tsv summary file
        genome_fasta: Path to the genome FASTA file
        genome_id: The genome ID to filter MGEs for
        output_fasta: Path to output MGE sequences FASTA
    """
    # Parse genome FASTA
    sequences = parse_fasta(genome_fasta)
    
    # Filter MGEs for this genome and extract sequences
    mge_records = []
    
    with open(mges_tsv, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['genome_id'] != genome_id:
                continue
            
            mge_id = row['ID']
            contig = row['contig']
            start = int(row['start'])
            end = int(row['end'])
            
            # Find matching contig in FASTA
            contig_seq = sequences.get(contig)
            if contig_seq is None:
                print(f"Warning: Contig {contig} not found in genome FASTA for MGE {mge_id}")
                continue
            
            # Extract subsequence (1-based coordinates, inclusive)
            # Convert to 0-based for Python slicing
            mge_seq = contig_seq[start - 1:end]
            
            if not mge_seq:
                print(f"Warning: Empty sequence for MGE {mge_id}")
                continue
            
            # Store MGE info in header for later coordinate conversion
            header = f"{mge_id} contig={contig} start={start} end={end}"
            mge_records.append((header, mge_seq))
    
    # Write output
    Path(output_fasta).parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_fasta, 'w') as f:
        for header, seq in mge_records:
            f.write(f">{header}\n")
            # Write sequence in 80-char lines
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')
    
    print(f"Extracted {len(mge_records)} MGE sequences for {genome_id}")


def main():
    parser = argparse.ArgumentParser(
        description='Extract MGE sequences from genome FASTA'
    )
    parser.add_argument(
        '--mges-tsv', required=True,
        help='Path to all_mges.tsv summary file'
    )
    parser.add_argument(
        '--genome-fasta', required=True,
        help='Path to genome FASTA file'
    )
    parser.add_argument(
        '--genome-id', required=True,
        help='Genome ID to filter MGEs for'
    )
    parser.add_argument(
        '--output-fasta', required=True,
        help='Output path for MGE sequences FASTA'
    )
    
    args = parser.parse_args()
    
    extract_mge_sequences(
        args.mges_tsv,
        args.genome_fasta,
        args.genome_id,
        args.output_fasta
    )


if __name__ == '__main__':
    main()

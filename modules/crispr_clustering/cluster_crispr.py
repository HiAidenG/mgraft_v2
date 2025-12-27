#!/usr/bin/env python3
"""
CRISPR Clustering Module

This module clusters CRISPR arrays based on:
1. Repeat sequences: 100% identity, 100% coverage
2. Spacer sequences: 90% identity, 90% coverage

Uses MMseqs2 for clustering with proper handling of reverse complement sequences.

Usage:
    python cluster_crispr.py --arrays <arrays.tsv> --spacers <spacers.tsv> \
        --output-dir <output_dir> [--threads <n>]
"""

import argparse
import subprocess
import sys
import os
import tempfile
from pathlib import Path
from collections import defaultdict
import shutil

import pandas as pd


def check_mmseqs_installed():
    """Check if MMseqs2 is installed and accessible."""
    try:
        result = subprocess.run(
            ["mmseqs", "version"],
            capture_output=True,
            text=True,
            check=True
        )
        print(f"Using MMseqs2 version: {result.stdout.strip()}")
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("ERROR: MMseqs2 is not installed or not in PATH", file=sys.stderr)
        print("Install with: conda install -c bioconda mmseqs2", file=sys.stderr)
        return False


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'a': 't', 't': 'a', 'g': 'c', 'c': 'g',
                  'N': 'N', 'n': 'n'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def canonicalize_sequence(seq: str) -> tuple[str, bool]:
    """
    Return the canonical form of a sequence (lexicographically smaller of seq and revcomp).
    
    Returns:
        tuple: (canonical_sequence, is_revcomp) where is_revcomp indicates if the
               canonical form is the reverse complement
    """
    seq_upper = seq.upper()
    revcomp = reverse_complement(seq_upper)
    if revcomp < seq_upper:
        return revcomp, True
    return seq_upper, False


def write_fasta(sequences: dict, output_file: Path) -> None:
    """
    Write sequences to FASTA file.
    
    Args:
        sequences: Dict mapping sequence ID to sequence string
        output_file: Path to output FASTA file
    """
    with open(output_file, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n{seq}\n")


def prepare_repeat_sequences(arrays_df: pd.DataFrame) -> tuple[dict, dict]:
    """
    Prepare repeat sequences for clustering.
    
    Canonicalizes sequences (uses lexicographically smaller of seq/revcomp)
    to ensure sequences that are reverse complements of each other cluster together.
    
    Args:
        arrays_df: DataFrame with CRISPR array data including repeat_sequence column
    
    Returns:
        tuple: (sequences dict, metadata dict with original orientation info)
    """
    sequences = {}
    metadata = {}
    
    for _, row in arrays_df.iterrows():
        array_id = row['ID']
        repeat_seq = row['repeat_sequence']
        
        if pd.isna(repeat_seq) or not repeat_seq:
            continue
        
        # Canonicalize the sequence
        canonical_seq, is_revcomp = canonicalize_sequence(repeat_seq)
        
        sequences[array_id] = canonical_seq
        metadata[array_id] = {
            'original_sequence': repeat_seq,
            'canonical_sequence': canonical_seq,
            'used_revcomp': is_revcomp,
            'genome_id': row['genome_id'],
            'contig': row['contig'],
            'start': row['start'],
            'end': row['end'],
            'strand': row['strand'],
            'n_spacers': row.get('n_spacers', 0),
            'array_type': row.get('array_type', 'unknown')
        }
    
    return sequences, metadata


def prepare_spacer_sequences(spacers_df: pd.DataFrame) -> tuple[dict, dict]:
    """
    Prepare spacer sequences for clustering.
    
    Canonicalizes sequences to handle reverse complement matching.
    
    Args:
        spacers_df: DataFrame with CRISPR spacer data
    
    Returns:
        tuple: (sequences dict, metadata dict)
    """
    sequences = {}
    metadata = {}
    
    for _, row in spacers_df.iterrows():
        spacer_id = row['ID']
        spacer_seq = row['sequence']
        
        if pd.isna(spacer_seq) or not spacer_seq:
            continue
        
        # Canonicalize the sequence
        canonical_seq, is_revcomp = canonicalize_sequence(spacer_seq)
        
        sequences[spacer_id] = canonical_seq
        metadata[spacer_id] = {
            'original_sequence': spacer_seq,
            'canonical_sequence': canonical_seq,
            'used_revcomp': is_revcomp,
            'parent_array': row['Parent'],
            'genome_id': row['genome_id'],
            'contig': row['contig'],
            'spacer_idx': row.get('spacer_idx', 0),
            'length': row.get('length', len(spacer_seq))
        }
    
    return sequences, metadata


def run_mmseqs_cluster(input_fasta: Path, output_prefix: Path, 
                       min_seq_id: float, coverage: float,
                       threads: int, tmp_dir: Path) -> Path:
    """
    Run MMseqs2 easy-linclust for sequence clustering.
    
    Uses linclust instead of cluster for better memory efficiency.
    Linclust is suitable for high-identity clustering (>90%).
    
    Args:
        input_fasta: Path to input FASTA file
        output_prefix: Prefix for output files
        min_seq_id: Minimum sequence identity (0-1)
        coverage: Minimum coverage (0-1)
        threads: Number of threads
        tmp_dir: Temporary directory for MMseqs2
    
    Returns:
        Path to the cluster TSV file
    """
    cluster_tsv = Path(f"{output_prefix}_cluster.tsv")
    
    cmd = [
        "mmseqs", "easy-linclust",
        str(input_fasta),
        str(output_prefix),
        str(tmp_dir),
        "--min-seq-id", str(min_seq_id),
        "-c", str(coverage),
        "--cov-mode", "0",  # Coverage of both query and target
        "--threads", str(threads),
        "--split-memory-limit", "2G",  # Limit memory usage
        "-v", "2"  # Verbosity
    ]
    
    print(f"Running MMseqs2 linclust: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    return cluster_tsv


def run_mmseqs_linclust(input_fasta: Path, output_prefix: Path,
                        min_seq_id: float, coverage: float,
                        threads: int, tmp_dir: Path) -> Path:
    """
    Run MMseqs2 easy-linclust for fast linear-time clustering.
    Better for large datasets with exact identity requirements.
    
    Args:
        input_fasta: Path to input FASTA file
        output_prefix: Prefix for output files
        min_seq_id: Minimum sequence identity (0-1)
        coverage: Minimum coverage (0-1)
        threads: Number of threads
        tmp_dir: Temporary directory for MMseqs2
    
    Returns:
        Path to the cluster TSV file
    """
    cluster_tsv = Path(f"{output_prefix}_cluster.tsv")
    
    cmd = [
        "mmseqs", "easy-linclust",
        str(input_fasta),
        str(output_prefix),
        str(tmp_dir),
        "--min-seq-id", str(min_seq_id),
        "-c", str(coverage),
        "--cov-mode", "0",  # Coverage of both query and target
        "--threads", str(threads),
        "--split-memory-limit", "2G",  # Limit memory usage
        "-v", "2"
    ]
    
    print(f"Running MMseqs2 linclust: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    return cluster_tsv


def parse_cluster_tsv(cluster_tsv: Path) -> dict:
    """
    Parse MMseqs2 cluster TSV output.
    
    The TSV has two columns: representative_id, member_id
    
    Args:
        cluster_tsv: Path to cluster TSV file
    
    Returns:
        Dict mapping member ID to cluster representative ID
    """
    member_to_rep = {}
    
    with open(cluster_tsv, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                rep_id, member_id = parts[0], parts[1]
                member_to_rep[member_id] = rep_id
    
    return member_to_rep


def assign_cluster_ids(member_to_rep: dict, prefix: str = "CL") -> tuple[dict, dict]:
    """
    Assign numeric cluster IDs based on representative sequences.
    
    Args:
        member_to_rep: Dict mapping member ID to representative ID
        prefix: Prefix for cluster IDs
    
    Returns:
        tuple: (member_to_cluster_id, cluster_info)
            - member_to_cluster_id: Dict mapping member ID to cluster ID
            - cluster_info: Dict mapping cluster ID to {representative, members}
    """
    # Group members by representative
    rep_to_members = defaultdict(list)
    for member, rep in member_to_rep.items():
        rep_to_members[rep].append(member)
    
    # Sort representatives by cluster size (descending) for consistent numbering
    sorted_reps = sorted(rep_to_members.keys(), 
                         key=lambda x: len(rep_to_members[x]), 
                         reverse=True)
    
    member_to_cluster = {}
    cluster_info = {}
    
    for i, rep in enumerate(sorted_reps, 1):
        cluster_id = f"{prefix}{i:06d}"
        members = rep_to_members[rep]
        
        for member in members:
            member_to_cluster[member] = cluster_id
        
        cluster_info[cluster_id] = {
            'representative': rep,
            'members': members,
            'size': len(members)
        }
    
    return member_to_cluster, cluster_info


def cluster_repeats(arrays_df: pd.DataFrame, output_dir: Path, 
                    threads: int = 4) -> tuple[pd.DataFrame, dict]:
    """
    Cluster CRISPR arrays by repeat sequence (100% identity, 100% coverage).
    
    Args:
        arrays_df: DataFrame with CRISPR array data
        output_dir: Output directory
        threads: Number of threads
    
    Returns:
        tuple: (updated arrays DataFrame with cluster IDs, cluster info dict)
    """
    print("\n=== Clustering CRISPR Repeats ===")
    print(f"Input: {len(arrays_df)} arrays")
    
    # Prepare sequences
    sequences, metadata = prepare_repeat_sequences(arrays_df)
    print(f"Valid repeat sequences: {len(sequences)}")
    
    if not sequences:
        print("WARNING: No valid repeat sequences found")
        arrays_df['repeat_cluster_id'] = None
        return arrays_df, {}
    
    # Write FASTA
    repeat_fasta = output_dir / "repeat_sequences.fasta"
    write_fasta(sequences, repeat_fasta)
    
    # Create temp directory for MMseqs2
    tmp_dir = output_dir / "mmseqs_tmp_repeats"
    tmp_dir.mkdir(exist_ok=True)
    
    # Run clustering with 100% identity, 100% coverage
    # Use linclust for exact matching (faster for high identity)
    output_prefix = output_dir / "repeat_clusters"
    
    try:
        cluster_tsv = run_mmseqs_linclust(
            input_fasta=repeat_fasta,
            output_prefix=output_prefix,
            min_seq_id=1.0,  # 100% identity
            coverage=1.0,     # 100% coverage
            threads=threads,
            tmp_dir=tmp_dir
        )
        
        # Parse results
        member_to_rep = parse_cluster_tsv(cluster_tsv)
        member_to_cluster, cluster_info = assign_cluster_ids(
            member_to_rep, prefix="RPT"
        )
        
        # Add cluster IDs to dataframe
        arrays_df['repeat_cluster_id'] = arrays_df['ID'].map(member_to_cluster)
        
        # Add metadata to cluster info
        for cluster_id, info in cluster_info.items():
            rep_id = info['representative']
            if rep_id in metadata:
                info['representative_sequence'] = metadata[rep_id]['canonical_sequence']
                info['sequence_length'] = len(metadata[rep_id]['canonical_sequence'])
        
        print(f"Repeat clusters: {len(cluster_info)}")
        
    finally:
        # Cleanup temp directory
        if tmp_dir.exists():
            shutil.rmtree(tmp_dir)
    
    return arrays_df, cluster_info


def cluster_spacers(spacers_df: pd.DataFrame, output_dir: Path,
                    threads: int = 4) -> tuple[pd.DataFrame, dict]:
    """
    Cluster CRISPR spacers (90% identity, 90% coverage).
    
    Args:
        spacers_df: DataFrame with CRISPR spacer data
        output_dir: Output directory
        threads: Number of threads
    
    Returns:
        tuple: (updated spacers DataFrame with cluster IDs, cluster info dict)
    """
    print("\n=== Clustering CRISPR Spacers ===")
    print(f"Input: {len(spacers_df)} spacers")
    
    # Prepare sequences
    sequences, metadata = prepare_spacer_sequences(spacers_df)
    print(f"Valid spacer sequences: {len(sequences)}")
    
    if not sequences:
        print("WARNING: No valid spacer sequences found")
        spacers_df['spacer_cluster_id'] = None
        return spacers_df, {}
    
    # Write FASTA
    spacer_fasta = output_dir / "spacer_sequences.fasta"
    write_fasta(sequences, spacer_fasta)
    
    # Create temp directory for MMseqs2
    tmp_dir = output_dir / "mmseqs_tmp_spacers"
    tmp_dir.mkdir(exist_ok=True)
    
    # Run clustering with 90% identity, 90% coverage
    output_prefix = output_dir / "spacer_clusters"
    
    try:
        cluster_tsv = run_mmseqs_cluster(
            input_fasta=spacer_fasta,
            output_prefix=output_prefix,
            min_seq_id=0.9,  # 90% identity
            coverage=0.9,     # 90% coverage
            threads=threads,
            tmp_dir=tmp_dir
        )
        
        # Parse results
        member_to_rep = parse_cluster_tsv(cluster_tsv)
        member_to_cluster, cluster_info = assign_cluster_ids(
            member_to_rep, prefix="SPC"
        )
        
        # Add cluster IDs to dataframe
        spacers_df['spacer_cluster_id'] = spacers_df['ID'].map(member_to_cluster)
        
        # Add metadata to cluster info
        for cluster_id, info in cluster_info.items():
            rep_id = info['representative']
            if rep_id in metadata:
                info['representative_sequence'] = metadata[rep_id]['canonical_sequence']
                info['sequence_length'] = len(metadata[rep_id]['canonical_sequence'])
                info['parent_arrays'] = list(set(
                    metadata[m]['parent_array'] for m in info['members'] 
                    if m in metadata
                ))
        
        print(f"Spacer clusters: {len(cluster_info)}")
        
    finally:
        # Cleanup temp directory
        if tmp_dir.exists():
            shutil.rmtree(tmp_dir)
    
    return spacers_df, cluster_info


def generate_summary_stats(arrays_df: pd.DataFrame, spacers_df: pd.DataFrame,
                           repeat_clusters: dict, spacer_clusters: dict,
                           output_dir: Path) -> None:
    """
    Generate summary statistics for clustering results.
    
    Args:
        arrays_df: DataFrame with array cluster assignments
        spacers_df: DataFrame with spacer cluster assignments
        repeat_clusters: Dict with repeat cluster info
        spacer_clusters: Dict with spacer cluster info
        output_dir: Output directory
    """
    summary_file = output_dir / "clustering_summary.txt"
    
    with open(summary_file, 'w') as f:
        f.write("CRISPR Clustering Summary\n")
        f.write("=" * 50 + "\n\n")
        
        # Repeat clustering stats
        f.write("Repeat Sequence Clustering (100% ID, 100% coverage)\n")
        f.write("-" * 50 + "\n")
        f.write(f"Total arrays: {len(arrays_df)}\n")
        f.write(f"Arrays with valid repeats: {arrays_df['repeat_cluster_id'].notna().sum()}\n")
        f.write(f"Number of repeat clusters: {len(repeat_clusters)}\n")
        
        if repeat_clusters:
            sizes = [c['size'] for c in repeat_clusters.values()]
            f.write(f"Largest cluster: {max(sizes)} arrays\n")
            f.write(f"Singleton clusters: {sum(1 for s in sizes if s == 1)}\n")
            f.write(f"Average cluster size: {sum(sizes)/len(sizes):.2f}\n")
        
        f.write("\n")
        
        # Spacer clustering stats
        f.write("Spacer Sequence Clustering (90% ID, 90% coverage)\n")
        f.write("-" * 50 + "\n")
        f.write(f"Total spacers: {len(spacers_df)}\n")
        f.write(f"Spacers with valid sequences: {spacers_df['spacer_cluster_id'].notna().sum()}\n")
        f.write(f"Number of spacer clusters: {len(spacer_clusters)}\n")
        
        if spacer_clusters:
            sizes = [c['size'] for c in spacer_clusters.values()]
            f.write(f"Largest cluster: {max(sizes)} spacers\n")
            f.write(f"Singleton clusters: {sum(1 for s in sizes if s == 1)}\n")
            f.write(f"Average cluster size: {sum(sizes)/len(sizes):.2f}\n")
            
            # Count spacers shared across genomes
            multi_genome_clusters = sum(
                1 for c in spacer_clusters.values() 
                if len(set(spacers_df[spacers_df['ID'].isin(c['members'])]['genome_id'])) > 1
            )
            f.write(f"Clusters spanning multiple genomes: {multi_genome_clusters}\n")
    
    print(f"\nSummary written to: {summary_file}")


def save_cluster_info(repeat_clusters: dict, spacer_clusters: dict,
                      output_dir: Path) -> None:
    """
    Save detailed cluster information to TSV files.
    
    Args:
        repeat_clusters: Dict with repeat cluster info
        spacer_clusters: Dict with spacer cluster info
        output_dir: Output directory
    """
    # Repeat clusters
    repeat_rows = []
    for cluster_id, info in repeat_clusters.items():
        repeat_rows.append({
            'cluster_id': cluster_id,
            'representative': info['representative'],
            'size': info['size'],
            'sequence': info.get('representative_sequence', ''),
            'sequence_length': info.get('sequence_length', 0),
            'members': ';'.join(info['members'])
        })
    
    repeat_df = pd.DataFrame(repeat_rows)
    repeat_df.to_csv(output_dir / "repeat_cluster_info.tsv", sep='\t', index=False)
    
    # Spacer clusters
    spacer_rows = []
    for cluster_id, info in spacer_clusters.items():
        spacer_rows.append({
            'cluster_id': cluster_id,
            'representative': info['representative'],
            'size': info['size'],
            'sequence': info.get('representative_sequence', ''),
            'sequence_length': info.get('sequence_length', 0),
            'n_parent_arrays': len(info.get('parent_arrays', [])),
            'parent_arrays': ';'.join(info.get('parent_arrays', [])),
            'members': ';'.join(info['members'])
        })
    
    spacer_df = pd.DataFrame(spacer_rows)
    spacer_df.to_csv(output_dir / "spacer_cluster_info.tsv", sep='\t', index=False)
    
    print(f"Cluster info saved to {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Cluster CRISPR arrays by repeat and spacer sequences using MMseqs2"
    )
    parser.add_argument(
        "--arrays", "-a",
        type=Path,
        required=True,
        help="Path to all_crispr_arrays.tsv"
    )
    parser.add_argument(
        "--spacers", "-s",
        type=Path,
        required=True,
        help="Path to all_crispr_spacers.tsv"
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=Path,
        required=True,
        help="Output directory for clustering results"
    )
    parser.add_argument(
        "--threads", "-t",
        type=int,
        default=4,
        help="Number of threads (default: 4)"
    )
    parser.add_argument(
        "--repeat-identity",
        type=float,
        default=1.0,
        help="Minimum identity for repeat clustering (default: 1.0)"
    )
    parser.add_argument(
        "--repeat-coverage",
        type=float,
        default=1.0,
        help="Minimum coverage for repeat clustering (default: 1.0)"
    )
    parser.add_argument(
        "--spacer-identity",
        type=float,
        default=0.9,
        help="Minimum identity for spacer clustering (default: 0.9)"
    )
    parser.add_argument(
        "--spacer-coverage",
        type=float,
        default=0.9,
        help="Minimum coverage for spacer clustering (default: 0.9)"
    )
    
    args = parser.parse_args()
    
    # Check MMseqs2 installation
    if not check_mmseqs_installed():
        sys.exit(1)
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    print(f"Loading arrays from: {args.arrays}")
    arrays_df = pd.read_csv(args.arrays, sep='\t')
    
    print(f"Loading spacers from: {args.spacers}")
    spacers_df = pd.read_csv(args.spacers, sep='\t')
    
    # Cluster repeats
    arrays_df, repeat_clusters = cluster_repeats(
        arrays_df, args.output_dir, args.threads
    )
    
    # Cluster spacers
    spacers_df, spacer_clusters = cluster_spacers(
        spacers_df, args.output_dir, args.threads
    )
    
    # Save updated dataframes
    arrays_output = args.output_dir / "crispr_arrays_clustered.tsv"
    arrays_df.to_csv(arrays_output, sep='\t', index=False)
    print(f"Clustered arrays saved to: {arrays_output}")
    
    spacers_output = args.output_dir / "crispr_spacers_clustered.tsv"
    spacers_df.to_csv(spacers_output, sep='\t', index=False)
    print(f"Clustered spacers saved to: {spacers_output}")
    
    # Save cluster info
    save_cluster_info(repeat_clusters, spacer_clusters, args.output_dir)
    
    # Generate summary
    generate_summary_stats(
        arrays_df, spacers_df, 
        repeat_clusters, spacer_clusters,
        args.output_dir
    )
    
    print("\n=== Clustering Complete ===")


if __name__ == "__main__":
    main()

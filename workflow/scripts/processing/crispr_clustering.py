#!/usr/bin/env python3
"""
CRISPR Spacer/Repeat Clustering

Clusters spacers and repeats across all genomes using MMseqs2.

Config parameters:
    - crispr.spacer_cluster_identity: Minimum identity for spacer clustering (default: 0.95)
    - crispr.repeat_cluster_identity: Minimum identity for repeat clustering (default: 0.90)
"""
import json
import pandas as pd
import subprocess
from pathlib import Path


def safe_read(path):
    """Read TSV file safely, returning empty DataFrame if missing."""
    if path and Path(path).exists() and Path(path).stat().st_size > 0:
        return pd.read_csv(path, sep='\t')
    return pd.DataFrame()


def cluster_sequences(sequences, seq_ids, work_dir, min_seq_id=0.9):
    """Cluster nucleotide sequences using MMseqs2.
    
    Args:
        sequences: List of nucleotide sequences
        seq_ids: List of sequence identifiers
        work_dir: Working directory for MMseqs2 temp files
        min_seq_id: Minimum sequence identity threshold
        
    Returns:
        Dict mapping sequence_id -> cluster_rep_id
    """
    if not sequences:
        return {}
    
    # Filter out invalid sequences
    valid_pairs = [(sid, seq) for sid, seq in zip(seq_ids, sequences) 
                   if seq and seq != 'NULL' and pd.notna(seq) and len(str(seq)) >= 10]
    
    if not valid_pairs:
        print("No valid sequences to cluster")
        return {}
    
    work_path = Path(work_dir)
    work_path.mkdir(parents=True, exist_ok=True)
    
    # Write FASTA
    fasta_path = work_path / 'seqs.fa'
    with open(fasta_path, 'w') as f:
        for seq_id, seq in valid_pairs:
            f.write(f'>{seq_id}\n{seq}\n')
    
    print(f"Clustering {len(valid_pairs)} sequences at {min_seq_id:.0%} identity...")
    
    # Run MMseqs2 with nucleotide mode (--dbtype 2)
    db_path = work_path / 'db'
    clust_path = work_path / 'clust'
    tsv_path = work_path / 'clusters.tsv'
    
    try:
        # Create nucleotide database
        subprocess.run(
            ['mmseqs', 'createdb', str(fasta_path), str(db_path), '--dbtype', '2'],
            check=True, capture_output=True, text=True
        )
        
        # Cluster
        subprocess.run(
            ['mmseqs', 'cluster', str(db_path), str(clust_path), str(work_path),
             '--min-seq-id', str(min_seq_id), '-c', '0.8', '--cov-mode', '1'],
            check=True, capture_output=True, text=True
        )
        
        # Convert to TSV
        subprocess.run(
            ['mmseqs', 'createtsv', str(db_path), str(db_path), str(clust_path), str(tsv_path)],
            check=True, capture_output=True, text=True
        )
        
        # Parse results
        clusters = {}
        with open(tsv_path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    rep, member = parts[0], parts[1]
                    clusters[member] = rep
        
        n_clusters = len(set(clusters.values()))
        print(f"Created {n_clusters} clusters from {len(clusters)} sequences")
        return clusters
        
    except subprocess.CalledProcessError as e:
        print(f"Clustering error: {e}")
        print(f"stderr: {e.stderr}")
        return {}
    except Exception as e:
        print(f"Unexpected error: {e}")
        return {}


def main(snakemake):
    """Main entry point for Snakemake."""
    # Read config
    crispr_config = snakemake.config.get("crispr", {})
    spacer_identity = crispr_config.get("spacer_cluster_identity", 0.95)
    repeat_identity = crispr_config.get("repeat_cluster_identity", 0.90)
    
    work_dir = snakemake.params.work_dir
    
    # Load data
    arrays_df = safe_read(snakemake.input.arrays)
    spacers_df = safe_read(snakemake.input.spacers)
    
    print(f"Loaded: {len(arrays_df)} arrays, {len(spacers_df)} spacers")
    print(f"Config: spacer_identity={spacer_identity}, repeat_identity={repeat_identity}")
    
    stats = {
        'config': {
            'spacer_cluster_identity': spacer_identity,
            'repeat_cluster_identity': repeat_identity,
            'coverage': 0.8,
            'cov_mode': 1,
        },
        'input_arrays': int(len(arrays_df)),
        'input_spacers': int(len(spacers_df)),
    }
    
    # Cluster repeats
    if not arrays_df.empty and 'repeat_sequence' in arrays_df.columns:
        repeat_clusters = cluster_sequences(
            arrays_df['repeat_sequence'].tolist(),
            arrays_df['ID'].tolist(),
            f"{work_dir}/repeats",
            min_seq_id=repeat_identity
        )
        arrays_df['repeat_cluster_id'] = arrays_df['ID'].map(repeat_clusters)
        
        if repeat_clusters:
            cluster_counts = pd.Series(repeat_clusters).value_counts()
            stats['repeat_clusters'] = int(len(cluster_counts))
            stats['repeat_singletons'] = int((cluster_counts == 1).sum())
            stats['repeat_largest_cluster'] = int(cluster_counts.max())
        else:
            stats['repeat_clusters'] = 0
            stats['repeat_singletons'] = 0
            stats['repeat_largest_cluster'] = 0
        print(f"Repeat clusters: {stats['repeat_clusters']} ({stats['repeat_singletons']} singletons, largest: {stats['repeat_largest_cluster']})")
    else:
        stats['repeat_clusters'] = 0
        stats['repeat_singletons'] = 0
        stats['repeat_largest_cluster'] = 0
    
    # Cluster spacers
    if not spacers_df.empty and 'sequence' in spacers_df.columns:
        spacer_clusters = cluster_sequences(
            spacers_df['sequence'].tolist(),
            spacers_df['ID'].tolist(),
            f"{work_dir}/spacers",
            min_seq_id=spacer_identity
        )
        spacers_df['spacer_cluster_id'] = spacers_df['ID'].map(spacer_clusters)
        
        if spacer_clusters:
            cluster_counts = pd.Series(spacer_clusters).value_counts()
            stats['spacer_clusters'] = int(len(cluster_counts))
            stats['spacer_singletons'] = int((cluster_counts == 1).sum())
            stats['spacer_largest_cluster'] = int(cluster_counts.max())
        else:
            stats['spacer_clusters'] = 0
            stats['spacer_singletons'] = 0
            stats['spacer_largest_cluster'] = 0
        print(f"Spacer clusters: {stats['spacer_clusters']} ({stats['spacer_singletons']} singletons, largest: {stats['spacer_largest_cluster']})")
    else:
        stats['spacer_clusters'] = 0
        stats['spacer_singletons'] = 0
        stats['spacer_largest_cluster'] = 0
    
    # Save outputs
    Path(snakemake.output.arrays).parent.mkdir(parents=True, exist_ok=True)
    arrays_df.to_csv(snakemake.output.arrays, sep='\t', index=False, na_rep='NULL')
    spacers_df.to_csv(snakemake.output.spacers, sep='\t', index=False, na_rep='NULL')
    
    with open(snakemake.output.stats, 'w') as f:
        json.dump(stats, f, indent=2)
    
    print(f"\nSaved to: {snakemake.output.arrays}, {snakemake.output.spacers}, {snakemake.output.stats}")


if __name__ == '__main__':
    main(snakemake)  # noqa: F821

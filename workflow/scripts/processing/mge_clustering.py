#!/usr/bin/env python3
"""
MGE Sequence Deduplication

Clusters MGEs based on whole nucleotide sequence similarity using vclust.

Algorithm:
1. Load MGE sequences from FASTA files (extracted by parse_mges.py)
2. Cluster by ANI threshold using vclust
3. Assign cluster IDs back to MGE TSV

Config params:
- mge.ani: default 0.95 (ANI threshold)
- mge.qcov: default 0.8 (query coverage)
- mge.rcov: default 0.8 (reference coverage)
"""

import json
import pandas as pd
import subprocess
from pathlib import Path
from Bio import SeqIO


def cluster_nucleotides_vclust(fasta_paths, work_dir, ani=0.95, qcov=0.8, rcov=0.8):
    """
    Cluster nucleotide sequences using vclust.
    
    Args:
        fasta_paths: List of paths to FASTA files
        work_dir: working directory
        ani: ANI threshold (--ani flag)
        qcov: query coverage threshold
        rcov: reference coverage threshold
    
    Returns:
        tuple: ({seq_id: cluster_rep_id}, stats_dict)
    """
    # Load all sequences
    seqs_dict = {}
    for fasta_path in fasta_paths:
        if Path(fasta_path).exists() and Path(fasta_path).stat().st_size > 0:
            for rec in SeqIO.parse(fasta_path, 'fasta'):
                seqs_dict[rec.id] = str(rec.seq)
    
    stats = {
        'config': {
            'ani_threshold': ani,
            'qcov_threshold': qcov,
            'rcov_threshold': rcov,
        },
        'n_clusters': 0,
    }
    
    if not seqs_dict:
        return {}, stats
    
    print(f"Loaded {len(seqs_dict)} MGE sequences for clustering")
    print(f"Clustering with ANI >= {ani}, qcov >= {qcov}, rcov >= {rcov}")
    
    work_path = Path(work_dir)
    work_path.mkdir(parents=True, exist_ok=True)
    
    # Write combined FASTA
    fasta_path = work_path / 'all_mges.fa'
    with open(fasta_path, 'w') as f:
        for seq_id, seq in seqs_dict.items():
            f.write(f'>{seq_id}\n{seq}\n')
    
    output_path = work_path / 'clusters.tsv'
    align_path = work_path / 'align.tsv'
    
    try:
        # Run vclust prefilter
        prefilter_path = work_path / 'prefilter.txt'
        result = subprocess.run(
            ['vclust', 'prefilter', '-i', str(fasta_path), '-o', str(prefilter_path)],
            capture_output=True, text=True, check=True
        )
        print(f"vclust prefilter: {result.stdout}")
        
        # Run vclust align to get ANI values
        result = subprocess.run(
            ['vclust', 'align', '-i', str(fasta_path), '-o', str(align_path),
             '--filter', str(prefilter_path)],
            capture_output=True, text=True, check=True
        )
        print(f"vclust align: {result.stdout}")
        
        # Run vclust cluster
        result = subprocess.run(
            ['vclust', 'cluster', '-i', str(align_path), '-o', str(output_path),
             '--ani', str(ani), '--qcov', str(qcov), '--tcov', str(rcov)],
            capture_output=True, text=True, check=True
        )
        print(f"vclust cluster: {result.stdout}")
        
        # Parse alignment results for stats
        ani_values = []
        qcov_values = []
        if align_path.exists():
            with open(align_path) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        try:
                            ani_values.append(float(parts[2]))
                            qcov_values.append(float(parts[3]))
                        except (ValueError, IndexError):
                            pass
        
        if ani_values:
            stats['ani_min'] = min(ani_values)
            stats['ani_max'] = max(ani_values)
            stats['ani_mean'] = sum(ani_values) / len(ani_values)
            stats['qcov_min'] = min(qcov_values)
            stats['qcov_max'] = max(qcov_values)
            stats['qcov_mean'] = sum(qcov_values) / len(qcov_values)
            stats['n_comparisons'] = len(ani_values)
        
        # Parse vclust cluster output
        clusters = {}
        if output_path.exists():
            with open(output_path) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        rep, member = parts[0], parts[1]
                        clusters[member] = rep
        
        stats['n_clusters'] = len(set(clusters.values()))
        
        print(f"\n=== Clustering Results ===")
        print(f"Total sequences: {len(seqs_dict)}")
        print(f"Clusters: {stats['n_clusters']}")
        if ani_values:
            print(f"ANI range: {stats['ani_min']:.3f} - {stats['ani_max']:.3f} (mean: {stats['ani_mean']:.3f})")
            print(f"Qcov range: {stats['qcov_min']:.3f} - {stats['qcov_max']:.3f} (mean: {stats['qcov_mean']:.3f})")
        
        return clusters, stats
        
    except subprocess.CalledProcessError as e:
        print(f"vclust error: {e}")
        print(f"stderr: {e.stderr}")
        return {}, stats
    except Exception as e:
        print(f"Unexpected error: {e}")
        return {}, stats


def main(snakemake):
    """Main entry point."""
    # Access inputs/outputs
    mges_tsvs = snakemake.input.mges
    mges_fastas = snakemake.input.mges_fasta
    output_mges = snakemake.output.mges
    output_stats = snakemake.output.stats
    work_dir = snakemake.params.work_dir
    
    # Clustering params from config
    mge_config = snakemake.config.get('mge', {})
    ANI = mge_config.get('ani', mge_config.get('identity', 0.95))
    QCOV = mge_config.get('qcov', 0.8)
    RCOV = mge_config.get('rcov', 0.8)
    
    print(f"MGE Deduplication with ANI={ANI}, qcov={QCOV}, rcov={RCOV}")
    
    # Load all MGE TSV data
    dfs = []
    for tsv_path in mges_tsvs:
        if Path(tsv_path).exists() and Path(tsv_path).stat().st_size > 0:
            df = pd.read_csv(tsv_path, sep='\t')
            dfs.append(df)
    
    all_mges = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
    print(f"Loaded {len(all_mges)} MGEs from {len(dfs)} samples")
    
    stats = {'input_mges': len(all_mges)}
    
    if all_mges.empty:
        print("No MGEs to cluster")
        all_mges['cluster_id'] = None
        all_mges['n_cluster_members'] = None
        stats['clusters'] = 0
        stats['clustered_mges'] = 0
    else:
        # Run vclust clustering on whole sequences
        clusters, cluster_stats = cluster_nucleotides_vclust(
            mges_fastas, work_dir, ani=ANI, qcov=QCOV, rcov=RCOV
        )
        stats.update(cluster_stats)
        
        # Map cluster IDs back to DataFrame
        all_mges['cluster_id'] = all_mges['ID'].map(clusters)
        
        # For sequences that weren't clustered (singletons), use their own ID
        mask = all_mges['cluster_id'].isna()
        all_mges.loc[mask, 'cluster_id'] = all_mges.loc[mask, 'ID']
        
        # Count cluster members
        cluster_counts = all_mges['cluster_id'].value_counts()
        all_mges['n_cluster_members'] = all_mges['cluster_id'].map(cluster_counts.to_dict())
        
        # Stats
        stats['n_singletons'] = int((cluster_counts == 1).sum())
        stats['largest_cluster_size'] = int(cluster_counts.max())
    
    print(f"\nCreated {stats.get('n_clusters', 0)} clusters ({stats.get('n_singletons', 0)} singletons, largest: {stats.get('largest_cluster_size', 0)})")
    
    # Save outputs
    Path(output_mges).parent.mkdir(parents=True, exist_ok=True)
    all_mges.to_csv(output_mges, sep='\t', index=False, na_rep='NULL')
    
    with open(output_stats, 'w') as f:
        json.dump(stats, f, indent=2)
    
    print(f"\nSaved {len(all_mges)} MGEs to: {output_mges}")
    print(f"Stats: {json.dumps(stats, indent=2)}")


if __name__ == '__main__':
    main(snakemake)  # noqa: F821

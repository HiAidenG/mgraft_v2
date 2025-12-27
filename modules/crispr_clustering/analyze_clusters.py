#!/usr/bin/env python3
"""
Analyze CRISPR Clustering Results

This script generates detailed analysis and visualizations of CRISPR 
clustering results, including:
- Shared spacer analysis across genomes
- Repeat type diversity
- Array similarity networks based on shared spacer clusters
"""

import argparse
from pathlib import Path
from collections import defaultdict
from itertools import combinations

import pandas as pd
import numpy as np


def load_clustering_results(output_dir: Path) -> tuple:
    """Load clustering result files."""
    arrays_df = pd.read_csv(output_dir / "crispr_arrays_clustered.tsv", sep='\t')
    spacers_df = pd.read_csv(output_dir / "crispr_spacers_clustered.tsv", sep='\t')
    repeat_info = pd.read_csv(output_dir / "repeat_cluster_info.tsv", sep='\t')
    spacer_info = pd.read_csv(output_dir / "spacer_cluster_info.tsv", sep='\t')
    
    return arrays_df, spacers_df, repeat_info, spacer_info


def analyze_shared_spacers(spacers_df: pd.DataFrame) -> pd.DataFrame:
    """
    Analyze spacers shared across different genomes.
    
    Returns DataFrame with spacer clusters found in multiple genomes.
    """
    # Group spacers by cluster
    cluster_genomes = spacers_df.groupby('spacer_cluster_id')['genome_id'].apply(
        lambda x: list(set(x))
    ).reset_index()
    cluster_genomes.columns = ['spacer_cluster_id', 'genomes']
    
    # Filter to clusters in multiple genomes
    cluster_genomes['n_genomes'] = cluster_genomes['genomes'].apply(len)
    shared = cluster_genomes[cluster_genomes['n_genomes'] > 1].copy()
    
    # Add spacer count
    spacer_counts = spacers_df.groupby('spacer_cluster_id').size()
    shared['n_spacers'] = shared['spacer_cluster_id'].map(spacer_counts)
    
    # Add representative sequence
    rep_seqs = spacers_df.groupby('spacer_cluster_id')['sequence'].first()
    shared['representative_sequence'] = shared['spacer_cluster_id'].map(rep_seqs)
    
    shared = shared.sort_values('n_genomes', ascending=False)
    
    return shared


def analyze_repeat_diversity(arrays_df: pd.DataFrame, 
                             repeat_info: pd.DataFrame) -> pd.DataFrame:
    """
    Analyze diversity of repeat sequences across genomes.
    
    Returns DataFrame with repeat clusters and their genome distribution.
    """
    # Map arrays to genomes
    cluster_genomes = arrays_df.groupby('repeat_cluster_id')['genome_id'].apply(
        lambda x: list(set(x))
    ).reset_index()
    cluster_genomes.columns = ['repeat_cluster_id', 'genomes']
    cluster_genomes['n_genomes'] = cluster_genomes['genomes'].apply(len)
    
    # Merge with repeat info
    result = repeat_info.merge(cluster_genomes, left_on='cluster_id', 
                                right_on='repeat_cluster_id', how='left')
    
    result = result.sort_values('n_genomes', ascending=False)
    
    return result


def compute_array_similarity_matrix(spacers_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute similarity matrix between arrays based on shared spacer clusters.
    
    Similarity = Jaccard index of spacer cluster sets
    """
    # Get spacer clusters for each array
    array_clusters = spacers_df.groupby('Parent')['spacer_cluster_id'].apply(set).to_dict()
    
    arrays = list(array_clusters.keys())
    n = len(arrays)
    
    # Compute pairwise Jaccard similarity
    similarity_data = []
    
    for i, j in combinations(range(n), 2):
        arr1, arr2 = arrays[i], arrays[j]
        set1, set2 = array_clusters[arr1], array_clusters[arr2]
        
        intersection = len(set1 & set2)
        union = len(set1 | set2)
        
        if union > 0:
            jaccard = intersection / union
            if jaccard > 0:  # Only store non-zero similarities
                similarity_data.append({
                    'array1': arr1,
                    'array2': arr2,
                    'shared_clusters': intersection,
                    'union_clusters': union,
                    'jaccard_similarity': jaccard
                })
    
    return pd.DataFrame(similarity_data)


def find_related_arrays(spacers_df: pd.DataFrame, 
                        min_shared: int = 3) -> pd.DataFrame:
    """
    Find arrays that share multiple spacer clusters.
    
    Args:
        spacers_df: Spacer dataframe with cluster assignments
        min_shared: Minimum number of shared spacer clusters
    
    Returns:
        DataFrame with array pairs and their shared spacers
    """
    # Get spacer clusters for each array
    array_clusters = spacers_df.groupby('Parent')['spacer_cluster_id'].apply(set).to_dict()
    
    # Also get genome info for each array
    array_genomes = spacers_df.groupby('Parent')['genome_id'].first().to_dict()
    
    related_pairs = []
    
    for (arr1, clusters1), (arr2, clusters2) in combinations(array_clusters.items(), 2):
        shared = clusters1 & clusters2
        if len(shared) >= min_shared:
            related_pairs.append({
                'array1': arr1,
                'genome1': array_genomes.get(arr1),
                'array2': arr2,
                'genome2': array_genomes.get(arr2),
                'shared_spacer_clusters': len(shared),
                'array1_total_clusters': len(clusters1),
                'array2_total_clusters': len(clusters2),
                'shared_cluster_ids': ';'.join(sorted(shared))
            })
    
    result = pd.DataFrame(related_pairs)
    if not result.empty:
        result = result.sort_values('shared_spacer_clusters', ascending=False)
    
    return result


def generate_network_file(similarity_df: pd.DataFrame, 
                          output_file: Path,
                          min_similarity: float = 0.1) -> None:
    """
    Generate network file for visualization (e.g., Cytoscape).
    
    Exports edge list in TSV format.
    """
    filtered = similarity_df[similarity_df['jaccard_similarity'] >= min_similarity]
    
    filtered.to_csv(output_file, sep='\t', index=False)
    print(f"Network file written: {output_file}")
    print(f"  Edges: {len(filtered)}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze CRISPR clustering results"
    )
    parser.add_argument(
        "--input-dir", "-i",
        type=Path,
        required=True,
        help="Directory with clustering results"
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=Path,
        help="Output directory (default: input-dir/analysis)"
    )
    parser.add_argument(
        "--min-shared-spacers",
        type=int,
        default=3,
        help="Minimum shared spacer clusters for related array detection"
    )
    parser.add_argument(
        "--min-similarity",
        type=float,
        default=0.1,
        help="Minimum Jaccard similarity for network edges"
    )
    
    args = parser.parse_args()
    
    if args.output_dir is None:
        args.output_dir = args.input_dir / "analysis"
    
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    print("Loading clustering results...")
    arrays_df, spacers_df, repeat_info, spacer_info = load_clustering_results(
        args.input_dir
    )
    
    # Shared spacer analysis
    print("\nAnalyzing shared spacers...")
    shared_spacers = analyze_shared_spacers(spacers_df)
    shared_spacers.to_csv(
        args.output_dir / "shared_spacers_across_genomes.tsv", 
        sep='\t', index=False
    )
    print(f"  Found {len(shared_spacers)} spacer clusters shared across genomes")
    
    # Repeat diversity analysis
    print("\nAnalyzing repeat diversity...")
    repeat_diversity = analyze_repeat_diversity(arrays_df, repeat_info)
    repeat_diversity.to_csv(
        args.output_dir / "repeat_diversity.tsv",
        sep='\t', index=False
    )
    print(f"  Found {len(repeat_diversity)} unique repeat clusters")
    
    # Related arrays analysis
    print("\nFinding related arrays...")
    related_arrays = find_related_arrays(spacers_df, args.min_shared_spacers)
    related_arrays.to_csv(
        args.output_dir / "related_arrays.tsv",
        sep='\t', index=False
    )
    print(f"  Found {len(related_arrays)} pairs of related arrays")
    
    # Similarity network
    print("\nComputing array similarity network...")
    similarity_matrix = compute_array_similarity_matrix(spacers_df)
    
    generate_network_file(
        similarity_matrix,
        args.output_dir / "array_similarity_network.tsv",
        args.min_similarity
    )
    
    # Summary report
    print("\n=== Analysis Summary ===")
    print(f"Total arrays: {len(arrays_df)}")
    print(f"Total spacers: {len(spacers_df)}")
    print(f"Unique repeat clusters: {arrays_df['repeat_cluster_id'].nunique()}")
    print(f"Unique spacer clusters: {spacers_df['spacer_cluster_id'].nunique()}")
    print(f"Spacer clusters in multiple genomes: {len(shared_spacers)}")
    print(f"Related array pairs (>={args.min_shared_spacers} shared): {len(related_arrays)}")
    
    print(f"\nResults saved to: {args.output_dir}")


if __name__ == "__main__":
    main()

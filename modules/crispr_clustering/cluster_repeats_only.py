#!/usr/bin/env python3
"""
Cluster CRISPR repeat sequences only.

This is a standalone script for Snakemake rule integration.
"""

import sys
from pathlib import Path

import pandas as pd

# Add parent module to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from modules.crispr_clustering.cluster_crispr import (
    cluster_repeats,
    check_mmseqs_installed,
    save_cluster_info
)


def main():
    # Get snakemake parameters
    arrays_file = Path(snakemake.input.arrays)
    output_dir = Path(snakemake.params.output_dir)
    threads = snakemake.threads
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Check MMseqs2
    if not check_mmseqs_installed():
        sys.exit(1)
    
    # Load and cluster
    arrays_df = pd.read_csv(arrays_file, sep='\t')
    arrays_df, repeat_clusters = cluster_repeats(arrays_df, output_dir, threads)
    
    # Save results
    save_cluster_info(repeat_clusters, {}, output_dir)
    
    print("Repeat clustering complete")


if __name__ == "__main__":
    main()

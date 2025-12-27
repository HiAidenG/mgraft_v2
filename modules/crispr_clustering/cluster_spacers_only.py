#!/usr/bin/env python3
"""
Cluster CRISPR spacer sequences only.

This is a standalone script for Snakemake rule integration.
"""

import sys
from pathlib import Path

import pandas as pd

# Add parent module to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from modules.crispr_clustering.cluster_crispr import (
    cluster_spacers,
    check_mmseqs_installed,
    save_cluster_info
)


def main():
    # Get snakemake parameters
    spacers_file = Path(snakemake.input.spacers)
    output_dir = Path(snakemake.params.output_dir)
    threads = snakemake.threads
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Check MMseqs2
    if not check_mmseqs_installed():
        sys.exit(1)
    
    # Load and cluster
    spacers_df = pd.read_csv(spacers_file, sep='\t')
    spacers_df, spacer_clusters = cluster_spacers(spacers_df, output_dir, threads)
    
    # Save results
    save_cluster_info({}, spacer_clusters, output_dir)
    
    print("Spacer clustering complete")


if __name__ == "__main__":
    main()

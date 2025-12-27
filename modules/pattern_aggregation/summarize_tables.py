#!/usr/bin/env python3
"""
Aggregate per-genome feature tables into a single summary table.

This script concatenates all per-genome TSV files of a given type into
a single combined summary table.
"""
import pandas as pd
from pathlib import Path


def main(snakemake):
    """Main entry point for Snakemake script."""
    input_files = snakemake.input
    output_file = snakemake.output[0]
    
    dfs = []
    for f in input_files:
        try:
            df = pd.read_csv(f, sep='\t')
            if not df.empty:
                dfs.append(df)
        except Exception as e:
            # Log but continue - some files might be empty or malformed
            print(f"Warning: Could not read {f}: {e}")
    
    if dfs:
        result = pd.concat(dfs, ignore_index=True)
        # Write with NULL for missing values
        result.to_csv(output_file, sep='\t', index=False, na_rep='NULL')
    else:
        # Create empty file with no content
        pd.DataFrame().to_csv(output_file, sep='\t', index=False)


if __name__ == '__main__':
    main(snakemake)

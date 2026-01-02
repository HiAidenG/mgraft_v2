#!/usr/bin/env python3
"""
CRISPR Array and Spacer Filtering

Filters CRISPR arrays and spacers based on quality metrics.

Purpose:
    Remove low-quality or artifact CRISPR arrays/spacers before downstream analysis.

Assumptions:
    - Spacers < 20 bp are likely annotation artifacts (incomplete spacers)
    - Spacers > 55 bp are likely merged repeats or misannotations
    - LCC (Linguistic Complexity) < 0.8 indicates repetitive/low-complexity sequence
"""
import json
import pandas as pd
from pathlib import Path
from Bio.SeqUtils.lcc import lcc_simp


def calculate_lcc(sequence: str) -> float:
    """Calculate linguistic complexity (LCC) for a sequence.
    
    LCC measures sequence complexity using k-mer diversity.
    Ranges from 0 (homopolymer) to ~2 (random sequence).
    """
    if not sequence or len(sequence) < 4:
        return 0.0
    try:
        return lcc_simp(sequence)
    except Exception:
        return 0.0


def safe_read(path):
    """Read TSV file safely, returning empty DataFrame if missing."""
    if path and Path(path).exists() and Path(path).stat().st_size > 0:
        return pd.read_csv(path, sep='\t')
    return pd.DataFrame()


def filter_arrays_by_lcc(arrays_df, threshold):
    """Filter arrays by LCC threshold.
    
    Returns:
        Tuple of (filtered_df, set of filtered array IDs)
    """
    if arrays_df.empty:
        return arrays_df, set()
    
    arrays_df = arrays_df.copy()
    arrays_df['lcc_score'] = arrays_df['repeat_sequence'].apply(
        lambda x: calculate_lcc(str(x)) if pd.notna(x) and x != 'NULL' else 0.0
    )
    
    pass_mask = arrays_df['lcc_score'] >= threshold
    filtered_ids = set(arrays_df.loc[~pass_mask, 'ID'].dropna())
    
    return arrays_df[pass_mask].copy(), filtered_ids


def filter_spacers_by_length(spacers_df, min_len, max_len):
    """Filter spacers by length.
    
    Returns:
        Tuple of (filtered_df, count removed)
    """
    if spacers_df.empty:
        return spacers_df, 0
    
    spacers_df = spacers_df.copy()
    
    if 'length' not in spacers_df.columns:
        spacers_df['length'] = spacers_df['sequence'].apply(
            lambda x: len(str(x)) if pd.notna(x) and x != 'NULL' else 0
        )
    
    mask = (spacers_df['length'] >= min_len) & (spacers_df['length'] <= max_len)
    removed = (~mask).sum()
    
    return spacers_df[mask].copy(), removed


def filter_spacers_by_lcc(spacers_df, threshold):
    """Filter spacers by LCC threshold.
    
    Returns:
        Tuple of (filtered_df, count removed)
    """
    if spacers_df.empty:
        return spacers_df, 0
    
    spacers_df = spacers_df.copy()
    spacers_df['lcc_score'] = spacers_df['sequence'].apply(
        lambda x: calculate_lcc(str(x)) if pd.notna(x) and x != 'NULL' else 0.0
    )
    
    mask = spacers_df['lcc_score'] >= threshold
    removed = (~mask).sum()
    
    return spacers_df[mask].copy(), removed


def main(snakemake):
    """Main entry point for Snakemake."""
    # Read config with defaults
    crispr_config = snakemake.config.get("crispr", {})
    spacer_min_length = crispr_config.get("spacer_min_length", 20)
    spacer_max_length = crispr_config.get("spacer_max_length", 55)
    lcc_threshold = crispr_config.get("lcc_threshold", 0.8)
    
    # Load input data
    arrays_df = safe_read(snakemake.input.arrays)
    spacers_df = safe_read(snakemake.input.spacers)
    
    print(f"Loaded {len(arrays_df)} arrays, {len(spacers_df)} spacers")
    print(f"Config: min_len={spacer_min_length}, max_len={spacer_max_length}, lcc={lcc_threshold}")
    
    stats = {
        'input_arrays': int(len(arrays_df)),
        'input_spacers': int(len(spacers_df)),
    }
    
    # Filter arrays by LCC
    arrays_df, filtered_array_ids = filter_arrays_by_lcc(arrays_df, lcc_threshold)
    stats['arrays_filtered_low_complexity'] = len(filtered_array_ids)
    print(f"Arrays: {stats['input_arrays']} -> {len(arrays_df)} (filtered {len(filtered_array_ids)} low-complexity)")
    
    # Remove spacers from filtered arrays
    if filtered_array_ids and not spacers_df.empty and 'Parent' in spacers_df.columns:
        spacer_keep_mask = ~spacers_df['Parent'].isin(filtered_array_ids)
        stats['spacers_removed_with_arrays'] = int((~spacer_keep_mask).sum())
        spacers_df = spacers_df[spacer_keep_mask].copy()
        print(f"Removed {stats['spacers_removed_with_arrays']} spacers from filtered arrays")
    else:
        stats['spacers_removed_with_arrays'] = 0
    
    # Filter spacers by length
    spacers_df, removed_length = filter_spacers_by_length(spacers_df, spacer_min_length, spacer_max_length)
    stats['spacers_filtered_length'] = int(removed_length)
    print(f"Filtered {removed_length} spacers by length")
    
    # Filter spacers by complexity
    spacers_df, removed_lcc = filter_spacers_by_lcc(spacers_df, lcc_threshold)
    stats['spacers_filtered_low_complexity'] = int(removed_lcc)
    print(f"Filtered {removed_lcc} low-complexity spacers")
    
    # Summary
    stats['output_arrays'] = int(len(arrays_df))
    stats['output_spacers'] = int(len(spacers_df))
    
    print(f"\n=== Filtering Summary ===")
    print(f"Arrays: {stats['input_arrays']} -> {stats['output_arrays']}")
    print(f"Spacers: {stats['input_spacers']} -> {stats['output_spacers']}")
    
    # Save outputs
    Path(snakemake.output.arrays).parent.mkdir(parents=True, exist_ok=True)
    arrays_df.to_csv(snakemake.output.arrays, sep='\t', index=False, na_rep='NULL')
    spacers_df.to_csv(snakemake.output.spacers, sep='\t', index=False, na_rep='NULL')
    
    with open(snakemake.output.stats, 'w') as f:
        json.dump(stats, f, indent=2)
    
    print(f"\nSaved outputs to {snakemake.output.arrays}, {snakemake.output.spacers}, {snakemake.output.stats}")


if __name__ == '__main__':
    main(snakemake)  # noqa: F821

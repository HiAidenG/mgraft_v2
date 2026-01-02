#!/usr/bin/env python3
"""
Array-Cas System Linking

Links CRISPR arrays to Cas systems based on genomic proximity.

Algorithm (see .agent-instructions/array_cas_linking.md):
1. Find candidate arrays within max_distance of any Cas gene
2. Select best (closest) array per system
3. Classify arrays as canonical/putative/orphan based on system wholeness
"""
import pandas as pd
from pathlib import Path
from collections import defaultdict


def safe_read(path):
    """Read TSV file safely, returning empty DataFrame if missing."""
    if path and Path(path).exists() and Path(path).stat().st_size > 0:
        return pd.read_csv(path, sep='\t')
    return pd.DataFrame()


def find_array_cas_candidates(arrays_df, cas_genes_df, max_distance):
    """Find candidate array-Cas links within distance threshold.
    
    Returns:
        List of candidate dicts with array_idx, array_id, sys_id, distance, quality
    """
    candidates = []
    
    if arrays_df.empty or cas_genes_df.empty:
        return candidates
    
    for idx, array_row in arrays_df.iterrows():
        same_contig = cas_genes_df[cas_genes_df['contig'] == array_row['contig']]
        if same_contig.empty:
            continue
        
        array_mid = (array_row['start'] + array_row['end']) / 2
        same_contig = same_contig.copy()
        same_contig['midpoint'] = (same_contig['start'] + same_contig['end']) / 2
        same_contig['dist'] = (same_contig['midpoint'] - array_mid).abs()
        
        nearest_gene = same_contig.loc[same_contig['dist'].idxmin()]
        min_dist = nearest_gene['dist']
        
        if min_dist <= max_distance:
            sys_id = nearest_gene['Parent']
            candidates.append({
                'array_idx': idx,
                'array_id': array_row['ID'],
                'sys_id': sys_id,
                'distance': min_dist,
                'quality': array_row.get('quality_score', 0) or 0
            })
    
    return candidates


def select_best_array_per_system(candidates):
    """Select best (closest) array for each Cas system.
    
    Tie-breaker: higher quality score wins.
    
    Returns:
        Dict of array_idx -> candidate dict for winning arrays
    """
    sys_candidates = defaultdict(list)
    for c in candidates:
        sys_candidates[c['sys_id']].append(c)
    
    winning_arrays = {}
    for sys_id, cands in sys_candidates.items():
        # Sort by distance (asc), then quality (desc)
        cands.sort(key=lambda x: (x['distance'], -x['quality']))
        winning_arrays[cands[0]['array_idx']] = cands[0]
    
    return winning_arrays


def classify_and_update_arrays(arrays_df, cas_systems_df, winning_arrays):
    """Update arrays with linking info and classification.
    
    Classifications:
        - canonical: linked to complete Cas system (wholeness=1.0)
        - putative: linked to incomplete Cas system
        - orphan: no Cas system within threshold
    """
    system_linked_arrays = defaultdict(list)
    
    for array_idx, link in winning_arrays.items():
        sys_id = link['sys_id']
        
        # Get system wholeness
        is_complete = False
        if not cas_systems_df.empty:
            sys_row = cas_systems_df[cas_systems_df['ID'] == sys_id]
            if not sys_row.empty:
                wholeness = float(sys_row.iloc[0].get('sys_wholeness', 0))
                is_complete = wholeness == 1.0
        
        arrays_df.at[array_idx, 'linked_cas_id'] = sys_id
        arrays_df.at[array_idx, 'distance_to_cas'] = int(link['distance'])
        arrays_df.at[array_idx, 'array_class'] = 'canonical' if is_complete else 'putative'
        system_linked_arrays[sys_id].append(link['array_id'])
    
    return arrays_df, system_linked_arrays


def update_cas_systems(cas_systems_df, system_linked_arrays):
    """Update Cas systems with linked array info."""
    if cas_systems_df.empty:
        return cas_systems_df
    
    for sys_id, arr_ids in system_linked_arrays.items():
        mask = cas_systems_df['ID'] == sys_id
        cas_systems_df.loc[mask, 'linked_array_ids'] = ','.join(arr_ids)
        cas_systems_df.loc[mask, 'n_linked_arrays'] = len(arr_ids)
    
    return cas_systems_df


def main(snakemake):
    """Main entry point for Snakemake."""
    # Read config
    crispr_config = snakemake.config.get("crispr", {})
    max_distance = crispr_config.get("max_distance_to_cas", 10000)
    
    # Load data
    arrays_df = safe_read(snakemake.input.arrays)
    spacers_df = safe_read(snakemake.input.spacers)
    cas_systems_df = safe_read(snakemake.input.cas_systems)
    cas_genes_df = safe_read(snakemake.input.cas_genes)
    
    print(f"Loaded: {len(arrays_df)} arrays, {len(cas_systems_df)} Cas systems")
    print(f"Config: max_distance={max_distance}")
    
    # Initialize columns
    if not arrays_df.empty:
        arrays_df['linked_cas_id'] = None
        arrays_df['array_class'] = 'orphan'
        arrays_df['distance_to_cas'] = None
    
    if not cas_systems_df.empty:
        cas_systems_df['linked_array_ids'] = None
        cas_systems_df['n_linked_arrays'] = 0
    
    # Find and select candidates
    candidates = find_array_cas_candidates(arrays_df, cas_genes_df, max_distance)
    winning_arrays = select_best_array_per_system(candidates)
    
    # Update dataframes
    arrays_df, system_linked_arrays = classify_and_update_arrays(
        arrays_df, cas_systems_df, winning_arrays
    )
    cas_systems_df = update_cas_systems(cas_systems_df, system_linked_arrays)
    
    print(f"Linked {len(winning_arrays)} arrays to Cas systems")
    
    # Save outputs
    Path(snakemake.output.arrays).parent.mkdir(parents=True, exist_ok=True)
    arrays_df.to_csv(snakemake.output.arrays, sep='\t', index=False, na_rep='NULL')
    cas_systems_df.to_csv(snakemake.output.cas_systems, sep='\t', index=False, na_rep='NULL')
    spacers_df.to_csv(snakemake.output.spacers, sep='\t', index=False, na_rep='NULL')
    
    print(f"Saved to: {snakemake.output.arrays}, {snakemake.output.cas_systems}, {snakemake.output.spacers}")


if __name__ == '__main__':
    main(snakemake)  # noqa: F821

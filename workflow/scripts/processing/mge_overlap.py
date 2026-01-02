#!/usr/bin/env python3
"""
MGE Overlap Assignment

Assigns `mge_id` to genes/systems/arrays that overlap with MGEs and 
calculates `mobility_score` for systems.

Logic:
- Genes: Assigned mge_id if gene coordinates overlap with any MGE
- Systems: 
  - mge_id = MGE containing majority of genes (or NULL if none)
  - mobility_score = fraction of system genes on any single MGE
- Arrays: Assigned mge_id if array coordinates overlap with any MGE
- Spacers: Inherit mge_id from parent array
"""

import pandas as pd
from pathlib import Path


def safe_read(path):
    """Read TSV file safely, returning empty DataFrame if missing."""
    if path and Path(path).exists() and Path(path).stat().st_size > 0:
        return pd.read_csv(path, sep='\t')
    return pd.DataFrame()


def feature_overlaps_mge(feature_row, mges_df):
    """Check if a feature overlaps with any MGE.
    
    Returns mge_id if feature overlaps, None otherwise.
    Overlap = feature midpoint within MGE boundaries on same contig.
    """
    if mges_df.empty:
        return None
    
    # Handle NULL/NaN coordinates
    if pd.isna(feature_row.get('start')) or pd.isna(feature_row.get('end')):
        return None
    
    feature_contig = feature_row.get('contig')
    feature_start = int(feature_row['start'])
    feature_end = int(feature_row['end'])
    
    # Filter MGEs on same contig
    contig_mges = mges_df[mges_df['contig'] == feature_contig]
    
    for _, mge in contig_mges.iterrows():
        mge_start = int(mge['start'])
        mge_end = int(mge['end'])
        
        # Check overlap (feature midpoint within MGE)
        feature_mid = (feature_start + feature_end) / 2
        if mge_start <= feature_mid <= mge_end:
            return mge['ID']
    
    return None


def assign_mge_to_features(features_df, mges_df):
    """Assign mge_id to each feature based on overlap."""
    if features_df.empty:
        features_df['mge_id'] = None
        return features_df
    
    features_df = features_df.copy()
    features_df['mge_id'] = features_df.apply(
        lambda row: feature_overlaps_mge(row, mges_df), axis=1
    )
    return features_df


def assign_mge_to_spacers(spacers_df, arrays_df):
    """Assign mge_id to spacers from their parent array."""
    if spacers_df.empty:
        spacers_df['mge_id'] = None
        return spacers_df
    
    spacers_df = spacers_df.copy()
    
    # Build array ID to mge_id map
    if arrays_df.empty or 'ID' not in arrays_df.columns:
        spacers_df['mge_id'] = None
        return spacers_df
    
    array_mge_map = dict(zip(arrays_df['ID'], arrays_df.get('mge_id', [None] * len(arrays_df))))
    
    # Assign mge_id from parent array
    spacers_df['mge_id'] = spacers_df['Parent'].map(array_mge_map)
    
    return spacers_df


def calculate_system_mobility(systems_df, genes_df):
    """Calculate mge_id and mobility_score for each system.
    
    - mge_id: MGE containing majority of genes (or NULL if none)
    - mobility_score: fraction of genes on the most common MGE
    """
    if systems_df.empty:
        systems_df['mge_id'] = None
        systems_df['mobility_score'] = None
        return systems_df
    
    systems_df = systems_df.copy()
    systems_df['mge_id'] = None
    systems_df['mobility_score'] = 0.0
    
    if genes_df.empty or 'Parent' not in genes_df.columns:
        return systems_df
    
    for idx, sys_row in systems_df.iterrows():
        sys_id = sys_row.get('ID')
        if pd.isna(sys_id):
            continue
        
        # Get genes for this system
        sys_genes = genes_df[genes_df['Parent'] == sys_id]
        if sys_genes.empty:
            continue
        
        n_genes = len(sys_genes)
        
        # Count genes per MGE
        mge_counts = sys_genes['mge_id'].dropna().value_counts()
        
        if len(mge_counts) == 0:
            # No genes on MGEs
            systems_df.at[idx, 'mge_id'] = None
            systems_df.at[idx, 'mobility_score'] = 0.0
        else:
            # Most common MGE
            top_mge = mge_counts.idxmax()
            top_count = mge_counts.max()
            
            systems_df.at[idx, 'mge_id'] = top_mge
            systems_df.at[idx, 'mobility_score'] = round(top_count / n_genes, 3)
    
    return systems_df


def main(snakemake):
    """Main entry point for Snakemake."""
    genome_id = snakemake.params.genome_id
    
    # Load MGEs
    mges_df = safe_read(snakemake.input.mges)
    if not mges_df.empty and 'genome_id' in mges_df.columns:
        mges_df = mges_df[mges_df['genome_id'] == genome_id].copy()
    
    print(f"Loaded {len(mges_df)} MGEs for {genome_id}")
    
    # Load features
    rm_systems_df = safe_read(snakemake.input.rm_systems)
    rm_genes_df = safe_read(snakemake.input.rm_genes)
    cas_systems_df = safe_read(snakemake.input.cas_systems)
    cas_genes_df = safe_read(snakemake.input.cas_genes)
    defense_systems_df = safe_read(snakemake.input.defense_systems)
    defense_genes_df = safe_read(snakemake.input.defense_genes)
    crispr_arrays_df = safe_read(snakemake.input.crispr_arrays)
    crispr_spacers_df = safe_read(snakemake.input.crispr_spacers)
    
    # Assign MGE IDs to genes
    rm_genes_df = assign_mge_to_features(rm_genes_df, mges_df)
    cas_genes_df = assign_mge_to_features(cas_genes_df, mges_df)
    defense_genes_df = assign_mge_to_features(defense_genes_df, mges_df)
    
    # Assign MGE IDs to arrays
    crispr_arrays_df = assign_mge_to_features(crispr_arrays_df, mges_df)
    
    # Assign MGE IDs to spacers from parent arrays
    crispr_spacers_df = assign_mge_to_spacers(crispr_spacers_df, crispr_arrays_df)
    
    # Calculate system mobility
    rm_systems_df = calculate_system_mobility(rm_systems_df, rm_genes_df)
    cas_systems_df = calculate_system_mobility(cas_systems_df, cas_genes_df)
    defense_systems_df = calculate_system_mobility(defense_systems_df, defense_genes_df)
    
    # Print stats
    print(f"RM genes on MGEs: {(rm_genes_df['mge_id'].notna()).sum()}/{len(rm_genes_df)}")
    print(f"Cas genes on MGEs: {(cas_genes_df['mge_id'].notna()).sum()}/{len(cas_genes_df)}")
    print(f"CRISPR arrays on MGEs: {(crispr_arrays_df['mge_id'].notna()).sum()}/{len(crispr_arrays_df)}")
    print(f"RM systems mobile: {(rm_systems_df['mobility_score'] > 0).sum()}/{len(rm_systems_df)}")
    print(f"Cas systems mobile: {(cas_systems_df['mobility_score'] > 0).sum()}/{len(cas_systems_df)}")
    
    # Save outputs
    out_dir = Path(snakemake.output.rm_systems).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    
    rm_systems_df.to_csv(snakemake.output.rm_systems, sep='\t', index=False, na_rep='NULL')
    rm_genes_df.to_csv(snakemake.output.rm_genes, sep='\t', index=False, na_rep='NULL')
    cas_systems_df.to_csv(snakemake.output.cas_systems, sep='\t', index=False, na_rep='NULL')
    cas_genes_df.to_csv(snakemake.output.cas_genes, sep='\t', index=False, na_rep='NULL')
    defense_systems_df.to_csv(snakemake.output.defense_systems, sep='\t', index=False, na_rep='NULL')
    defense_genes_df.to_csv(snakemake.output.defense_genes, sep='\t', index=False, na_rep='NULL')
    crispr_arrays_df.to_csv(snakemake.output.crispr_arrays, sep='\t', index=False, na_rep='NULL')
    crispr_spacers_df.to_csv(snakemake.output.crispr_spacers, sep='\t', index=False, na_rep='NULL')
    
    print(f"Saved outputs to {out_dir}")


if __name__ == '__main__':
    main(snakemake)  # noqa: F821

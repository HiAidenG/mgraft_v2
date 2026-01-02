#!/usr/bin/env python3
"""
Build master GFF3 file from feature-type-specific TSVs.

Reads from features/ directory:
- mges.tsv
- {cas,rm,defense,antidefense}_systems.tsv
- {cas,rm,defense,antidefense}_genes.tsv
- crispr_arrays.tsv
- crispr_spacers.tsv
- virulence_genes.tsv

All TSVs have final IDs already assigned, making GFF building straightforward.
This module dynamically reads all columns from input TSVs and writes them as
GFF attributes, avoiding hardcoded field mappings.
"""

import pandas as pd
from pathlib import Path


# Standard GFF columns that should not be written as attributes
GFF_STANDARD_COLS = {'contig', 'seqid', 'start', 'end', 'strand', 'score', 'phase', 'source', 'feature_type'}

# Debug/internal columns that should not be written to GFF attributes
DEBUG_COLS = {'genome_id', 'feature_type', 'tool_id'}


def format_gff_attributes(attr_dict):
    """Format attributes dict as GFF3-compliant attribute string."""
    # Remove None and NULL values
    clean_attrs = {k: v for k, v in attr_dict.items() 
                   if v is not None and v != 'NULL' and str(v) != 'nan' and not (isinstance(v, float) and pd.isna(v))}
    # Format as key=value pairs
    return ';'.join([f"{k}={v}" for k, v in clean_attrs.items()])


def write_gff_feature(f, seqid, source, feature_type, start, end, score, strand, phase, attributes):
    """Write a single GFF3 feature line."""
    score_str = '.' if score is None or pd.isna(score) else str(score)
    strand_str = '.' if strand is None or pd.isna(strand) else str(strand)
    phase_str = '.' if phase is None or pd.isna(phase) else str(phase)
    attr_str = format_gff_attributes(attributes)
    
    f.write(f"{seqid}\t{source}\t{feature_type}\t{start}\t{end}\t{score_str}\t{strand_str}\t{phase_str}\t{attr_str}\n")


def load_tsv(path):
    """Load TSV file, returning empty DataFrame if missing or empty."""
    if path is None or not Path(path).exists() or Path(path).stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep='\t')


def get_source_for_feature(feature_type):
    """Get appropriate source field for each feature type."""
    defensefinder_types = [
        'defense_system', 'defense_gene',
        'cas_system', 'cas_gene',
        'rm_system', 'rm_gene',
        'antidefense_system', 'antidefense_gene'
    ]
    if feature_type in defensefinder_types:
        return 'DefenseFinder'
    elif feature_type == 'mge':
        return 'proMGEflow'
    elif feature_type in ['CRISPR_array', 'CRISPR_spacer']:
        return 'CRISPRDetect'
    elif feature_type == 'virulence_gene':
        return 'VFDB'
    return 'mgraft_v2'


def build_features_from_tsv(df, feature_type):
    """
    Build GFF features from any TSV by dynamically reading all columns.
    
    Standard GFF columns (contig, start, end, strand, score, phase) are used
    for the main GFF fields. All other columns become GFF attributes.
    """
    if df.empty:
        return []
    
    features = []
    source = get_source_for_feature(feature_type)
    
    for _, row in df.iterrows():
        # Build attributes from all non-standard columns (exclude GFF fields and debug columns)
        attrs = {}
        excluded_cols = GFF_STANDARD_COLS | DEBUG_COLS
        for col in df.columns:
            if col.lower() not in excluded_cols and col not in excluded_cols:
                val = row.get(col)
                if val is not None and not (isinstance(val, float) and pd.isna(val)) and str(val) != 'nan' and val != 'NULL':
                    attrs[col] = val
        
        features.append({
            'contig': row.get('contig') or row.get('seqid'),
            'start': row['start'],
            'end': row['end'],
            'source': source,
            'feature_type': feature_type,
            'score': row.get('score'),
            'strand': row.get('strand', '.'),
            'phase': row.get('phase'),
            'ID': row.get('ID'),
            'Parent': row.get('Parent'),
            'attributes': attrs
        })
    
    return features


def sort_features_hierarchical(features):
    """Sort features by coordinates while preserving parent-child relationships.
    
    Sort order: mge > system > array > gene = spacer
    Within each group, sort by (contig, start, end).
    Systems/arrays are followed immediately by their children (genes/spacers).
    
    Hierarchy:
    - MGEs are top-level standalone features
    - Systems (cas, rm, defense, antidefense) can be children of MGEs
    - Arrays can be children of systems (via Parent attribute)
    - Genes are children of their respective systems
    - Spacers are children of arrays
    
    Each feature is written exactly once, positioned after its parent.
    """
    from collections import defaultdict
    
    # Define feature type precedence (lower = higher precedence)
    type_precedence = {
        'mge': 0,
        'rm_system': 1, 'defense_system': 1, 'antidefense_system': 1, 'cas_system': 1,
        'CRISPR_array': 2,
        'rm_gene': 3, 'defense_gene': 3, 'antidefense_gene': 3, 'cas_gene': 3,
        'CRISPR_spacer': 3,
        'virulence_gene': 4
    }
    
    # Types that can have children
    parent_types = {'rm_system', 'defense_system', 'antidefense_system', 'cas_system', 'CRISPR_array'}
    
    # Build parent-child map and categorize features
    all_features_by_id = {}  # ID -> feature
    children = defaultdict(list)  # parent_ID -> [children]
    top_level_features = []  # Features with no parent or parent not in our feature set
    
    # First pass: index all features by ID
    for feat in features:
        feat_id = feat.get('ID')
        if feat_id:
            all_features_by_id[feat_id] = feat
    
    # Second pass: categorize as top-level or child
    for feat in features:
        feat_id = feat.get('ID')
        parent_id = feat.get('Parent')
        
        # Check if this feature has a parent that exists in our feature set
        if parent_id and parent_id in all_features_by_id:
            children[parent_id].append(feat)
        else:
            # Top-level feature (no parent, or parent not in our set like external MGE reference)
            top_level_features.append(feat)
    
    # Sort key for coordinates
    def sort_key(f):
        return (f['contig'], f['start'], f['end'], type_precedence.get(f['feature_type'], 99))
    
    # Recursive function to add a feature and its descendants
    def add_with_descendants(feat, result_list, visited):
        feat_id = feat.get('ID')
        
        # Avoid duplicates
        if feat_id in visited:
            return
        visited.add(feat_id)
        
        result_list.append(feat)
        
        # Add children sorted by position
        if feat_id in children:
            for child in sorted(children[feat_id], key=sort_key):
                add_with_descendants(child, result_list, visited)
    
    # Build result by processing top-level features in order
    result = []
    visited = set()
    
    # Group top-level features by contig, then sort
    by_contig = defaultdict(list)
    for feat in top_level_features:
        by_contig[feat['contig']].append(feat)
    
    for contig in sorted(by_contig.keys()):
        contig_features = sorted(by_contig[contig], key=sort_key)
        for feat in contig_features:
            add_with_descendants(feat, result, visited)
    
    return result


def build_master_gff(genome_id, mges_tsv, cas_systems_tsv, cas_genes_tsv,
                     rm_systems_tsv, rm_genes_tsv, defense_systems_tsv, defense_genes_tsv,
                     antidefense_systems_tsv, antidefense_genes_tsv,
                     arrays_tsv, spacers_tsv, virulence_genes_tsv, output_gff):
    """Build master GFF3 file from feature-type-specific TSVs.
    
    Dynamically reads all columns from input TSVs and writes them as GFF attributes.
    No hardcoded field mappings - upstream TSVs define the schema.
    """
    features = []
    
    # Define TSV -> feature type mappings
    tsv_feature_map = [
        (mges_tsv, 'mge'),
        (cas_systems_tsv, 'cas_system'),
        (cas_genes_tsv, 'cas_gene'),
        (rm_systems_tsv, 'rm_system'),
        (rm_genes_tsv, 'rm_gene'),
        (defense_systems_tsv, 'defense_system'),
        (defense_genes_tsv, 'defense_gene'),
        (antidefense_systems_tsv, 'antidefense_system'),
        (antidefense_genes_tsv, 'antidefense_gene'),
        (arrays_tsv, 'CRISPR_array'),
        (spacers_tsv, 'CRISPR_spacer'),
        (virulence_genes_tsv, 'virulence_gene'),
    ]
    
    # Load and build features from each TSV
    for tsv_path, feature_type in tsv_feature_map:
        df = load_tsv(tsv_path)
        if not df.empty:
            features.extend(build_features_from_tsv(df, feature_type))
    
    # Sort features hierarchically
    features = sort_features_hierarchical(features)
    
    # Write GFF3
    with open(output_gff, 'w') as f:
        f.write("##gff-version 3\n")
        f.write("##mgraft-schema-version 1.2.0\n")
        f.write(f"##genome-id {genome_id}\n")
        
        for feat in features:
            write_gff_feature(
                f,
                seqid=feat['contig'],
                source=feat['source'],
                feature_type=feat['feature_type'],
                start=feat['start'],
                end=feat['end'],
                score=feat['score'],
                strand=feat['strand'],
                phase=feat['phase'],
                attributes=feat['attributes']
            )


def main(snakemake):
    """Main entry point for Snakemake."""
    import shutil
    
    # mges_tsv is optional - may not be provided if no MGEs in samplesheet
    mges_tsv = getattr(snakemake.input, 'mges', None)
    if mges_tsv == []:
        mges_tsv = None
    
    build_master_gff(
        genome_id=snakemake.params.genome_id,
        mges_tsv=mges_tsv,
        cas_systems_tsv=snakemake.input.cas_systems,
        cas_genes_tsv=snakemake.input.cas_genes,
        rm_systems_tsv=snakemake.input.rm_systems,
        rm_genes_tsv=snakemake.input.rm_genes,
        defense_systems_tsv=snakemake.input.defense_systems,
        defense_genes_tsv=snakemake.input.defense_genes,
        antidefense_systems_tsv=snakemake.input.antidefense_systems,
        antidefense_genes_tsv=snakemake.input.antidefense_genes,
        arrays_tsv=snakemake.input.crispr_arrays,
        spacers_tsv=snakemake.input.crispr_spacers,
        virulence_genes_tsv=snakemake.input.virulence_genes,
        output_gff=snakemake.output.gff
    )
    
    # Copy feature-specific TSVs to final output directory
    out_tsv_dir = Path(snakemake.params.out_tsv_dir)
    out_tsv_dir.mkdir(parents=True, exist_ok=True)
    
    # Define which TSVs to copy and their output names
    genome_id = snakemake.params.genome_id
    tsvs_to_copy = [
        (snakemake.input.cas_systems, 'cas_systems.tsv'),
        (snakemake.input.cas_genes, 'cas_genes.tsv'),
        (snakemake.input.rm_systems, 'rm_systems.tsv'),
        (snakemake.input.rm_genes, 'rm_genes.tsv'),
        (snakemake.input.defense_systems, 'defense_systems.tsv'),
        (snakemake.input.defense_genes, 'defense_genes.tsv'),
        (snakemake.input.antidefense_systems, 'antidefense_systems.tsv'),
        (snakemake.input.antidefense_genes, 'antidefense_genes.tsv'),
        (snakemake.input.virulence_genes, 'virulence_genes.tsv'),
    ]
    
    # For arrays and spacers, filter to this genome from the clustered file
    for src_path, out_name in tsvs_to_copy:
        if src_path and Path(src_path).exists():
            shutil.copy(src_path, out_tsv_dir / out_name)
    
    # Handle arrays and spacers (filter from clustered files)
    arrays_clustered = snakemake.input.crispr_arrays
    spacers_clustered = snakemake.input.crispr_spacers
    
    if arrays_clustered and Path(arrays_clustered).exists():
        df = pd.read_csv(arrays_clustered, sep='\t')
        df_genome = df[df['genome_id'] == genome_id]
        df_genome.to_csv(out_tsv_dir / 'crispr_arrays.tsv', sep='\t', index=False, na_rep='NULL')
    
    if spacers_clustered and Path(spacers_clustered).exists():
        df = pd.read_csv(spacers_clustered, sep='\t')
        df_genome = df[df['genome_id'] == genome_id]
        df_genome.to_csv(out_tsv_dir / 'crispr_spacers.tsv', sep='\t', index=False, na_rep='NULL')
    
    # Create marker file
    tsv_marker = Path(snakemake.output.tsv_marker)
    tsv_marker.parent.mkdir(parents=True, exist_ok=True)
    tsv_marker.touch()
    
    print(f"Built GFF: {snakemake.output.gff}")
    print(f"Wrote TSVs to: {out_tsv_dir}")


if __name__ == '__main__':
    main(snakemake)  # noqa: F821

#!/usr/bin/env python3
"""
Parse DefenseFinder outputs to extract defense features with deterministic ID assignment.

Outputs feature-type-specific TSVs with final IDs ready for GFF building:
- rm_systems.tsv, rm_genes.tsv
- defense_systems.tsv, defense_genes.tsv  
- antidefense_systems.tsv, antidefense_genes.tsv
- cas_systems.tsv, cas_genes.tsv
"""
import pandas as pd
from pathlib import Path
from Bio import SeqIO


def parse_faa_coords(faa_path):
    """Parse FAA file to get gene coordinates.
    
    Args:
        faa_path: Path to FAA file
        
    Returns:
        dict: gene_id -> {'contig': str, 'start': int, 'end': int, 'strand': str}
    """
    coords = {}
    for record in SeqIO.parse(faa_path, 'fasta'):
        header = record.description
        # Header format: >genome_contig_pos # start # end # strand # ID=...
        parts = header.split(' # ')
        if len(parts) >= 4:
            gene_id = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            strand = '+' if parts[3] == '1' else '-'
            contig = '_'.join(gene_id.split('_')[:-1])  # Remove the last _pos part
            coords[gene_id] = {
                'contig': contig,
                'start': start,
                'end': end,
                'strand': strand
            }
    return coords


def system_feature_type(sys_type: str, activity: str) -> str:
    """Map DefenseFinder system type/activity to schema feature_type."""
    if sys_type == 'Cas':
        return 'cas_system'
    if sys_type == 'RM':
        return 'rm_system'
    if activity == 'Antidefense' or (sys_type and sys_type.startswith('Anti')):
        return 'antidefense_system'
    return 'defense_system'


def gene_feature_type(gene_type: str, activity: str) -> str:
    """Map DefenseFinder gene type/activity to schema feature_type."""
    if gene_type == 'Cas':
        return 'cas_gene'
    if gene_type == 'RM':
        return 'rm_gene'
    if activity == 'Antidefense' or (gene_type and gene_type.startswith('Anti')):
        return 'antidefense_gene'
    return 'defense_gene'


def normalize_rm_type(subtype_val: str):
    """Normalize RM type from DefenseFinder subtype to canonical schema format.
    
    DefenseFinder subtypes: RM_Type_I, RM_Type_II, RM_Type_IIG_2, RM_Type_III, RM_Type_IV_1
    Schema types: Type I, Type II, Type IIG, Type III, Type IV
    
    Args:
        subtype_val: DefenseFinder subtype string
        
    Returns:
        Normalized type string (e.g., 'Type II') or None
    """
    import re
    
    if not subtype_val or not isinstance(subtype_val, str):
        return None
    
    # Strip RM_ prefix if present
    if subtype_val.startswith('RM_'):
        subtype_val = subtype_val[3:]
    
    # Strip subtype suffix (e.g., _1, _2) - match underscore followed by digits at end
    subtype_val = re.sub(r'_\d+$', '', subtype_val)
    
    # Replace underscores with spaces (Type_I -> Type I)
    subtype_val = subtype_val.replace('_', ' ')
    
    return subtype_val


def classify_rm_gene_subunit(gene_name, system_type):
    """
    Classify RM gene by subunit type based on gene_name pattern.
    
    Args:
        gene_name: DefenseFinder gene_name field
        system_type: Normalized RM type (Type I, Type II, Type IIG, Type III, Type IV)
    
    Returns:
        One of: 'S', 'R', 'M', 'RM', or the original gene_name if unclassified
    """
    if not gene_name or not isinstance(gene_name, str):
        return gene_name
    
    gene_name_upper = gene_name.upper()
    
    if system_type == 'Type I':
        if '_S_' in gene_name_upper:
            return 'S'
        elif 'MTASE' in gene_name_upper:
            return 'M'
        elif 'REASE' in gene_name_upper:
            return 'R'
    
    elif system_type == 'Type II':
        if 'REASE' in gene_name_upper:
            return 'R'
        elif 'MTASE' in gene_name_upper:
            return 'M'
    
    elif system_type == 'Type IIG':
        # All Type IIG genes are bifunctional RM
        return 'RM'
    
    elif system_type == 'Type III':
        if 'REASE' in gene_name_upper:
            return 'R'
        elif 'MTASE' in gene_name_upper:
            return 'M'
    
    elif system_type == 'Type IV':
        if 'REASE' in gene_name_upper:
            return 'R'
    
    return gene_name  # Return original if no match


def parse_cas_class_and_subtype(subtype_val: str):
    """Extract cas class/subtype from DefenseFinder subtype string.
    
    Examples:
        'CAS_Class2-Subtype-II-C' -> ('Class 2', 'Type II-C')
        'CAS-Class1-Subtype-I-C' -> ('Class 1', 'Type I-C')
    """
    import re
    
    cas_class = None
    cas_subtype = None
    
    if not subtype_val or not isinstance(subtype_val, str):
        return cas_class, cas_subtype
    
    if 'CAS' in subtype_val.upper():
        # Extract class number (e.g., Class1, Class2)
        class_match = re.search(r'Class(\d+)', subtype_val, re.IGNORECASE)
        if class_match:
            cas_class = f"Class {class_match.group(1)}"
        
        # Extract subtype (everything after "Subtype-")
        subtype_match = re.search(r'Subtype[_-](.+)', subtype_val, re.IGNORECASE)
        if subtype_match:
            cas_subtype = f"Type {subtype_match.group(1)}"
    
    return cas_class, cas_subtype


def get_type_code(feature_type: str) -> str:
    """Get ID type code for a feature type."""
    codes = {
        'cas_system': 'CAS',
        'rm_system': 'RMS',
        'defense_system': 'DFS',
        'antidefense_system': 'ADS',
    }
    return codes.get(feature_type, 'UNK')


def parse_defensefinder(genes_tsv, systems_tsv, hmmer_tsv, faa_path, genome_id):
    """Parse DefenseFinder outputs into schema-aligned features with deterministic IDs.
    
    Args:
        genes_tsv: Path to DefenseFinder genes output
        systems_tsv: Path to DefenseFinder systems output
        hmmer_tsv: Path to DefenseFinder HMMER output
        faa_path: Path to protein FASTA for coordinate extraction
        genome_id: Genome identifier for deterministic ID assignment
        
    Returns:
        dict of DataFrames keyed by feature_type
    """
    coords = parse_faa_coords(faa_path)

    genes_input_df = pd.read_csv(genes_tsv, sep='\t') if Path(genes_tsv).exists() else pd.DataFrame()
    
    # Build sys_wholeness map from genes
    sys_wholeness_map = {}
    if not genes_input_df.empty and 'sys_wholeness' in genes_input_df.columns:
        sys_wholeness_map = genes_input_df.groupby('sys_id')['sys_wholeness'].max().to_dict()

    # Collect all systems and build sys_id to type map for RM gene classification
    all_systems = []
    sys_type_map = {}  # Maps sys_id -> normalized RM type (for RM systems only)
    
    if Path(systems_tsv).exists():
        systems_input_df = pd.read_csv(systems_tsv, sep='\t')
        for _, row in systems_input_df.iterrows():
            sys_id = row['sys_id']
            sys_type = row['type']
            subtype = row['subtype']
            activity = row['activity']
            sys_beg = row['sys_beg']
            sys_end = row['sys_end']
            sys_feat_type = system_feature_type(sys_type, activity)
            cas_class, cas_subtype = parse_cas_class_and_subtype(subtype) if sys_feat_type == 'cas_system' else (None, None)

            beg_coords = coords.get(sys_beg, {})
            end_coords = coords.get(sys_end, {})
            if not (beg_coords and end_coords):
                continue

            contig = beg_coords['contig']
            start = min(beg_coords['start'], end_coords['start'])
            end = max(beg_coords['end'], end_coords['end'])
            strand = beg_coords['strand']

            # Normalize RM type if this is an RM system
            normalized_type = normalize_rm_type(subtype) if sys_feat_type == 'rm_system' else subtype
            
            # Store normalized type for RM gene classification
            if sys_feat_type == 'rm_system':
                sys_type_map[sys_id] = normalized_type
            
            all_systems.append({
                'genome_id': genome_id,
                'feature_type': sys_feat_type,
                'tool_id': sys_id,
                'ID': None,
                'contig': contig,
                'start': start,
                'end': end,
                'strand': strand,
                'score': row.get('sys_score', None),
                'type': normalized_type if sys_feat_type != 'cas_system' else None,
                'sys_wholeness': sys_wholeness_map.get(sys_id),
                'n_genes': row.get('genes_count'),
                'cas_class': cas_class,
                'cas_subtype': cas_subtype,
            })

    # Collect all genes
    all_genes = []
    if not genes_input_df.empty:
        for _, row in genes_input_df.iterrows():
            hit_id = row['hit_id']
            gene_coords = coords.get(hit_id, {})
            if not gene_coords:
                continue

            raw_gene_name = row['gene_name']
            sys_id = row['sys_id']
            gene_type = row['type']
            activity = row['activity']
            feature_type = gene_feature_type(gene_type, activity)
            
            # Classify RM gene subunits (S, R, M, RM) based on system type
            if feature_type == 'rm_gene' and sys_id in sys_type_map:
                gene_name = classify_rm_gene_subunit(raw_gene_name, sys_type_map[sys_id])
            else:
                gene_name = raw_gene_name

            all_genes.append({
                'genome_id': genome_id,
                'feature_type': feature_type,
                'tool_id': hit_id,
                'tool_id_parent': sys_id,
                'ID': None,
                'Parent': None,
                'contig': gene_coords['contig'],
                'start': gene_coords['start'],
                'end': gene_coords['end'],
                'strand': gene_coords['strand'],
                'score': row.get('hit_score', None),
                'gene_name': gene_name,
            })

    systems_df = pd.DataFrame(all_systems)
    genes_df = pd.DataFrame(all_genes)
    
    # Assign deterministic IDs to systems (grouped by feature_type)
    tool_id_to_det_id = {}
    
    for feat_type in ['cas_system', 'rm_system', 'defense_system', 'antidefense_system']:
        type_mask = systems_df['feature_type'] == feat_type if not systems_df.empty else pd.Series(dtype=bool)
        if systems_df.empty or not type_mask.any():
            continue
        
        type_systems = systems_df[type_mask].sort_values(['contig', 'start'])
        type_code = get_type_code(feat_type)
        
        for i, (idx, row) in enumerate(type_systems.iterrows()):
            det_id = f"{genome_id}__{type_code}{i+1:04d}"
            tool_id_to_det_id[row['tool_id']] = det_id
            systems_df.at[idx, 'ID'] = det_id
    
    # Assign gene IDs based on parent system
    if not genes_df.empty:
        genes_df['ID'] = genes_df['ID'].astype('object')
        genes_df['Parent'] = genes_df['Parent'].astype('object')
        
        for tool_id, det_id in tool_id_to_det_id.items():
            sys_genes = genes_df[genes_df['tool_id_parent'] == tool_id].sort_values(['contig', 'start'])
            
            for gene_idx, (idx, gene_row) in enumerate(sys_genes.iterrows()):
                gene_det_id = f"{det_id}_g{gene_idx+1:02d}"
                genes_df.at[idx, 'ID'] = gene_det_id
                genes_df.at[idx, 'Parent'] = det_id
    
    # Split into separate DataFrames by feature type
    result = {}
    
    # Systems
    for feat_type in ['cas_system', 'rm_system', 'defense_system', 'antidefense_system']:
        if systems_df.empty:
            if feat_type == 'cas_system':
                result[feat_type] = get_empty_cas_system_df()
            else:
                result[feat_type] = get_empty_system_df()
            continue
            
        type_df = systems_df[systems_df['feature_type'] == feat_type].copy()
        base_cols = ['genome_id', 'feature_type', 'ID', 'contig', 'start', 'end', 'strand', 
                     'score', 'type', 'n_genes', 'sys_wholeness', 'tool_id']
        if feat_type == 'cas_system':
            base_cols.extend(['cas_class', 'cas_subtype'])
        
        cols = [c for c in base_cols if c in type_df.columns]
        result[feat_type] = type_df[cols].reset_index(drop=True)
    
    # Genes
    for feat_type in ['cas_gene', 'rm_gene', 'defense_gene', 'antidefense_gene']:
        if genes_df.empty:
            result[feat_type] = get_empty_gene_df()
            continue
            
        type_df = genes_df[genes_df['feature_type'] == feat_type].copy()
        base_cols = ['genome_id', 'feature_type', 'ID', 'Parent', 'contig', 'start', 'end', 
                     'strand', 'score', 'gene_name', 'tool_id']
        cols = [c for c in base_cols if c in type_df.columns]
        result[feat_type] = type_df[cols].reset_index(drop=True)
    
    return result


def get_empty_system_df():
    """Return empty DataFrame with system schema."""
    return pd.DataFrame(columns=[
        'genome_id', 'feature_type', 'ID', 'contig', 'start', 'end', 'strand',
        'score', 'type', 'n_genes', 'sys_wholeness', 'tool_id'
    ])


def get_empty_cas_system_df():
    """Return empty DataFrame with cas_system schema."""
    return pd.DataFrame(columns=[
        'genome_id', 'feature_type', 'ID', 'contig', 'start', 'end', 'strand',
        'score', 'type', 'n_genes', 'sys_wholeness', 'tool_id',
        'cas_class', 'cas_subtype'
    ])


def get_empty_gene_df():
    """Return empty DataFrame with gene schema."""
    return pd.DataFrame(columns=[
        'genome_id', 'feature_type', 'ID', 'Parent', 'contig', 'start', 'end',
        'strand', 'score', 'gene_name', 'tool_id'
    ])


def main(snakemake):
    """Main function for Snakemake script."""
    genes_tsv = snakemake.input.genes_tsv
    systems_tsv = snakemake.input.systems_tsv
    hmmer_tsv = snakemake.input.get('hmmer_tsv', None)
    faa_path = snakemake.input.faa
    genome_id = snakemake.params.genome_id
    
    # Parse DefenseFinder outputs
    features = parse_defensefinder(genes_tsv, systems_tsv, hmmer_tsv, faa_path, genome_id)
    
    # Write output TSVs - systems
    cas_systems = features.get('cas_system', get_empty_cas_system_df())
    rm_systems = features.get('rm_system', get_empty_system_df())
    defense_systems = features.get('defense_system', get_empty_system_df())
    antidefense_systems = features.get('antidefense_system', get_empty_system_df())
    
    cas_systems.to_csv(snakemake.output.cas_systems_tsv, sep='\t', index=False, na_rep='NULL')
    rm_systems.to_csv(snakemake.output.rm_systems_tsv, sep='\t', index=False, na_rep='NULL')
    defense_systems.to_csv(snakemake.output.defense_systems_tsv, sep='\t', index=False, na_rep='NULL')
    antidefense_systems.to_csv(snakemake.output.antidefense_systems_tsv, sep='\t', index=False, na_rep='NULL')
    
    # Write output TSVs - genes
    cas_genes = features.get('cas_gene', get_empty_gene_df())
    rm_genes = features.get('rm_gene', get_empty_gene_df())
    defense_genes = features.get('defense_gene', get_empty_gene_df())
    antidefense_genes = features.get('antidefense_gene', get_empty_gene_df())
    
    cas_genes.to_csv(snakemake.output.cas_genes_tsv, sep='\t', index=False, na_rep='NULL')
    rm_genes.to_csv(snakemake.output.rm_genes_tsv, sep='\t', index=False, na_rep='NULL')
    defense_genes.to_csv(snakemake.output.defense_genes_tsv, sep='\t', index=False, na_rep='NULL')
    antidefense_genes.to_csv(snakemake.output.antidefense_genes_tsv, sep='\t', index=False, na_rep='NULL')


if __name__ == '__main__':
    main(snakemake)

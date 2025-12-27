#!/usr/bin/env python3
"""
Generate per-genome summary statistics from all feature TSVs.

Aggregates counts, type breakdowns, quality metrics, mobility metrics, and derived
statistics for each genome in the dataset. Outputs a single genome_summary.tsv with
one row per genome.

Input files per genome (from features/ directory):
- cas_systems.linked.tsv
- crispr_arrays.linked.tsv
- crispr_spacers.linked.tsv
- defense_systems.tsv
- rm_systems.rebase.tsv
- virulence_genes.tsv
- antidefense_systems.tsv
- mges.tsv
"""
import pandas as pd
import subprocess
from pathlib import Path
from typing import Dict, Any, Optional


def safe_read_tsv(path: str) -> pd.DataFrame:
    """Read TSV file, returning empty DataFrame if file doesn't exist or is empty."""
    try:
        if Path(path).exists() and Path(path).stat().st_size > 0:
            return pd.read_csv(path, sep='\t')
    except Exception:
        pass
    return pd.DataFrame()


def count_by_wholeness(df: pd.DataFrame, complete_col: str = 'sys_wholeness') -> tuple:
    """Count complete (wholeness=1) and partial (wholeness<1) systems."""
    if df.empty or complete_col not in df.columns:
        return 0, 0
    complete = (df[complete_col] == 1.0).sum()
    partial = (df[complete_col] < 1.0).sum()
    return int(complete), int(partial)


def count_on_mge(df: pd.DataFrame, mobility_col: str = 'mobility') -> int:
    """Count features with mobility == 1.0 (fully on MGE)."""
    if df.empty or mobility_col not in df.columns:
        return 0
    return int((df[mobility_col] == 1.0).sum())


def count_with_mge_id(df: pd.DataFrame, mge_col: str = 'mge_id') -> int:
    """Count features with a non-null mge_id."""
    if df.empty or mge_col not in df.columns:
        return 0
    return int(df[mge_col].notna().sum())


def format_type_breakdown(df: pd.DataFrame, type_col: str) -> str:
    """Format type counts as semicolon-delimited string (e.g., 'I-E:2;II-A:1')."""
    if df.empty or type_col not in df.columns:
        return 'NULL'
    counts = df[type_col].value_counts()
    if counts.empty:
        return 'NULL'
    return ';'.join(f"{t}:{c}" for t, c in counts.items())


def safe_mean(series: pd.Series) -> Optional[float]:
    """Calculate mean, returning None if empty or all NaN."""
    if series.empty:
        return None
    valid = series.dropna()
    if valid.empty:
        return None
    return round(float(valid.mean()), 4)


def get_genome_stats(genome_fasta: str) -> Dict[str, Any]:
    """
    Get genome length and GC content using seqkit stats.
    
    Args:
        genome_fasta: Path to genome FASTA file
    
    Returns:
        Dictionary with 'genome_length' and 'gc_content' keys
    """
    result = {'genome_length': None, 'gc_content': None}
    
    if not genome_fasta or not Path(genome_fasta).exists():
        return result
    
    try:
        # Run seqkit stats to get genome length and GC content
        cmd = ['seqkit', 'stats', '-T', genome_fasta]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        if proc.returncode == 0:
            lines = proc.stdout.strip().split('\n')
            if len(lines) >= 2:
                # Header: file  format  type  num_seqs  sum_len  min_len  avg_len  max_len
                # Some versions include sum_gap and N50
                header = lines[0].split('\t')
                values = lines[1].split('\t')
                
                # Create dict from header/values
                stats = dict(zip(header, values))
                
                # Extract sum_len as genome_length
                if 'sum_len' in stats:
                    result['genome_length'] = int(stats['sum_len'].replace(',', ''))
        
        # Run seqkit fx2tab to get GC content
        cmd_gc = ['seqkit', 'fx2tab', '-n', '-g', genome_fasta]
        proc_gc = subprocess.run(cmd_gc, capture_output=True, text=True, timeout=60)
        
        if proc_gc.returncode == 0:
            lines = proc_gc.stdout.strip().split('\n')
            gc_values = []
            for line in lines:
                parts = line.split('\t')
                if len(parts) >= 2:
                    try:
                        gc_values.append(float(parts[1]))
                    except ValueError:
                        pass
            
            if gc_values:
                # Weighted average would be better, but simple mean is close enough
                result['gc_content'] = round(sum(gc_values) / len(gc_values), 2)
    
    except (subprocess.TimeoutExpired, FileNotFoundError, Exception):
        # seqkit not available or other error - return None values
        pass
    
    return result


def summarize_genome(sample: str, genomes_dir: str, genome_fasta: str = None) -> Dict[str, Any]:
    """
    Compute summary statistics for a single genome.
    
    Args:
        sample: Genome identifier
        genomes_dir: Base directory containing per-genome outputs
    
    Returns:
        Dictionary of summary statistics
    """
    features_dir = Path(genomes_dir) / sample / 'features'
    
    # Load all feature TSVs
    cas_systems = safe_read_tsv(features_dir / 'cas_systems.linked.tsv')
    crispr_arrays = safe_read_tsv(features_dir / 'crispr_arrays.linked.tsv')
    crispr_spacers = safe_read_tsv(features_dir / 'crispr_spacers.linked.tsv')
    defense_systems = safe_read_tsv(features_dir / 'defense_systems.tsv')
    rm_systems = safe_read_tsv(features_dir / 'rm_systems.rebase.tsv')
    virulence_genes = safe_read_tsv(features_dir / 'virulence_genes.tsv')
    antidefense_systems = safe_read_tsv(features_dir / 'antidefense_systems.tsv')
    mges = safe_read_tsv(features_dir / 'mges.tsv')
    
    # Initialize result with genome_id
    result: Dict[str, Any] = {'genome_id': sample}
    
    # ======================== Genome Stats ========================
    genome_stats = get_genome_stats(genome_fasta)
    result['genome_length'] = genome_stats['genome_length']
    result['gc_content'] = genome_stats['gc_content']
    
    # ======================== Basic Counts ========================
    # Cas systems (complete vs partial)
    cas_complete, cas_partial = count_by_wholeness(cas_systems)
    result['n_cas_systems_complete'] = cas_complete
    result['n_cas_systems_partial'] = cas_partial
    
    # CRISPR arrays and spacers
    result['n_crispr_arrays'] = len(crispr_arrays)
    result['n_spacers'] = len(crispr_spacers)
    
    # RM systems (complete vs partial)
    rm_complete, rm_partial = count_by_wholeness(rm_systems)
    result['n_rm_systems_complete'] = rm_complete
    result['n_rm_systems_partial'] = rm_partial
    
    # Defense systems (complete vs partial)
    def_complete, def_partial = count_by_wholeness(defense_systems)
    result['n_defense_systems_complete'] = def_complete
    result['n_defense_systems_partial'] = def_partial
    
    # Virulence genes
    result['n_virulence_genes'] = len(virulence_genes)
    
    # Antidefense systems (complete vs partial)
    antidef_complete, antidef_partial = count_by_wholeness(antidefense_systems)
    result['n_antidefense_systems_complete'] = antidef_complete
    result['n_antidefense_systems_partial'] = antidef_partial
    
    # MGEs
    result['n_mges'] = len(mges)
    
    # ======================== Type Breakdowns ========================
    result['cas_subtypes'] = format_type_breakdown(cas_systems, 'cas_subtype')
    result['rm_types'] = format_type_breakdown(rm_systems, 'type')
    result['defense_types'] = format_type_breakdown(defense_systems, 'type')
    result['antidefense_types'] = format_type_breakdown(antidefense_systems, 'type')
    result['mge_types'] = format_type_breakdown(mges, 'type')
    
    # ======================== Cas System Metrics ========================
    # Array types (canonical, orphan, putative)
    if not crispr_arrays.empty and 'array_type' in crispr_arrays.columns:
        result['n_canonical_arrays'] = int((crispr_arrays['array_type'] == 'canonical').sum())
        result['n_putative_arrays'] = int((crispr_arrays['array_type'] == 'putative').sum())
    else:
        result['n_canonical_arrays'] = 0
        result['n_putative_arrays'] = 0
    
    if not crispr_arrays.empty and 'is_orphan' in crispr_arrays.columns:
        # is_orphan could be bool or string
        if crispr_arrays['is_orphan'].dtype == bool:
            result['n_orphan_arrays'] = int(crispr_arrays['is_orphan'].sum())
        else:
            # Handle string 'True'/'true' values
            orphan_mask = crispr_arrays['is_orphan'].astype(str).str.lower() == 'true'
            result['n_orphan_arrays'] = int(orphan_mask.sum())
    else:
        result['n_orphan_arrays'] = 0
    
    # Orphan vs linked cas systems
    if not cas_systems.empty and 'array_ID' in cas_systems.columns:
        result['n_linked_cas'] = int(cas_systems['array_ID'].notna().sum())
        result['n_orphan_cas'] = int(cas_systems['array_ID'].isna().sum())
    else:
        result['n_orphan_cas'] = len(cas_systems)
        result['n_linked_cas'] = 0
    
    # Mean distance to cas (for linked arrays)
    if not crispr_arrays.empty and 'distance_to_cas' in crispr_arrays.columns:
        linked_arrays = crispr_arrays[crispr_arrays['distance_to_cas'].notna()]
        result['mean_cas_array_distance'] = safe_mean(linked_arrays['distance_to_cas'])
    else:
        result['mean_cas_array_distance'] = None
    
    # Mean quality scores (linked vs orphan)
    if not crispr_arrays.empty and 'quality_score' in crispr_arrays.columns:
        if 'is_orphan' in crispr_arrays.columns:
            # Handle both bool and string types for is_orphan
            if crispr_arrays['is_orphan'].dtype == bool:
                linked = crispr_arrays[~crispr_arrays['is_orphan']]
                orphan = crispr_arrays[crispr_arrays['is_orphan']]
            else:
                linked = crispr_arrays[crispr_arrays['is_orphan'].astype(str).str.lower() != 'true']
                orphan = crispr_arrays[crispr_arrays['is_orphan'].astype(str).str.lower() == 'true']
            result['mean_linked_array_quality'] = safe_mean(linked['quality_score'])
            result['mean_orphan_array_quality'] = safe_mean(orphan['quality_score'])
        else:
            result['mean_linked_array_quality'] = safe_mean(crispr_arrays['quality_score'])
            result['mean_orphan_array_quality'] = None
    else:
        result['mean_linked_array_quality'] = None
        result['mean_orphan_array_quality'] = None
    
    # Arrays with unknown strand
    if not crispr_arrays.empty and 'strand' in crispr_arrays.columns:
        result['num_arrays_no_strand'] = int((crispr_arrays['strand'] == '.').sum())
    else:
        result['num_arrays_no_strand'] = 0
    
    # ======================== Mobility Metrics ========================
    result['n_cas_on_mge'] = count_on_mge(cas_systems)
    result['n_rm_on_mge'] = count_on_mge(rm_systems)
    result['n_defense_on_mge'] = count_on_mge(defense_systems)
    result['n_antidefense_on_mge'] = count_on_mge(antidefense_systems)
    result['n_virulence_genes_on_mge'] = count_with_mge_id(virulence_genes)
    
    # ======================== Spacer Metrics ========================
    if not crispr_spacers.empty:
        if 'length' in crispr_spacers.columns:
            result['avg_spacer_length'] = safe_mean(crispr_spacers['length'])
        else:
            result['avg_spacer_length'] = None
        
        if 'entropy' in crispr_spacers.columns:
            result['avg_spacer_entropy'] = safe_mean(crispr_spacers['entropy'])
            result['n_spacers_low_entropy'] = int((crispr_spacers['entropy'] <= 1.0).sum())
        else:
            result['avg_spacer_entropy'] = None
            result['n_spacers_low_entropy'] = 0
    else:
        result['avg_spacer_length'] = None
        result['avg_spacer_entropy'] = None
        result['n_spacers_low_entropy'] = 0
    
    # ======================== MGE Metrics ========================
    if not mges.empty and 'start' in mges.columns and 'end' in mges.columns:
        result['total_mge_length'] = int((mges['end'] - mges['start']).sum())
    else:
        result['total_mge_length'] = 0
    
    # ======================== RM Quality Metrics ========================
    if not rm_systems.empty:
        # Check for REBASE_recseq column (from REBASE alignment)
        if 'REBASE_recseq' in rm_systems.columns:
            # Has recseq if not null and not 'NULL' string
            has_recseq = rm_systems['REBASE_recseq'].notna() & (rm_systems['REBASE_recseq'] != 'NULL') & (rm_systems['REBASE_recseq'] != '')
            
            # Count annotated systems by type
            if 'type' in rm_systems.columns:
                for rm_type in ['Type I', 'Type II', 'Type III', 'Type IIG']:
                    # Match exact type or type with space (e.g., "Type II" matches "Type II" but not "Type IIG")
                    if rm_type == 'Type II':
                        # Special handling: match "Type II" but not "Type IIG"
                        type_mask = (rm_systems['type'] == 'Type II')
                    else:
                        type_mask = rm_systems['type'].str.contains(rm_type, case=False, na=False, regex=False)
                    
                    col_name = f"n_annotated_{rm_type.lower().replace(' ', '_')}"
                    result[col_name] = int((type_mask & has_recseq).sum())
            else:
                result['n_annotated_type_i'] = 0
                result['n_annotated_type_ii'] = 0
                result['n_annotated_type_iii'] = 0
                result['n_annotated_type_iig'] = 0
            
            # Fraction annotated
            total_rm = len(rm_systems)
            n_annotated = has_recseq.sum()
            result['fraction_annotated_rm'] = round(n_annotated / total_rm, 4) if total_rm > 0 else None
        else:
            result['n_annotated_type_i'] = 0
            result['n_annotated_type_ii'] = 0
            result['n_annotated_type_iii'] = 0
            result['n_annotated_type_iig'] = 0
            result['fraction_annotated_rm'] = None
    else:
        result['n_annotated_type_i'] = 0
        result['n_annotated_type_ii'] = 0
        result['n_annotated_type_iii'] = 0
        result['n_annotated_type_iig'] = 0
        result['fraction_annotated_rm'] = None
    
    # ======================== Derived Metrics ========================
    # Defense density (systems per Mb)
    total_defense = (cas_complete + cas_partial + rm_complete + rm_partial + 
                     def_complete + def_partial)
    result['total_defense_systems'] = total_defense
    
    # Calculate defense density if genome length is available
    if result.get('genome_length') and result['genome_length'] > 0:
        result['defense_density'] = round(total_defense / (result['genome_length'] / 1_000_000), 4)
        result['mge_density'] = round(len(mges) / (result['genome_length'] / 1_000_000), 4)
    else:
        result['defense_density'] = None
        result['mge_density'] = None
    
    return result


def main(snakemake):
    """Main entry point for Snakemake script."""
    samples = snakemake.params.samples
    genomes_dir = snakemake.params.genomes_dir
    genome_fastas = snakemake.params.get('genome_fastas', {})
    output_path = snakemake.output[0]
    
    # Process each genome
    rows = []
    for sample in samples:
        genome_fasta = genome_fastas.get(sample) if genome_fastas else None
        row = summarize_genome(sample, genomes_dir, genome_fasta)
        rows.append(row)
    
    # Create DataFrame with consistent column order
    df = pd.DataFrame(rows)
    
    # Define column order (matching plan)
    column_order = [
        # Genome stats
        'genome_id', 'genome_length', 'gc_content',
        # Basic counts
        'n_cas_systems_complete', 'n_cas_systems_partial',
        'n_crispr_arrays', 'n_spacers',
        'n_rm_systems_complete', 'n_rm_systems_partial',
        'n_defense_systems_complete', 'n_defense_systems_partial',
        'n_virulence_genes',
        'n_antidefense_systems_complete', 'n_antidefense_systems_partial',
        'n_mges',
        # Type breakdowns
        'cas_subtypes', 'rm_types', 'defense_types', 'antidefense_types', 'mge_types',
        # Cas system metrics
        'n_canonical_arrays', 'n_orphan_arrays', 'n_putative_arrays',
        'n_orphan_cas', 'n_linked_cas',
        'mean_cas_array_distance', 'mean_linked_array_quality', 'mean_orphan_array_quality',
        'num_arrays_no_strand',
        # Mobility metrics
        'n_cas_on_mge', 'n_rm_on_mge', 'n_defense_on_mge', 'n_antidefense_on_mge',
        'n_virulence_genes_on_mge',
        # Spacer metrics
        'avg_spacer_length', 'avg_spacer_entropy', 'n_spacers_low_entropy',
        # MGE metrics
        'total_mge_length',
        # RM quality metrics
        'n_annotated_type_i', 'n_annotated_type_ii', 'n_annotated_type_iii', 'n_annotated_type_iig',
        'fraction_annotated_rm',
        # Derived metrics
        'total_defense_systems', 'defense_density', 'mge_density',
    ]
    
    # Reorder columns, keeping only those that exist
    existing_cols = [c for c in column_order if c in df.columns]
    df = df[existing_cols]
    
    # Write output
    df.to_csv(output_path, sep='\t', index=False, na_rep='NULL')


if __name__ == '__main__':
    main(snakemake)

#!/usr/bin/env python3
"""
Genome Summary Statistics - Script version

Aggregates statistics across all genomes in the dataset.
Creates all expected output files, even if empty.
"""
import gzip
import pandas as pd
from pathlib import Path
from collections import defaultdict


def safe_read_tsv(path):
    """Read TSV file, returning empty DataFrame if file doesn't exist or is empty."""
    if path and Path(path).exists() and Path(path).stat().st_size > 0:
        try:
            return pd.read_csv(path, sep='\t')
        except:
            pass
    return pd.DataFrame()


def main(snakemake):
    """Main function for Snakemake script."""
    # Access outputs/params
    output_summary = snakemake.output.summary
    output_cas = snakemake.output.cas
    output_crispr_arrays = snakemake.output.crispr_arrays
    output_crispr_spacers = snakemake.output.crispr_spacers
    output_rm_systems = snakemake.output.rm_systems
    
    samples = snakemake.params.samples
    input_dirs = snakemake.params.input_dirs
    genome_fastas = snakemake.params.get('genome_fastas', {})
    
    # Aggregators for tables
    all_tables = defaultdict(list)
    table_names = ['cas_systems', 'rm_systems', 'defense_systems', 'crispr_arrays', 'crispr_spacers', 'virulence_genes']
    
    genome_rows = []
    
    print(f"Processing {len(samples)} genomes...")
    
    for sample in samples:
        base_dir = Path(input_dirs.get(sample, ''))
        row = {'genome_id': sample}
        
        # Load features
        features = {}
        for name in table_names:
            df = safe_read_tsv(base_dir / f"{name}.tsv")
            features[name] = df
            
            # Accumulate for master tables
            if not df.empty:
                df_copy = df.copy()
                if 'genome_id' not in df_copy.columns:
                    df_copy.insert(0, 'genome_id', sample)
                all_tables[name].append(df_copy)
        
        # Basic counts
        row['n_cas_systems'] = len(features['cas_systems'])
        row['n_rm_systems'] = len(features['rm_systems'])
        row['n_defense_systems'] = len(features['defense_systems'])
        row['n_crispr_arrays'] = len(features['crispr_arrays'])
        row['n_spacers'] = len(features['crispr_spacers'])
        row['n_virulence_genes'] = len(features['virulence_genes'])
        
        # CRISPR metrics
        arr = features['crispr_arrays']
        if not arr.empty and 'array_class' in arr.columns:
            row['n_canonical_arrays'] = (arr['array_class'] == 'canonical').sum()
            row['n_orphan_arrays'] = (arr['array_class'] == 'orphan').sum()
        else:
            row['n_canonical_arrays'] = 0
            row['n_orphan_arrays'] = len(arr)
        
        genome_rows.append(row)
    
    # Write Genome Summary
    summary_df = pd.DataFrame(genome_rows)
    Path(output_summary).parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(output_summary, sep='\t', index=False, na_rep='NULL')
    print(f"Written summary for {len(summary_df)} genomes to {output_summary}")
    
    # Write Aggregated Tables - ALWAYS create output files even if empty
    summary_dir = Path(output_summary).parent
    
    # Map of output file names to their corresponding table names
    output_map = {
        'cas_systems': output_cas,
        'crispr_arrays': output_crispr_arrays,
        'crispr_spacers': output_crispr_spacers,
        'rm_systems': output_rm_systems,
    }
    
    for name, outfile in output_map.items():
        dfs = all_tables.get(name, [])
        if dfs:
            combined = pd.concat(dfs, ignore_index=True)
            combined.to_csv(outfile, sep='\t', index=False, na_rep='NULL', compression='gzip')
            print(f"Written aggregated {name} ({len(combined)} rows) to {outfile}")
        else:
            # Create empty gzipped file with empty content
            with gzip.open(outfile, 'wt') as f:
                f.write("")
            print(f"Created empty {name} output: {outfile}")
    
    # Also write defense_systems aggregated file (not a declared output, but helpful)
    if all_tables['defense_systems']:
        combined = pd.concat(all_tables['defense_systems'], ignore_index=True)
        outfile = summary_dir / "all_defense_systems.tsv.gz"
        combined.to_csv(outfile, sep='\t', index=False, na_rep='NULL', compression='gzip')
        print(f"Written aggregated defense_systems ({len(combined)} rows) to {outfile}")


if __name__ == '__main__':
    main(snakemake)  # noqa: F821

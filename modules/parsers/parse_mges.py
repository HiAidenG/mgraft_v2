#!/usr/bin/env python3
"""
Parse MGE GFF from proMGEflow to extract MGE metadata and sequences.
"""
import pandas as pd
from pathlib import Path
from Bio import SeqIO


def parse_mge_gff(gff_path, genome_id):
    """Parse MGE GFF file and extract metadata.
    
    Args:
        gff_path: Path to MGE GFF file
        genome_id: Genome identifier for deterministic ID assignment
        
    Returns:
        pd.DataFrame with columns: ID, type, contig, start, end, length, n_genes, is_nested,
                                   genome_type, recombinase
    """
    if not gff_path or gff_path == '' or not Path(gff_path).exists():
        # Return empty dataframe with correct schema
        return pd.DataFrame(columns=[
            'genome_id', 'feature_type', 'tool_id', 'ID', 'type', 'contig', 'start', 'end', 
            'strand', 'length', 'n_genes', 'is_nested', 'genome_type', 'recombinase'
        ])
    
    mges = []
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            contig_id, _, feature_type, start, end, _, strand = fields[0], fields[1], fields[2], fields[3], fields[4], fields[5], fields[6]
            attributes = fields[8]
            
            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value
            
            # Extract MGE information from mobile_genetic_element features
            if feature_type == 'mobile_genetic_element':
                mge_id = attr_dict.get('ID', f"{contig_id}_{start}_{end}")

                mge_attr = attr_dict.get('mge', '')
                mge_type_parts = [p.split(':')[0].strip().lower() for p in mge_attr.split(',') if ':' in p]
                mge_type_parts = [p for p in mge_type_parts if p]

                def short_form(t):
                    if t == 'mobility_island':
                        return 'mi'
                    if t == 'conjugative_element':
                        return 'ce'
                    return t

                mge_types = [short_form(t) for t in mge_type_parts] if mge_type_parts else ['unknown']
                mge_type = ','.join(mge_types)

                nested_flag_attr = attr_dict.get('mge_type', '').lower()
                is_nested_attr = attr_dict.get('is_nested', '').lower() == 'true'
                is_nested = nested_flag_attr == 'nested' or is_nested_attr or len(mge_types) > 1

                length = int(end) - int(start) + 1
                n_genes = int(attr_dict.get('n_genes', 0))
                
                # Extract genome_type attribute (e.g., 'ACC', 'CHR')
                genome_type_attr = attr_dict.get('genome_type', None)
                
                # Extract recombinase from mgeR attribute (e.g., 'huh_y1:1,phage_integrase:1' -> 'huh_y1,phage_integrase')
                mge_r = attr_dict.get('mgeR', '')
                if mge_r:
                    # Parse mgeR format: 'recombinase_type:count,recombinase_type:count,...'
                    recombinase_parts = [p.split(':')[0].strip() for p in mge_r.split(',') if ':' in p]
                    recombinase = ','.join(recombinase_parts) if recombinase_parts else None
                else:
                    recombinase = None

                mges.append({
                    'genome_id': genome_id,
                    'feature_type': 'mge',
                    'tool_id': mge_id,  # Keep original ID for intermediate tracking
                    'ID': None,  # Will be assigned deterministically later
                    'type': mge_type,
                    'contig': contig_id,
                    'start': int(start),
                    'end': int(end),
                    'strand': strand if strand in ['+', '-'] else '.',
                    'length': length,
                    'n_genes': n_genes,
                    'is_nested': is_nested,
                    'genome_type': genome_type_attr,
                    'recombinase': recombinase
                })
    
    df = pd.DataFrame(mges)
    
    # Assign deterministic IDs: sort by contig, start
    if not df.empty:
        df = df.sort_values(['contig', 'start'])
        df['ID'] = [f"{genome_id}__MGE{i+1:04d}" for i in range(len(df))]
    else:
        # Ensure ID column exists even for empty dataframes
        df['ID'] = pd.Series(dtype='object')
    
    return df


def extract_mge_sequences(mges_df, genome_fasta, output_fasta):
    """Extract MGE sequences from genome FASTA based on coordinates.
    
    Args:
        mges_df: DataFrame with MGE coordinates (contig, start, end, strand, ID)
        genome_fasta: Path to genome FASTA file
        output_fasta: Path to output FASTA file for MGE sequences
    """
    if mges_df.empty:
        # Write empty FASTA file
        Path(output_fasta).write_text('')
        return
    
    # Load genome sequences
    genome_seqs = {rec.id: rec for rec in SeqIO.parse(genome_fasta, 'fasta')}
    
    with open(output_fasta, 'w') as f:
        for _, row in mges_df.iterrows():
            contig = row['contig']
            if contig not in genome_seqs:
                print(f"WARNING: Contig {contig} not found in genome FASTA for {row['ID']}")
                continue
            
            # Extract sequence (convert to 0-based)
            start = int(row['start']) - 1
            end = int(row['end'])
            seq = genome_seqs[contig].seq[start:end]
            
            # Reverse complement if on minus strand
            if row.get('strand') == '-':
                seq = seq.reverse_complement()
            
            # Build header with key attributes
            attrs = f"type={row.get('type', '')};length={len(seq)}"
            f.write(f">{row['ID']} {attrs}\n")
            f.write(f"{seq}\n")


def main(snakemake):
    """Main function for Snakemake script."""
    mge_gff = snakemake.input.mge_gff
    genome_fasta = snakemake.input.genome_fasta
    genome_id = snakemake.params.genome_id
    output_tsv = snakemake.output.mges_tsv
    output_fasta = snakemake.output.mges_fasta
    
    # Parse MGE GFF
    mges_df = parse_mge_gff(mge_gff, genome_id)
    
    # Write TSV output with NULL for empty values
    mges_df.to_csv(output_tsv, sep='\t', index=False, na_rep='NULL')
    
    # Extract and write MGE sequences
    extract_mge_sequences(mges_df, genome_fasta, output_fasta)


if __name__ == '__main__':
    main(snakemake)

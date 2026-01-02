#!/usr/bin/env python3
"""
Parse MGE GFF from proMGEflow to extract MGE metadata and recombinase sequences.

Input: mge_islands.gff3 (ONLY this file - other proMGE outputs are intermediates)
Output: Standardized MGE TSV with recombinase nucleotide sequences for clustering
"""
import pandas as pd
from pathlib import Path
from Bio import SeqIO


def parse_mge_gff(gff_path, ffn_path, genome_id):
    """Parse mge_islands.gff3 and extract MGE metadata with recombinase sequences.
    
    Args:
        gff_path: Path to mge_islands.gff3 file
        ffn_path: Path to Prodigal FFN file for nucleotide sequence lookup
        genome_id: Genome identifier for deterministic ID assignment
        
    Returns:
        pd.DataFrame with schema columns
    """
    if not gff_path or gff_path == '' or not Path(gff_path).exists():
        return get_empty_mge_df()
    
    # Load FFN sequences for recombinase lookup
    ffn_seqs = {}
    if ffn_path and Path(ffn_path).exists():
        for rec in SeqIO.parse(ffn_path, 'fasta'):
            # Prodigal FFN headers: >contig_genenum # start # end # strand # ...
            ffn_seqs[rec.id] = str(rec.seq)
    
    # First pass: collect MGE records
    mge_records = {}  # mge_id -> record dict
    # Second pass: collect recombinase genes per MGE
    mge_recombinases = {}  # mge_id -> [(gene_id, recombinase_type, seq), ...]
    
    # Track ID matching statistics
    match_log = []
    
    with open(gff_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            contig_id = fields[0]
            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6] if fields[6] in ['+', '-'] else '.'
            attributes = fields[8]
            
            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value
            
            # Process MGE features
            if feature_type == 'mobile_genetic_element':
                mge_id = attr_dict.get('ID', f"{contig_id}_{start}_{end}")
                
                # Parse mge_type (nested vs non-nested)
                mge_type_attr = attr_dict.get('mge_type', 'non-nested').lower()
                mge_type = 'nested' if mge_type_attr == 'nested' else 'non-nested'
                
                # Parse mge class (is_tn, ce, phage, mi)
                mge_attr = attr_dict.get('mge', '')
                if mge_attr:
                    mge_classes = []
                    for part in mge_attr.split(','):
                        if ':' in part:
                            cls = part.split(':')[0].strip().lower()
                            # Normalize class names
                            if cls in ('is_tn', 'transposon', 'insertion_sequence'):
                                mge_classes.append('is_tn')
                            elif cls in ('ce', 'conjugative_element'):
                                mge_classes.append('ce')
                            elif cls in ('phage', 'prophage'):
                                mge_classes.append('phage')
                            elif cls in ('mi', 'integron', 'mobile_integron'):
                                mge_classes.append('mi')
                            else:
                                mge_classes.append(cls)
                    mge_class = ','.join(sorted(set(mge_classes))) if mge_classes else 'unknown'
                else:
                    mge_class = 'unknown'
                
                # Parse recombinase types from mgeR attribute
                recombinase_type = attr_dict.get('mgeR', None)
                
                mge_records[mge_id] = {
                    'genome_id': genome_id,
                    'feature_type': 'MGE',
                    'tool_id': mge_id,
                    'ID': None,
                    'contig': contig_id,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'mge_type': mge_type,
                    'mge_class': mge_class,
                    'size': int(attr_dict.get('size', end - start + 1)),
                    'n_genes': int(attr_dict.get('n_genes', 0)),
                    'recombinase_type': recombinase_type,
                    'recombinase_ids': None,
                    'recombinase_seqs': None
                }
                mge_recombinases[mge_id] = []
            
            # Process gene features to find recombinases
            elif feature_type == 'gene':
                parent_mge = attr_dict.get('Parent', '')
                recombinase_type = attr_dict.get('recombinase', None)
                
                if parent_mge and recombinase_type:
                    gene_id = attr_dict.get('ID', '')
                    
                    # Try to find nucleotide sequence in FFN
                    # Gene ID format varies - try exact match first
                    seq = ffn_seqs.get(gene_id, None)
                    match_type = 'exact' if seq else None
                    matched_ffn_id = gene_id if seq else None
                    
                    # If not found, try matching by coordinate pattern
                    if not seq and gene_id:
                        # Some FFN IDs use underscore pattern: contig_genenum
                        for ffn_id, ffn_seq in ffn_seqs.items():
                            if gene_id.endswith(ffn_id.split('_')[-1]) or ffn_id.endswith(gene_id.split('_')[-1]):
                                seq = ffn_seq
                                match_type = 'fuzzy'
                                matched_ffn_id = ffn_id
                                break
                    
                    # Log match result for debugging
                    if match_type:
                        match_log.append({
                            'gene_id': gene_id,
                            'ffn_id': matched_ffn_id,
                            'match_type': match_type,
                            'has_seq': seq is not None
                        })
                    else:
                        match_log.append({
                            'gene_id': gene_id,
                            'ffn_id': None,
                            'match_type': 'none',
                            'has_seq': False
                        })
                    
                    if parent_mge not in mge_recombinases:
                        mge_recombinases[parent_mge] = []
                    
                    mge_recombinases[parent_mge].append({
                        'gene_id': gene_id,
                        'recombinase_type': recombinase_type,
                        'seq': seq
                    })
    
    # Merge recombinase info into MGE records
    for mge_id, rec in mge_records.items():
        recomb_list = mge_recombinases.get(mge_id, [])
        if recomb_list:
            rec['recombinase_ids'] = ','.join(r['gene_id'] for r in recomb_list if r['gene_id'])
    
    # Print ID matching statistics
    if match_log:
        exact_count = sum(1 for m in match_log if m['match_type'] == 'exact')
        fuzzy_count = sum(1 for m in match_log if m['match_type'] == 'fuzzy')
        none_count = sum(1 for m in match_log if m['match_type'] == 'none')
        print(f"Recombinase ID matching stats for {genome_id}: "
              f"{exact_count} exact, {fuzzy_count} fuzzy, {none_count} unmatched")
        
        # Print details for fuzzy matches (for debugging)
        for m in match_log:
            if m['match_type'] == 'fuzzy':
                print(f"  FUZZY: {m['gene_id']} -> {m['ffn_id']}")
            elif m['match_type'] == 'none':
                print(f"  UNMATCHED: {m['gene_id']}")
    
    df = pd.DataFrame(list(mge_records.values()))
    
    # Assign deterministic IDs
    if not df.empty:
        df = df.sort_values(['contig', 'start'])
        df['ID'] = [f"{genome_id}__MGE{i+1:04d}" for i in range(len(df))]
    
    # Ensure correct column order
    columns = [
        'genome_id', 'feature_type', 'ID', 'tool_id', 'contig', 'start', 'end', 'strand',
        'mge_type', 'mge_class', 'size', 'n_genes', 'recombinase_type',
        'recombinase_ids'
    ]
    
    for col in columns:
        if col not in df.columns:
            df[col] = None
    
    return df[columns]


def get_empty_mge_df():
    """Return empty DataFrame with MGE schema."""
    return pd.DataFrame(columns=[
        'genome_id', 'feature_type', 'ID', 'tool_id', 'contig', 'start', 'end', 'strand',
        'mge_type', 'mge_class', 'size', 'n_genes', 'recombinase_type',
        'recombinase_ids'
    ])


def extract_mge_sequences(mges_df, genome_fasta, output_fasta):
    """Extract MGE sequences from genome FASTA based on coordinates.
    
    Args:
        mges_df: DataFrame with MGE coordinates (contig, start, end, strand, ID)
        genome_fasta: Path to genome FASTA file
        output_fasta: Path to output FASTA file for MGE sequences
    """
    if mges_df.empty:
        Path(output_fasta).write_text('')
        return
    
    genome_seqs = {rec.id: rec for rec in SeqIO.parse(genome_fasta, 'fasta')}
    
    with open(output_fasta, 'w') as f:
        for _, row in mges_df.iterrows():
            contig = row['contig']
            if contig not in genome_seqs:
                print(f"WARNING: Contig {contig} not found in genome FASTA for {row['ID']}")
                continue
            
            start = int(row['start']) - 1
            end = int(row['end'])
            seq = genome_seqs[contig].seq[start:end]
            
            if row.get('strand') == '-':
                seq = seq.reverse_complement()
            
            attrs = f"type={row.get('mge_class', '')};size={len(seq)}"
            f.write(f">{row['ID']} {attrs}\n")
            f.write(f"{seq}\n")


def main(snakemake):
    """Main function for Snakemake script."""
    mge_gff = snakemake.input.mge_gff
    ffn_path = snakemake.input.get('ffn', None)
    genome_fasta = snakemake.input.genome
    genome_id = snakemake.params.genome_id
    output_tsv = snakemake.output.mges_tsv
    output_fasta = snakemake.output.mges_fasta
    
    # Parse MGE GFF with FFN for recombinase sequences
    mges_df = parse_mge_gff(mge_gff, ffn_path, genome_id)
    
    # Write TSV output
    mges_df.to_csv(output_tsv, sep='\t', index=False, na_rep='NULL')
    
    # Extract MGE island sequences
    extract_mge_sequences(mges_df, genome_fasta, output_fasta)
    
    # Stats
    n_with_recomb = (mges_df['recombinase_ids'].notna()).sum() if not mges_df.empty else 0
    print(f"Parsed {len(mges_df)} MGEs, {n_with_recomb} with recombinase genes")


if __name__ == '__main__':
    main(snakemake)

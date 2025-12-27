#!/usr/bin/env python3
"""
Map spacers and RM recognition sites to MGE sequences.
"""
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq


def extract_mge_sequences(genome_fasta, mges_df):
    """Extract MGE sequences from genome."""
    sequences = {rec.id: str(rec.seq) for rec in SeqIO.parse(genome_fasta, 'fasta')}
    
    mge_sequences = {}
    for _, mge in mges_df.iterrows():
        contig_seq = sequences.get(mge['contig_id'])
        if contig_seq:
            mge_seq = contig_seq[mge['start']-1:mge['end']]  # 0-based slicing
            mge_sequences[mge['mge_id']] = {
                'sequence': mge_seq,
                'contig_id': mge['contig_id'],
                'mge_start': mge['start'],
                'mge_end': mge['end']
            }
    
    return mge_sequences


def search_mge_sequences(mge_sequences, queries, complement=True):
    """Search for query sequences in MGE sequences."""
    hits = []
    
    for query_id, query_seq, query_type in queries:
        query_upper = query_seq.upper()
        query_rc = str(Seq(query_upper).reverse_complement())
        
        for mge_id, mge_data in mge_sequences.items():
            mge_seq_upper = mge_data['sequence'].upper()
            
            # Search forward strand
            pos = 0
            while True:
                pos = mge_seq_upper.find(query_upper, pos)
                if pos == -1:
                    break
                
                # Convert to genomic coordinates
                genomic_start = mge_data['mge_start'] + pos
                genomic_end = genomic_start + len(query_upper) - 1
                
                hits.append({
                    'query_id': query_id,
                    'query_sequence': query_seq,
                    'query_type': query_type,
                    'contig_id': mge_data['contig_id'],
                    'start': genomic_start,
                    'end': genomic_end,
                    'strand': '+',
                    'score': 100.0,
                    'mge_id': mge_id
                })
                pos += 1
            
            # Search reverse complement
            if complement:
                pos = 0
                while True:
                    pos = mge_seq_upper.find(query_rc, pos)
                    if pos == -1:
                        break
                    
                    genomic_start = mge_data['mge_start'] + pos
                    genomic_end = genomic_start + len(query_rc) - 1
                    
                    hits.append({
                        'query_id': query_id,
                        'query_sequence': query_seq,
                        'query_type': query_type,
                        'contig_id': mge_data['contig_id'],
                        'start': genomic_start,
                        'end': genomic_end,
                        'strand': '-',
                        'score': 100.0,
                        'mge_id': mge_id
                    })
                    pos += 1
    
    return hits


def create_gff_line(hit, run_id):
    """Create GFF3 format line from hit with MGE information."""
    query_type = 'protospacer' if hit['query_type'] == 'spacer' else 'restriction_site'
    
    attributes = (f"run_id={run_id};query_id={hit['query_id']};"
                 f"query_sequence={hit['query_sequence']};mge_id={hit['mge_id']}")
    
    return '\t'.join([
        hit['contig_id'],
        'mgraft_v2',
        query_type,
        str(hit['start']),
        str(hit['end']),
        str(hit['score']),
        hit['strand'],
        '.',
        attributes
    ])


def main(snakemake):  # noqa: F821
    """Main function for Snakemake script."""
    genome = snakemake.input.genome
    mges_tsv = snakemake.input.mges_tsv
    spacers_tsv = snakemake.input.spacers
    alignments_tsv = snakemake.input.alignments
    mobile_targets = snakemake.output.mobile_targets
    
    run_id = snakemake.params.run_id
    
    # Load MGEs
    try:
        mges_df = pd.read_csv(mges_tsv, sep='\t')
    except pd.errors.EmptyDataError:
        mges_df = pd.DataFrame()
    
    if mges_df.empty:
        # Create empty output
        Path(mobile_targets).touch()
        return
    
    # Extract MGE sequences
    mge_sequences = extract_mge_sequences(genome, mges_df)
    
    # Load queries
    queries = []
    
    # Add spacers
    if Path(spacers_tsv).exists():
        try:
            spacers_df = pd.read_csv(spacers_tsv, sep='\t')
            for _, spacer in spacers_df.iterrows():
                if spacer['entropy'] > 1.0:
                    queries.append((f"spacer|{spacer['spacer_id']}", spacer['sequence'], 'spacer'))
        except pd.errors.EmptyDataError:
            pass
    
    # Add RM recognition sites from alignments
    # Note: This is a placeholder - REBASE alignments don't directly give recognition sites
    # In practice, you would need to parse REBASE metadata or use another approach
    # For now, we skip RM site mapping until recognition sites can be properly extracted
    if Path(alignments_tsv).exists():
        try:
            alignments_df = pd.read_csv(alignments_tsv, sep='\t')
            # Future: Parse recognition sites from REBASE metadata
            # rebase_alignments = alignments_df[alignments_df['database'] == 'REBASE']
            pass
        except pd.errors.EmptyDataError:
            pass
    
    if not queries:
        Path(mobile_targets).touch()
        return
    
    # Search MGE sequences
    hits = search_mge_sequences(mge_sequences, queries, complement=True)
    
    # Write GFF
    with open(mobile_targets, 'w') as out:
        out.write("##gff-version 3\n")
        for hit in hits:
            line = create_gff_line(hit, run_id)
            out.write(line + '\n')


if __name__ == '__main__':
    main(snakemake)  # noqa: F821

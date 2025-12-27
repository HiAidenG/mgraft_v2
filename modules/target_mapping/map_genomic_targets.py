#!/usr/bin/env python3
"""
Map spacers and RM recognition sites to genomic sequences.
"""
import pandas as pd
from pathlib import Path
from Bio import SeqIO


def search_exact_matches(genome_fasta, queries, complement=True):
    """Search for exact matches of query sequences in genome.
    
    Returns list of dicts with: query_id, contig_id, start, end, strand, score
    """
    from Bio.Seq import Seq
    
    hits = []
    sequences = {rec.id: str(rec.seq) for rec in SeqIO.parse(genome_fasta, 'fasta')}
    
    for query_id, query_seq in queries:
        query_upper = query_seq.upper()
        query_rc = str(Seq(query_upper).reverse_complement())
        
        for contig_id, genome_seq in sequences.items():
            genome_upper = genome_seq.upper()
            
            # Search forward strand
            pos = 0
            while True:
                pos = genome_upper.find(query_upper, pos)
                if pos == -1:
                    break
                hits.append({
                    'query_id': query_id,
                    'query_sequence': query_seq,
                    'contig_id': contig_id,
                    'start': pos + 1,  # 1-based
                    'end': pos + len(query_upper),
                    'strand': '+',
                    'score': 100.0
                })
                pos += 1
            
            # Search reverse complement
            if complement:
                pos = 0
                while True:
                    pos = genome_upper.find(query_rc, pos)
                    if pos == -1:
                        break
                    hits.append({
                        'query_id': query_id,
                        'query_sequence': query_seq,
                        'contig_id': contig_id,
                        'start': pos + 1,
                        'end': pos + len(query_rc),
                        'strand': '-',
                        'score': 100.0
                    })
                    pos += 1
    
    return hits


def create_gff_line(hit, run_id, query_type):
    """Create GFF3 format line from hit."""
    attributes = f"run_id={run_id};query_id={hit['query_id']};query_sequence={hit['query_sequence']}"
    
    return '\t'.join([
        hit['contig_id'],
        'mgraft_v2',
        'protospacer' if query_type == 'spacer' else 'restriction_site',
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
    spacers_tsv = snakemake.input.spacers
    alignments_tsv = snakemake.input.alignments
    genomic_targets = snakemake.output.genomic_targets
    
    run_id = snakemake.params.run_id
    backend = snakemake.params.backend
    
    # Load queries
    queries = []
    
    # Add spacers
    if Path(spacers_tsv).exists():
        try:
            spacers_df = pd.read_csv(spacers_tsv, sep='\t')
            for _, spacer in spacers_df.iterrows():
                if spacer['entropy'] > 1.0:  # Filter low-complexity spacers
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
        # Create empty output
        Path(genomic_targets).touch()
        return
    
    # Search for matches
    if backend == 'exact':
        query_list = [(qid, seq) for qid, seq, _ in queries]
        hits = search_exact_matches(genome, query_list, complement=True)
    else:
        raise ValueError(f"Unsupported backend: {backend}")
    
    # Write GFF
    with open(genomic_targets, 'w') as out:
        out.write("##gff-version 3\n")
        for hit in hits:
            # Determine query type
            query_type = 'spacer' if hit['query_id'].startswith('spacer|') else 'rm'
            line = create_gff_line(hit, run_id, query_type)
            out.write(line + '\n')


if __name__ == '__main__':
    main(snakemake)  # noqa: F821

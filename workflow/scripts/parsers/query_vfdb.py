#!/usr/bin/env python3
"""
Query VFDB database with all proteins.
Aligns entire proteome against virulence factor database.
"""

import pandas as pd
import subprocess
import sys
import re
from pathlib import Path
from Bio import SeqIO


def run_diamond(query_fasta, db_path, output_tsv, threads=4, min_identity=90.0,
                min_query_cov=80.0, min_subject_cov=80.0, evalue=1e-5):
    """Run Diamond BLASTP alignment."""
    # Check if query exists and not empty
    if not Path(query_fasta).exists() or Path(query_fasta).stat().st_size == 0:
        Path(output_tsv).touch()
        return
    
    cmd = [
        'diamond', 'blastp',
        '--db', db_path,
        '--query', query_fasta,
        '--out', output_tsv,
        '--threads', str(threads),
        '--id', str(min_identity),
        '--query-cover', str(min_query_cov),
        '--subject-cover', str(min_subject_cov),
        '--evalue', str(evalue),
        '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen',
                   'qcovhsp', 'scovhsp', 'evalue', 'bitscore', 'stitle',
        '--max-target-seqs', '1',  # Keep only best hit
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Diamond error: {result.stderr}", file=sys.stderr)
        raise RuntimeError(f"Diamond failed with return code {result.returncode}")


def parse_vfdb_header(stitle):
    """
    Parse VFDB header to extract metadata.
    Format: VFG037170(gb|WP_001081754) (plc1) phospholipase C [Phospholipase C (VF0470) - Exotoxin (VFC0235)] [Acinetobacter baumannii 1656-2]
    
    Returns: dict with VF_id, gene_name, product, vf_category, organism
    """
    result = {
        'VF_id': None,
        'VF_gene': None,
        'VF_product': None,
        'VF_category': None,
        'VF_organism': None,
        'VF_accession': None
    }
    
    if not stitle or pd.isna(stitle):
        return result
    
    # Extract VF ID (e.g., VFG037170)
    vf_match = re.search(r'^(VFG\d+)', stitle)
    if vf_match:
        result['VF_id'] = vf_match.group(1)
    
    # Extract accession (e.g., gb|WP_001081754)
    acc_match = re.search(r'\(([^)]+)\)', stitle)
    if acc_match:
        result['VF_accession'] = acc_match.group(1)
    
    # Extract gene name (in parentheses after accession)
    gene_match = re.search(r'\)\s+\(([^)]+)\)', stitle)
    if gene_match:
        result['VF_gene'] = gene_match.group(1)
    
    # Extract product description (between gene and first bracket)
    prod_match = re.search(r'\)\s+\([^)]+\)\s+([^\[]+)', stitle)
    if prod_match:
        result['VF_product'] = prod_match.group(1).strip()
    
    # Extract VF category (e.g., "Exotoxin (VFC0235)")
    cat_match = re.search(r'\[([^\]]+)\s+-\s+([^\]]+)\]', stitle)
    if cat_match:
        result['VF_category'] = f"{cat_match.group(1)} - {cat_match.group(2)}"
    
    # Extract organism (last bracketed item)
    org_match = re.search(r'\[([^\]]+)\]$', stitle)
    if org_match:
        result['VF_organism'] = org_match.group(1).strip()
    
    return result


def parse_orf_coords(orf_id, protein_fasta):
    """
    Parse ORF ID to extract actual contig, start, end, strand from protein FASTA.
    Prodigal format: >ID # start # end # strand # ...
    
    Args:
        orf_id: e.g. 'MGYG000001432_1_878'
        protein_fasta: Path to protein FASTA file
        
    Returns:
        dict: {'contig': str, 'start': int, 'end': int, 'strand': str} or None
    """
    # Try to read coordinates from FASTA header
    try:
        for record in SeqIO.parse(protein_fasta, 'fasta'):
            if record.id == orf_id:
                # Parse Prodigal header format: >ID # start # end # strand
                desc_parts = record.description.split('#')
                if len(desc_parts) >= 4:
                    start = int(desc_parts[1].strip())
                    end = int(desc_parts[2].strip())
                    strand_val = int(desc_parts[3].strip())
                    strand = '+' if strand_val == 1 else '-'
                    
                    # Extract contig from ID (format: genome_contig_gene)
                    contig_prefix = '_'.join(orf_id.rsplit('_', 1)[:-1])
                    
                    return {'contig': contig_prefix, 'start': start, 'end': end, 'strand': strand}
    except Exception as e:
        print(f"Warning: Could not parse coordinates for {orf_id}: {e}", file=sys.stderr)
    
    # Fallback: use gene number
    parts = orf_id.rsplit('_', 1)
    if len(parts) == 2:
        try:
            contig_prefix = parts[0]
            gene_num = int(parts[1])
            return {'contig': contig_prefix, 'start': gene_num, 'end': gene_num, 'strand': '+'}
        except ValueError:
            pass
    return None


def process_alignments(alignments_tsv, genome_id, mges_tsv, protein_fasta):
    """Process diamond alignments and parse VFDB headers.
    
    Args:
        alignments_tsv: Path to diamond output
        genome_id: Genome identifier
        mges_tsv: Path to MGEs TSV for mge_id assignment
        protein_fasta: Path to protein FASTA for coordinate extraction
    
    Returns:
        DataFrame with virulence gene features
    """
    # Check if file exists and not empty
    if not Path(alignments_tsv).exists() or Path(alignments_tsv).stat().st_size == 0:
        # Return empty dataframe with schema columns
        return pd.DataFrame(columns=[
            'genome_id', 'feature_type', 'ID', 'contig', 'start', 'end', 'strand',
            'mge_id', 'VF_id', 'VF_gene', 'VF_product', 'VF_category',
            'VF_organism', 'VF_accession', 'identity', 'query_coverage',
            'target_coverage', 'alignment_evalue', 'tool_id'
        ])
    
    # Load MGEs for mge_id assignment
    mges_df = pd.DataFrame()
    if mges_tsv and Path(mges_tsv).exists():
        mges_df = pd.read_csv(mges_tsv, sep='\t')
    
    # Read diamond output
    aln_df = pd.read_csv(
        alignments_tsv,
        sep='\t',
        header=None,
        names=['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen',
               'qcovhsp', 'scovhsp', 'evalue', 'bitscore', 'stitle']
    )
    
    # Parse VFDB headers
    vfdb_info = aln_df['stitle'].apply(parse_vfdb_header).apply(pd.Series)
    aln_df = pd.concat([aln_df, vfdb_info], axis=1)
    
    # Parse ORF coordinates (now includes actual coordinates and strand)
    coord_info = aln_df['qseqid'].apply(lambda x: parse_orf_coords(x, protein_fasta)).apply(pd.Series)
    aln_df = pd.concat([aln_df, coord_info], axis=1)
    
    # Sort by contig and start position for deterministic ID assignment
    aln_df = aln_df.sort_values(['contig', 'start', 'end']).reset_index(drop=True)
    
    # Assign deterministic IDs
    aln_df['ID'] = [f"{genome_id}__VFG{i+1:04d}" for i in range(len(aln_df))]
    
    # Add metadata columns
    aln_df['genome_id'] = genome_id
    aln_df['feature_type'] = 'virulence_gene'
    aln_df['tool_id'] = aln_df['qseqid']
    aln_df['mge_id'] = None
    
    # Assign mge_id based on coordinate overlap
    if not mges_df.empty:
        aln_df['mge_id'] = aln_df['mge_id'].astype('object')
        for idx, row in aln_df.iterrows():
            for _, mge_row in mges_df.iterrows():
                if (row['contig'] == mge_row['contig'] and
                    row['start'] >= mge_row['start'] and
                    row['end'] <= mge_row['end']):
                    aln_df.at[idx, 'mge_id'] = mge_row['ID']
                    break
    
    # Rename and select columns
    aln_df = aln_df.rename(columns={
        'pident': 'identity',
        'qcovhsp': 'query_coverage',
        'scovhsp': 'target_coverage',
        'evalue': 'alignment_evalue'
    })
    
    # Select output columns matching schema
    output_cols = [
        'genome_id', 'feature_type', 'ID', 'contig', 'start', 'end', 'strand',
        'mge_id', 'VF_id', 'VF_gene', 'VF_product', 'VF_category',
        'VF_organism', 'VF_accession', 'identity', 'query_coverage',
        'target_coverage', 'alignment_evalue', 'tool_id'
    ]
    
    return aln_df[output_cols]


def main(snakemake):
    """Main entry point."""
    import tempfile
    import json
    from Bio import SeqIO as SeqIOLocal
    
    # Create temp file for diamond output
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as tmp_aln:
        tmp_aln_path = tmp_aln.name
    
    try:
        run_diamond(
            query_fasta=snakemake.input.proteins,
            db_path=snakemake.input.vfdb_db,
            output_tsv=tmp_aln_path,
            threads=snakemake.threads,
            min_identity=snakemake.params.min_identity,
            min_query_cov=snakemake.params.min_query_cov,
            min_subject_cov=snakemake.params.min_subject_cov,
            evalue=snakemake.params.evalue
        )
        
        results_df = process_alignments(
            alignments_tsv=tmp_aln_path,
            genome_id=snakemake.params.genome_id,
            mges_tsv=snakemake.input.get('mges_tsv', None),
            protein_fasta=snakemake.input.proteins
        )
        
        # Write diamond hits
        import shutil
        shutil.copy(tmp_aln_path, snakemake.output.hits)
        
        # Write final output
        Path(snakemake.output.parsed).parent.mkdir(parents=True, exist_ok=True)
        results_df.to_csv(snakemake.output.parsed, sep='\t', index=False, na_rep='NULL')
        
        # Generate stats
        proteins_queried = sum(1 for _ in SeqIOLocal.parse(snakemake.input.proteins, 'fasta')) if Path(snakemake.input.proteins).exists() else 0
        hits_found = len(results_df)
        
        if hits_found > 0:
            min_identity = float(results_df['identity'].min())
            avg_identity = float(results_df['identity'].mean())
            max_evalue = float(results_df['alignment_evalue'].max())
        else:
            min_identity = 0.0
            avg_identity = 0.0
            max_evalue = 0.0
        
        stats = {
            'proteins_queried': proteins_queried,
            'hits_found': hits_found,
            'min_identity': round(min_identity, 2),
            'avg_identity': round(avg_identity, 2),
            'max_evalue': max_evalue
        }
        
        with open(snakemake.output.stats, 'w') as f:
            json.dump(stats, f, indent=2)
    
    finally:
        # Cleanup temp file
        Path(tmp_aln_path).unlink(missing_ok=True)


if __name__ == '__main__':
    main(snakemake)


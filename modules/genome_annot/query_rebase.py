#!/usr/bin/env python3
"""
Query REBASE database with RM gene proteins.
Extracts RM genes, aligns to REBASE, merges annotations into rm_genes.tsv and rm_systems.tsv.
"""

import pandas as pd
import subprocess
import sys
import re
from pathlib import Path
from Bio import SeqIO
import tempfile


def extract_rm_proteins_by_type(rm_genes_tsv, rm_systems_tsv, faa_path, output_dir):
    """
    Extract RM gene protein sequences grouped by system type.
    Includes genes from ALL systems regardless of sys_wholeness.
    Type IV systems are skipped entirely (no recognition sequence).
    
    Returns:
        dict: {rm_type: (fasta_path, genes_df)} for each RM type present
    """
    # Read RM genes and systems
    if not Path(rm_genes_tsv).exists() or Path(rm_genes_tsv).stat().st_size == 0:
        return {}
    
    rm_genes_df = pd.read_csv(rm_genes_tsv, sep='\t')
    if rm_genes_df.empty:
        return {}
    
    # Read systems to get type info
    if not Path(rm_systems_tsv).exists() or Path(rm_systems_tsv).stat().st_size == 0:
        return {}
    
    rm_systems_df = pd.read_csv(rm_systems_tsv, sep='\t')
    if rm_systems_df.empty:
        return {}
    
    # Filter out Type IV systems (no recognition sequence)
    non_type_iv_systems = rm_systems_df[rm_systems_df['type'] != 'Type IV'].copy()
    
    if non_type_iv_systems.empty:
        return {}
    
    # Merge genes with system types (via Parent relationship)
    rm_genes_df = rm_genes_df.merge(
        non_type_iv_systems[['ID', 'type']],
        left_on='Parent',
        right_on='ID',
        how='inner',  # Only keep genes from non-Type IV systems
        suffixes=('', '_system')
    )
    rm_genes_df.rename(columns={'type': 'system_type'}, inplace=True)
    rm_genes_df.drop(columns=['ID_system'], inplace=True, errors='ignore')
    
    # Group by system type
    type_groups = {}
    for rm_type, type_genes in rm_genes_df.groupby('system_type', dropna=True):
        # Get gene IDs for this type
        gene_ids = set(type_genes['tool_id'].dropna())
        
        # Create type-specific FASTA
        type_fasta = Path(output_dir) / f"rm_{rm_type.replace(' ', '_').lower()}.faa"
        with open(type_fasta, 'w') as out:
            for record in SeqIO.parse(faa_path, 'fasta'):
                if record.id in gene_ids:
                    SeqIO.write(record, out, 'fasta')
        
        # Store path and genes dataframe
        type_groups[rm_type] = (str(type_fasta), type_genes)
    
    return type_groups


def run_diamond(query_fasta, db_path, output_tsv, threads=4, min_identity=100.0,
                min_query_cov=90.0, min_subject_cov=90.0, evalue=1e-10):
    """Run Diamond BLASTP alignment."""
    # Check if query is empty
    if Path(query_fasta).stat().st_size == 0:
        # Create empty output
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
        '--sensitive'
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Diamond error: {result.stderr}", file=sys.stderr)
        raise RuntimeError(f"Diamond failed with return code {result.returncode}")


def parse_rebase_header(stitle):
    """
    Parse REBASE header to extract metadata.
    Format: REBASE|M.Aac9709I RecSeq:GATC Org:Aggregatibacter actinomycetemcomitans NCTC9709 ProteinId:SSY84327.1 UniProt:NA
    
    Returns: dict with accession, organism, recseq, protein_id
    """
    result = {
        'REBASE_accession': None,
        'REBASE_organism': None,
        'REBASE_recseq': None,
        'REBASE_protein_id': None
    }
    
    if not stitle or pd.isna(stitle):
        return result
    
    # Extract accession (e.g., M.Aac9709I)
    acc_match = re.search(r'REBASE\|([^\s]+)', stitle)
    if acc_match:
        result['REBASE_accession'] = acc_match.group(1)
    
    # Extract recognition sequence
    recseq_match = re.search(r'RecSeq:([^\s]+)', stitle)
    if recseq_match:
        result['REBASE_recseq'] = recseq_match.group(1)
    
    # Extract organism (between "Org:" and "ProteinId:")
    org_match = re.search(r'Org:(.+?)\s+ProteinId:', stitle)
    if org_match:
        result['REBASE_organism'] = org_match.group(1).strip()
    
    # Extract protein ID
    prot_match = re.search(r'ProteinId:([^\s]+)', stitle)
    if prot_match:
        result['REBASE_protein_id'] = prot_match.group(1)
    
    return result


def merge_alignments(rm_genes_df, alignments_tsv):
    """Merge diamond alignments with RM gene metadata."""
    # Read alignments if file exists and not empty
    if not Path(alignments_tsv).exists() or Path(alignments_tsv).stat().st_size == 0:
        # No alignments - add NULL columns
        rm_genes_df = rm_genes_df.copy()
        rm_genes_df['REBASE_accession'] = None
        rm_genes_df['REBASE_protein_id'] = None
        rm_genes_df['identity'] = None
        rm_genes_df['query_coverage'] = None
        rm_genes_df['target_coverage'] = None
        rm_genes_df['alignment_evalue'] = None
        rm_genes_df['REBASE_organism'] = None
        rm_genes_df['REBASE_recseq'] = None
        return rm_genes_df
    
    # Read diamond output
    aln_df = pd.read_csv(
        alignments_tsv,
        sep='\t',
        header=None,
        names=['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen',
               'qcovhsp', 'scovhsp', 'evalue', 'bitscore', 'stitle']
    )
    
    # Parse REBASE headers
    rebase_info = aln_df['stitle'].apply(parse_rebase_header).apply(pd.Series)
    aln_df = pd.concat([aln_df, rebase_info], axis=1)
    
    # Rename columns to match schema
    aln_df = aln_df.rename(columns={
        'qseqid': 'tool_id',
        'pident': 'identity',
        'qcovhsp': 'query_coverage',
        'scovhsp': 'target_coverage',
        'evalue': 'alignment_evalue'
    })
    
    # Select relevant columns
    aln_cols = ['tool_id', 'identity', 'query_coverage', 'target_coverage',
                'alignment_evalue', 'REBASE_accession', 'REBASE_organism', 'REBASE_recseq',
                'REBASE_protein_id']
    aln_df = aln_df[aln_cols]
    
    # Merge with rm_genes (left join to keep all RM genes)
    merged = rm_genes_df.merge(aln_df, on='tool_id', how='left')
    
    return merged


def _get_gene_recseq(gene_row):
    """Extract recseq from gene row, returning None if missing/NULL."""
    recseq = gene_row.get('REBASE_recseq')
    if recseq and pd.notna(recseq) and recseq != 'NULL':
        return recseq
    return None


def validate_type_i_recseq(sys_genes):
    """
    Validate Type I system recseq.
    
    Type I requires R, M, and S subunits present, but only S determines recseq.
    - If no S genes: missing_annotation
    - If S genes present but none have recseq: missing_annotation
    - If multiple S genes with different recseqs: conflicting
    - If S gene(s) have single unique recseq AND R and M present: annotated
    
    Returns: (recseq, recseq_note)
    """
    # Classify genes using pre-classified gene_name (S, R, M)
    s_genes = []
    has_r = False
    has_m = False
    
    for _, gene_row in sys_genes.iterrows():
        gene_name = gene_row.get('gene_name', '')
        if gene_name == 'S':
            s_genes.append(gene_row)
        elif gene_name == 'R':
            has_r = True
        elif gene_name == 'M':
            has_m = True
    
    # Must have R and M subunits
    if not has_r or not has_m:
        return None, 'missing_annotation'
    
    # Must have at least one S gene
    if not s_genes:
        return None, 'missing_annotation'
    
    # Collect recseqs from S genes
    s_recseqs = []
    for gene_row in s_genes:
        recseq = _get_gene_recseq(gene_row)
        if recseq:
            s_recseqs.append(recseq)
    
    if not s_recseqs:
        return None, 'missing_annotation'
    
    unique_recseqs = set(s_recseqs)
    if len(unique_recseqs) == 1:
        return s_recseqs[0], 'annotated'
    else:
        return None, 'conflicting'


def validate_type_ii_iii_recseq(sys_genes, system_type):
    """
    Validate Type II or Type III system recseq.
    
    Requires at least one cognate R+M pair with matching recseqs.
    - If no R or no M genes: missing_annotation
    - If R and M genes present but no matching recseq pair: conflicting (if both have recseqs) or missing_annotation
    - If at least one R+M pair has matching recseq: annotated
    
    Returns: (recseq, recseq_note)
    """
    # Classify genes using pre-classified gene_name (R, M)
    r_genes = []
    m_genes = []
    
    for _, gene_row in sys_genes.iterrows():
        gene_name = gene_row.get('gene_name', '')
        if gene_name == 'R':
            r_genes.append(gene_row)
        elif gene_name == 'M':
            m_genes.append(gene_row)
    
    # Must have at least one R and one M gene
    if not r_genes or not m_genes:
        return None, 'missing_annotation'
    
    # Collect recseqs from R and M genes
    r_recseqs = set()
    m_recseqs = set()
    
    for gene_row in r_genes:
        recseq = _get_gene_recseq(gene_row)
        if recseq:
            r_recseqs.add(recseq)
    
    for gene_row in m_genes:
        recseq = _get_gene_recseq(gene_row)
        if recseq:
            m_recseqs.add(recseq)
    
    # Find matching recseqs between R and M
    matching_recseqs = r_recseqs & m_recseqs
    
    if matching_recseqs:
        # At least one valid R+M pair exists
        if len(matching_recseqs) == 1:
            return matching_recseqs.pop(), 'annotated'
        else:
            # Multiple different matching pairs - conflicting
            return None, 'conflicting'
    
    # No matching pairs
    if r_recseqs and m_recseqs:
        # Both have recseqs but they don't match
        return None, 'conflicting'
    else:
        # At least one side is missing recseqs
        return None, 'missing_annotation'


def validate_type_iig_recseq(sys_genes):
    """
    Validate Type IIG system recseq.
    
    Type IIG has bifunctional RM genes. All genes are treated as RM.
    - If no genes have recseq: missing_annotation
    - If multiple genes with different recseqs: conflicting
    - If gene(s) have single unique recseq: annotated
    
    Returns: (recseq, recseq_note)
    """
    recseqs = []
    
    for _, gene_row in sys_genes.iterrows():
        recseq = _get_gene_recseq(gene_row)
        if recseq:
            recseqs.append(recseq)
    
    if not recseqs:
        return None, 'missing_annotation'
    
    unique_recseqs = set(recseqs)
    if len(unique_recseqs) == 1:
        return recseqs[0], 'annotated'
    else:
        return None, 'conflicting'


def propagate_recseq_to_systems(rm_genes_df, rm_systems_df):
    """Propagate REBASE_recseq from rm_genes to rm_systems using type-specific logic.
    
    Type-specific rules:
    - Type I: Requires R, M, S subunits. Only S subunit recseq is used.
    - Type II/III: Requires at least one cognate R+M pair with matching recseqs.
    - Type IIG: Bifunctional RM genes - all must have same recseq.
    - Type IV: No recognition sequence (no_recognition_site).
    
    Note: sys_wholeness is NOT checked - annotation depends only on required genes being present.
    """
    if rm_systems_df.empty:
        return rm_systems_df
    
    rm_systems_df = rm_systems_df.copy()
    rm_systems_df['REBASE_recseq'] = None
    rm_systems_df['recseq_note'] = None
    
    for idx, sys_row in rm_systems_df.iterrows():
        sys_id = sys_row['ID']
        sys_type = sys_row.get('type', '')
        
        # Type IV systems have no recognition sequence - check first before gene lookup
        if sys_type == 'Type IV':
            rm_systems_df.at[idx, 'recseq_note'] = 'no_recognition_site'
            continue
        
        # Get genes for this system
        sys_genes = rm_genes_df[rm_genes_df['Parent'] == sys_id]
        
        if sys_genes.empty:
            rm_systems_df.at[idx, 'recseq_note'] = 'missing_annotation'
            continue
        
        # Dispatch to type-specific validation
        if sys_type == 'Type I':
            recseq, note = validate_type_i_recseq(sys_genes)
        elif sys_type in ('Type II', 'Type III'):
            recseq, note = validate_type_ii_iii_recseq(sys_genes, sys_type)
        elif sys_type == 'Type IIG':
            recseq, note = validate_type_iig_recseq(sys_genes)
        else:
            # Unknown type - fallback to missing
            recseq, note = None, 'missing_annotation'
        
        rm_systems_df.at[idx, 'REBASE_recseq'] = recseq
        rm_systems_df.at[idx, 'recseq_note'] = note
    
    return rm_systems_df


def get_alignment_params(config, rm_type):
    """
    Get alignment parameters for a specific RM type.
    Falls back to global defaults if type-specific config is missing.
    
    Args:
        config: query_rebase config dict
        rm_type: RM system type (e.g., 'Type I', 'Type II', etc.)
    
    Returns:
        dict with threads, min_identity, min_query_cov, min_subject_cov, evalue
    """
    # Map schema type names to config keys
    type_map = {
        'Type I': 'type_I',
        'Type II': 'type_II',
        'Type IIG': 'type_IIG',
        'Type III': 'type_III',
        'Type IV': 'type_IV'
    }
    
    # Global defaults
    defaults = {
        'threads': config.get('threads', 4),
        'min_identity': config.get('min_identity', 100.0),
        'min_query_cov': config.get('min_query_coverage', 90.0),
        'min_subject_cov': config.get('min_subject_coverage', 90.0),
        'evalue': config.get('evalue', 1e-10),
        'skip_alignment': False
    }
    
    # Get type-specific config if available
    config_key = type_map.get(rm_type)
    if config_key and config_key in config:
        type_config = config[config_key]
        return {
            'threads': defaults['threads'],  # threads always from global
            'min_identity': type_config.get('min_identity', defaults['min_identity']),
            'min_query_cov': type_config.get('min_query_coverage', defaults['min_query_cov']),
            'min_subject_cov': type_config.get('min_subject_coverage', defaults['min_subject_cov']),
            'evalue': type_config.get('evalue', defaults['evalue']),
            'skip_alignment': type_config.get('skip_alignment', False)
        }
    
    return defaults


def main(snakemake):
    """Main entry point for Snakemake."""
    # Inputs
    rm_genes_tsv = snakemake.input.rm_genes_tsv
    rm_systems_tsv = snakemake.input.rm_systems_tsv
    faa_path = snakemake.input.proteins
    rebase_db = snakemake.input.rebase_db
    
    # Outputs
    rm_genes_out = snakemake.output.rm_genes_tsv
    rm_systems_out = snakemake.output.rm_systems_tsv
    
    # Get full config dict
    config = snakemake.params.get('query_rebase_config', {})
    
    # Load rm_systems
    rm_systems_df = pd.read_csv(rm_systems_tsv, sep='\t') if Path(rm_systems_tsv).exists() else pd.DataFrame()
    
    # Create temp directory for type-specific FASTAs
    temp_dir = tempfile.mkdtemp(prefix='rebase_')
    temp_files = []
    
    try:
        # Step 1: Extract RM proteins grouped by type
        type_groups = extract_rm_proteins_by_type(rm_genes_tsv, rm_systems_tsv, faa_path, temp_dir)
        
        if not type_groups:
            # No RM genes - write empty outputs with schema columns
            empty_genes = pd.DataFrame(columns=[
                'genome_id', 'feature_type', 'ID', 'Parent', 'contig', 'start', 'end',
                'strand', 'score', 'gene_name', 'mge_id', 'tool_id',
                'REBASE_accession', 'REBASE_protein_id', 'identity', 'query_coverage',
                'target_coverage', 'alignment_evalue', 'REBASE_organism', 'REBASE_recseq',
                'system_type'
            ])
            empty_genes.to_csv(rm_genes_out, sep='\t', index=False, na_rep='NULL')
            
            # Write systems with recseq columns
            if not rm_systems_df.empty:
                rm_systems_df['REBASE_recseq'] = None
                rm_systems_df['recseq_note'] = 'incomplete_system'
            rm_systems_df.to_csv(rm_systems_out, sep='\t', index=False, na_rep='NULL')
            return
        
        # Step 2: Run Diamond alignment for each type with type-specific parameters
        all_genes_dfs = []
        
        for rm_type, (type_fasta, type_genes_df) in type_groups.items():
            temp_files.append(type_fasta)
            
            # Get type-specific alignment parameters
            params = get_alignment_params(config, rm_type)
            
            # Skip alignment for Type IV if configured
            if params['skip_alignment']:
                print(f"Skipping alignment for {rm_type} (skip_alignment=true)", file=sys.stderr)
                all_genes_dfs.append(type_genes_df)
                continue
            
            # Run Diamond for this type
            tmp_aln = tempfile.NamedTemporaryFile(mode='w', suffix=f'_{rm_type.replace(" ", "_")}.tsv', delete=False)
            tmp_aln_path = tmp_aln.name
            tmp_aln.close()
            temp_files.append(tmp_aln_path)
            
            print(f"Aligning {rm_type} genes with params: identity={params['min_identity']}, "
                  f"query_cov={params['min_query_cov']}, subject_cov={params['min_subject_cov']}, "
                  f"evalue={params['evalue']}", file=sys.stderr)
            
            run_diamond(
                type_fasta, rebase_db, tmp_aln_path,
                threads=params['threads'],
                min_identity=params['min_identity'],
                min_query_cov=params['min_query_cov'],
                min_subject_cov=params['min_subject_cov'],
                evalue=params['evalue']
            )
            
            # Merge alignments for this type
            type_genes_df = merge_alignments(type_genes_df, tmp_aln_path)
            all_genes_dfs.append(type_genes_df)
        
        # Step 3: Combine all type-specific results
        rm_genes_df = pd.concat(all_genes_dfs, ignore_index=True)
        
        # Drop the system_type column added during merge (keep original structure)
        rm_genes_df.drop(columns=['system_type'], inplace=True, errors='ignore')
        
        # Step 4: Propagate recseq to rm_systems
        rm_systems_df = propagate_recseq_to_systems(rm_genes_df, rm_systems_df)
        
        # Step 5: Write outputs
        rm_genes_df.to_csv(rm_genes_out, sep='\t', index=False, na_rep='NULL')
        rm_systems_df.to_csv(rm_systems_out, sep='\t', index=False, na_rep='NULL')
        
    finally:
        # Cleanup temp files
        for temp_file in temp_files:
            Path(temp_file).unlink(missing_ok=True)
        # Remove temp directory
        import shutil
        shutil.rmtree(temp_dir, ignore_errors=True)


if __name__ == '__main__':
    main(snakemake)

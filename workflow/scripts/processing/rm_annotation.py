#!/usr/bin/env python3
"""
RM System REBASE Annotation

Annotates RM systems and genes with REBASE hits including recognition sequences.

Recognition Sequence Logic by RM Type:
    - Type I: S subunit determines recseq. Multiple S genes with different sequences → conflicting.
    - Type II / Type III: Requires R+M pair with matching recseq. Different sequences → conflicting.
    - Type IIG: All bifunctional RM genes must share same recseq. Different → conflicting.
    - Type IV: No recognition sequence. Always `no_recognition_site`.

recseq_note Values:
    - annotated: Required genes have matching REBASE recognition sequences
    - conflicting: Required genes have recognition sites but they don't match
    - missing_annotation: Required genes are missing or lack REBASE annotations
    - no_recognition_site: Type IV system
"""
import re
import pandas as pd
from pathlib import Path


def safe_read(path):
    """Read TSV file safely, returning empty DataFrame if missing."""
    if path and Path(path).exists() and Path(path).stat().st_size > 0:
        return pd.read_csv(path, sep='\t')
    return pd.DataFrame()


def parse_rebase_stitle(stitle):
    """Parse REBASE stitle to extract recognition sequence.
    
    REBASE stitle format contains tab-separated key:value pairs like:
    REBASE:Msa17ORFC2P  EnzType:putative Type I  RecSeq:AACNNNNNNGTGC  GenBank:XXX
    """
    result = {'enzyme': None, 'recseq': None}
    
    if not stitle or pd.isna(stitle):
        return result
    
    stitle = str(stitle)
    
    # Extract enzyme name from REBASE:XXX
    enzyme_match = re.search(r'REBASE:([^\t]+)', stitle)
    if enzyme_match:
        result['enzyme'] = enzyme_match.group(1).strip()
    
    # Extract recognition sequence from RecSeq:XXX
    recseq_match = re.search(r'RecSeq:([ACGTNMRWSYKVHDB/^]+)', stitle, re.IGNORECASE)
    if recseq_match:
        result['recseq'] = recseq_match.group(1).upper()
    
    return result


def get_gene_subunit(gene_name):
    """Classify gene subunit from gene_name."""
    if not gene_name or pd.isna(gene_name):
        return None
    gene_name = str(gene_name).upper()
    if gene_name in ('S', 'R', 'M', 'RM'):
        return gene_name
    return None


def load_rebase_hits(path):
    """Load REBASE hits from Diamond output.
    
    Diamond outfmt 6 with stitle containing extra tab-separated fields.
    Parse manually to handle variable number of stitle fields.
    """
    rows = []
    if not Path(path).exists() or Path(path).stat().st_size == 0:
        return pd.DataFrame()
    
    with open(path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 12:
                row = {
                    'qseqid': parts[0],
                    'sseqid': parts[1],
                    'pident': float(parts[2]),
                    'length': int(parts[3]),
                    'mismatch': int(parts[4]),
                    'gapopen': int(parts[5]),
                    'qstart': int(parts[6]),
                    'qend': int(parts[7]),
                    'sstart': int(parts[8]),
                    'send': int(parts[9]),
                    'evalue': float(parts[10]),
                    'bitscore': float(parts[11]),
                    'stitle': '\t'.join(parts[12:])  # Join remaining fields
                }
                rows.append(row)
    
    return pd.DataFrame(rows)


def annotate_genes_with_rebase(rm_genes_df, rebase_df):
    """Add REBASE annotations to RM genes using prodigal_id column."""
    if rm_genes_df.empty:
        return rm_genes_df
    
    rm_genes_df = rm_genes_df.copy()
    rm_genes_df['REBASE_hit'] = None
    rm_genes_df['REBASE_pident'] = None
    rm_genes_df['REBASE_evalue'] = None
    rm_genes_df['REBASE_recseq'] = None
    
    if rebase_df.empty or 'prodigal_id' not in rm_genes_df.columns:
        return rm_genes_df
    
    for prodigal_id in rm_genes_df['prodigal_id'].dropna().unique():
        hits = rebase_df[rebase_df['qseqid'] == prodigal_id]
        if not hits.empty:
            best = hits.loc[hits['pident'].idxmax()]
            rebase_info = parse_rebase_stitle(best['stitle'])
            
            mask = rm_genes_df['prodigal_id'] == prodigal_id
            rm_genes_df.loc[mask, 'REBASE_hit'] = best['sseqid']
            rm_genes_df.loc[mask, 'REBASE_pident'] = best['pident']
            rm_genes_df.loc[mask, 'REBASE_evalue'] = best['evalue']
            rm_genes_df.loc[mask, 'REBASE_recseq'] = rebase_info['recseq']
    
    return rm_genes_df


def annotate_system_recseq(sys_row, sys_genes_df):
    """Annotate system with recognition sequence based on RM type logic.
    
    Returns: (recseq, recseq_note)
    """
    sys_type = str(sys_row.get('type', '')).strip()
    
    # Type IV: No recognition sequence
    if 'IV' in sys_type:
        return (None, 'no_recognition_site')
    
    if sys_genes_df.empty:
        return (None, 'missing_annotation')
    
    # Get subunit info
    sys_genes_df = sys_genes_df.copy()
    sys_genes_df['subunit'] = sys_genes_df['gene_name'].apply(get_gene_subunit)
    
    # Type I: S subunit determines recseq
    if 'Type I' in sys_type and 'IIG' not in sys_type:
        s_genes = sys_genes_df[sys_genes_df['subunit'] == 'S']
        if s_genes.empty:
            return (None, 'missing_annotation')
        
        s_recseqs = s_genes['REBASE_recseq'].dropna().unique()
        if len(s_recseqs) == 0:
            return (None, 'missing_annotation')
        elif len(s_recseqs) == 1:
            return (s_recseqs[0], 'annotated')
        else:
            return (None, 'conflicting')
    
    # Type IIG: All bifunctional RM genes must match
    if 'IIG' in sys_type:
        rm_recseqs = sys_genes_df['REBASE_recseq'].dropna().unique()
        if len(rm_recseqs) == 0:
            return (None, 'missing_annotation')
        elif len(rm_recseqs) == 1:
            return (rm_recseqs[0], 'annotated')
        else:
            return (None, 'conflicting')
    
    # Type II and Type III: R+M pairs with matching recseq
    if 'II' in sys_type or 'III' in sys_type:
        r_genes = sys_genes_df[sys_genes_df['subunit'] == 'R']
        m_genes = sys_genes_df[sys_genes_df['subunit'] == 'M']
        
        r_recseqs = set(r_genes['REBASE_recseq'].dropna())
        m_recseqs = set(m_genes['REBASE_recseq'].dropna())
        
        # Find matching recseqs
        matching = r_recseqs & m_recseqs
        
        if len(matching) == 0:
            if r_recseqs and m_recseqs:
                return (None, 'conflicting')
            all_recseqs = r_recseqs | m_recseqs
            if len(all_recseqs) == 0:
                return (None, 'missing_annotation')
            elif len(all_recseqs) == 1:
                return (list(all_recseqs)[0], 'annotated')
            else:
                return (None, 'conflicting')
        elif len(matching) == 1:
            return (list(matching)[0], 'annotated')
        else:
            return (None, 'conflicting')
    
    # Default: check all genes
    all_recseqs = sys_genes_df['REBASE_recseq'].dropna().unique()
    if len(all_recseqs) == 0:
        return (None, 'missing_annotation')
    elif len(all_recseqs) == 1:
        return (all_recseqs[0], 'annotated')
    else:
        return (None, 'conflicting')


def annotate_systems(rm_systems_df, rm_genes_df):
    """Annotate all systems with recognition sequences."""
    if rm_systems_df.empty:
        return rm_systems_df
    
    rm_systems_df = rm_systems_df.copy()
    rm_systems_df['recseq'] = None
    rm_systems_df['recseq_note'] = None
    
    for idx, sys_row in rm_systems_df.iterrows():
        sys_id = sys_row.get('ID')
        if pd.isna(sys_id):
            continue
        
        # Get genes for this system
        if not rm_genes_df.empty and 'Parent' in rm_genes_df.columns:
            sys_genes = rm_genes_df[rm_genes_df['Parent'] == sys_id]
        else:
            sys_genes = pd.DataFrame()
        
        recseq, recseq_note = annotate_system_recseq(sys_row, sys_genes)
        rm_systems_df.at[idx, 'recseq'] = recseq
        rm_systems_df.at[idx, 'recseq_note'] = recseq_note
    
    return rm_systems_df


def main(snakemake):
    """Main entry point for Snakemake."""
    # Load data
    rm_systems_df = safe_read(snakemake.input.rm_systems)
    rm_genes_df = safe_read(snakemake.input.rm_genes)
    rebase_df = load_rebase_hits(snakemake.input.rebase_hits)
    
    print(f"Loaded: {len(rm_systems_df)} RM systems, {len(rm_genes_df)} genes, {len(rebase_df)} REBASE hits")
    
    # Annotate genes with REBASE
    rm_genes_df = annotate_genes_with_rebase(rm_genes_df, rebase_df)
    annotated_count = (rm_genes_df['REBASE_recseq'].notna()).sum() if not rm_genes_df.empty else 0
    print(f"Annotated {annotated_count} genes with REBASE recseq")
    
    # Annotate systems
    rm_systems_df = annotate_systems(rm_systems_df, rm_genes_df)
    
    # Summary
    if not rm_systems_df.empty:
        print("\nSystem annotation summary:")
        print(rm_systems_df['recseq_note'].value_counts())
    
    # Save outputs
    Path(snakemake.output.rm_systems).parent.mkdir(parents=True, exist_ok=True)
    rm_systems_df.to_csv(snakemake.output.rm_systems, sep='\t', index=False, na_rep='NULL')
    rm_genes_df.to_csv(snakemake.output.rm_genes, sep='\t', index=False, na_rep='NULL')
    
    print(f"\nSaved to: {snakemake.output.rm_systems}, {snakemake.output.rm_genes}")


if __name__ == '__main__':
    main(snakemake)  # noqa: F821

#!/usr/bin/env python3
"""
Build summary matrices for target mapping results.

Produces:
- crispr_vs_mge_matrix.tsv: CRISPR spacer patterns vs unique MGE IDs
  - Rows: "{pattern_id} - {array_ids}" for spacers
  - Columns: Unique MGE IDs (e.g., MGYG000001345__MGE0011)
  - Cells: Hit counts

- rm_vs_mge_matrix.tsv: RM systems vs unique MGE IDs
  - Rows: "{rm_type} - {recseq}" for RM systems
  - Columns: Unique MGE IDs
  - Cells: Hit counts

- crispr_vs_genome_matrix.tsv: CRISPR spacer patterns vs target genomes
  - Rows: "{pattern_id} - {array_ids}" for spacers
  - Columns: Target genome IDs
  - Cells: Hit counts

- rm_vs_genome_matrix.tsv: RM systems vs target genomes
  - Rows: "{rm_type} - {recseq}" for RM systems
  - Columns: Target genome IDs
  - Cells: Hit counts
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple


def parse_gff_attributes(attrs_str: str) -> Dict[str, str]:
    """Parse GFF3 attributes string into dict."""
    attrs = {}
    for pair in attrs_str.split(';'):
        if '=' in pair:
            key, value = pair.split('=', 1)
            # URL-decode
            value = value.replace('%2C', ',').replace('%3D', '=').replace('%3B', ';').replace('%25', '%')
            attrs[key] = value
    return attrs


def load_mobile_hits(genomes: List[str], output_dir: str) -> List[Dict]:
    """
    Load all mobile hits from per-genome GFF files.
    
    Returns:
        List of hit dicts with target_genome, feature_type, mge_id, and defense identifier info
    """
    all_hits = []
    
    for genome_id in genomes:
        mobile_gff = Path(output_dir) / genome_id / "targets" / "mobile_hits.gff"
        
        if not mobile_gff.exists():
            continue
        
        with open(mobile_gff, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                feature_type = fields[2]
                attrs = parse_gff_attributes(fields[8])
                
                mge_id = attrs.get('within_mge', '')
                if not mge_id:
                    continue  # Should not happen for mobile_hits.gff
                
                hit = {
                    'target_genome': genome_id,
                    'feature_type': feature_type,
                    'mge_id': mge_id,
                    'pattern_id': attrs.get('pattern_id', ''),
                    # For RM systems
                    'rm_type': attrs.get('rm_type', ''),
                    'recseq': attrs.get('recseq', ''),
                    'source_systems': attrs.get('source_systems', ''),
                    # For spacers
                    'cas_type': attrs.get('cas_type', ''),
                    'array_id': attrs.get('array_id', ''),
                    'source_genomes': attrs.get('source_genomes', '')
                }
                all_hits.append(hit)
    
    return all_hits


def get_defense_identifier(hit: Dict) -> str:
    """
    Get unique defense system identifier for a hit.
    
    For RM: "{rm_type} - {recseq}" (e.g., "Type II - GATC")
    For spacers: "{pattern_id} - {array_ids}" where array_ids is comma-separated
                 (e.g., "spacer_878134ac - MGYG000001547__CRA0004")
    
    Returns:
        Unique defense identifier string
    """
    if hit['feature_type'] == 'restriction_site':
        rm_type = hit.get('rm_type', 'unknown')
        recseq = hit.get('recseq', 'unknown')
        if rm_type and recseq:
            return f"{rm_type} - {recseq}"
        elif rm_type:
            return rm_type
        else:
            return 'unknown_rm'
    
    elif hit['feature_type'] == 'protospacer':
        pattern_id = hit.get('pattern_id', '')
        array_id = hit.get('array_id', '')  # Already comma-separated if multiple arrays
        
        if pattern_id and array_id:
            return f"{pattern_id} - {array_id}"
        elif pattern_id:
            return pattern_id
        elif array_id:
            return f"unknown - {array_id}"
        else:
            return 'unknown_spacer'
    
    return 'unknown'


def build_defense_vs_mge_matrix(
    hits: List[Dict],
    output_tsv: str,
    feature_type_filter: str = None,
    matrix_name: str = "defense"
) -> None:
    """
    Build matrix of unique defense systems vs unique MGE IDs.
    
    Rows: Unique defense identifiers (rm_type-recseq or spacer-array_id)
    Columns: Unique MGE IDs
    Cells: Hit counts
    
    Args:
        hits: List of hit dicts
        output_tsv: Output TSV path
        feature_type_filter: If provided, only include hits of this feature type
                            ('protospacer' or 'restriction_site')
        matrix_name: Name for logging
    """
    # Count hits by defense identifier and MGE ID
    counts = defaultdict(lambda: defaultdict(int))
    
    for hit in hits:
        # Filter by feature type if specified
        if feature_type_filter and hit.get('feature_type') != feature_type_filter:
            continue
        
        mge_id = hit.get('mge_id', '')
        if not mge_id:
            continue
        
        defense_id = get_defense_identifier(hit)
        counts[defense_id][mge_id] += 1
    
    if not counts:
        Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
        with open(output_tsv, 'w') as f:
            f.write("defense_system\n")
        print(f"No mobile hits found for {matrix_name} vs MGE matrix")
        return
    
    # Get all unique MGE IDs (sorted)
    all_mge_ids = sorted(set(
        mge_id for defense_counts in counts.values() 
        for mge_id in defense_counts
    ))
    
    # Write matrix
    Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(output_tsv, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        
        # Header: defense_system + all MGE IDs
        writer.writerow(['defense_system'] + all_mge_ids)
        
        # Rows: sorted defense identifiers
        for defense_id in sorted(counts.keys()):
            row = [defense_id]
            for mge_id in all_mge_ids:
                row.append(counts[defense_id][mge_id])
            writer.writerow(row)
    
    print(f"Built {matrix_name} vs MGE matrix: {len(counts)} systems x {len(all_mge_ids)} MGEs")


def build_defense_vs_genome_matrix(
    hits: List[Dict],
    all_genomes: List[str],
    output_tsv: str,
    feature_type_filter: str = None,
    matrix_name: str = "defense"
) -> None:
    """
    Build matrix of unique defense systems vs target genomes.
    
    Rows: Unique defense identifiers (rm_type-recseq or spacer-array_id)
    Columns: Target genome IDs
    Cells: Hit counts
    
    Args:
        hits: List of hit dicts
        all_genomes: List of all genome IDs for consistent columns
        output_tsv: Output TSV path
        feature_type_filter: If provided, only include hits of this feature type
                            ('protospacer' or 'restriction_site')
        matrix_name: Name for logging
    """
    # Count hits by defense identifier and target genome
    counts = defaultdict(lambda: defaultdict(int))
    
    for hit in hits:
        # Filter by feature type if specified
        if feature_type_filter and hit.get('feature_type') != feature_type_filter:
            continue
        
        target_genome = hit.get('target_genome', '')
        if not target_genome:
            continue
        
        defense_id = get_defense_identifier(hit)
        counts[defense_id][target_genome] += 1
    
    if not counts:
        Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
        with open(output_tsv, 'w') as f:
            f.write("defense_system\n")
        print(f"No hits found for {matrix_name} vs genome matrix")
        return
    
    # Use all genomes for consistent columns (sorted)
    genome_list = sorted(all_genomes)
    
    # Write matrix
    Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(output_tsv, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        
        # Header: defense_system + all genome IDs
        writer.writerow(['defense_system'] + genome_list)
        
        # Rows: sorted defense identifiers
        for defense_id in sorted(counts.keys()):
            row = [defense_id]
            for genome_id in genome_list:
                row.append(counts[defense_id][genome_id])
            writer.writerow(row)
    
    print(f"Built {matrix_name} vs genome matrix: {len(counts)} systems x {len(genome_list)} genomes")


def build_matrices(
    genomes: List[str],
    output_dir: str,
    crispr_vs_mge_output: str,
    rm_vs_mge_output: str,
    crispr_vs_genome_output: str,
    rm_vs_genome_output: str
) -> None:
    """
    Build all summary matrices from mobile hits GFFs.
    """
    # Load all mobile hits
    hits = load_mobile_hits(genomes, output_dir)
    
    print(f"Loaded {len(hits)} mobile hits from {len(genomes)} genomes")
    
    # Count by feature type
    spacer_hits = sum(1 for h in hits if h['feature_type'] == 'protospacer')
    rm_hits = sum(1 for h in hits if h['feature_type'] == 'restriction_site')
    print(f"  - Spacer (protospacer) hits: {spacer_hits}")
    print(f"  - RM (restriction_site) hits: {rm_hits}")
    
    # Build separate matrices for CRISPR and RM vs MGE
    build_defense_vs_mge_matrix(
        hits, crispr_vs_mge_output,
        feature_type_filter='protospacer',
        matrix_name='CRISPR'
    )
    
    build_defense_vs_mge_matrix(
        hits, rm_vs_mge_output,
        feature_type_filter='restriction_site',
        matrix_name='RM'
    )
    
    # Build separate matrices for CRISPR and RM vs genome
    build_defense_vs_genome_matrix(
        hits, genomes, crispr_vs_genome_output,
        feature_type_filter='protospacer',
        matrix_name='CRISPR'
    )
    
    build_defense_vs_genome_matrix(
        hits, genomes, rm_vs_genome_output,
        feature_type_filter='restriction_site',
        matrix_name='RM'
    )


def main():
    parser = argparse.ArgumentParser(
        description='Build target mapping summary matrices'
    )
    parser.add_argument(
        '--genomes', required=True, nargs='+',
        help='List of genome IDs'
    )
    parser.add_argument(
        '--output-dir', required=True,
        help='Output directory containing per-genome targets'
    )
    parser.add_argument(
        '--crispr-vs-mge-output', required=True,
        help='Output path for CRISPR vs MGE matrix'
    )
    parser.add_argument(
        '--rm-vs-mge-output', required=True,
        help='Output path for RM vs MGE matrix'
    )
    parser.add_argument(
        '--crispr-vs-genome-output', required=True,
        help='Output path for CRISPR vs genome matrix'
    )
    parser.add_argument(
        '--rm-vs-genome-output', required=True,
        help='Output path for RM vs genome matrix'
    )
    
    args = parser.parse_args()
    
    build_matrices(
        args.genomes,
        args.output_dir,
        args.crispr_vs_mge_output,
        args.rm_vs_mge_output,
        args.crispr_vs_genome_output,
        args.rm_vs_genome_output
    )


if __name__ == '__main__':
    main()

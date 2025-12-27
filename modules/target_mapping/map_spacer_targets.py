#!/usr/bin/env python3
"""
Map CRISPR spacer patterns to target genomes using BLAST.

Outputs GFF3 directly with hits classified as mobile (within MGE) or genomic.
Suspicious hits (spacers hitting within their own array) are written separately.
"""

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path
from typing import List, Dict, Tuple, Optional


def parse_fasta(fasta_path: str) -> Dict[str, str]:
    """
    Parse a FASTA file into a dictionary of sequences.
    
    Args:
        fasta_path: Path to the FASTA file
        
    Returns:
        dict: sequence_id -> sequence (uppercase)
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq).upper()
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq).upper()
    
    return sequences


def load_spacer_patterns(patterns_tsv: str) -> Dict[str, Dict]:
    """
    Load spacer patterns from metadata TSV file.
    
    Args:
        patterns_tsv: Path to patterns_spacers.tsv
        
    Returns:
        Dict: pattern_id -> metadata dict
    """
    patterns = {}
    
    if not Path(patterns_tsv).exists():
        return patterns
    
    with open(patterns_tsv, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # Parse source coordinates for self-hit depletion
            # Format: "genome:contig:start-end,genome:contig:start-end,..."
            source_coords = []
            if row.get('source_coords'):
                for coord_str in row['source_coords'].split(','):
                    if ':' in coord_str and '-' in coord_str:
                        try:
                            parts = coord_str.split(':')
                            genome_id = parts[0]
                            contig = parts[1]
                            start_end = parts[2].split('-')
                            start = int(start_end[0])
                            end = int(start_end[1])
                            source_coords.append((genome_id, contig, start, end))
                        except (ValueError, IndexError):
                            continue
            
            # Parse source array IDs for suspicious hit detection
            source_array_ids = []
            if row.get('source_array_ids'):
                source_array_ids = [a.strip() for a in row['source_array_ids'].split(',')]
            
            patterns[row['pattern_id']] = {
                'pattern_id': row['pattern_id'],
                'sequence': row['sequence'].upper(),
                'cas_type': row.get('source_cas_types', ''),
                'array_id': row.get('source_array_ids', ''),
                'source_genome_ids': row.get('source_genome_ids', ''),
                'source_coords': source_coords,
                'source_array_ids': source_array_ids
            }
    
    return patterns


def load_mge_bounds(mges_tsv: str, genome_id: str) -> Dict[str, List[Tuple[int, int, str]]]:
    """
    Load MGE bounds for a genome, organized by contig.
    
    Args:
        mges_tsv: Path to all_mges.tsv
        genome_id: Genome ID to filter for
        
    Returns:
        Dict: contig -> list of (start, end, mge_id) tuples
    """
    bounds = {}
    
    if not Path(mges_tsv).exists():
        return bounds
    
    with open(mges_tsv, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['genome_id'] != genome_id:
                continue
            
            contig = row['contig']
            start = int(row['start'])
            end = int(row['end'])
            mge_id = row['ID']
            
            if contig not in bounds:
                bounds[contig] = []
            bounds[contig].append((start, end, mge_id))
    
    # Sort by start position
    for contig in bounds:
        bounds[contig].sort(key=lambda x: x[0])
    
    return bounds


def load_crispr_array_bounds(arrays_tsv: str, genome_id: str) -> Dict[str, Tuple[str, int, int]]:
    """
    Load CRISPR array bounds for a genome, keyed by array ID.
    
    Args:
        arrays_tsv: Path to all_crispr_arrays.tsv
        genome_id: Genome ID to filter for
        
    Returns:
        Dict: array_id -> (contig, start, end)
    """
    bounds = {}
    
    if not Path(arrays_tsv).exists():
        return bounds
    
    with open(arrays_tsv, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['genome_id'] != genome_id:
                continue
            
            array_id = row['ID']
            contig = row['contig']
            start = int(row['start'])
            end = int(row['end'])
            
            bounds[array_id] = (contig, start, end)
    
    return bounds


def find_containing_mge(
    contig: str,
    start: int,
    end: int,
    mge_bounds: Dict[str, List[Tuple[int, int, str]]]
) -> str:
    """
    Check if a hit is 100% contained within any MGE.
    
    Returns:
        MGE ID if contained, empty string otherwise
    """
    if contig not in mge_bounds:
        return ''
    
    for mge_start, mge_end, mge_id in mge_bounds[contig]:
        if start >= mge_start and end <= mge_end:
            return mge_id
    
    return ''


def is_suspicious_array_hit(
    hit_genome: str,
    hit_contig: str,
    hit_start: int,
    hit_end: int,
    source_array_ids: List[str],
    array_bounds: Dict[str, Tuple[str, int, int]]
) -> Optional[str]:
    """
    Check if a hit lands within any of the source arrays for this spacer.
    
    This detects potential false positives from tandem repeats being
    miscalled as CRISPR arrays - a spacer shouldn't have protospacer
    hits within its own array.
    
    Args:
        hit_genome: Genome ID where hit was found
        hit_contig: Contig where hit was found
        hit_start: Hit start position
        hit_end: Hit end position
        source_array_ids: List of array IDs this spacer belongs to
        array_bounds: Dict of array_id -> (contig, start, end) for this genome
    
    Returns:
        Array ID if hit is within a source array, None otherwise
    """
    for array_id in source_array_ids:
        if array_id not in array_bounds:
            continue
        
        arr_contig, arr_start, arr_end = array_bounds[array_id]
        
        # Must be same contig
        if hit_contig != arr_contig:
            continue
        
        # Check if hit is within array bounds
        if hit_start >= arr_start and hit_end <= arr_end:
            return array_id
    
    return None


def is_self_hit(
    hit_genome: str,
    hit_contig: str,
    hit_start: int,
    hit_end: int,
    source_coords: List[Tuple[str, str, int, int]]
) -> bool:
    """
    Check if a hit overlaps with any source spacer location.
    
    This filters out protospacer hits that map back to the original
    spacer sequence in a CRISPR array.
    
    Args:
        hit_genome: Genome ID where hit was found
        hit_contig: Contig where hit was found
        hit_start: Hit start position
        hit_end: Hit end position
        source_coords: List of (genome_id, contig, start, end) tuples
                       representing original spacer locations
    
    Returns:
        True if hit overlaps with any source spacer, False otherwise
    """
    for src_genome, src_contig, src_start, src_end in source_coords:
        # Must be same genome and contig
        if hit_genome != src_genome or hit_contig != src_contig:
            continue
        
        # Check for overlap (any overlap counts as self-hit)
        if hit_start <= src_end and hit_end >= src_start:
            return True
    
    return False


def run_blastn(
    query_fasta: str,
    subject_fasta: str,
    min_identity: float = 90.0,
    min_coverage: float = 100.0,
    max_mismatches: int = 2,
    evalue: float = 10,
    num_threads: int = 4
) -> List[Dict]:
    """
    Run blastn for spacer mapping and return parsed hits.
    
    Args:
        query_fasta: Path to spacer patterns FASTA
        subject_fasta: Path to target genome FASTA
        min_identity: Minimum percent identity of aligned portion
        min_coverage: Minimum query coverage (percentage of spacer aligned, default 100%)
        max_mismatches: Maximum allowed mismatches
        evalue: E-value threshold
        num_threads: Number of threads for BLAST
        
    Returns:
        List of hit dictionaries
    """
    # Check inputs
    if not Path(query_fasta).exists():
        print(f"Query file not found: {query_fasta}")
        return []
    
    if not Path(subject_fasta).exists() or Path(subject_fasta).stat().st_size == 0:
        print(f"Subject file not found or empty: {subject_fasta}")
        return []
    
    # blastn output format - include alignment length and query length
    outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
    
    cmd = [
        'blastn',
        '-task', 'blastn-short',
        '-query', query_fasta,
        '-subject', subject_fasta,
        '-outfmt', outfmt,
        '-perc_identity', str(min_identity),
        '-qcov_hsp_perc', str(min_coverage),  # Enforce query coverage at BLAST level
        '-evalue', str(evalue),
        '-num_threads', str(num_threads),
        '-dust', 'no',
        '-word_size', '7',
        '-ungapped',
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        raw_output = result.stdout
    except subprocess.CalledProcessError as e:
        print(f"BLAST failed: {e.stderr}")
        return []
    except FileNotFoundError:
        print("blastn not found. Please ensure BLAST+ is installed.")
        return []
    
    # Parse BLAST output
    # Note: BLAST's -qcov_hsp_perc already filters by coverage, so we only
    # need to apply mismatch filter here
    hits = []
    
    for line in raw_output.strip().split('\n'):
        if not line:
            continue
        
        fields = line.split('\t')
        if len(fields) < 14:
            continue
        
        qseqid = fields[0]
        sseqid = fields[1]
        pident = float(fields[2])  # Identity of aligned portion
        aln_length = int(fields[3])  # Alignment length
        mismatch = int(fields[4])
        sstart = int(fields[8])
        send = int(fields[9])
        evalue_hit = float(fields[10])
        qlen = int(fields[12])   # Full query (spacer) length
        
        # Calculate query coverage for reporting
        query_coverage = (aln_length / qlen) * 100 if qlen > 0 else 0
        
        # Filter by mismatches
        if mismatch > max_mismatches:
            continue
        
        # Determine strand and coordinates
        if sstart < send:
            strand = '+'
            start = sstart
            end = send
        else:
            strand = '-'
            start = send
            end = sstart
        
        hits.append({
            'pattern_id': qseqid,
            'contig': sseqid,
            'start': start,
            'end': end,
            'strand': strand,
            'identity': round(pident, 3),  # Identity of aligned portion
            'coverage': round(query_coverage, 1),  # Fraction of spacer aligned
            'aln_length': aln_length,
            'qlen': qlen,
            'mismatches': mismatch,
            'evalue': evalue_hit
        })
    
    return hits


def format_gff_attributes(attrs: Dict) -> str:
    """Format attributes as GFF3 attribute string."""
    parts = []
    for key, value in attrs.items():
        if value:
            value = str(value).replace('%', '%25').replace(';', '%3B').replace('=', '%3D').replace(',', '%2C')
            parts.append(f"{key}={value}")
    return ';'.join(parts)


def write_gff(
    hits: List[Dict],
    pattern_metadata: Dict[str, Dict],
    genome_id: str,
    output_gff: str,
    feature_type: str = "protospacer"
) -> None:
    """
    Write hits to GFF3 file with metadata.
    
    Args:
        hits: List of hit dicts with pattern_id, contig, start, end, strand, mge_id
        pattern_metadata: Dict of pattern_id -> metadata
        genome_id: Target genome ID
        output_gff: Output GFF path
        feature_type: GFF feature type (default: protospacer)
    """
    Path(output_gff).parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_gff, 'w') as f:
        f.write("##gff-version 3\n")
        
        for idx, hit in enumerate(hits, 1):
            pattern_id = hit['pattern_id']
            meta = pattern_metadata.get(pattern_id, {})
            
            attrs = {
                'ID': f"{genome_id}_spacer_hit_{idx:06d}",
                'pattern_id': pattern_id,
                'cas_type': meta.get('cas_type', ''),
                'array_id': meta.get('array_id', ''),
                'source_genomes': meta.get('source_genome_ids', ''),
                'identity': hit.get('identity', ''),
                'coverage': hit.get('coverage', ''),
                'aln_length': hit.get('aln_length', ''),
                'spacer_length': hit.get('qlen', ''),
                'mismatches': hit.get('mismatches', ''),
                'evalue': hit.get('evalue', '')
            }
            
            if hit.get('mge_id'):
                attrs['within_mge'] = hit['mge_id']
            
            # For suspicious hits, note which array the hit is within
            if hit.get('suspicious_array'):
                attrs['suspicious_array'] = hit['suspicious_array']
            
            attrs_str = format_gff_attributes(attrs)
            
            line = f"{hit['contig']}\ttarget_mapping\t{feature_type}\t{hit['start']}\t{hit['end']}\t.\t{hit['strand']}\t.\t{attrs_str}"
            f.write(line + '\n')


def map_spacer_targets(
    genome_id: str,
    patterns_fasta: str,
    patterns_tsv: str,
    genome_fasta: str,
    mges_tsv: str,
    arrays_tsv: str,
    output_gff: str,
    suspicious_gff: str,
    min_identity: float = 90.0,
    min_coverage: float = 80.0,
    max_mismatches: int = 2,
    num_threads: int = 4
) -> None:
    """
    Map spacer patterns to genome using BLAST and output GFF.
    
    Args:
        genome_id: Target genome ID
        patterns_fasta: Path to patterns_spacers.fa
        patterns_tsv: Path to patterns_spacers.tsv with metadata
        genome_fasta: Path to genome FASTA
        mges_tsv: Path to all_mges.tsv for mobility classification
        arrays_tsv: Path to all_crispr_arrays.tsv for suspicious hit detection
        output_gff: Output GFF path for valid hits
        suspicious_gff: Output GFF path for suspicious hits (within source array)
        min_identity: Minimum percent identity of aligned portion
        min_coverage: Minimum query coverage (percentage of spacer aligned)
        max_mismatches: Maximum allowed mismatches
        num_threads: Number of threads for BLAST
    """
    # Load pattern metadata
    pattern_metadata = load_spacer_patterns(patterns_tsv)
    if not pattern_metadata:
        Path(output_gff).parent.mkdir(parents=True, exist_ok=True)
        with open(output_gff, 'w') as f:
            f.write("##gff-version 3\n")
        with open(suspicious_gff, 'w') as f:
            f.write("##gff-version 3\n")
        print("No spacer patterns to search")
        return
    
    # Check genome FASTA
    if not Path(genome_fasta).exists():
        Path(output_gff).parent.mkdir(parents=True, exist_ok=True)
        with open(output_gff, 'w') as f:
            f.write("##gff-version 3\n")
        with open(suspicious_gff, 'w') as f:
            f.write("##gff-version 3\n")
        print(f"Genome FASTA not found: {genome_fasta}")
        return
    
    # Load MGE bounds for mobility classification
    mge_bounds = load_mge_bounds(mges_tsv, genome_id)
    
    # Load CRISPR array bounds for suspicious hit detection
    array_bounds = load_crispr_array_bounds(arrays_tsv, genome_id)
    
    # Run BLAST with coverage filter
    raw_hits = run_blastn(
        patterns_fasta,
        genome_fasta,
        min_identity=min_identity,
        min_coverage=min_coverage,
        max_mismatches=max_mismatches,
        num_threads=num_threads
    )
    
    # Filter out self-hits (protospacers mapping back to original spacer location)
    # Detect suspicious hits (hits within source array but not at source position)
    # Classify remaining hits as mobile or genomic
    valid_hits = []
    suspicious_hits = []
    self_hit_count = 0
    
    for hit in raw_hits:
        pattern_id = hit['pattern_id']
        meta = pattern_metadata.get(pattern_id, {})
        source_coords = meta.get('source_coords', [])
        source_array_ids = meta.get('source_array_ids', [])
        
        # Check if this hit maps back to the original spacer location (exact self-hit)
        if is_self_hit(genome_id, hit['contig'], hit['start'], hit['end'], source_coords):
            self_hit_count += 1
            continue
        
        # Check if this hit lands within a source array (suspicious - possible tandem repeat)
        suspicious_array = is_suspicious_array_hit(
            genome_id, hit['contig'], hit['start'], hit['end'],
            source_array_ids, array_bounds
        )
        
        if suspicious_array:
            hit['suspicious_array'] = suspicious_array
            # Also classify MGE for context
            mge_id = find_containing_mge(hit['contig'], hit['start'], hit['end'], mge_bounds)
            hit['mge_id'] = mge_id
            suspicious_hits.append(hit)
            continue
        
        # Valid hit - classify as mobile or genomic
        mge_id = find_containing_mge(hit['contig'], hit['start'], hit['end'], mge_bounds)
        hit['mge_id'] = mge_id
        valid_hits.append(hit)
    
    # Write valid hits GFF
    write_gff(valid_hits, pattern_metadata, genome_id, output_gff)
    
    # Write suspicious hits GFF
    write_gff(suspicious_hits, pattern_metadata, genome_id, suspicious_gff, 
              feature_type="suspicious_protospacer")
    
    # Count mobile vs genomic
    mobile_count = sum(1 for h in valid_hits if h['mge_id'])
    genomic_count = len(valid_hits) - mobile_count
    
    print(f"Found {len(raw_hits)} raw spacer hits for {genome_id}")
    print(f"  - Self-hits filtered: {self_hit_count}")
    print(f"  - Suspicious (within source array): {len(suspicious_hits)}")
    print(f"  - Valid hits: {len(valid_hits)}")
    print(f"  - Mobile (within MGE): {mobile_count}")
    print(f"  - Genomic (outside MGE): {genomic_count}")


def main():
    parser = argparse.ArgumentParser(
        description='Map CRISPR spacer patterns to genome and output GFF'
    )
    parser.add_argument(
        '--genome-id', required=True,
        help='Target genome ID'
    )
    parser.add_argument(
        '--patterns-fasta', required=True,
        help='Path to patterns_spacers.fa'
    )
    parser.add_argument(
        '--patterns-tsv', required=True,
        help='Path to patterns_spacers.tsv with metadata'
    )
    parser.add_argument(
        '--genome-fasta', required=True,
        help='Path to genome FASTA'
    )
    parser.add_argument(
        '--mges-tsv', required=True,
        help='Path to all_mges.tsv for mobility classification'
    )
    parser.add_argument(
        '--arrays-tsv', required=True,
        help='Path to all_crispr_arrays.tsv for suspicious hit detection'
    )
    parser.add_argument(
        '--output-gff', required=True,
        help='Output GFF path for valid hits'
    )
    parser.add_argument(
        '--suspicious-gff', required=True,
        help='Output GFF path for suspicious hits (within source array)'
    )
    parser.add_argument(
        '--min-identity', type=float, default=90.0,
        help='Minimum percent identity of aligned portion (default: 90)'
    )
    parser.add_argument(
        '--min-coverage', type=float, default=100.0,
        help='Minimum query coverage - percentage of spacer aligned (default: 100)'
    )
    parser.add_argument(
        '--max-mismatches', type=int, default=2,
        help='Maximum allowed mismatches (default: 2)'
    )
    parser.add_argument(
        '--threads', type=int, default=4,
        help='Number of threads (default: 4)'
    )
    
    args = parser.parse_args()
    
    map_spacer_targets(
        args.genome_id,
        args.patterns_fasta,
        args.patterns_tsv,
        args.genome_fasta,
        args.mges_tsv,
        args.arrays_tsv,
        args.output_gff,
        args.suspicious_gff,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        max_mismatches=args.max_mismatches,
        num_threads=args.threads
    )


if __name__ == '__main__':
    main()

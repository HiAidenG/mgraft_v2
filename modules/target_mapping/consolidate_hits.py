#!/usr/bin/env python3
"""
Consolidate spacer and RM target mapping GFF outputs into merged files.

Merges the two input GFFs (spacer_hits.gff and rm_hits.gff) and produces:
- mobile_hits.gff: Hits within MGE bounds (have within_mge attribute)
- genomic_hits.gff: Hits outside MGE bounds
"""

import argparse
from pathlib import Path
from typing import List, Tuple


def parse_gff_line(line: str) -> Tuple[str, dict]:
    """
    Parse a GFF3 line into fields and attributes.
    
    Args:
        line: GFF3 line
        
    Returns:
        Tuple of (line, attributes dict)
    """
    if line.startswith('#') or not line.strip():
        return line, {}
    
    fields = line.strip().split('\t')
    if len(fields) < 9:
        return line, {}
    
    # Parse attributes
    attrs = {}
    for pair in fields[8].split(';'):
        if '=' in pair:
            key, value = pair.split('=', 1)
            # URL-decode
            value = value.replace('%2C', ',').replace('%3D', '=').replace('%3B', ';').replace('%25', '%')
            attrs[key] = value
    
    return line, attrs


def read_gff(gff_path: str) -> List[Tuple[str, dict]]:
    """
    Read a GFF file and return lines with parsed attributes.
    
    Args:
        gff_path: Path to GFF file
        
    Returns:
        List of (line, attributes) tuples
    """
    lines = []
    
    if not Path(gff_path).exists():
        return lines
    
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            if not line.strip():
                continue
            parsed_line, attrs = parse_gff_line(line)
            if attrs:  # Only include data lines
                lines.append((parsed_line, attrs))
    
    return lines


def consolidate_hits(
    spacer_gff: str,
    rm_gff: str,
    output_mobile_gff: str,
    output_genomic_gff: str
) -> None:
    """
    Merge spacer and RM GFF files into mobile and genomic outputs.
    
    Args:
        spacer_gff: Path to spacer hits GFF
        rm_gff: Path to RM hits GFF
        output_mobile_gff: Output path for mobile hits
        output_genomic_gff: Output path for genomic hits
    """
    # Create output directories
    for path in [output_mobile_gff, output_genomic_gff]:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
    
    # Read both GFFs
    spacer_lines = read_gff(spacer_gff)
    rm_lines = read_gff(rm_gff)
    
    # Combine all lines
    all_lines = spacer_lines + rm_lines
    
    # Separate mobile and genomic
    mobile_lines = []
    genomic_lines = []
    
    for line, attrs in all_lines:
        if attrs.get('within_mge'):
            mobile_lines.append(line)
        else:
            genomic_lines.append(line)
    
    # Write outputs
    gff_header = "##gff-version 3\n"
    
    with open(output_mobile_gff, 'w') as f:
        f.write(gff_header)
        for line in mobile_lines:
            f.write(line if line.endswith('\n') else line + '\n')
    
    with open(output_genomic_gff, 'w') as f:
        f.write(gff_header)
        for line in genomic_lines:
            f.write(line if line.endswith('\n') else line + '\n')
    
    print(f"Consolidated {len(all_lines)} hits")
    print(f"  - Mobile (within MGE): {len(mobile_lines)}")
    print(f"  - Genomic (outside MGE): {len(genomic_lines)}")


def main():
    parser = argparse.ArgumentParser(
        description='Consolidate spacer and RM GFF outputs'
    )
    parser.add_argument(
        '--spacer-gff', required=True,
        help='Path to spacer hits GFF'
    )
    parser.add_argument(
        '--rm-gff', required=True,
        help='Path to RM hits GFF'
    )
    parser.add_argument(
        '--output-mobile-gff', required=True,
        help='Output path for mobile hits GFF'
    )
    parser.add_argument(
        '--output-genomic-gff', required=True,
        help='Output path for genomic hits GFF'
    )
    
    args = parser.parse_args()
    
    consolidate_hits(
        args.spacer_gff,
        args.rm_gff,
        args.output_mobile_gff,
        args.output_genomic_gff
    )


if __name__ == '__main__':
    main()

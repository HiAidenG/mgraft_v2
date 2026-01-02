#!/usr/bin/env python3
import pandas as pd
import sys
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Split tabular data by genome_id")
    parser.add_argument("--input", required=True, help="Input TSV file with genome_id column")
    parser.add_argument("--output-dir", required=True, help="Base output directory pattern (e.g., results/{}/features)")
    parser.add_argument("--filename", required=True, help="Output filename (e.g., crispr_spacers.clustered.tsv)")
    args = parser.parse_args()

    if not Path(args.input).exists():
        print(f"Input {args.input} not found")
        sys.exit(1)

    # Read efficiently
    try:
        df = pd.read_csv(args.input, sep='\t')
    except Exception as e:
        print(f"Error reading {args.input}: {e}")
        sys.exit(1)
        
    if 'genome_id' not in df.columns:
        print("genome_id column missing")
        sys.exit(1)
        
    # Group and write
    for genome_id, group in df.groupby('genome_id'):
        # Construct path: output_dir.format(genome_id) / filename
        out_path = Path(args.output_dir.replace("{}", str(genome_id))) / args.filename
        out_path.parent.mkdir(parents=True, exist_ok=True)
        group.to_csv(out_path, sep='\t', index=False, na_rep='NULL')

if __name__ == "__main__":
    main()

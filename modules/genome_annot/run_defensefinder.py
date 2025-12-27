#!/usr/bin/env python3
"""
Run DefenseFinder and rename outputs to standard names.

DefenseFinder outputs files with the input filename prefix (e.g., MGYG000001345_defense_finder_genes.tsv).
We rename them to standard names (defense_finder_genes.tsv) for consistent downstream access.
"""
import subprocess
import shutil
from pathlib import Path


def main(snakemake): 
    """Main function for Snakemake script."""
    protein_fasta = snakemake.input.protein_fasta
    
    output_genes = Path(snakemake.output.defense_genes_tsv)
    output_systems = Path(snakemake.output.defense_systems_tsv)
    output_hmmer = Path(snakemake.output.defense_hmmer_tsv)
    
    # DefenseFinder writes directly to output directory
    output_dir = output_genes.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Clean any existing DefenseFinder outputs to avoid stale files
    for pattern in ["*_defense_finder_*.tsv", "defense_finder_*.tsv"]:
        for f in output_dir.glob(pattern):
            f.unlink()

    # Run DefenseFinder with -a flag to get all systems including antidefense
    cmd = [
        "defense-finder", "run",
        "-a",
        "--out-dir", str(output_dir),
        "--workers", str(snakemake.threads),
        protein_fasta
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("DefenseFinder STDOUT:", result.stdout)
        print("DefenseFinder STDERR:", result.stderr)
        result.check_returncode()
    
    # Find and rename output files from prefixed names to standard names
    genes_tsv = list(output_dir.glob("*_defense_finder_genes.tsv"))
    systems_tsv = list(output_dir.glob("*_defense_finder_systems.tsv"))
    hmmer_tsv = list(output_dir.glob("*_defense_finder_hmmer.tsv"))
    
    if not genes_tsv:
        # No output - create empty files
        output_genes.touch()
        output_systems.touch()
        output_hmmer.touch()
        return
    
    # Rename outputs to standard names (move, not copy)
    genes_tsv[0].rename(output_genes)
    
    if systems_tsv:
        systems_tsv[0].rename(output_systems)
    else:
        output_systems.touch()
        
    if hmmer_tsv:
        hmmer_tsv[0].rename(output_hmmer)
    else:
        output_hmmer.touch()


if __name__ == '__main__':
    main(snakemake)  # noqa: F821

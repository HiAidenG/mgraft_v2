#!/usr/bin/env python3
"""
Run MinCED to identify CRISPR arrays and parse output.
"""
import subprocess
from pathlib import Path


def run_minced(genome_fasta, output_gff, min_nr):
    """Run MinCED to detect CRISPR arrays."""
    cmd = [
        "minced",
        "-gffFull",
        "-spacers",
        "-minNR", str(min_nr),
        str(genome_fasta),
        str(output_gff)
    ]
    subprocess.run(cmd, check=True, capture_output=True)


def main(snakemake):  # noqa: F821
    """Main function for Snakemake script."""
    genome_fasta = snakemake.input.genome_fasta
    minced_gff = snakemake.output.minced_gff
    min_nr = snakemake.params.min_nr
    minced_spacers_fa = snakemake.output.minced_spacers_fa
    
    # Run MinCED to the expected output location
    run_minced(genome_fasta, minced_gff, min_nr)

    # Rename spacers output file
    spacers_output = str(minced_gff).rsplit('.', 1)[0] + "_spacers.fasta"
    if Path(spacers_output).exists():
        subprocess.run(["mv", spacers_output, minced_spacers_fa], check=True)
    else:
        # Touch empty file if no spacers found
        Path(minced_spacers_fa).touch()
    
    
if __name__ == '__main__':
    main(snakemake)  # noqa: F821

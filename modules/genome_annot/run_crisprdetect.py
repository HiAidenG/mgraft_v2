#!/usr/bin/env python3
"""
Run CRISPRDetect 3.0 to identify CRISPR arrays with strand information.

CRISPRDetect determines strand orientation via repeat degeneracy analysis
and RNA secondary structure prediction, providing reliable strand calls.
"""
import subprocess
import shutil
from pathlib import Path


def run_crisprdetect(
    genome_fasta: str,
    output_prefix: str,
    crisprdetect_path: str,
    quality_cutoff: int = 0,
    min_repeats: int = 3,
    threads: int = 4
):
    """
    Run CRISPRDetect3 to detect CRISPR arrays.
    
    Args:
        genome_fasta: Path to input genome FASTA
        output_prefix: Prefix for output files (will create .gff, .spacers.fa, etc.)
        crisprdetect_path: Path to CRISPRDetect3 executable
        quality_cutoff: Minimum array quality score (0 for max sensitivity)
        min_repeats: Minimum number of repeats required
        threads: Number of parallel threads
    
    Outputs created by CRISPRDetect:
        {output_prefix}.txt        - Main text output
        {output_prefix}.txt.gff    - GFF3 annotations with strand
        {output_prefix}.txt.spacers.fa - Spacer sequences
        {output_prefix}.txt.fp     - Filtered predictions (false positives)
    """
    cmd = [
        "perl",
        str(crisprdetect_path),
        "-f", str(genome_fasta),
        "-o", str(output_prefix),
        "-array_quality_score_cutoff", str(quality_cutoff),
        "-minimum_no_of_repeats", str(min_repeats),
        "-T", str(threads),
        # Skip cas gene annotation since DefenseFinder handles that
        # Note: -annotate_cas_genes is 0 by default
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # CRISPRDetect may write to stderr even on success
    if result.returncode != 0:
        raise RuntimeError(
            f"CRISPRDetect failed with return code {result.returncode}\n"
            f"stdout: {result.stdout}\n"
            f"stderr: {result.stderr}"
        )
    
    return result


def main(snakemake):  # noqa: F821
    """Main function for Snakemake script."""
    genome_fasta = snakemake.input.genome_fasta
    output_gff = snakemake.output.crisprdetect_gff
    output_spacers = snakemake.output.crisprdetect_spacers_fa
    output_txt = snakemake.output.crisprdetect_txt
    
    # Get config params
    quality_cutoff = snakemake.params.get('quality_cutoff', 0)
    min_repeats = snakemake.params.get('min_repeats', 3)
    threads = snakemake.threads
    
    # CRISPRDetect path (relative to workflow dir)
    crisprdetect_path = snakemake.params.crisprdetect_path
    
    # Output directory setup
    output_dir = Path(output_txt).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # CRISPRDetect expects output prefix, adds extensions itself
    # We use output_txt as the base, CRISPRDetect will create:
    #   {output_txt}.gff
    #   {output_txt}.spacers.fa
    output_prefix = str(output_txt)
    
    # Run CRISPRDetect
    run_crisprdetect(
        genome_fasta=genome_fasta,
        output_prefix=output_prefix,
        crisprdetect_path=crisprdetect_path,
        quality_cutoff=quality_cutoff,
        min_repeats=min_repeats,
        threads=threads
    )
    
    # CRISPRDetect creates files with .txt extension appended
    # {prefix}.gff, {prefix}.spacers.fa
    cd_gff = Path(f"{output_prefix}.gff")
    cd_spacers = Path(f"{output_prefix}.spacers.fa")
    
    # Move/rename to expected output locations
    if cd_gff.exists():
        shutil.move(str(cd_gff), str(output_gff))
    else:
        # Create empty GFF if no CRISPRs found
        Path(output_gff).touch()
    
    if cd_spacers.exists():
        shutil.move(str(cd_spacers), str(output_spacers))
    else:
        # Create empty spacers file if none found
        Path(output_spacers).touch()
    
    # Ensure txt output exists (even if empty)
    if not Path(output_txt).exists():
        Path(output_txt).touch()


if __name__ == '__main__':
    main(snakemake)  # noqa: F821

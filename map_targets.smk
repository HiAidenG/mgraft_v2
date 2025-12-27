"""
Target Mapping Workflow - Snakefile_targets

Independent workflow that maps defense system patterns (CRISPR spacers, RM sites)
against input genomes to identify potential targets.

Requires summary outputs from the main mgraft_v2 workflow:
- all_crispr_spacers.tsv
- all_rm_systems.tsv
- all_cas_systems.tsv
- all_mges.tsv

Usage:
    snakemake -s workflows/Snakefile_targets --configfile target_mapping_config.yml --use-conda

Config requires:
- samplesheet: Path to original samplesheet with genome paths
- summary_dir: Path to directory containing summary TSVs from main workflow
- output_dir: Output directory for target mapping results
- target_mapping:
    min_spacer_identity: Minimum % identity for spacer BLAST (default: 90)
    min_spacer_coverage: Minimum query coverage for spacer BLAST (default: 100)
    max_mismatches: Maximum mismatches allowed (default: 2)
    threads: Number of threads for BLAST (default: 4)
"""

import os
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(workflow.basedir).parent))
from workflows.common import load_samplesheet, genome_for


# ============================= Config & Setup =============================

configfile: "target_mapping_config.yml"

# Load samplesheet
SAMPLESHEET = load_samplesheet(config["samplesheet"])
SAMPLES = list(SAMPLESHEET.keys())

# Paths
SUMMARY_DIR = config["summary_dir"]
OUTPUT_DIR = config["output_dir"]

# Target mapping parameters
TARGET_CONFIG = config.get("target_mapping", {})
MIN_SPACER_IDENTITY = TARGET_CONFIG.get("min_spacer_identity", 90)
MIN_SPACER_COVERAGE = TARGET_CONFIG.get("min_spacer_coverage", 100)
MAX_MISMATCHES = TARGET_CONFIG.get("max_mismatches", 2)
THREADS = TARGET_CONFIG.get("threads", 4)

# Module paths
MODULES_DIR = Path(workflow.basedir).parent / "modules"
PATTERN_AGG_SCRIPT = MODULES_DIR / "pattern_aggregation" / "aggregate_patterns.py"
MAP_SPACER_SCRIPT = MODULES_DIR / "target_mapping" / "map_spacer_targets.py"
MAP_RM_SCRIPT = MODULES_DIR / "target_mapping" / "map_rm_targets.py"
CONSOLIDATE_SCRIPT = MODULES_DIR / "target_mapping" / "consolidate_hits.py"
MATRICES_SCRIPT = MODULES_DIR / "target_mapping" / "build_target_matrices.py"

# Conda environment
TARGET_MAPPING_ENV = str(Path(workflow.basedir).parent / "envs" / "target_mapping.yaml")


# ============================= Helper Functions =============================

def get_genome_fasta(wildcards):
    """Get genome FASTA path for a sample."""
    return genome_for(wildcards.sample, SAMPLESHEET)


# ============================= Target Rules =============================

rule all:
    """Main target: produce all outputs."""
    input:
        # Aggregated patterns
        f"{OUTPUT_DIR}/patterns/patterns_spacers.fa",
        f"{OUTPUT_DIR}/patterns/patterns_rm.fa",
        # Per-genome target GFF outputs (mobile and genomic only)
        expand(f"{OUTPUT_DIR}/{{sample}}/targets/mobile_hits.gff", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/targets/genomic_hits.gff", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/targets/suspicious_hits.gff", sample=SAMPLES),
        # Summary matrices (separate for CRISPR and RM)
        f"{OUTPUT_DIR}/summary/crispr_vs_mge_matrix.tsv",
        f"{OUTPUT_DIR}/summary/rm_vs_mge_matrix.tsv",
        f"{OUTPUT_DIR}/summary/crispr_vs_genome_matrix.tsv",
        f"{OUTPUT_DIR}/summary/rm_vs_genome_matrix.tsv"


# ============================= Pattern Aggregation =============================

rule aggregate_spacer_patterns:
    """Deduplicate CRISPR spacers across all genomes."""
    input:
        spacers_tsv = f"{SUMMARY_DIR}/all_crispr_spacers.tsv",
        cas_systems_tsv = f"{SUMMARY_DIR}/all_cas_systems.tsv"
    output:
        fasta = f"{OUTPUT_DIR}/patterns/patterns_spacers.fa",
        tsv = f"{OUTPUT_DIR}/patterns/patterns_spacers.tsv"
    log:
        f"{OUTPUT_DIR}/logs/aggregate_spacer_patterns.log"
    conda:
        TARGET_MAPPING_ENV
    shell:
        """
        python {PATTERN_AGG_SCRIPT} spacers \
            --spacers-tsv {input.spacers_tsv} \
            --cas-systems-tsv {input.cas_systems_tsv} \
            --output-fasta {output.fasta} \
            --output-tsv {output.tsv} \
            > {log} 2>&1
        """


rule aggregate_rm_patterns:
    """Deduplicate RM recognition sites across all genomes."""
    input:
        rm_systems_tsv = f"{SUMMARY_DIR}/all_rm_systems.tsv"
    output:
        fasta = f"{OUTPUT_DIR}/patterns/patterns_rm.fa",
        tsv = f"{OUTPUT_DIR}/patterns/patterns_rm.tsv"
    log:
        f"{OUTPUT_DIR}/logs/aggregate_rm_patterns.log"
    conda:
        TARGET_MAPPING_ENV
    shell:
        """
        python {PATTERN_AGG_SCRIPT} rm \
            --rm-systems-tsv {input.rm_systems_tsv} \
            --output-fasta {output.fasta} \
            --output-tsv {output.tsv} \
            > {log} 2>&1
        """


# ============================= Spacer Target Mapping =============================

rule map_spacer_targets:
    """Map spacer patterns to genome using BLAST and output temp GFF."""
    input:
        patterns_fasta = f"{OUTPUT_DIR}/patterns/patterns_spacers.fa",
        patterns_tsv = f"{OUTPUT_DIR}/patterns/patterns_spacers.tsv",
        genome_fasta = get_genome_fasta,
        mges_tsv = f"{SUMMARY_DIR}/all_mges.tsv",
        arrays_tsv = f"{SUMMARY_DIR}/all_crispr_arrays.tsv"
    output:
        gff = temp(f"{OUTPUT_DIR}/{{sample}}/targets/.spacer_hits.gff"),
        suspicious_gff = f"{OUTPUT_DIR}/{{sample}}/targets/suspicious_hits.gff"
    params:
        genome_id = lambda wildcards: wildcards.sample,
        min_identity = MIN_SPACER_IDENTITY,
        min_coverage = MIN_SPACER_COVERAGE,
        max_mismatches = MAX_MISMATCHES
    threads: THREADS
    log:
        f"{OUTPUT_DIR}/logs/{{sample}}/map_spacer_targets.log"
    conda:
        TARGET_MAPPING_ENV
    shell:
        """
        python {MAP_SPACER_SCRIPT} \
            --genome-id {params.genome_id} \
            --patterns-fasta {input.patterns_fasta} \
            --patterns-tsv {input.patterns_tsv} \
            --genome-fasta {input.genome_fasta} \
            --mges-tsv {input.mges_tsv} \
            --arrays-tsv {input.arrays_tsv} \
            --output-gff {output.gff} \
            --suspicious-gff {output.suspicious_gff} \
            --min-identity {params.min_identity} \
            --min-coverage {params.min_coverage} \
            --max-mismatches {params.max_mismatches} \
            --threads {threads} \
            > {log} 2>&1
        """


# ============================= RM Target Mapping =============================

rule map_rm_targets:
    """Map RM recognition sites to genome and output temp GFF."""
    input:
        patterns_tsv = f"{OUTPUT_DIR}/patterns/patterns_rm.tsv",
        genome_fasta = get_genome_fasta,
        mges_tsv = f"{SUMMARY_DIR}/all_mges.tsv"
    output:
        gff = temp(f"{OUTPUT_DIR}/{{sample}}/targets/.rm_hits.gff")
    params:
        genome_id = lambda wildcards: wildcards.sample
    log:
        f"{OUTPUT_DIR}/logs/{{sample}}/map_rm_targets.log"
    conda:
        TARGET_MAPPING_ENV
    shell:
        """
        python {MAP_RM_SCRIPT} \
            --genome-id {params.genome_id} \
            --patterns-tsv {input.patterns_tsv} \
            --genome-fasta {input.genome_fasta} \
            --mges-tsv {input.mges_tsv} \
            --output-gff {output.gff} \
            > {log} 2>&1
        """


# ============================= Hit Consolidation =============================

rule consolidate_hits:
    """Merge spacer and RM temp GFFs into mobile and genomic outputs."""
    input:
        spacer_gff = f"{OUTPUT_DIR}/{{sample}}/targets/.spacer_hits.gff",
        rm_gff = f"{OUTPUT_DIR}/{{sample}}/targets/.rm_hits.gff"
    output:
        mobile_gff = f"{OUTPUT_DIR}/{{sample}}/targets/mobile_hits.gff",
        genomic_gff = f"{OUTPUT_DIR}/{{sample}}/targets/genomic_hits.gff"
    log:
        f"{OUTPUT_DIR}/logs/{{sample}}/consolidate_hits.log"
    conda:
        TARGET_MAPPING_ENV
    shell:
        """
        python {CONSOLIDATE_SCRIPT} \
            --spacer-gff {input.spacer_gff} \
            --rm-gff {input.rm_gff} \
            --output-mobile-gff {output.mobile_gff} \
            --output-genomic-gff {output.genomic_gff} \
            > {log} 2>&1
        """


# ============================= Summary Matrices =============================

rule build_target_matrices:
    """Build cross-genome summary matrices."""
    input:
        # Ensure all per-genome mobile GFFs are done first
        mobile_gffs = expand(f"{OUTPUT_DIR}/{{sample}}/targets/mobile_hits.gff", sample=SAMPLES)
    output:
        crispr_vs_mge = f"{OUTPUT_DIR}/summary/crispr_vs_mge_matrix.tsv",
        rm_vs_mge = f"{OUTPUT_DIR}/summary/rm_vs_mge_matrix.tsv",
        crispr_vs_genome = f"{OUTPUT_DIR}/summary/crispr_vs_genome_matrix.tsv",
        rm_vs_genome = f"{OUTPUT_DIR}/summary/rm_vs_genome_matrix.tsv"
    params:
        genomes = SAMPLES,
        output_dir = OUTPUT_DIR
    log:
        f"{OUTPUT_DIR}/logs/build_target_matrices.log"
    conda:
        TARGET_MAPPING_ENV
    run:
        import subprocess
        genome_args = " ".join(params.genomes)
        cmd = f"""
        python {MATRICES_SCRIPT} \
            --genomes {genome_args} \
            --output-dir {params.output_dir} \
            --crispr-vs-mge-output {output.crispr_vs_mge} \
            --rm-vs-mge-output {output.rm_vs_mge} \
            --crispr-vs-genome-output {output.crispr_vs_genome} \
            --rm-vs-genome-output {output.rm_vs_genome} \
            > {log} 2>&1
        """
        shell(cmd)

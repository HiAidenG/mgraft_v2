# mgraft_v2 main workflow
# Processes genomes to create a comprehensive defense network database
import sys
from pathlib import Path

# Add workflows directory to path for shared imports
workflow_dir = Path(workflow.basedir)
sys.path.insert(0, str(workflow_dir))

from common import load_samplesheet, proteins_for, genome_for, gene_fasta_for, mges_for, get_run_id

# ======================== Configuration ========================
OUTPUT_DIR = config["output_dir"]
SUMMARY_DIR = f"{OUTPUT_DIR}/summary"
GENOMES_DIR = f"{OUTPUT_DIR}/genomes"
LOGS_DIR = config.get("logs_dir", f"{OUTPUT_DIR}/logs")

# Load samples (genome IDs) from samplesheet
samplesheet = load_samplesheet(config["samplesheet"])
SAMPLES = tuple(samplesheet.keys())

# Run ID for tracking
RUN_ID = get_run_id()

# ======================== Target Rules ========================
# Summary table configuration: maps table name to per-genome source file pattern
# Note: mges pattern is determined dynamically based on clustering config
SUMMARY_TABLE_SOURCES = {
    "cas_systems": "cas_systems.linked.tsv",
    "crispr_arrays": "crispr_arrays.linked.tsv",
    "crispr_spacers": "crispr_spacers.linked.tsv",
    "defense_systems": "defense_systems.tsv",
    "rm_systems": "rm_systems.rebase.tsv",
    "virulence_genes": "virulence_genes.tsv",
    "mges": "mges.clustered.tsv" if config.get("mge_clustering", {}).get("enabled", False) else "mges.tsv",
}

def get_summary_outputs():
    """Get summary table outputs based on config."""
    outputs = []
    
    # Feature summary tables
    enabled_tables = config.get("summary", {}).get("tables", list(SUMMARY_TABLE_SOURCES.keys()))
    outputs.extend([f"{SUMMARY_DIR}/all_{table}.tsv" for table in enabled_tables if table in SUMMARY_TABLE_SOURCES])
    
    # Genome summary table (enabled by default)
    if config.get("summary", {}).get("genome_summary", True):
        outputs.append(f"{SUMMARY_DIR}/genome_summary.tsv")
    
    return outputs

def get_clustering_outputs():
    """Get MGE clustering outputs if enabled in config."""
    if config.get("mge_clustering", {}).get("enabled", False):
        return [
            f"{SUMMARY_DIR}/mge_cluster_stats.tsv",
            f"{SUMMARY_DIR}/mge_cluster_assignments.tsv"
        ]
    return []

all_inputs = (
    # Master GFF (final output that depends on all features)
    expand(f"{GENOMES_DIR}/{{sample}}/all_features.gff", sample=SAMPLES) +
    # CRISPR feature TSVs (retained for downstream analysis)
    expand(f"{GENOMES_DIR}/{{sample}}/features/crispr_arrays.linked.tsv", sample=SAMPLES) +
    expand(f"{GENOMES_DIR}/{{sample}}/features/crispr_spacers.linked.tsv", sample=SAMPLES) +
    # Summary tables (configurable)
    get_summary_outputs() +
    # MGE clustering outputs (if enabled)
    get_clustering_outputs()
)

rule all:
    input: all_inputs


# ======================== Genome Annotation Rules ========================



rule parse_mges:
    """Parse MGE GFF to standardized TSV and extract sequences."""
    conda: "../envs/python_base.yaml"
    input:
        mge_gff = lambda wc: mges_for(wc.sample, samplesheet),
        genome_fasta = lambda wc: genome_for(wc.sample, samplesheet)
    output:
        mges_tsv = f"{GENOMES_DIR}/{{sample}}/features/mges.tsv",
        mges_fasta = f"{GENOMES_DIR}/{{sample}}/features/mges.fasta"
    params:
        genome_id = "{sample}"
    script:
        "../modules/parsers/parse_mges.py"


rule link_arrays_to_cas:
    """Link CRISPR arrays to Cas systems based on proximity with strand-aware processing."""
    conda: "../envs/python_base.yaml"
    input:
        cas_systems_tsv = f"{GENOMES_DIR}/{{sample}}/features/cas_systems.tsv",
        cas_genes_tsv = f"{GENOMES_DIR}/{{sample}}/features/cas_genes.tsv",
        crispr_arrays_tsv = f"{GENOMES_DIR}/{{sample}}/features/crispr_arrays.tsv",
        crispr_spacers_tsv = f"{GENOMES_DIR}/{{sample}}/features/crispr_spacers.tsv",
        genome_fasta = lambda wc: genome_for(wc.sample, samplesheet)
    output:
        crispr_arrays_tsv = f"{GENOMES_DIR}/{{sample}}/features/crispr_arrays.linked.tsv",
        cas_systems_tsv = f"{GENOMES_DIR}/{{sample}}/features/cas_systems.linked.tsv",
        crispr_spacers_tsv = f"{GENOMES_DIR}/{{sample}}/features/crispr_spacers.linked.tsv"
    params:
        max_distance = config.get("assign_arrays", {}).get("max_distance_to_cas", 10000),
        revcomp_if_cas_opposite = config.get("assign_arrays", {}).get("revcomp_if_cas_opposite", True)
    script:
        "../modules/genome_annot/link_arrays_to_cas.py"


rule query_rebase:
    """Align RM system proteins to REBASE database and merge annotations."""
    conda: "../envs/diamond.yaml"
    input:
        rm_genes_tsv = f"{GENOMES_DIR}/{{sample}}/features/rm_genes.tsv",
        rm_systems_tsv = f"{GENOMES_DIR}/{{sample}}/features/rm_systems.tsv",
        proteins = lambda wc: proteins_for(wc.sample, samplesheet),
        rebase_db = config["rebase_db"]
    output:
        rm_genes_tsv = f"{GENOMES_DIR}/{{sample}}/features/rm_genes.rebase.tsv",
        rm_systems_tsv = f"{GENOMES_DIR}/{{sample}}/features/rm_systems.rebase.tsv"
    params:
        query_rebase_config = config.get("query_rebase", {})
    threads: config.get("query_rebase", {}).get("threads", 4)
    resources:
        mem_mb=16000,
        runtime=60
    script:
        "../modules/genome_annot/query_rebase.py"




rule build_master_gff:
    """Build master GFF3 file from all feature-type-specific TSVs."""
    conda: "../envs/python_base.yaml"
    input:
        mges_tsv = get_mges_tsv_for_gff,
        cas_systems_tsv = f"{GENOMES_DIR}/{{sample}}/features/cas_systems.linked.tsv",
        cas_genes_tsv = f"{GENOMES_DIR}/{{sample}}/features/cas_genes.tsv",
        rm_systems_tsv = f"{GENOMES_DIR}/{{sample}}/features/rm_systems.rebase.tsv",
        rm_genes_tsv = f"{GENOMES_DIR}/{{sample}}/features/rm_genes.rebase.tsv",
        defense_systems_tsv = f"{GENOMES_DIR}/{{sample}}/features/defense_systems.tsv",
        defense_genes_tsv = f"{GENOMES_DIR}/{{sample}}/features/defense_genes.tsv",
        antidefense_systems_tsv = f"{GENOMES_DIR}/{{sample}}/features/antidefense_systems.tsv",
        antidefense_genes_tsv = f"{GENOMES_DIR}/{{sample}}/features/antidefense_genes.tsv",
        crispr_arrays_tsv = f"{GENOMES_DIR}/{{sample}}/features/crispr_arrays.linked.tsv",
        crispr_spacers_tsv = f"{GENOMES_DIR}/{{sample}}/features/crispr_spacers.linked.tsv",
        virulence_genes_tsv = f"{GENOMES_DIR}/{{sample}}/features/virulence_genes.tsv"
    output:
        gff = f"{GENOMES_DIR}/{{sample}}/all_features.gff"
    params:
        genome_id = "{sample}"
    script:
        "../modules/parsers/build_master_gff.py"


# ======================== Summary Rules ========================
def get_source_files_for_summary(wildcards):
    """Get per-genome source files for a summary table."""
    source_pattern = SUMMARY_TABLE_SOURCES[wildcards.table]
    return expand(f"{GENOMES_DIR}/{{sample}}/features/{source_pattern}", sample=SAMPLES)

rule summarize_tables:
    """Aggregate per-genome feature tables into summary tables.
    
    Configure which tables to generate via config:
        summary:
          tables:
            - cas_systems
            - crispr_arrays
            - crispr_spacers
            - defense_systems
            - rm_systems
            - virulence_genes
            - mges
    
    By default, all tables are generated.
    """
    conda: "../envs/python_base.yaml"
    input:
        get_source_files_for_summary
    output:
        f"{SUMMARY_DIR}/all_{{table}}.tsv"
    wildcard_constraints:
        table = "|".join(SUMMARY_TABLE_SOURCES.keys())
    script:
        "../modules/pattern_aggregation/summarize_tables.py"


rule summarize_genomes:
    """Generate per-genome summary statistics table.
    
    Aggregates counts, type breakdowns, quality metrics, and mobility metrics
    for each genome. Enable/disable via config:
        summary:
          genome_summary: true  # default
    """
    conda: "../envs/genome_processing.yaml"
    input:
        cas = expand(f"{GENOMES_DIR}/{{sample}}/features/cas_systems.linked.tsv", sample=SAMPLES),
        arrays = expand(f"{GENOMES_DIR}/{{sample}}/features/crispr_arrays.linked.tsv", sample=SAMPLES),
        spacers = expand(f"{GENOMES_DIR}/{{sample}}/features/crispr_spacers.linked.tsv", sample=SAMPLES),
        defense = expand(f"{GENOMES_DIR}/{{sample}}/features/defense_systems.tsv", sample=SAMPLES),
        rm = expand(f"{GENOMES_DIR}/{{sample}}/features/rm_systems.rebase.tsv", sample=SAMPLES),
        vf = expand(f"{GENOMES_DIR}/{{sample}}/features/virulence_genes.tsv", sample=SAMPLES),
        antidefense = expand(f"{GENOMES_DIR}/{{sample}}/features/antidefense_systems.tsv", sample=SAMPLES),
        mges = expand(f"{GENOMES_DIR}/{{sample}}/features/mges.clustered.tsv", sample=SAMPLES) if config.get("mge_clustering", {}).get("enabled", False) else expand(f"{GENOMES_DIR}/{{sample}}/features/mges.tsv", sample=SAMPLES),
        genomes = [genome_for(sample, samplesheet) for sample in SAMPLES],
    output:
        f"{SUMMARY_DIR}/genome_summary.tsv"
    params:
        samples = SAMPLES,
        genomes_dir = GENOMES_DIR,
        genome_fastas = {sample: genome_for(sample, samplesheet) for sample in SAMPLES}
    script:
        "../modules/pattern_aggregation/summarize_genomes.py"

# ======================== MGE Clustering Rules ========================
# These rules cluster MGEs for deduplication purposes using vclust

rule combine_mge_sequences:
    """Concatenate per-genome MGE FASTAs for clustering."""
    input:
        mges_fastas = expand(f"{GENOMES_DIR}/{{sample}}/features/mges.fasta", sample=SAMPLES)
    output:
        combined_fasta = f"{SUMMARY_DIR}/clustering/mge_sequences.fasta"
    shell:
        "cat {input.mges_fastas} > {output.combined_fasta}"


rule run_vclust_clustering:
    """Run vclust to cluster MGE sequences by ANI."""
    conda: "../envs/vclust.yaml"
    input:
        combined_fasta = f"{SUMMARY_DIR}/clustering/mge_sequences.fasta"
    output:
        prefilter = f"{SUMMARY_DIR}/clustering/vclust_prefilter.txt",
        ani = f"{SUMMARY_DIR}/clustering/vclust_ani.tsv",
        ids = f"{SUMMARY_DIR}/clustering/vclust_ani.ids.tsv",
        clusters = f"{SUMMARY_DIR}/clustering/vclust_clusters.tsv"
    params:
        prefilter_identity = config.get("mge_clustering", {}).get("prefilter", {}).get("identity", 0.95),
        align_outfmt = config.get("mge_clustering", {}).get("align", {}).get("outfmt", "lite"),
        cluster_algorithm = config.get("mge_clustering", {}).get("cluster", {}).get("algorithm", "leiden"),
        cluster_metric = config.get("mge_clustering", {}).get("cluster", {}).get("metric", "ani"),
        cluster_ani = config.get("mge_clustering", {}).get("cluster", {}).get("ani_threshold", 0.95),
        cluster_qcov = config.get("mge_clustering", {}).get("cluster", {}).get("query_coverage", 0.85),
        cluster_rcov = config.get("mge_clustering", {}).get("cluster", {}).get("ref_coverage", 0.85)
    threads: config.get("mge_clustering", {}).get("threads", 8)
    resources:
        mem_mb=32000,
        runtime=240,
        cpus_per_task=config.get("mge_clustering", {}).get("threads", 8)
    script:
        "../modules/mge_clustering/run_vclust.py"


rule assign_mge_cluster_ids:
    """Assign cluster IDs to MGE features in per-genome TSVs."""
    conda: "../envs/python_base.yaml"
    input:
        clusters = f"{SUMMARY_DIR}/clustering/vclust_clusters.tsv",
        mges_tsvs = expand(f"{GENOMES_DIR}/{{sample}}/features/mges.tsv", sample=SAMPLES)
    output:
        mges_tsvs = expand(f"{GENOMES_DIR}/{{sample}}/features/mges.clustered.tsv", sample=SAMPLES),
        cluster_stats = f"{SUMMARY_DIR}/mge_cluster_stats.tsv",
        cluster_assignments = f"{SUMMARY_DIR}/mge_cluster_assignments.tsv"
    script:
        "../modules/mge_clustering/assign_cluster_ids.py"

# DefenseFinder workflow
# Run DefenseFinder to identify defense systems in genomes
import sys
from pathlib import Path

# Add workflows directory to path for shared imports
workflow_dir = Path(workflow.basedir)
sys.path.insert(0, str(workflow_dir))

from common import load_samplesheet, proteins_for

# ======================== Configuration ========================
OUTPUT_DIR = config["output_dir"].rstrip("/")
LOGS_DIR = config.get("logs_dir", f"{OUTPUT_DIR}/logs").rstrip("/")

# Load samples (genome IDs) from samplesheet
samplesheet = load_samplesheet(config["samplesheet"])
SAMPLES = tuple(samplesheet.keys())

# ======================== Target Rules ========================
rule all:
    input:
        expand(f"{OUTPUT_DIR}/{{sample}}/cas_systems.tsv", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/cas_genes.tsv", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/rm_systems.tsv", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/rm_genes.tsv", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/defense_systems.tsv", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/defense_genes.tsv", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/antidefense_systems.tsv", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/antidefense_genes.tsv", sample=SAMPLES)


# ======================== DefenseFinder Rule ========================
rule run_defensefinder:
    """Run DefenseFinder to identify defense systems."""
    conda: "../envs/defensefinder.yaml"
    input:
        protein_fasta = lambda wc: proteins_for(wc.sample, samplesheet)
    output:
        systems_tsv = f"{OUTPUT_DIR}/{{sample}}/defensefinder/defense_finder_systems.tsv",
        genes_tsv = f"{OUTPUT_DIR}/{{sample}}/defensefinder/defense_finder_genes.tsv",
        hmmer_tsv = f"{OUTPUT_DIR}/{{sample}}/defensefinder/defense_finder_hmmer.tsv"
    params:
        outdir = f"{OUTPUT_DIR}/{{sample}}/defensefinder"
    threads: config.get("defensefinder", {}).get("threads", 4)
    resources:
        mem_mb = config.get("defensefinder", {}).get("mem_mb", 4000),
        runtime = config.get("defensefinder", {}).get("runtime", 60),
        cpus_per_task = config.get("defensefinder", {}).get("threads", 4)
    log:
        f"{LOGS_DIR}/defensefinder/{{sample}}.log"
    shell:
        """
        defense-finder run -a --out-dir {params.outdir} --workers {threads} {input.protein_fasta} 2> {log}
        # Ensure outputs exist
        touch {output.systems_tsv}
        touch {output.genes_tsv}
        touch {output.hmmer_tsv}
        """


# ======================== Parse DefenseFinder Rule ========================
rule parse_defensefinder:
    """Parse DefenseFinder output to extract defense systems and genes."""
    conda: "../envs/python_base.yaml"
    input:
        systems_tsv = rules.run_defensefinder.output.systems_tsv,
        genes_tsv = rules.run_defensefinder.output.genes_tsv,
        hmmer_tsv = rules.run_defensefinder.output.hmmer_tsv,
        faa = lambda wc: proteins_for(wc.sample, samplesheet)
    output:
        cas_systems_tsv = f"{OUTPUT_DIR}/{{sample}}/cas_systems.tsv",
        cas_genes_tsv = f"{OUTPUT_DIR}/{{sample}}/cas_genes.tsv",
        rm_systems_tsv = f"{OUTPUT_DIR}/{{sample}}/rm_systems.tsv",
        rm_genes_tsv = f"{OUTPUT_DIR}/{{sample}}/rm_genes.tsv",
        defense_systems_tsv = f"{OUTPUT_DIR}/{{sample}}/defense_systems.tsv",
        defense_genes_tsv = f"{OUTPUT_DIR}/{{sample}}/defense_genes.tsv",
        antidefense_systems_tsv = f"{OUTPUT_DIR}/{{sample}}/antidefense_systems.tsv",
        antidefense_genes_tsv = f"{OUTPUT_DIR}/{{sample}}/antidefense_genes.tsv"
    params:
        genome_id = lambda wc: wc.sample
    log:
        f"{LOGS_DIR}/parse_defensefinder/{{sample}}.log"
    script:
        "../modules/parsers/parse_defensefinder.py"

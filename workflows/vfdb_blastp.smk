# VFDB BLASTP workflow
# Run Diamond BLASTP against VFDB database to identify virulence genes
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
        expand(f"{OUTPUT_DIR}/{{sample}}/vfdb_blastp/vfdb_hits.tsv", sample=SAMPLES)


# ======================== VFDB BLASTP Rule ========================
rule vfdb_blastp:
    """Run Diamond BLASTP against VFDB database to identify virulence genes."""
    conda: "../envs/diamond.yaml"
    input:
        proteins = lambda wc: proteins_for(wc.sample, samplesheet),
        vfdb_db = config["vfdb_db"]
    output:
        vfdb_hits = f"{OUTPUT_DIR}/{{sample}}/vfdb_blastp/vfdb_hits.tsv"
    params:
        outdir = f"{OUTPUT_DIR}/{{sample}}/vfdb_blastp",
        min_identity = config.get("vfdb_blastp", {}).get("min_identity", 95.0),
        min_query_cov = config.get("vfdb_blastp", {}).get("min_query_coverage", 85.0),
        min_subject_cov = config.get("vfdb_blastp", {}).get("min_subject_coverage", 85.0),
        evalue = config.get("vfdb_blastp", {}).get("evalue", 1e-5)
    threads: config.get("vfdb_blastp", {}).get("threads", 4)
    resources:
        mem_mb = config.get("vfdb_blastp", {}).get("mem_mb", 16000),
        runtime = config.get("vfdb_blastp", {}).get("runtime", 60),
        cpus_per_task = config.get("vfdb_blastp", {}).get("threads", 4)
    log:
        f"{LOGS_DIR}/vfdb_blastp/{{sample}}.log"
    shell:
        """
        mkdir -p {params.outdir}
        diamond blastp \
            --db {input.vfdb_db} \
            --query {input.proteins} \
            --out {output.vfdb_hits} \
            --threads {threads} \
            --id {params.min_identity} \
            --query-cover {params.min_query_cov} \
            --subject-cover {params.min_subject_cov} \
            --evalue {params.evalue} \
            --outfmt 6 qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue bitscore stitle \
            --max-target-seqs 1 \
            --sensitive 2> {log}
        # Ensure output exists
        touch {output.vfdb_hits}
        """

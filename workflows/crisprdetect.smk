# CRISPRDetect workflow
# Run CRISPRDetect 3.0 to identify CRISPR arrays in genomes
import sys
from pathlib import Path

# Add workflows directory to path for shared imports
workflow_dir = Path(workflow.basedir)
sys.path.insert(0, str(workflow_dir))

from common import load_samplesheet, genome_for

# ======================== Configuration ========================
OUTPUT_DIR = config["output_dir"].rstrip("/")
LOGS_DIR = config.get("logs_dir", f"{OUTPUT_DIR}/logs").rstrip("/")

# Load samples (genome IDs) from samplesheet
samplesheet = load_samplesheet(config["samplesheet"])
SAMPLES = tuple(samplesheet.keys())

# CRISPRDetect path
CRISPRDETECT_PATH = str(Path(workflow.basedir).parent / "tools" / "CRISPRDetect" / "CRISPRDetect3")

# ======================== Target Rules ========================
rule all:
    input:
        expand(f"{OUTPUT_DIR}/{{sample}}/crispr_arrays.tsv", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/{{sample}}/crispr_spacers.tsv", sample=SAMPLES)


# ======================== CRISPRDetect Rule ========================
rule run_crisprdetect:
    """Run CRISPRDetect 3.0 to identify CRISPR arrays with strand information."""
    conda: "../envs/crisprdetect.yaml"
    input:
        genome_fasta = lambda wc: genome_for(wc.sample, samplesheet)
    output:
        txt = f"{OUTPUT_DIR}/{{sample}}/crisprdetect/output.txt",
        gff = f"{OUTPUT_DIR}/{{sample}}/crisprdetect/output.gff"
    params:
        outdir = f"{OUTPUT_DIR}/{{sample}}/crisprdetect",
        output_prefix = f"{OUTPUT_DIR}/{{sample}}/crisprdetect/output",
        quality_cutoff = config.get("crisprdetect", {}).get("quality_cutoff", 0),
        min_repeats = config.get("crisprdetect", {}).get("min_repeats", 3),
        crisprdetect_path = CRISPRDETECT_PATH
    threads: config.get("crisprdetect", {}).get("threads", 4)
    resources:
        mem_mb = config.get("crisprdetect", {}).get("mem_mb", 8000),
        runtime = config.get("crisprdetect", {}).get("runtime", 60),
        cpus_per_task = config.get("crisprdetect", {}).get("threads", 4)
    log:
        f"{LOGS_DIR}/crisprdetect/{{sample}}.log"
    shell:
        """
        mkdir -p {params.outdir}
        perl {params.crisprdetect_path} \
            -f {input.genome_fasta} \
            -o {params.output_prefix} \
            -array_quality_score_cutoff {params.quality_cutoff} \
            -minimum_no_of_repeats {params.min_repeats} \
            -T {threads} 2> {log}
        # Ensure outputs exist
        touch {output.txt}
        touch {output.gff}
        """


# ======================== Parse CRISPRDetect Rule ========================
rule parse_crisprdetect:
    """Parse CRISPRDetect output to extract CRISPR arrays and spacers in TSV format."""
    conda: "../envs/python_base.yaml"
    input:
        crisprdetect_txt = rules.run_crisprdetect.output.txt,
        crisprdetect_gff = rules.run_crisprdetect.output.gff
    output:
        crispr_arrays_tsv = f"{OUTPUT_DIR}/{{sample}}/crispr_arrays.tsv",
        crispr_spacers_tsv = f"{OUTPUT_DIR}/{{sample}}/crispr_spacers.tsv"
    params:
        genome_id = lambda wc: wc.sample
    log:
        f"{LOGS_DIR}/parse_crisprdetect/{{sample}}.log"
    script:
        "../modules/parsers/parse_crisprdetect.py"

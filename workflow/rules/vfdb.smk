# VFDB workflow
# Runs Diamond BLASTP against VFDB database for virulence factor annotation

RETAIN_INTERMEDIATES = config.get("retain_intermediates", False)

rule query_vfdb:
    """Align proteins to VFDB database for virulence factor annotation."""
    input:
        proteins = lambda wc: proteins_for(wc.sample, SAMPLESHEET, OUTPUT_DIR),
        vfdb_db = config["vfdb_db"]
    output:
        hits = OUTPUT_DIR + "/{sample}/intermediates/vfdb-hits-raw.tsv" if RETAIN_INTERMEDIATES else temp(OUTPUT_DIR + "/{sample}/vfdb-hits-raw.tsv"),
        parsed = OUTPUT_DIR + "/{sample}/intermediates/parsed/virulence_genes.tsv" if RETAIN_INTERMEDIATES else temp(OUTPUT_DIR + "/{sample}/parsed/virulence_genes.tsv"),
        stats = OUTPUT_DIR + "/{sample}/vfdb-alignment-stats.json"
    params:
        min_identity = config.get("vfdb", {}).get("min_identity", 95.0),
        min_query_cov = config.get("vfdb", {}).get("min_query_coverage", 85.0),
        min_subject_cov = config.get("vfdb", {}).get("min_subject_coverage", 85.0),
        evalue = config.get("vfdb", {}).get("evalue", 1e-5),
        genome_id = "{sample}"
    log:
        OUTPUT_DIR + "/logs/{sample}/query_vfdb.log"
    conda:
        "../envs/diamond.yaml"
    threads: config.get("vfdb", {}).get("threads", 4)
    resources:
        mem_mb = 16000,
        runtime = 60
    script:
        "../scripts/parsers/query_vfdb.py"

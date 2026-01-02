# Processing rules - CRISPR filtering, linking, clustering, MGE overlap, GFF generation
import os

BASE_ENV = config.get("base_env", "../envs/base.yaml")
RETAIN_INTERMEDIATES = config.get("retain_intermediates", False)

# Simplified path helper - always use intermediates/, mark as temp if not retaining
def temp_if_not_retained(path):
    """Wrap path in temp() if not retaining intermediates."""
    if RETAIN_INTERMEDIATES:
        return path
    return temp(path)

rule crispr_filter:
    """Filter CRISPR arrays and spacers by quality metrics."""
    input:
        arrays = OUTPUT_DIR + "/{sample}/intermediates/parsed/crispr_arrays.tsv",
        spacers = OUTPUT_DIR + "/{sample}/intermediates/parsed/crispr_spacers.tsv"
    output:
        arrays = temp_if_not_retained(OUTPUT_DIR + "/{sample}/intermediates/crispr-filtering/crispr_arrays.tsv"),
        spacers = temp_if_not_retained(OUTPUT_DIR + "/{sample}/intermediates/crispr-filtering/crispr_spacers.tsv"),
        stats = OUTPUT_DIR + "/{sample}/crispr-filtering-stats.json"
    log:
        OUTPUT_DIR + "/logs/{sample}/crispr_filtering.log"
    threads: 1
    conda: BASE_ENV
    script:
        "../scripts/processing/crispr_filtering.py"

rule array_link:
    """Link CRISPR arrays to Cas systems based on proximity."""
    input:
        cas_systems = OUTPUT_DIR + "/{sample}/intermediates/parsed/cas_systems.tsv",
        cas_genes = OUTPUT_DIR + "/{sample}/intermediates/parsed/cas_genes.tsv",
        arrays = OUTPUT_DIR + "/{sample}/intermediates/crispr-filtering/crispr_arrays.tsv",
        spacers = OUTPUT_DIR + "/{sample}/intermediates/crispr-filtering/crispr_spacers.tsv"
    output:
        cas_systems = temp_if_not_retained(OUTPUT_DIR + "/{sample}/intermediates/linked-annotations/cas_systems.linked.tsv"),
        arrays = temp_if_not_retained(OUTPUT_DIR + "/{sample}/intermediates/linked-annotations/crispr_arrays.linked.tsv"),
        spacers = temp_if_not_retained(OUTPUT_DIR + "/{sample}/intermediates/linked-annotations/crispr_spacers.linked.tsv")
    log:
        OUTPUT_DIR + "/logs/{sample}/array_cas_linking.log"
    conda: BASE_ENV
    script:
        "../scripts/processing/array_cas_linking.py"

rule rm_annotate:
    """Annotate RM systems with REBASE recognition sequences."""
    input:
        rm_systems = OUTPUT_DIR + "/{sample}/intermediates/parsed/rm_systems.tsv",
        rm_genes = OUTPUT_DIR + "/{sample}/intermediates/parsed/rm_genes.tsv",
        rebase_hits = OUTPUT_DIR + "/{sample}/intermediates/rebase-hits-raw.tsv"
    output:
        rm_systems = temp_if_not_retained(OUTPUT_DIR + "/{sample}/intermediates/linked-annotations/rm_systems.rebase.tsv"),
        rm_genes = temp_if_not_retained(OUTPUT_DIR + "/{sample}/intermediates/linked-annotations/rm_genes.rebase.tsv")
    log:
        OUTPUT_DIR + "/logs/{sample}/rm_annotation.log"
    conda: BASE_ENV
    script:
        "../scripts/processing/rm_annotation.py"

# Clustering requires aggregation of ALL samples
rule concat_for_clustering:
    input:
        arrays = expand(OUTPUT_DIR + "/{sample}/intermediates/linked-annotations/crispr_arrays.linked.tsv", sample=SAMPLES),
        spacers = expand(OUTPUT_DIR + "/{sample}/intermediates/linked-annotations/crispr_spacers.linked.tsv", sample=SAMPLES)
    output:
        arrays = temp(OUTPUT_DIR + "/clustering/all_arrays.raw.tsv"),
        spacers = temp(OUTPUT_DIR + "/clustering/all_spacers.raw.tsv")
    shell:
        """
        mkdir -p $(dirname {output.arrays})
        awk 'FNR==1 && NR!=1{{next;}}{{print}}' {input.arrays} > {output.arrays}
        awk 'FNR==1 && NR!=1{{next;}}{{print}}' {input.spacers} > {output.spacers}
        """

rule crispr_cluster:
    """Cluster CRISPR spacers and repeats across all samples."""
    input:
        arrays = OUTPUT_DIR + "/clustering/all_arrays.raw.tsv",
        spacers = OUTPUT_DIR + "/clustering/all_spacers.raw.tsv"
    output:
        arrays = temp(OUTPUT_DIR + "/clustering/all_arrays.clustered.tsv"),
        spacers = temp(OUTPUT_DIR + "/clustering/all_spacers.clustered.tsv"),
        stats = OUTPUT_DIR + "/crispr-clustering-stats.json"
    params:
        work_dir = OUTPUT_DIR + "/clustering/mmseqs_tmp"
    log:
        OUTPUT_DIR + "/logs/crispr_clustering.log"
    threads: 8
    conda: BASE_ENV
    script:
        "../scripts/processing/crispr_clustering.py"

# MGE Clustering (if MGE annotations available)
if config.get("run_mge"):
    rule mge_cluster:
        """Deduplicate MGEs based on whole sequence similarity using vclust."""
        input:
            mges = expand(OUTPUT_DIR + "/{sample}/final/mges.tsv", sample=SAMPLES),
            mges_fasta = expand(OUTPUT_DIR + "/{sample}/intermediates/parsed/mges.fasta", sample=SAMPLES)
        output:
            mges = OUTPUT_DIR + "/clustering/all_mges.clustered.tsv",
            stats = OUTPUT_DIR + "/mge-clustering-stats.json"
        params:
            work_dir = OUTPUT_DIR + "/clustering/mge_vclust_tmp"
        log:
            OUTPUT_DIR + "/logs/mge_clustering.log"
        threads: config.get("mge", {}).get("threads", 8)
        conda: "../envs/vclust.yaml"
        script:
            "../scripts/processing/mge_clustering.py"

    # Split clustered CRISPR data back to per-sample files in final/
    rule split_clustered_crispr:
        """Split cross-sample clustered CRISPR data back to per-sample TSVs in final/ for MGE overlap analysis."""
        input:
            arrays = OUTPUT_DIR + "/clustering/all_arrays.clustered.tsv",
            spacers = OUTPUT_DIR + "/clustering/all_spacers.clustered.tsv"
        output:
            arrays = OUTPUT_DIR + "/{sample}/final/crispr_arrays.tsv",
            spacers = OUTPUT_DIR + "/{sample}/final/crispr_spacers.tsv"
        params:
            genome_id = "{sample}"
        conda: BASE_ENV
        shell:
            """
            awk -F'\\t' -v genome="{params.genome_id}" 'NR==1 || $1==genome' {input.arrays} > {output.arrays}
            awk -F'\\t' -v genome="{params.genome_id}" 'NR==1 || $1==genome' {input.spacers} > {output.spacers}
            """

    # Copy other feature TSVs to final/ directory for MGE overlap
    rule prepare_features_for_mge:
        """Copy linked annotation TSVs to final/ directory before MGE overlap analysis."""
        input:
            cas_systems = OUTPUT_DIR + "/{sample}/intermediates/linked-annotations/cas_systems.linked.tsv",
            rm_systems = OUTPUT_DIR + "/{sample}/intermediates/linked-annotations/rm_systems.rebase.tsv",
            rm_genes = OUTPUT_DIR + "/{sample}/intermediates/linked-annotations/rm_genes.rebase.tsv",
            cas_genes = OUTPUT_DIR + "/{sample}/intermediates/parsed/cas_genes.tsv",
            defense_systems = OUTPUT_DIR + "/{sample}/intermediates/parsed/defense_systems.tsv",
            defense_genes = OUTPUT_DIR + "/{sample}/intermediates/parsed/defense_genes.tsv"
        output:
            cas_systems = OUTPUT_DIR + "/{sample}/final/cas_systems.tsv",
            rm_systems = OUTPUT_DIR + "/{sample}/final/rm_systems.tsv",
            rm_genes = OUTPUT_DIR + "/{sample}/final/rm_genes.tsv",
            cas_genes = OUTPUT_DIR + "/{sample}/final/cas_genes.tsv",
            defense_systems = OUTPUT_DIR + "/{sample}/final/defense_systems.tsv",
            defense_genes = OUTPUT_DIR + "/{sample}/final/defense_genes.tsv"
        shell:
            """
            mkdir -p $(dirname {output.cas_systems})
            cp {input.cas_systems} {output.cas_systems}
            cp {input.rm_systems} {output.rm_systems}
            cp {input.rm_genes} {output.rm_genes}
            cp {input.cas_genes} {output.cas_genes}
            cp {input.defense_systems} {output.defense_systems}
            cp {input.defense_genes} {output.defense_genes}
            """

    rule mge_overlap:
        """Assign mge_id to genes/arrays and calculate mobility_score for systems."""
        input:
            mges = OUTPUT_DIR + "/clustering/all_mges.clustered.tsv",
            rm_systems = OUTPUT_DIR + "/{sample}/final/rm_systems.tsv",
            rm_genes = OUTPUT_DIR + "/{sample}/final/rm_genes.tsv",
            cas_systems = OUTPUT_DIR + "/{sample}/final/cas_systems.tsv",
            cas_genes = OUTPUT_DIR + "/{sample}/final/cas_genes.tsv",
            defense_systems = OUTPUT_DIR + "/{sample}/final/defense_systems.tsv",
            defense_genes = OUTPUT_DIR + "/{sample}/final/defense_genes.tsv",
            crispr_arrays = OUTPUT_DIR + "/{sample}/final/crispr_arrays.tsv",
            crispr_spacers = OUTPUT_DIR + "/{sample}/final/crispr_spacers.tsv"
        output:
            rm_systems = temp_if_not_retained(OUTPUT_DIR + "/{sample}/intermediates/mge-overlap/rm_systems.mge.tsv"),
            rm_genes = temp_if_not_retained(OUTPUT_DIR + "/{sample}/intermediates/mge-overlap/rm_genes.mge.tsv"),
            cas_systems = temp_if_not_retained(OUTPUT_DIR + "/{sample}/intermediates/mge-overlap/cas_systems.mge.tsv"),
            cas_genes = temp_if_not_retained(OUTPUT_DIR + "/{sample}/intermediates/mge-overlap/cas_genes.mge.tsv"),
            defense_systems = temp_if_not_retained(OUTPUT_DIR + "/{sample}/intermediates/mge-overlap/defense_systems.mge.tsv"),
            defense_genes = temp_if_not_retained(OUTPUT_DIR + "/{sample}/intermediates/mge-overlap/defense_genes.mge.tsv"),
            crispr_arrays = temp_if_not_retained(OUTPUT_DIR + "/{sample}/intermediates/mge-overlap/crispr_arrays.mge.tsv"),
            crispr_spacers = temp_if_not_retained(OUTPUT_DIR + "/{sample}/intermediates/mge-overlap/crispr_spacers.mge.tsv")
        params:
            genome_id = "{sample}"
        log:
            OUTPUT_DIR + "/logs/{sample}/mge_overlap.log"
        conda: BASE_ENV
        script:
            "../scripts/processing/mge_overlap.py"

# Path helper for inputs that may come from MGE overlap or upstream
def _feature_input_path(sample, feature_base, upstream_subdir, upstream_basename):
    """Get input path for a feature type."""
    if config.get("run_mge"):
        return OUTPUT_DIR + f"/{sample}/intermediates/mge-overlap/{feature_base}.mge.tsv"
    else:
        return OUTPUT_DIR + f"/{sample}/intermediates/{upstream_subdir}/{upstream_basename}.tsv"

rule build_gff:
    input:
        mges = (OUTPUT_DIR + "/clustering/all_mges.clustered.tsv") if config.get("run_mge") else [],
        cas_systems = lambda wc: _feature_input_path(wc.sample, "cas_systems", "linked-annotations", "cas_systems.linked"),
        rm_systems = lambda wc: _feature_input_path(wc.sample, "rm_systems", "linked-annotations", "rm_systems.rebase"),
        rm_genes = lambda wc: _feature_input_path(wc.sample, "rm_genes", "linked-annotations", "rm_genes.rebase"),
        cas_genes = lambda wc: _feature_input_path(wc.sample, "cas_genes", "parsed", "cas_genes"),
        defense_systems = lambda wc: _feature_input_path(wc.sample, "defense_systems", "parsed", "defense_systems"),
        defense_genes = lambda wc: _feature_input_path(wc.sample, "defense_genes", "parsed", "defense_genes"),
        antidefense_systems = OUTPUT_DIR + "/{sample}/intermediates/parsed/antidefense_systems.tsv",
        antidefense_genes = OUTPUT_DIR + "/{sample}/intermediates/parsed/antidefense_genes.tsv",
        virulence_genes = OUTPUT_DIR + "/{sample}/intermediates/parsed/virulence_genes.tsv",
        crispr_arrays = (OUTPUT_DIR + "/clustering/all_arrays.clustered.tsv") if not config.get("run_mge") else (lambda wc: _feature_input_path(wc.sample, "crispr_arrays", "", "")),
        crispr_spacers = (OUTPUT_DIR + "/clustering/all_spacers.clustered.tsv") if not config.get("run_mge") else (lambda wc: _feature_input_path(wc.sample, "crispr_spacers", "", ""))
    output:
        gff = OUTPUT_DIR + "/{sample}/all_features.gff",
        tsv_marker = OUTPUT_DIR + "/{sample}/final/.tsvs_written"
    params:
        genome_id = "{sample}",
        out_tsv_dir = OUTPUT_DIR + "/{sample}/final"
    log:
        OUTPUT_DIR + "/logs/{sample}/build_master_gff.log"
    conda: BASE_ENV
    script:
        "../scripts/reporting/build_master_gff.py"

rule genome_summary:
    input:
        markers = expand(OUTPUT_DIR + "/{sample}/final/.tsvs_written", sample=SAMPLES)
    output:
        summary = OUTPUT_DIR + "/summary/genome_summary.tsv",
        cas = OUTPUT_DIR + "/summary/all_cas_systems.tsv.gz",
        crispr_arrays = OUTPUT_DIR + "/summary/all_crispr_arrays.tsv.gz",
        crispr_spacers = OUTPUT_DIR + "/summary/all_crispr_spacers.tsv.gz",
        rm_systems = OUTPUT_DIR + "/summary/all_rm_systems.tsv.gz"
    params:
        samples = SAMPLES,
        input_dirs = {s: OUTPUT_DIR + "/" + s + "/final" for s in SAMPLES},
        genome_fastas = {s: genome_for(s, SAMPLESHEET) for s in SAMPLES}
    log:
        OUTPUT_DIR + "/logs/genome_summary.log"
    conda: BASE_ENV
    script:
        "../scripts/reporting/genome_summary.py"

# Clean up clustering directory after all samples are processed
if not RETAIN_INTERMEDIATES:
    onsuccess:
        import shutil
        clustering_dir = os.path.join(config.get("output_dir", "results"), "clustering")
        if os.path.exists(clustering_dir):
            shutil.rmtree(clustering_dir)

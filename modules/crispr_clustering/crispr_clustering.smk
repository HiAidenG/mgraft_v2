"""
Snakemake rules for CRISPR clustering.

Clusters CRISPR arrays by repeat and spacer sequences using MMseqs2.
"""


rule cluster_crispr_all:
    """Run complete CRISPR clustering pipeline."""
    input:
        arrays = "summary/all_crispr_arrays.tsv",
        spacers = "summary/all_crispr_spacers.tsv"
    output:
        arrays_clustered = "summary/crispr_arrays_clustered.tsv",
        spacers_clustered = "summary/crispr_spacers_clustered.tsv",
        repeat_info = "summary/crispr_clustering/repeat_cluster_info.tsv",
        spacer_info = "summary/crispr_clustering/spacer_cluster_info.tsv",
        summary = "summary/crispr_clustering/clustering_summary.txt"
    params:
        output_dir = "summary/crispr_clustering",
        repeat_identity = config.get("crispr_repeat_identity", 1.0),
        repeat_coverage = config.get("crispr_repeat_coverage", 1.0),
        spacer_identity = config.get("crispr_spacer_identity", 0.9),
        spacer_coverage = config.get("crispr_spacer_coverage", 0.9)
    threads: config.get("crispr_clustering_threads", 4)
    conda: "../envs/mmseqs.yaml"
    log: "logs/crispr_clustering.log"
    shell:
        """
        python modules/crispr_clustering/cluster_crispr.py \
            --arrays {input.arrays} \
            --spacers {input.spacers} \
            --output-dir {params.output_dir} \
            --threads {threads} \
            --repeat-identity {params.repeat_identity} \
            --repeat-coverage {params.repeat_coverage} \
            --spacer-identity {params.spacer_identity} \
            --spacer-coverage {params.spacer_coverage} \
            2>&1 | tee {log}
        
        # Copy main output files to summary directory
        cp {params.output_dir}/crispr_arrays_clustered.tsv {output.arrays_clustered}
        cp {params.output_dir}/crispr_spacers_clustered.tsv {output.spacers_clustered}
        """


rule cluster_crispr_repeats_only:
    """Cluster only repeat sequences (faster, for testing)."""
    input:
        arrays = "summary/all_crispr_arrays.tsv"
    output:
        repeat_fasta = "summary/crispr_clustering/repeat_sequences.fasta",
        repeat_clusters = "summary/crispr_clustering/repeat_cluster_info.tsv"
    params:
        output_dir = "summary/crispr_clustering",
        identity = config.get("crispr_repeat_identity", 1.0),
        coverage = config.get("crispr_repeat_coverage", 1.0)
    threads: config.get("crispr_clustering_threads", 4)
    conda: "../envs/mmseqs.yaml"
    log: "logs/crispr_repeat_clustering.log"
    script:
        "../modules/crispr_clustering/cluster_repeats_only.py"


rule cluster_crispr_spacers_only:
    """Cluster only spacer sequences."""
    input:
        spacers = "summary/all_crispr_spacers.tsv"
    output:
        spacer_fasta = "summary/crispr_clustering/spacer_sequences.fasta",
        spacer_clusters = "summary/crispr_clustering/spacer_cluster_info.tsv"
    params:
        output_dir = "summary/crispr_clustering",
        identity = config.get("crispr_spacer_identity", 0.9),
        coverage = config.get("crispr_spacer_coverage", 0.9)
    threads: config.get("crispr_clustering_threads", 4)
    conda: "../envs/mmseqs.yaml"
    log: "logs/crispr_spacer_clustering.log"
    script:
        "../modules/crispr_clustering/cluster_spacers_only.py"

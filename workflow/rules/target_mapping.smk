# Target mapping workflow
# Maps CRISPR spacers and RM recognition sites to target genomes
# Only runs when MGE annotations are provided

if config.get("run_mge"):
    rule aggregate_patterns:
        """Deduplicate CRISPR spacers and RM sites across genomes."""
        input:
            spacers = OUTPUT_DIR + "/summary/all_crispr_spacers.tsv.gz",
            cas_systems = OUTPUT_DIR + "/summary/all_cas_systems.tsv.gz",
            rm_systems = OUTPUT_DIR + "/summary/all_rm_systems.tsv.gz"
        output:
            spacer_fasta = OUTPUT_DIR + "/patterns/patterns_spacers.fa",
            spacer_tsv = OUTPUT_DIR + "/patterns/patterns_spacers.tsv",
            rm_fasta = OUTPUT_DIR + "/patterns/patterns_rm.fa",
            rm_tsv = OUTPUT_DIR + "/patterns/patterns_rm.tsv"
        log:
            OUTPUT_DIR + "/logs/aggregate_patterns.log"
        conda:
            "../envs/base.yaml"
        script:
            "../scripts/targets/aggregate_patterns.py"


    rule map_spacer_targets:
        """Map spacer patterns to genome using BLAST."""
        input:
            patterns_fasta = OUTPUT_DIR + "/patterns/patterns_spacers.fa",
            patterns_tsv = OUTPUT_DIR + "/patterns/patterns_spacers.tsv",
            genome = lambda wc: genome_for(wc.sample, SAMPLESHEET),
            mges = OUTPUT_DIR + "/summary/all_mges.tsv.gz",
            arrays = OUTPUT_DIR + "/summary/all_crispr_arrays.tsv.gz"
        output:
            gff = OUTPUT_DIR + "/targets/{sample}/spacer_hits.gff",
            suspicious = OUTPUT_DIR + "/targets/{sample}/suspicious_hits.gff"
        params:
            genome_id = "{sample}",
            min_identity = config.get("target_mapping", {}).get("min_spacer_identity", 90),
            min_coverage = config.get("target_mapping", {}).get("min_spacer_coverage", 100),
            max_mismatches = config.get("target_mapping", {}).get("max_mismatches", 2)
        log:
            OUTPUT_DIR + "/logs/{sample}/map_spacer_targets.log"
        conda:
            "../envs/base.yaml"
        threads: config.get("target_mapping", {}).get("threads", 4)
        script:
            "../scripts/targets/map_spacer_targets.py"


    rule map_rm_targets:
        """Map RM recognition sites to genome."""
        input:
            patterns_tsv = OUTPUT_DIR + "/patterns/patterns_rm.tsv",
            genome = lambda wc: genome_for(wc.sample, SAMPLESHEET),
            mges = OUTPUT_DIR + "/summary/all_mges.tsv.gz"
        output:
            gff = OUTPUT_DIR + "/targets/{sample}/rm_hits.gff"
        params:
            genome_id = "{sample}"
        log:
            OUTPUT_DIR + "/logs/{sample}/map_rm_targets.log"
        conda:
            "../envs/base.yaml"
        script:
            "../scripts/targets/map_rm_targets.py"


    rule consolidate_hits:
        """Merge spacer and RM hits into mobile/genomic outputs."""
        input:
            spacer_gff = OUTPUT_DIR + "/targets/{sample}/spacer_hits.gff",
            rm_gff = OUTPUT_DIR + "/targets/{sample}/rm_hits.gff"
        output:
            mobile = OUTPUT_DIR + "/targets/{sample}/mobile_hits.gff",
            genomic = OUTPUT_DIR + "/targets/{sample}/genomic_hits.gff"
        log:
            OUTPUT_DIR + "/logs/{sample}/consolidate_hits.log"
        conda:
            "../envs/base.yaml"
        script:
            "../scripts/targets/consolidate_hits.py"


    rule build_target_matrices:
        """Build cross-genome summary matrices."""
        input:
            expand(OUTPUT_DIR + "/targets/{sample}/mobile_hits.gff", sample=SAMPLES)
        output:
            crispr_vs_mge = OUTPUT_DIR + "/targets/matrices/crispr_vs_mge.tsv",
            rm_vs_mge = OUTPUT_DIR + "/targets/matrices/rm_vs_mge.tsv",
            crispr_vs_genome = OUTPUT_DIR + "/targets/matrices/crispr_vs_genome.tsv",
            rm_vs_genome = OUTPUT_DIR + "/targets/matrices/rm_vs_genome.tsv"
        params:
            samples = SAMPLES,
            output_dir = OUTPUT_DIR
        log:
            OUTPUT_DIR + "/logs/build_target_matrices.log"
        conda:
            "../envs/base.yaml"
        script:
            "../scripts/targets/build_target_matrices.py"

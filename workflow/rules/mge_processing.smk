# MGE processing workflow
# Parses MGE GFF files and clusters via recombinase nucleotide sequences

RETAIN_INTERMEDIATES = config.get("retain_intermediates", False)

if config.get("run_mge"):
    rule parse_mges:
        """Parse mge_islands.gff3 to standardized TSV and extract sequences."""
        input:
            mge_gff = lambda wc: mges_for(wc.sample, SAMPLESHEET),
            genome = lambda wc: genome_for(wc.sample, SAMPLESHEET),
            # FFN for recombinase nucleotide sequences
            ffn = lambda wc: ffn_for(wc.sample, SAMPLESHEET, OUTPUT_DIR)
        output:
            mges_tsv = OUTPUT_DIR + "/{sample}/final/mges.tsv",
            mges_fasta = (OUTPUT_DIR + "/{sample}/intermediates/parsed/mges.fasta") if RETAIN_INTERMEDIATES else temp(OUTPUT_DIR + "/{sample}/parsed/mges.fasta")
        params:
            genome_id = "{sample}"
        log:
            OUTPUT_DIR + "/logs/{sample}/parse_mges.log"
        conda:
            "../envs/base.yaml"
        script:
            "../scripts/parsers/parse_mges.py"


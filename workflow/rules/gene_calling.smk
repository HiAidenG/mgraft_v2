# Gene calling with Pyrodigal
# Only runs when protein FASTAs are not provided in samplesheet


rule pyrodigal:
    """Run Pyrodigal to predict genes from genome FASTA."""
    input:
        genome = lambda wc: genome_for(wc.sample, SAMPLESHEET)
    output:
        proteins = OUTPUT_DIR + "/{sample}/tools-raw-out/pyrodigal/{sample}.faa",
        genes = OUTPUT_DIR + "/{sample}/tools-raw-out/pyrodigal/{sample}.ffn",
        gff = OUTPUT_DIR + "/{sample}/tools-raw-out/pyrodigal/{sample}.gff"
    log:
        OUTPUT_DIR + "/logs/{sample}/pyrodigal.log"
    conda:
        "../envs/pyrodigal.yaml"
    threads: 4
    resources:
        mem_mb = 4000,
        runtime = 30
    shell:
        """
        mkdir -p $(dirname {output.proteins})
        pyrodigal -i {input.genome} \
            -j {threads} \
            -a {output.proteins} \
            -d {output.genes} \
            -o {output.gff} \
            -f gff 2> {log}
        """

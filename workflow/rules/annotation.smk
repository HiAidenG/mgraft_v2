# Annotation workflow
# Runs primary annotation tools and parsers

# Environments
DEFENSEFINDER_ENV = config.get("defensefinder", {}).get("conda_env", "../envs/defensefinder.yaml")
CRISPRDETECT_ENV = config.get("crispr", {}).get("conda_env", "../envs/crisprdetect.yaml")
DIAMOND_ENV = config.get("diamond", {}).get("conda_env", "../envs/diamond.yaml")

# Config option for retaining intermediates
RETAIN_INTERMEDIATES = config.get("retain_intermediates", False)

# Helper for intermediate paths
def intermediate_dir(sample):
    if RETAIN_INTERMEDIATES:
        return OUTPUT_DIR + "/" + sample + "/intermediates"
    else:
        return OUTPUT_DIR + "/" + sample

# --- DefenseFinder ---

rule defensefinder:
    """Run DefenseFinder to identify defense systems."""
    input:
        proteins = lambda wc: proteins_for(wc.sample, SAMPLESHEET, OUTPUT_DIR)
    output:
        systems = OUTPUT_DIR + "/{sample}/tools-raw-out/defensefinder/defense_finder_systems.tsv",
        genes = OUTPUT_DIR + "/{sample}/tools-raw-out/defensefinder/defense_finder_genes.tsv",
        hmmer = OUTPUT_DIR + "/{sample}/tools-raw-out/defensefinder/defense_finder_hmmer.tsv"
    params:
        outdir = directory(OUTPUT_DIR + "/{sample}/tools-raw-out/defensefinder")
    log:
        OUTPUT_DIR + "/logs/{sample}/defensefinder.log"
    conda:
        DEFENSEFINDER_ENV
    threads: config.get("defensefinder", {}).get("threads", 4)
    resources:
        mem_mb = config.get("defensefinder", {}).get("mem_mb", 4000),
        runtime = config.get("defensefinder", {}).get("runtime", 60)
    shell:
        """
        defense-finder run -a --out-dir {params.outdir} --workers {threads} {input.proteins} 2> {log}
        
        # DefenseFinder outputs files with sample prefix, rename to generic
        sample=$(basename {input.proteins} .faa)
        for file in {params.outdir}/*_systems.tsv; do
            if [ -f "$file" ]; then mv "$file" {output.systems}; fi
        done
        for file in {params.outdir}/*_genes.tsv; do
             if [ -f "$file" ]; then mv "$file" {output.genes}; fi
        done
        for file in {params.outdir}/*_hmmer.tsv; do
             if [ -f "$file" ]; then mv "$file" {output.hmmer}; fi
        done
        """

rule parse_defensefinder:
    """Parse DefenseFinder output into normalized TSVs."""
    input:
        systems = OUTPUT_DIR + "/{sample}/tools-raw-out/defensefinder/defense_finder_systems.tsv",
        genes = OUTPUT_DIR + "/{sample}/tools-raw-out/defensefinder/defense_finder_genes.tsv",
        hmmer = OUTPUT_DIR + "/{sample}/tools-raw-out/defensefinder/defense_finder_hmmer.tsv",
        proteins = lambda wc: proteins_for(wc.sample, SAMPLESHEET, OUTPUT_DIR)
    output:
        cas_systems = OUTPUT_DIR + "/{sample}/intermediates/parsed/cas_systems.tsv" if RETAIN_INTERMEDIATES else temp(OUTPUT_DIR + "/{sample}/parsed/cas_systems.tsv"),
        cas_genes = OUTPUT_DIR + "/{sample}/intermediates/parsed/cas_genes.tsv" if RETAIN_INTERMEDIATES else temp(OUTPUT_DIR + "/{sample}/parsed/cas_genes.tsv"),
        rm_systems = OUTPUT_DIR + "/{sample}/intermediates/parsed/rm_systems.tsv" if RETAIN_INTERMEDIATES else temp(OUTPUT_DIR + "/{sample}/parsed/rm_systems.tsv"),
        rm_genes = OUTPUT_DIR + "/{sample}/intermediates/parsed/rm_genes.tsv" if RETAIN_INTERMEDIATES else temp(OUTPUT_DIR + "/{sample}/parsed/rm_genes.tsv"),
        defense_systems = OUTPUT_DIR + "/{sample}/intermediates/parsed/defense_systems.tsv" if RETAIN_INTERMEDIATES else temp(OUTPUT_DIR + "/{sample}/parsed/defense_systems.tsv"),
        defense_genes = OUTPUT_DIR + "/{sample}/intermediates/parsed/defense_genes.tsv" if RETAIN_INTERMEDIATES else temp(OUTPUT_DIR + "/{sample}/parsed/defense_genes.tsv"),
        antidefense_systems = OUTPUT_DIR + "/{sample}/intermediates/parsed/antidefense_systems.tsv" if RETAIN_INTERMEDIATES else temp(OUTPUT_DIR + "/{sample}/parsed/antidefense_systems.tsv"),
        antidefense_genes = OUTPUT_DIR + "/{sample}/intermediates/parsed/antidefense_genes.tsv" if RETAIN_INTERMEDIATES else temp(OUTPUT_DIR + "/{sample}/parsed/antidefense_genes.tsv")
    params:
        genome_id = "{sample}"
    log:
        OUTPUT_DIR + "/logs/{sample}/parse_defensefinder.log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/parsers/parse_defensefinder.py"

# --- CRISPRDetect ---

CRISPRDETECT_PATH = config.get("crispr", {}).get("crisprdetect_path", "workflow/resources/crisprdetect/CRISPRDetect_3.0/CRISPRDetect3")
CRISPRDETECT_BATCH_SIZE = config.get("crispr", {}).get("batch_size", 200)

rule crisprdetect:
    """Run CRISPRDetect with optional batching for fragmented genomes."""
    input:
        genome = lambda wc: genome_for(wc.sample, SAMPLESHEET)
    output:
        txt = OUTPUT_DIR + "/{sample}/tools-raw-out/crisprdetect/output",
        gff = OUTPUT_DIR + "/{sample}/tools-raw-out/crisprdetect/output.gff"
    params:
        outdir = directory(OUTPUT_DIR + "/{sample}/tools-raw-out/crisprdetect"),
        output_prefix = OUTPUT_DIR + "/{sample}/tools-raw-out/crisprdetect/output",
        quality_cutoff = config.get("crispr", {}).get("quality_cutoff", 0),
        min_repeats = config.get("crispr", {}).get("min_repeats", 2),
        crisprdetect_path = CRISPRDETECT_PATH,
        batch_size = CRISPRDETECT_BATCH_SIZE
    log:
        OUTPUT_DIR + "/logs/{sample}/crisprdetect.log"
    conda:
        CRISPRDETECT_ENV
    threads: config.get("crispr", {}).get("threads", 4)
    shell:
        """
        mkdir -p {params.outdir}
        
        # Count contigs
        n_contigs=$(grep -c "^>" {input.genome} || echo 0)
        
        if [ "{params.batch_size}" -eq 0 ] || [ "$n_contigs" -le "{params.batch_size}" ]; then
            # No batching needed
            perl {params.crisprdetect_path} \\
                -f {input.genome} \\
                -o {params.output_prefix} \\
                -array_quality_score_cutoff {params.quality_cutoff} \\
                -minimum_no_of_repeats {params.min_repeats} \\
                -T {threads} 2> {log}
        else
            # Batch processing
            echo "Processing $n_contigs contigs in batches of {params.batch_size}" >> {log}
            
            batch_dir="{params.outdir}/batches"
            mkdir -p "$batch_dir"
            
            # Split genome into batch FASTAs using awk
            awk -v batch_size={params.batch_size} -v outdir="$batch_dir" '
                BEGIN {{ batch=1; count=0; outfile=outdir"/batch_"batch".fa" }}
                /^>/ {{
                    count++
                    if (count > batch_size) {{
                        close(outfile)
                        batch++
                        count=1
                        outfile=outdir"/batch_"batch".fa"
                    }}
                }}
                {{ print >> outfile }}
            ' {input.genome}
            
            # Process each batch
            batch_count=0
            batch_success=0
            batch_failed=0
            
            for batch_fa in "$batch_dir"/batch_*.fa; do
                batch_name=$(basename "$batch_fa" .fa)
                batch_out="$batch_dir/$batch_name"
                batch_count=$((batch_count + 1))
                
                echo "Processing $batch_name..." >> {log}
                # Note: CRISPRDetect writes main report to 'output_prefix' (no extension)
                # and GFF to 'output_prefix.gff'
                perl {params.crisprdetect_path} \\
                    -f "$batch_fa" \\
                    -o "$batch_out" \\
                    -array_quality_score_cutoff {params.quality_cutoff} \\
                    -minimum_no_of_repeats {params.min_repeats} \\
                    -T {threads} 2>> {log}
                batch_exit=$?
                
                if [ $batch_exit -ne 0 ]; then
                    echo "WARNING: CRISPRDetect failed for $batch_name (exit code $batch_exit)" >> {log}
                    batch_failed=$((batch_failed + 1))
                else
                    batch_success=$((batch_success + 1))
                fi
            done
            
            echo "Batch processing complete: $batch_success/$batch_count succeeded, $batch_failed failed" >> {log}
            
            # Merge outputs
            echo "Merging batch outputs..." >> {log}
            
            # Merge txt outputs: find files that are NOT .fa or .gff (main reports have no extension)
            > {output.txt}
            merged_count=0
            for batch_file in "$batch_dir"/batch_*; do
                # Skip .fa and .gff files - only process main report files (no extension)
                if [[ "$batch_file" =~ \.(fa|gff)$ ]]; then
                    continue
                fi
                if [ -f "$batch_file" ] && [ -s "$batch_file" ]; then
                    echo "Merging main report: $(basename $batch_file)" >> {log}
                    cat "$batch_file" >> {output.txt}
                    merged_count=$((merged_count + 1))
                fi
            done
            echo "Merged $merged_count main report files" >> {log}
            
            # Merge GFF outputs (combine headers + features)
            echo "##gff-version 3" > {output.gff}
            gff_features=0
            for gff in "$batch_dir"/batch_*.gff; do
                if [ -f "$gff" ] && [ -s "$gff" ]; then
                    feature_count=$(grep -v "^#" "$gff" 2>/dev/null | wc -l || echo 0)
                    gff_features=$((gff_features + feature_count))
                    grep -v "^#" "$gff" >> {output.gff} 2>/dev/null || true
                fi
            done
            echo "Merged $gff_features GFF features total" >> {log}
            
            # Validation: warn if output is empty but batches ran
            if [ ! -s {output.txt} ] && [ $batch_success -gt 0 ]; then
                echo "WARNING: Main output is empty despite $batch_success successful batches" >> {log}
            fi
            
            # Cleanup batch files
            rm -rf "$batch_dir"
        fi
        
        # Ensure outputs exist
        touch {output.txt} {output.gff}
        """

rule parse_crisprdetect:
    """Parse CRISPRDetect output."""
    input:
        txt = OUTPUT_DIR + "/{sample}/tools-raw-out/crisprdetect/output",
        gff = OUTPUT_DIR + "/{sample}/tools-raw-out/crisprdetect/output.gff"
    output:
        arrays = OUTPUT_DIR + "/{sample}/intermediates/parsed/crispr_arrays.tsv" if RETAIN_INTERMEDIATES else temp(OUTPUT_DIR + "/{sample}/parsed/crispr_arrays.tsv"),
        spacers = OUTPUT_DIR + "/{sample}/intermediates/parsed/crispr_spacers.tsv" if RETAIN_INTERMEDIATES else temp(OUTPUT_DIR + "/{sample}/parsed/crispr_spacers.tsv")
    params:
        genome_id = "{sample}"
    log:
        OUTPUT_DIR + "/logs/{sample}/parse_crisprdetect.log"
    conda:
        "../envs/base.yaml"
    script:
        "../scripts/parsers/parse_crisprdetect.py"

# --- REBASE Alignment ---

rule extract_rm_proteins:
    """Extract proteins for RM genes for REBASE alignment."""
    input:
        rm_genes = (OUTPUT_DIR + "/{sample}/intermediates/parsed/rm_genes.tsv") if RETAIN_INTERMEDIATES else (OUTPUT_DIR + "/{sample}/parsed/rm_genes.tsv"),
        proteins = lambda wc: proteins_for(wc.sample, SAMPLESHEET, OUTPUT_DIR)
    output:
        faa = temp(OUTPUT_DIR + "/{sample}/tmp/rm_genes.faa")
    conda:
        "../envs/base.yaml"
    shell:
        """
        mkdir -p $(dirname {output.faa})
        # Extract prodigal_id column (column 11) which contains the protein IDs
        awk -F'\\t' 'NR>1 && $11!="NULL" && $11!="" {{print $11}}' {input.rm_genes} > {output.faa}.ids
        if [ -s {output.faa}.ids ]; then
            seqkit grep -f {output.faa}.ids {input.proteins} -o {output.faa}
        else
            touch {output.faa}
        fi
        rm -f {output.faa}.ids
        """

rule diamond_rebase:
    """Align RM proteins to REBASE and generate QC stats."""
    input:
        query = OUTPUT_DIR + "/{sample}/tmp/rm_genes.faa",
        db = config["rebase_db"]
    output:
        tsv = OUTPUT_DIR + "/{sample}/intermediates/rebase-hits-raw.tsv" if RETAIN_INTERMEDIATES else temp(OUTPUT_DIR + "/{sample}/rebase-hits-raw.tsv"),
        stats = OUTPUT_DIR + "/{sample}/rebase-alignment-stats.json"
    conda:
        DIAMOND_ENV
    threads: config.get("rebase", {}).get("threads", 4)
    params:
        evalue = config.get("rebase", {}).get("evalue", 1e-10),
        id = config.get("rebase", {}).get("identity", 80),
        cov = config.get("rebase", {}).get("coverage", 80)
    shell:
        """
        if [ -s {input.query} ]; then
            diamond blastp \
                --db {input.db} \
                --query {input.query} \
                --out {output.tsv} \
                --threads {threads} \
                --evalue {params.evalue} \
                --id {params.id} \
                --query-cover {params.cov} \
                --subject-cover {params.cov} \
                --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
                --max-target-seqs 1 \
                --sensitive
            
            # Generate QC stats
            genes_queried=$(grep -c ">" {input.query} || echo 0)
            hits_found=$(wc -l < {output.tsv} || echo 0)
            if [ "$hits_found" -gt 0 ]; then
                min_identity=$(awk -F'\\t' 'NR==1{{min=$3}} $3<min{{min=$3}} END{{print min}}' {output.tsv})
                avg_identity=$(awk -F'\\t' '{{sum+=$3; n++}} END{{print sum/n}}' {output.tsv})
                max_evalue=$(awk -F'\\t' 'NR==1{{max=$11}} $11>max{{max=$11}} END{{print max}}' {output.tsv})
            else
                min_identity=0
                avg_identity=0
                max_evalue=0
            fi
            echo "{{\\"genes_queried\\": $genes_queried, \\"hits_found\\": $hits_found, \\"min_identity\\": $min_identity, \\"avg_identity\\": $avg_identity, \\"max_evalue\\": $max_evalue}}" > {output.stats}
        else
            touch {output.tsv}
            echo "{{\\"genes_queried\\": 0, \\"hits_found\\": 0, \\"min_identity\\": 0, \\"avg_identity\\": 0, \\"max_evalue\\": 0}}" > {output.stats}
        fi
        """

# mgraft_v2 Updates

## Type-Specific RM Recognition Sequence Assignment

Updated the REBASE recognition sequence propagation logic to use type-specific validation rules instead of requiring all genes to have matching target sequences.

### Logic by RM Type

- **Type I**: Requires R, M, and S subunits present. Only the S subunit determines the recognition sequence. Multiple S genes in the same system with different target sequences results in `conflicting`.

- **Type II / Type III**: Requires at least one cognate R+M pair with matching recognition sequences. Multiple valid pairs with the same sequence is acceptable. However, if multiple pairs have different sequences, the system is marked as `conflicting`. Genes with missing REBASE annotations are permitted as long as at least one valid pair exists.

- **Type IIG**: All genes are bifunctional RM. Annotated if all genes in the system share the same recognition sequence; `conflicting` if they differ. Genes with missing REBASE annotations are permitted as long as at least one gene has an annotation.

- **Type IV**: No recognition sequence exists for this type. Always assigned `no_recognition_site`.

### Changes to recseq_note Values

- `annotated`: Required genes have matching REBASE recognition sequences
- `conflicting`: Required genes have recognition sites but they don't match
- `missing_annotation`: Required genes are missing or lack REBASE annotations
- `no_recognition_site`: Type IV systems

This is technically less stringent than before, which required all genes within a DefenseFinder system to have the same recognition sequence. However, it better reflects biological reality and avoids penalizing systems with accessory components.

## CRISPR Calling and Array Assignment
minCED has been replaced with CRISPRDetect for calling CRISPR arrays. CRISPRDetect offers improved accuracy and additional features such as array orientation and confidence scoring. I'm unsure if this improved accuracy comes at the cost of lower sensitivity - have to look into whether anyone has benchmarked this. 

CRISPRDetect is run without cas calling via BLASTp since we already have cas genes from DefenseFinder. Instead, mGRAFT manually cross-references arrays with nearby cas genes - assigning arrays to cas systems based on a configurable distance thresholds (`max_distance_to_cas`, default 10kb). mGRAFT will search for the nearest cas gene within this distance ON BOTH STRANDS. If found, the array is assigned to that cas system; otherwise, it is marked `is_orphan: true`. Since CRISPRDetect does not have cas information when determining directionality mGRAFT will, by default, reverse-complement arrays that are on the opposite strand of their assigned cas system, if applicable. This is configurable via the `revcomp_if_cas_opposite` flag in the config. 

If an array (spacers and repeats) is reverse-complemented by mGRAFT, the `is_revcomp` attribute is set to true in the output GFF3 and TSV files for that array.

If CRIPRDetect cannot determine the strand of an array (i.e. strand is `.`), mGRAFT will attempt to infer the correct orientation by locating the repeat sequence in the genome using `seqkit locate`. The strand is determined based on where the repeat sequence is found. If this inferred strand is opposite to that of the assigned cas system, and revcomp is enabled, the array will be reverse-complemented accordingly.

In addition to the `is_orphan` attribute, arrays are annotated as either `canonical` or `putative` based on whether they are assigned to a cas system and the completeness of that system:
 - `canonical`: Assigned to a complete cas system (all required cas genes present)
 - `putative`: Not assigned to a cas system, or assigned to an incomplete cas system

 mGRAFT will only ever assign one array per cas system. If multiple arrays are found within the distance threshold of the same cas system, the nearest array is assigned and the others are marked as orphans.
 For tie breakers (i.e. two arrays equidistant from the same cas system), the array on the same strand as the cas system is preferred. If both are on the same strand, the one with the higher quality score (as reported by CRISPRDetect) is chosen.


## Defense System Mobility Annotation
Added `mobility` metric to defense_system features, calculated as the proportion of genes in the system located on the same MGE. This is defined as the largest possible proportion of genes in the system that share the same `mge_id`. For example, if a defense system has 5 genes, with 3 located on MGE_A and 2 on MGE_B, the mobility would be 0.6 (3/5). If all genes are on the same MGE, mobility is 1.0; if none are on an MGE, mobility is 0.0.

## MGE Clustering
- Detected MGEs are clustered for the purpose of deduplication. Each MGE is assigned a `mge_cluster_id` based on clustering results. Clustering is performed using vclust:
- prefilter with at least 95% identity over the shorter sequence: `vclust prefilter -i mge_sequences.fasta -o vclust_mge_prefilter.txt -id 0.95`
- align filtered pairs with vclust lz-ani: `vclust align -i mge_sequences.fasta -o vclust_mge_ani.tsv --outfmt lite --filter vclust_mge_prefilter.txt --threads 64`
- cluster with leiden algorithm: `vclust cluster -i vclust_mge_ani.tsv --qcov 0.85 --rcov 0.85 --ids mge_sequences.ids.tsv -o vclust_mge_clusters.tsv --metric ani --method leiden --ani 0.95 --min_cluster_size 2`

Defense system and CRISPR array annotations are still handled on a per-MGE basis, but deduplicated MGEs are used for downstream analyses and plots.

Notes:
- Will need to check the concordence of defense calls, crispr arrays, etc.. within clusters.

#
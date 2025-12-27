# CRISPR Clustering Module

This module clusters CRISPR arrays based on repeat and spacer sequences using MMseqs2.

## Clustering Parameters

### Repeat Clustering
- **Identity**: 100% (exact match)
- **Coverage**: 100% (full length)
- **Method**: MMseqs2 linclust (linear-time, optimized for high identity)

### Spacer Clustering
- **Identity**: 90%
- **Coverage**: 90%
- **Method**: MMseqs2 easy-cluster (greedy set cover)

## Handling Reverse Complements

Both repeat and spacer sequences are **canonicalized** before clustering:
- For each sequence, we compute both the forward and reverse complement
- The lexicographically smaller sequence is used as the canonical form
- This ensures that sequences that are reverse complements of each other cluster together

## Usage

### Command Line

```bash
python modules/crispr_clustering/cluster_crispr.py \
    --arrays summary/all_crispr_arrays.tsv \
    --spacers summary/all_crispr_spacers.tsv \
    --output-dir summary/crispr_clustering \
    --threads 4
```

### With Custom Parameters

```bash
python modules/crispr_clustering/cluster_crispr.py \
    --arrays summary/all_crispr_arrays.tsv \
    --spacers summary/all_crispr_spacers.tsv \
    --output-dir summary/crispr_clustering \
    --threads 8 \
    --repeat-identity 1.0 \
    --repeat-coverage 1.0 \
    --spacer-identity 0.9 \
    --spacer-coverage 0.9
```

### Snakemake Integration

Add to your workflow:

```python
include: "modules/crispr_clustering/crispr_clustering.smk"
```

Then run:

```bash
snakemake cluster_crispr_all --cores 4
```

## Output Files

| File | Description |
|------|-------------|
| `crispr_arrays_clustered.tsv` | Original arrays with `repeat_cluster_id` column |
| `crispr_spacers_clustered.tsv` | Original spacers with `spacer_cluster_id` column |
| `repeat_cluster_info.tsv` | Repeat cluster details (representative, size, members) |
| `spacer_cluster_info.tsv` | Spacer cluster details |
| `repeat_sequences.fasta` | Canonicalized repeat sequences |
| `spacer_sequences.fasta` | Canonicalized spacer sequences |
| `clustering_summary.txt` | Summary statistics |

## Analysis

Run additional analysis on clustering results:

```bash
python modules/crispr_clustering/analyze_clusters.py \
    --input-dir summary/crispr_clustering \
    --output-dir summary/crispr_clustering/analysis \
    --min-shared-spacers 3
```

This generates:
- `shared_spacers_across_genomes.tsv` - Spacer clusters found in multiple genomes
- `repeat_diversity.tsv` - Repeat cluster distribution
- `related_arrays.tsv` - Array pairs sharing spacer clusters
- `array_similarity_network.tsv` - Network file for visualization

## Cluster ID Format

- **Repeat clusters**: `RPT000001`, `RPT000002`, ...
- **Spacer clusters**: `SPC000001`, `SPC000002`, ...

Clusters are numbered by size (largest cluster = lowest number).

## Dependencies

- Python 3.10+
- pandas
- numpy
- MMseqs2 (install via conda: `conda install -c bioconda mmseqs2`)

## Conda Environment

```bash
conda env create -f envs/mmseqs.yaml
conda activate mgraft_v2_mmseqs
```

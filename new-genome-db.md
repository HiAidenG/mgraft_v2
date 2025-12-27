# UHGG problems
UHGG is being replaced in the new analysis. The reason for this is that UHGG contains genomes solely based on metadata annotations, which are often incorrect. For example, isolate Neisseria gonorrhoeae genomes from 'rectum' samples are not true human-gut associated genomes, but are included in UHGG because of their metadata source.

The solution to this is to analyze the prevalence of prokaryotic genomes across samples. This ensures that only genomes present in a sufficiently high number of samples will be included - thereby ruling out contaminants or genomes included solely based on metadata. Moreover, the presence of a MAG for a given species increases confidence that the species is truly present in the human gut, as MAGs are assembled directly from the community. Isolates, meanwhile, may have their source misannotated.

This is exactly what's done in this paper: https://link.springer.com/article/10.1186/s40168-021-01114-w. 
So, I'll be using this going forward. Still going to apply the same minimum 50% completion and maximum 5% contamination filters as before. 
Notably, Neisseria gonorrhoeae is not present in this set.

## Genome Selection
- Inferred pg cluster representatives for all humgut genomes based on taxonomy 
- of the 31,226 genomes in HumGut, 17,345 matched directly, 1057 matched a species in pg (not necessarily the exact strain), 12,823 did not map at all.
- of the 18,402 genomes that mapped to pg clusters, 16,493 had >= 2 representatives in pg clusters. 

Further filtered this to only include 2 genomes per humgut cluster - prioritizing an isolate and a MAG, where available. Then prioritizing highest completeness and lowest contamination.
This resulted in a final set of 3594 genomes, corresponding to 2392 humgut95 clusters. These mapped to 1233 pg clusters. While both are clustered at estimated 95% ANI based on MASH distance, this discrepancy may be due to the fact that HumGut95 clustered based on complete linkage, while pg clustered on single linkage.

Moreover, the corresponding pg cluster was only inferred based on taxonomy, not actual genome distances. (???) But for species complexes, it should be a complex in both humgut and pg? 
In some cases, these may actually map to different clusters in pg when considering actual genome distances. However, in cases where there are highly promiscuous pg clusters and relatively few pg isolates for that cluster, it's likely to cause problematic MGE calls. 
- pg: ANI_pg4_8165, tax: Phascolarctobacterium_A succinatutens, maps to 51 humgut clusters, 3 representatives in pg - remove
There are others but not as extreme. For the rest we just have to check after mge calling to exclude problematic complexes.

So now: 3530 genomes corresponding to 2391 humgut95 clusters and 1232 pg clusters. 





microbeatlas


# Type I Systems
MGYG000002349	rm_gene	MGYG000002349__RMS0003_g01	MGYG000002349__RMS0003	MGYG000002349_64	2	937	-	62.8	S	MGYG000002349__MGE0014	MGYG000002349_64_1	NULL	NULL	NULL	NULL	NULL	NULL	NULL	NULL
MGYG000002349	rm_gene	MGYG000002349__RMS0003_g02	MGYG000002349__RMS0003	MGYG000002349_64	937	2400	-	590.4	M	MGYG000002349__MGE0014	MGYG000002349_64_2	100.0	99.8	100.0	0.0	M.Sag515CII	Streptococcus agalactiae 515	GACNNNNNNNCTT	PRJNA594846:GRB
MGYG000002349	rm_gene	MGYG000002349__RMS0003_g03	MGYG000002349__RMS0003	MGYG000002349_64	2413	4737	-	767.3	R	MGYG000002349__MGE0014	MGYG000002349_64_3	100.0	99.9	100.0	0.0	Sag515CIIP	Streptococcus agalactiae 515	GACNNNNNNNCTT	PRJNA594846:GRB

Missing recognition site for the S subunit, but known for R and M subunits.
I think in this case the ORF is probably truncated or incomplete, leading to missing annotation in REBASE - HMM hit score for S is quite low (62.8) compared to others. 

No, this is just a diverged S unit, can't know recognition site without experimental validation.
Coverage is 100% on query and target, so not truncated. Identity is only 85.6%. 

# Thresholding
MAGs can be assemblages of multiple strains present in the same sample, leading to potential variation within the assembly. This can result in defense system genes having mutations that affect their annotation.
Should we consider going to like 95% identity and 90% coverage for defense system gene annotation to account for this potential variation?

# Spacer problems
- Some spacers are very short (e.g. 4nt) others are very long (~140 nt). Should we set min/max length thresholds for spacers to be included in clustering?
- Filter out low entropy?



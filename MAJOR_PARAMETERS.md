# Major parameters and their provenance

These entries distinguish supplied-code settings from manuscript statements. They are not a substitute for the original run logs.

| Pipeline | Parameters recorded in the supplied files or manuscript |
|---|---|
| Study inputs | UC, BD, MDD, SCZ; manuscript cases 12366/41917/170756/52017; controls 33609/371549/329443/75889 |
| GWAS preparation | PLACO and scRNA preparation set maf=0.01; MetaXcan preparation removes chr6 positions strictly between 25 and 35 Mb and leaves the loader's MAF default; inspect downstream filters separately |
| LDSC/HDL | UC as target; matching European references; HDL reference sample size in supplied code 335265; retain actual per-input sample-size fields |
| LAVA | phenos UC/BD/MDD/SCZ, target UC, hg19 locus file; supplied loop hard-coded to 2495 intervals; univariate threshold inherited from installed LAVA defaults unless specified |
| QTLEnrich | Separate eQTL and sQTL runs, GTEx v8; inspect calls for permutations and matched-null resources. Manuscript states BH across 392 tests at q<0.005, but supporting tables contain a conflicting family definition |
| gsMap | Supplied quick-mode script and reference inputs; manuscript describes E16.5 mouse atlas and Bonferroni across 25 organs; supplementary tables also mention BH |
| eCAVIAR | Manuscript: r2>0.1 region plus 50 kb, up to 2 causal variants, CLPP>0.01; verify source code and coordinate convention against final logs |
| fastENLOC | Manuscript: DAP-G fine mapping, same window convention, RCP>0.1; historical DAP-G version not recorded |
| PLACO | UC against BD/MDD/SCZ; manuscript threshold 5e-8; check actual wrapper defaults and FUMA task configurations |
| condFDR/conjFDR | UC as target; three comparators; bidirectional conditional estimates; manuscript conjFDR<0.05, MAF>=0.005, 500 random prunings; supplied wrappers also contain explicit region exclusions |
| scPagwas | Supplied code: 24 cores, RNA assay, cell_type identities, Genes_by_pathway_kegg, singlecell=TRUE, celltype=TRUE, iters_singlecell=100; derive maf=min(eaf,1-eaf); block_annotation and chrom_ld resources; no explicit RNG seed in the supplied script |
| Candidate genes | KIAA1109, MED24, RHOA, HSD17B8; equal donor weighting for expression; >=20 cells/donor-cell type and >=3 donors/condition for DE; raw counts; filterByExpr; TMM; voom; condition-only model; robust empirical Bayes; moderated-t 95% CI; cohort-wise BH over estimable candidate comparisons |

Keep exact database files, annotation versions, external executables and random seeds with future run records. Unspecified or conflicting settings have not been invented or silently reconciled here.

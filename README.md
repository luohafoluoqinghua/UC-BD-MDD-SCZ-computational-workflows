# UC and psychiatric disorders computational materials

Code and processed inputs accompanying **Integrated Genetic and Single-Cell Analyses Reveal the Shared Genetic Architecture Between Ulcerative Colitis and Three Psychiatric Disorders**.

Traits: ulcerative colitis (UC), bipolar disorder (BD), major depressive disorder (MDD), and schizophrenia (SCZ). Pairwise comparisons use UC as the target and BD, MDD, and SCZ as comparators.

## Contents and scope

This version (0.1.0, 2026-10-21) contains 31 supplied R scripts with descriptive filenames, a software-version inventory, parameter documentation, and a checksum manifest for four processed GWAS tables and three Seurat objects. The scripts were configured for these traits and renamed during revision preparation. This is a documented code/input deposit, not a claim that every script reproduces the final manuscript without further configuration. The computational pipelines were not executed for this deposit.

| Location | Analysis |
|---|---|
| `MetaXcan/` | GWAS preparation, Manhattan/QQ plots, condFDR/conjFDR; the retained legacy folder name does not mean a MetaXcan TWAS is implemented here |
| `PLACO/` | GWAS preparation, LDSC, HDL, PLACO |
| `scRNA/` | Coordinate conversion, QTLEnrich, FUMA inputs, eCAVIAR, fastENLOC, gsMap; this folder does not contain Seurat preprocessing |
| `01_run_lava_local_genetic_correlation.R` | LAVA local genetic correlation and result aggregation |
| `02_run_scpagwas_celltype_mapping.R` | Supplied UC–GSE214695 scPagwas execution script |
| `03_run_candidate_gene_expression_pseudobulk.R` | Four candidate genes in two UC cohorts; donor-level expression and pseudobulk |

## Input data

The separate `UC_BD_MDD_SCZ_processed_inputs.zip` contains the original seven input files without changing their contents. See `data_manifest.csv` for sizes, source accessions and SHA-256 hashes. The data archive is deposited separately from GitHub. Its publication location is recorded in `DATA_AVAILABILITY.md` when assigned.

GWAS sources: UC ebi-a-GCST004133; BD ieu-b-5110; MDD ieu-b-102; SCZ ieu-b-5102 (IEU OpenGWAS). Single-cell sources: GEO GSE182270, GSE214695 and GSE254569. Upstream studies retain their original attribution and access/use conditions.

The provided GWAS tables have columns `chrom,pos,rsid,effect_allele,other_allele,eaf,beta,se,pval,N`. They are directly formatted for the supplied scPagwas loader; they are not interchangeable with every Fast2TWAS input. Do not infer that every manuscript QC filter was already applied solely from these filenames. Preserve row-specific N: for example, the inspected SCZ first row has N=126849, whereas the sample-index scripts contain manuscript-level case/control totals.

## Execution and original directory layout

1. Inspect `script_name_mapping.csv`. It maps each descriptive name to its original filename and stage-directory name.
2. The legacy wrappers expect separate stage folders, such as `06.pre_data/` and `07.export_index/`, under each analysis project, alongside `03.reference/` and `04.rawdata/`. Place each renamed script in the original stage directory given in the mapping, or adapt its paths consistently. Numeric script names alone do not create those stage directories.
3. The three GWAS preparation scripts expect four downloaded `.vcf.gz` files named after the OpenGWAS accessions in `../04.rawdata/`; these raw VCFs are not included. They output UC.txt, BD.txt, MDD.txt and SCZ.txt in Fast2TWAS format. Reference genotypes, LD scores, GTEx QTLs, liftover chains, gsMap atlas/resources and external executables must be supplied separately.
4. Follow numeric order within each module. Several wrappers use RStudio's active-document directory; others retain server-specific `setwd()` paths. Update these to the deployment directory before execution. Filename changes in this release do not change statistical algorithms or resolve these pre-existing path assumptions.
5. For candidate genes, from the script directory run:

```text
Rscript 03_run_candidate_gene_expression_pseudobulk.R /path/GSE182270.rds /path/GSE214695.rds
```

The candidate script accepts two positional input paths, produces donor-level summaries, and writes full-gene pseudobulk counts, sample metadata, design matrices, retained genes and fitted objects under `work/scrna/` during execution. Those generated model objects are not present in this deposit.

The supplied scPagwas script is currently configured for UC with GSE214695. For each other run set GWAS_FILE, SEURAT_FILE, OUTPUT_PREFIX, the working directory and FAST_TEMP_DIR consistently. Four GWAS × three references gives 12 distinct combinations; the deposit does not claim this one hard-coded script is a complete recorded history of all 12 runs.

## Versions, parameters and limitations

See `software_versions.csv`, `candidate_gene_sessionInfo.txt`, `r_packages.txt`, `MAJOR_PARAMETERS.md`, and `REPRODUCIBILITY_STATUS.md`. A recorded candidate-gene session is distinguished from versions merely reported in the manuscript and from unknown versions. Do not substitute current software releases for unrecorded historical versions.

All R files passed syntax parsing during packaging. No GWAS, QTL, single-cell or model-fitting computation was rerun. Missing original analysis records and discrepancies are itemized rather than represented as completed checks.

## Attribution and reuse

Original author and Fast2TWAS template notices are preserved. Public access to this deposit does not grant a new blanket license over third-party code or source datasets. Consult `THIRD_PARTY_NOTICES.md` and the original sources before redistribution or reuse.

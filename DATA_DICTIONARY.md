# Input data dictionary

GWAS text files are tab-delimited with one variant per row. `chrom`: chromosome; `pos`: genomic position (the manuscript describes hg19 inputs; provenance must be checked before combining with hg38 QTL resources); `rsid`: variant identifier; `effect_allele` and `other_allele`: allele orientation; `eaf`: effect-allele frequency; `beta`: supplied effect estimate; `se`: its standard error; `pval`: supplied association P value; `N`: supplied row-level sample size. Preserve source effect scales and alleles; do not assume per-variant N equals a manuscript cohort total.

RDS files are serialized Seurat objects. The candidate-gene script requires an RNA counts assay, `orig.ident` donor/sample identifiers, `Type` values Control/Treat, and `cell_type` annotations. The brain object is used for scPagwas rather than this UC-versus-control pseudobulk script. See `single_cell_object_inventory.csv` for inspected dimensions and metadata column names when available.

Checksums identify the inspected local input files. The seven data files themselves are not included in this public code release; their original filenames are recorded to connect the scripts with input provenance.

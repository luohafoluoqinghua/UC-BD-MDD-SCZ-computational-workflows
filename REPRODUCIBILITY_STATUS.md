# Reproducibility status

## Included

- 31 supplied scripts renamed by function; study labels configured for UC, BD, MDD and SCZ.
- Provenance records, data dictionaries and byte-level checksums for four inspected local GWAS inputs and three single-cell objects. The seven data files themselves are not included in the public release.
- Actual candidate-gene R session record, reported software-version inventory, major-parameter notes and source mapping.

## Still to recover or reconcile

- Original run logs/commands and exact historical environments for most genetic workflows; Fast2TWAS version and external binary versions.
- GNOVA execution code; original single-cell QC/integration/annotation code; brain-reference sampling code and selected-cell records; complete twelve-run scPagwas configuration/history.
- Construction of the gene-by-eight evidence matrix and LOMO code/results are not in the supplied desktop directory and are not included in this release.
- Fine-mapping PIP/signal-cluster files, per-locus LD/Z inputs, matched nulls, complete tested LAVA/QTL result families, FUMA tasks and figure-source outputs are not in the supplied desktop directory.
- Existing path assumptions require deployment configuration. In LA.R (renamed), the save directory and later results directory differ; the script writes a tab-separated .xlsx-named file but later reads a .csv filename. The original aggregation/BH code also differs from the manuscript's within-pair BH and MHC-reporting rule. These pre-existing issues are retained and documented, not represented as a tested final pipeline.
- QTLEnrich BH-family definitions and gsMap correction methods differ between manuscript and supplementary tables. Seurat/Harmony dimensions and hg19/hg38 labels require reconciliation.
- The scPagwas input tables are a particular processed format. Matching variants across disease pairs and QC provenance must be documented before describing them as fully harmonized inputs for all other pipelines.

Packaging validation covers readable files, syntax, consistent naming and checksums. It does not establish numerical reproduction of the manuscript. A later version should close these gaps and bind final run outputs to the manuscript figures and tables.

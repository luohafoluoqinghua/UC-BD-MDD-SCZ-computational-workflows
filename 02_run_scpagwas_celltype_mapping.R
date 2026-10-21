# C drive: original inputs and fast temporary FBM matrices (deleted per pathway).
# D drive: working directory, final results and figures.
# Run in a fresh R session / Rscript process. Moving this script does not change these absolute paths.
dir.create("D:/scPaGWAS/UC_GSE214695_fast_run", recursive = TRUE, showWarnings = FALSE)
setwd("D:/scPaGWAS/UC_GSE214695_fast_run")

#清除一切变量
rm(list=ls())
#释放不再使用的内存
gc()
#显式路径
getwd()
#显示R包安装位置
.libPaths()
# 加载必要的包
library(scPagwas)
library(RColorBrewer)
library(Seurat)
library(ggpubr)
library(tidyverse)
library(rhdf5)
library(ggplot2)
library(grDevices)
library(stats)
library(FactoMineR)
library(scales)
library(reshape2)
library(ggdendro)
library(grImport2)
library(gridExtra)
library(grid)
library(sisal)

# 设置参数
GWAS_FILE <- "C:/Users/1/Desktop/scPaGWAS/UC.txt"
SEURAT_FILE <- "C:/Users/1/Desktop/scPaGWAS/GSE214695.rds"
OUTPUT_PREFIX <- "UC_GSE214695"
OUTPUT_DIR <- "scPagwas_results_run"
N_CORES <- 24L  # 9950X（16核32线程）；本机回归小规模测试采用24线程
# Read the supplied cell-type annotations; scPagwas groups by Seurat identities.
stopifnot(file.exists(GWAS_FILE), file.exists(SEURAT_FILE))
seurat_input <- readRDS(SEURAT_FILE)
if (!inherits(seurat_input, "Seurat") || !("RNA" %in% names(seurat_input@assays))) {
  stop("GSE214695.rds must be a Seurat object containing the RNA assay.")
}
if (!("cell_type" %in% colnames(seurat_input@meta.data))) {
  stop("Missing cell_type annotation in GSE214695.rds.")
}
cell_labels <- as.character(seurat_input$cell_type)
if (anyNA(cell_labels) || any(!nzchar(trimws(cell_labels)))) {
  stop("Missing or empty cell_type labels.")
}
DefaultAssay(seurat_input) <- "RNA"
Idents(seurat_input) <- "cell_type"
CELLTYPE_ORDER <- levels(droplevels(Idents(seurat_input)))
message("Input: UC GWAS + GSE214695.rds; cells: ", ncol(seurat_input))
print(table(Idents(seurat_input)))

# UC.txt provides effect-allele frequency; derive minor-allele frequency.
gwas_header <- names(data.table::fread(GWAS_FILE, nrows = 0L))
required_gwas_cols <- c("chrom", "pos", "rsid", "beta", "se")
missing_gwas_cols <- setdiff(required_gwas_cols, gwas_header)
if (length(missing_gwas_cols)) {
  stop("Missing GWAS columns: ", paste(missing_gwas_cols, collapse = ", "))
}
frequency_col <- if ("maf" %in% gwas_header) "maf" else "eaf"
if (!(frequency_col %in% gwas_header)) stop("GWAS requires maf or eaf.")
gwas_input <- as.data.frame(data.table::fread(
  GWAS_FILE, select = c(required_gwas_cols, frequency_col), nThread = N_CORES
))
for (nm in c("pos", "beta", "se", frequency_col)) {
  gwas_input[[nm]] <- as.numeric(gwas_input[[nm]])
}
freq <- gwas_input[[frequency_col]]
if (any(!is.finite(freq) | freq < 0 | freq > 1)) {
  stop("Invalid GWAS allele frequencies; expected finite values in [0, 1].")
}
gwas_input$maf <- pmin(freq, 1 - freq)
if (frequency_col == "eaf") gwas_input$eaf <- NULL
if (any(!is.finite(gwas_input$beta)) ||
    any(!is.finite(gwas_input$se) | gwas_input$se <= 0) ||
    any(!is.finite(gwas_input$pos) | gwas_input$pos <= 0)) {
  stop("Invalid GWAS beta, se or pos values.")
}
rm(freq, cell_labels)
message("GWAS variants loaded: ", nrow(gwas_input))
options(bigstatsr.ncores.max = N_CORES)
options(bigstatsr.check.parallel.blas = FALSE)
# Windows: release and remove each pathway's temporary matrix after regression.
dir.create(file.path(OUTPUT_DIR, "temp"), recursive = TRUE, showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
# Session-only wrapper: release mapped matrices before removing each pathway file.
install_scpagwas_cleanup <- function(allowed_temp_dir, fast_temp_dir) {
  allowed_temp_dir <- normalizePath(allowed_temp_dir, winslash = '/', mustWork = TRUE)
  fast_temp_dir <- normalizePath(fast_temp_dir, winslash = '/', mustWork = TRUE)
  message('FBM temporary matrices: ', fast_temp_dir)
  ns <- asNamespace('scPagwas')
  original <- get('scParameter_regression', envir = ns)
  if (isTRUE(attr(original, 'windows_ssd_cleanup_fixed'))) return(invisible(NULL))
  fixed <- function(Pagwas_x, Pagwas_y, noise_per_snp, Rns,
                    n.cores = 1, backingpath) {
    parent <- normalizePath(backingpath, winslash = '/', mustWork = TRUE)
    if (!identical(tolower(parent), tolower(allowed_temp_dir)) ||
        length(Rns) != 1L || !grepl('^[A-Za-z0-9_-]+$', Rns)) {
      stop('Unexpected temporary matrix location; no cleanup attempted.')
    }
    targets <- file.path(fast_temp_dir, paste0(Rns, c('.bk', '.rds')))
    if (any(file.exists(targets))) stop('Temporary matrix name already exists: ', Rns)
    on.exit({
      # The original function has returned, so its FBM can now be finalized.
      invisible(gc(full = TRUE))
      present <- targets[file.exists(targets)]
      if (length(present)) unlink(present, recursive = FALSE)
      if (any(file.exists(targets))) {
        stop('Temporary matrix could not be removed; stopping to prevent disk accumulation: ', Rns)
      }
    }, add = TRUE)
    original(Pagwas_x = Pagwas_x, Pagwas_y = Pagwas_y,
             noise_per_snp = noise_per_snp, Rns = Rns,
             n.cores = n.cores, backingpath = fast_temp_dir)
  }
  attr(fixed, 'windows_ssd_cleanup_fixed') <- TRUE
  assignInNamespace('scParameter_regression', fixed, ns = 'scPagwas')
  invisible(NULL)
}

FAST_TEMP_DIR <- "C:/Users/1/Desktop/scPaGWAS/scPagwas_temp/UC_GSE214695_fast_run"
dir.create(FAST_TEMP_DIR, recursive = TRUE, showWarnings = FALSE)
install_scpagwas_cleanup(file.path(OUTPUT_DIR, "temp"), FAST_TEMP_DIR)
message("Working directory: ", getwd(), "; R temporary directory: ", tempdir())

# 运行scPagwas主分析
analysis_start <- Sys.time()
message("开始主分析；回归线程数：", N_CORES)
scPagwas_result <- scPagwas_main(
  Pagwas = NULL,
  gwas_data = gwas_input,
  Single_data = seurat_input,
  output.prefix = OUTPUT_PREFIX, 
  output.dirs = OUTPUT_DIR,
  Pathway_list = Genes_by_pathway_kegg,
  n.cores = N_CORES,
  assay = "RNA",
  singlecell = TRUE,
  iters_singlecell = 100,
  celltype = TRUE, 
  block_annotation = block_annotation,
  chrom_ld = chrom_ld
)

message("scPagwas主分析耗时（分钟）：",
        round(as.numeric(difftime(Sys.time(), analysis_start, units = "mins")), 2))

# 保存结果
save(scPagwas_result, file = file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "_scPagwas_final.RData")), compress = FALSE)
rm(gwas_input, seurat_input)
invisible(gc())

# 结果仍在内存中，后续直接使用；无需再次load。

# 检查数据结构
print(scPagwas_result)
names(scPagwas_result)
all_celltypes <- unique(scPagwas_result$cell_type)  # 或 seurat_obj$seurat_clusters
print(all_celltypes)
CELLTYPE_ORDER <- CELLTYPE_ORDER[CELLTYPE_ORDER %in% as.character(all_celltypes)]
# ==================== 可视化分析 ====================

# 1. 细胞类型比例的条形图（阳性/阴性细胞）
plot_bar_positie_nagtive(
  seurat_obj = scPagwas_result,
  var_ident = "cell_type",
  var_group = "positiveCells",
  vec_group_colors = c("#E8D0B3", "#7EB5A6"),
  do_plot = TRUE
)

# 2. UMAP可视化
scPagwas_Visualization(
  Single_data = scPagwas_result,
  p_thre = 0.05,
  FigureType = "umap",
  width = 7,
  height = 7,
  lowColor = "white",
  highColor = "red",
  output.dirs = "figures",
  size = 0.5,
  do_plot = TRUE
)

# 3. 细胞类型在阳性/阴性细胞中的分布
plot_bar_positie_nagtive(
  seurat_obj = scPagwas_result,
  var_ident = "positiveCells",
  var_group = "cell_type",
  p_thre = 0.01,
  vec_group_colors = NULL,
  f_color = colorRampPalette(brewer.pal(n = 10, name = "RdYlBu")),
  do_plot = TRUE
)

# 4. 遗传相关性散点图
heritability_cor_scatterplot(
  gene_heri_cor = scPagwas_result@misc$PCC,
  topn_genes_label = 10,
  color_low = "#035397",
  color_high = "#F32424",
  color_mid = "white",
  text_size = 2,
  do_plot = TRUE,
  max.overlaps = 20,
  width = 7,
  height = 7
)

# 5. Bootstrap结果可视化
Bootstrap_estimate_Plot(
  bootstrap_results = scPagwas_result@misc$bootstrap_results,
  figurenames = "figures/bootstrap_estimate.pdf",
  width = 9,
  height = 7,
  do_plot = TRUE
)

# 6. 通路贡献度点图（使用内置函数）
source(system.file("extdata", "plot_scpathway_contri_dot.R", package = "scPagwas"))

plot_scpathway_dot(
  Pagwas = scPagwas_result,
  Celltype = CELLTYPE_ORDER,
  topn_path_celltype = 5,
  filter_p = 0.05,
  max_logp = 15,
  display_max_sizes = FALSE,
  size_var = "logrankPvalue",
  col_var = "proportion",
  shape.scale = 8,
  cols.use = c("lightgrey", "#E45826"),
  dend_x_var = "logrankPvalue",
  dist_method = "euclidean",
  hclust_method = "ward.D",
  do_plot = TRUE,
  width = 7,
  height = 7
)

# ====================分割线=========================
# 感谢您浏览科研牛马成长日记课程
# 如有任何关于技术咨询的疑问请您务必联系我！(⑅˃◡˂⑅)
# 如果您有关于大课程咨询和建议请毫不犹豫的联系我！(⑅˃◡˂⑅)
# 1.邮箱：xiechengyong@yeah.net
# 2.b站：https://space.bilibili.com/19106653
# 3.微信号：chen9477777
# 4.淘宝商店ID：科研牛马工作室
# 5.小红书：小红薯677FE4D7
# ====================分割线==========================
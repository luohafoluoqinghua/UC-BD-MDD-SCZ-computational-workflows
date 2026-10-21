options(warn=1)
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(limma)
})
outdir <- "."
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2) stop("Usage: Rscript Reproduce_analysis.R discovery.rds validation.rds")
paths <- setNames(normalizePath(args,mustWork=TRUE),c("GSE182270","GSE214695"))
script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
if(length(script_arg)) setwd(dirname(normalizePath(sub("^--file=", "", script_arg[1]))))
dir.create("work/scrna", recursive=TRUE, showWarnings=FALSE)


targets <- c("KIAA1109","MED24","RHOA","HSD17B8")
common_types <- c("B cells","T cells","Plasma cells","Monocytes",
                  "Endothelial cells","Fibroblasts","Epithelial cells")
all_types <- c(common_types,"Macrophages","Mast cells")
min_cells <- 20L
min_donors <- 3L
all_dots <- all_donor_dots <- all_de <- all_coverage <- all_qc <- list()
for (ds in names(paths)) {
  cat("Processing",ds,"\n")
  obj <- readRDS(paths[[ds]])
  meta <- obj[[]]
  counts <- GetAssayData(obj, assay="RNA", slot="counts")
  stopifnot(identical(colnames(counts),rownames(meta)),
            all(targets %in% rownames(counts)),
            all(counts@x >= 0), all(abs(counts@x-round(counts@x)) < 1e-8),
            all(meta$Type %in% c("Control","Treat")))
  meta$condition <- factor(ifelse(meta$Type=="Control","Control","UC"),
                           levels=c("Control","UC"))
  meta$donor <- as.character(meta$orig.ident)
  stopifnot(all(tapply(as.character(meta$condition),meta$donor,
                      function(z)length(unique(z)))==1))
  libs <- Matrix::colSums(counts)
  stopifnot(all(libs>0))
  tar <- as.matrix(counts[targets,,drop=FALSE])
  norm <- sweep(tar,2,libs/1e4,"/")
  donor_records <- list()
  coverage <- list()
  for (ct in all_types) {
    for (dn in unique(meta$donor)) {
      ix <- which(meta$donor==dn & meta$cell_type==ct)
      nc <- length(ix)
      cond <- as.character(meta$condition[match(dn,meta$donor)])
      coverage[[length(coverage)+1]] <- data.frame(dataset=ds,donor=dn,
          condition=cond,cell_type=ct,n_cells=nc,eligible_pseudobulk=nc>=min_cells)
      if(nc==0) next
      donor_records[[length(donor_records)+1]] <- data.frame(dataset=ds,
          donor=dn,condition=cond,cell_type=ct,gene=targets,n_cells=nc,
          mean_cp10k=rowMeans(norm[,ix,drop=FALSE]),
          pct_expressing=100*rowMeans(tar[,ix,drop=FALSE]>0))
    }
  }
  donor_records <- do.call(rbind,donor_records)
  coverage <- do.call(rbind,coverage)
  all_coverage[[ds]] <- coverage
  all_donor_dots[[ds]] <- donor_records
  dots <- aggregate(cbind(mean_cp10k,pct_expressing) ~ dataset+condition+cell_type+gene,
                    data=donor_records,FUN=mean)
  dots$n_donors <- as.integer(table(interaction(donor_records$condition,
      donor_records$cell_type,donor_records$gene,drop=TRUE))[
        as.character(interaction(dots$condition,dots$cell_type,dots$gene,drop=TRUE))])
  # Explicit key matching avoids differences in factor ordering.
  keys <- paste(donor_records$condition,donor_records$cell_type,donor_records$gene,sep="|")
  dots$n_donors <- as.integer(table(keys)[paste(dots$condition,dots$cell_type,dots$gene,sep="|")])
  dots$log1p_mean_cp10k <- log1p(dots$mean_cp10k)
  all_dots[[ds]] <- dots
  result <- list()
  for (ct in common_types) {
    covct <- coverage[coverage$cell_type==ct & coverage$eligible_pseudobulk,]
    nn <- table(factor(covct$condition,levels=c("Control","UC")))
    res <- data.frame(dataset=ds,cell_type=ct,gene=targets,
      n_control=as.integer(nn["Control"]),n_uc=as.integer(nn["UC"]),
      log2FC=NA_real_,SE=NA_real_,CI_low=NA_real_,CI_high=NA_real_,
      df_total=NA_real_,t=NA_real_,PValue=NA_real_,FDR=NA_real_,
      status="Insufficient donors", stringsAsFactors=FALSE)
    if(all(nn>=min_donors)) {
      dn <- covct$donor
      ix <- which(meta$cell_type==ct & meta$donor %in% dn)
      indicator <- sparseMatrix(i=seq_along(ix),j=match(meta$donor[ix],dn),
          x=1,dims=c(length(ix),length(dn)))
      pb <- as.matrix(counts[,ix,drop=FALSE] %*% indicator)
      colnames(pb) <- dn
      group <- factor(covct$condition,levels=c("Control","UC"))
      design <- model.matrix(~group)
      y <- DGEList(pb,group=group)
      keep <- filterByExpr(y,design=design)
      y <- y[keep,,keep.lib.sizes=FALSE]
      y <- calcNormFactors(y,method="TMM")
      v <- voom(y,design,plot=FALSE)
      fit <- eBayes(lmFit(v,design),robust=TRUE)
      saveRDS(list(counts=pb,sample_metadata=covct,design=design,
                   retained_genes=rownames(y),fit=fit,voom=v),
              file=paste0("work/scrna/",ds,"_",gsub(" ","_",ct),"_pseudobulk.rds"))
      for (g in targets) {
        k <- match(g,res$gene)
        if(!g %in% rownames(fit$coefficients)) {res$status[k]<-"Low expression"; next}
        j <- match(g,rownames(fit$coefficients))
        se <- fit$stdev.unscaled[j,2]*sqrt(fit$s2.post[j])
        df <- fit$df.total[j]
        beta <- fit$coefficients[j,2]
        res[k,c("log2FC","SE","CI_low","CI_high","df_total","t","PValue")] <-
          list(beta,se,beta-qt(.975,df)*se,beta+qt(.975,df)*se,df,
               fit$t[j,2],fit$p.value[j,2])
        res$status[k] <- "Estimated"
      }
      cat(ct,"donors",nn,"retained genes",sum(keep),"\n")
    } else cat(ct,"NOT ESTIMABLE donors",nn,"\n")
    result[[ct]] <- res
  }
  result <- do.call(rbind,result)
  result$FDR <- p.adjust(result$PValue,method="BH",n=sum(is.finite(result$PValue)))
  all_de[[ds]] <- result
  all_qc[[ds]] <- data.frame(dataset=ds,n_genes=nrow(counts),n_cells=ncol(counts),
      n_control=length(unique(meta$donor[meta$condition=="Control"])),
      n_uc=length(unique(meta$donor[meta$condition=="UC"])),
      min_cells_per_pseudobulk=min_cells,min_donors_per_group=min_donors,
      tests=sum(result$status=="Estimated"))
  rm(obj,counts,tar,norm);gc()
}
dots <- do.call(rbind,all_dots)
donor_dots <- do.call(rbind,all_donor_dots)
de <- do.call(rbind,all_de)
coverage <- do.call(rbind,all_coverage)
qc <- do.call(rbind,all_qc)
# The localization score is standardized within gene over all observed
# dataset-condition-cell-type combinations; values are not truncated.
dots$gene_z <- ave(dots$log1p_mean_cp10k,dots$gene,
    FUN=function(x)if(sd(x)>0)(x-mean(x))/sd(x) else rep(0,length(x)))
write.csv(dots,file.path(outdir,"Expression_localization_source.csv"),row.names=FALSE)
write.csv(donor_dots,file.path(outdir,"Donor_expression_source.csv"),row.names=FALSE)
write.csv(de,file.path(outdir,"Pseudobulk_DE_source.csv"),row.names=FALSE)
write.csv(coverage,file.path(outdir,"Donor_cell_counts.csv"),row.names=FALSE)
write.csv(qc,file.path(outdir,"Cohort_QC.csv"),row.names=FALSE)
saveRDS(list(dots=dots,donor_dots=donor_dots,de=de,coverage=coverage,qc=qc,
             genes=targets,common_types=common_types,all_types=all_types),
        "work/scrna/analysis_results.rds")
capture.output(sessionInfo(),file=file.path(outdir,"R_sessionInfo.txt"))
cat("\nSignificant targeted tests:\n")
print(de[!is.na(de$FDR)&de$FDR<.05,c("dataset","cell_type","gene","log2FC","CI_low","CI_high","FDR")],row.names=FALSE)
cat("\nQC:\n");print(qc,row.names=FALSE)
cat("\nRange of plotted log expression:\n");print(range(dots$log1p_mean_cp10k))

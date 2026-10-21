# ====================分割线=========================
# 感谢您使用Fast2TWAS
# 如有任何关于该包安装和使用方面的疑问请您务必联系我！(⑅˃◡˂⑅)
# 如果您有关于该包的任何意见和建议请毫不犹豫的联系我！(⑅˃◡˂⑅)
# 1.邮箱：xiechengyong@yeah.net
# 2.b站：https://space.bilibili.com/19106653
# 3.微信号：chen9477777（科研咨询）；shengxindazhaxie（开发维护）
# 4.公众号：科研牛马成长日记
# 5.淘宝商店ID：科研牛马工作室
# 6.小红书：小红薯677FE4D7
# ====================分割线==========================


#设置当前工作路径
#Rsrudio打开该文件，根据R代码所在位置，设置所在文件夹为工作路径
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

setwd("/home/data/t230406/scRNA/14.pre_eCAVIAR/")
#清除一切变量
rm(list=ls())
#释放不再使用的内存
gc()
#显式路径
getwd()
#显示R包安装位置
.libPaths()

#下载以及加载包
if (!requireNamespace("tidyr", quietly = TRUE))install.packages("tidyr")
library(Fast2TWAS)
library(data.table)

#参数设置
inputFile = c("UC", "BD", "MDD", "SCZ")  # 本稿研究对象
#指定lead SNP的上下游范围
wind_size = 50000
#基因组参考人群
ref_pop="g1000_eur"    #c("g1000_eur","g1000_afr","g1000_amr","g1000_eas","g1000_sas")
#是否交换GWAS的effect_allele和other_allele(即对beta取负值)然后重新生成GWAS_variant列
trans=FALSE
#读入来自GTEx v8的组织列表文件
gtex_tissue_names <- fread("../03.reference/QTLEnrich/gtex_tissue_names.txt",header = F)$V1

################################################################################
#############循环构建eCAVIAR输入文件，该步骤耗时较长，注意掌握时间##############
################################################################################
write(x = "", file = "14.pre_eCAVIAR.log", append = F)
for (trait in inputFile) {
  #列出GenomicRiskLoci.txt文件位置
  GenomicRiskLoci_dir <- list.files(
    path = paste0("../13.fuma/",trait,"/"),
    pattern = "^FUMA.*",
    full.names = TRUE,
    include.dirs = TRUE
  ) %>% .[file.info(.)$isdir]
  #读入FUMA风险基因座结果
  GenomicRiskLoci=fread(paste0(GenomicRiskLoci_dir,"/GenomicRiskLoci.txt"))
  #将GenomicRiskLoci数据框中CHR列的23或小写x转换为大写X
  GenomicRiskLoci$chr <- gsub("^23$|^x$", "X", as.character(GenomicRiskLoci$chr))
  
  #读入GWAS数据
  gwas_data <- fread(paste0("../07.hg19Tohg38/",trait,"_hg38_qc.txt"))
  #将gwas_data数据框中CHR列的23或小写x转换为大写X
  gwas_data$CHR <- gsub("^23$|^x$", "X", as.character(gwas_data$CHR))
  for (snp in GenomicRiskLoci$rsID) {
    #显示运行信息  
    i=which(snp == GenomicRiskLoci$rsID)
    message("\n当前性状：",trait,"，当前第",i,"个GenomicRiskLoci,LeadSNPs：",GenomicRiskLoci[i,]$rsID)
    
    #筛选GWAS数据lead SNP上下游50kb范围内的SNPs
    filter <- subset(gwas_data,SNP==snp)[1]
    chrom <- filter$CHR
    pos <- filter$BP
    endbp <- pos+wind_size
    startbp <- pos-wind_size
    markers <- subset(gwas_data,(CHR==chrom & BP <= endbp & BP >= startbp))
    head(markers)  
    
    #是否翻转GWAS的effect_allele和other_allele
    if (trans) {
      #计算z值
      markers=transform(markers,z_gwas=-beta/se)
      head(markers)
      
      #翻转
      if ("eaf" %in% colnames(markers)) {
        markers$eaf <- 1-markers$eaf
      }
      
      #添加GWAS_variant
      markers$GWAS_variant <- paste0("chr",markers$CHR,
                                     "_",markers$BP,
                                     "_",markers$effect_allele,
                                     "_",markers$other_allele,
                                     "_b38")
    }else{
      #计算z值
      markers=transform(markers,z_gwas=beta/se)
      head(markers)
      
      #添加GWAS_variant
      markers$GWAS_variant <- paste0("chr",markers$CHR,
                                     "_",markers$BP,
                                     "_",markers$other_allele,
                                     "_",markers$effect_allele,
                                     "_b38")
    }
    
    #eQTL和sQTL循环
    for (qtl_type in c("eQTL","sQTL")) {
      #gtex_tissue循环
      for (tissue in gtex_tissue_names) {  
        #显示运行信息  
        message("\n当前性状：",trait,"，当前QTL类型：",qtl_type,"，当前组织：",tissue)
        
        #读入QTL数据
        if (qtl_type=="eQTL") {
          qtl_data_signif_variant_file <- paste0("../03.reference/QTLEnrich/GTEx_Analysis_v8_",qtl_type,"/",tissue,".v8.signif_variant_gene_pairs.txt.gz")
        }else if(qtl_type=="sQTL"){
          qtl_data_signif_variant_file <- paste0("../03.reference/QTLEnrich/GTEx_Analysis_v8_",qtl_type,"/",tissue,".v8.sqtl_signifpairs.txt.gz")
        }
        
        qtl_data_signif_variant <- fread(qtl_data_signif_variant_file)
        
        # 将 variant_id 按 "_" 分割为 4 列（默认列名）
        qtl_data_signif_variant_split <- qtl_data_signif_variant %>%
          tidyr::separate(
            col = variant_id,
            into = c("chr", "pos", "ref", "alt","build"),
            sep = "_",
            remove = FALSE  # 保留原始 variant_id 列
          )
        
        #筛选QTL数据lead SNP上下游50kb范围内的SNPs，并去除重复的SNP
        markers_QTL <- subset(qtl_data_signif_variant_split,(chr==paste0("chr",chrom) & pos <= endbp & pos >= startbp)) %>% dplyr::distinct(variant_id,.keep_all = TRUE)
        head(markers_QTL) 
        
        #合并markers_QTL和markers
        markers_QTL_gwas <- merge(markers_QTL,markers,by.x = "variant_id", by.y = "GWAS_variant")
        head(markers_QTL_gwas)
        #判断合并后的数据是否nrow(markers_QTL_gwas)>0
        if (nrow(markers_QTL_gwas)==0) {
          line=paste0("当前性状：",trait,"，当前lead SNP：",snp,"，当前QTL类型：",qtl_type,"，当前组织：",tissue,"，无交集SNP")
          message(line)
          write(x = line, file = "14.pre_eCAVIAR.log", append = T)
          next
        }
        
        #构建corr：生成 LD 矩阵
        corr= try(ieugwasr::ld_matrix_local(markers_QTL_gwas$SNP, 
                                            bfile=paste0("../03.reference/magma/g1000/",ref_pop), 
                                            # plink_bin=friendlypostGWASsoft::plink_binary(),
                                            plink_bin=plinkbinr::get_plink_exe(),
                                            with_alleles = F),silent = T)
        if ("try-error" %in% class(corr)) {
          line=paste0("当前性状：",trait,"，当前lead SNP：",snp,"，当前QTL类型：",qtl_type,"，当前组织：",tissue,"，生成 LD 矩阵失败")
          message(line)
          write(x = line, file = "14.pre_eCAVIAR.log", append = T)
          next
        }
        
        #求交集SNP,重新构建markers
        final_SNP=intersect(markers_QTL_gwas$SNP, colnames(corr))
        
        #had at least 5 e/sVariants within the GWAS locus LD interval
        if (length(final_SNP)<5) {
          line=paste0("当前性状：",trait,"，当前lead SNP：",snp,"，当前QTL类型：",qtl_type,"，当前组织：",tissue,"，交集SNP小于5个")
          message(line)
          write(x = line, file = "14.pre_eCAVIAR.log", append = T)
          next
        }
        
        markers_QTL_gwas_final_SNP <- markers_QTL_gwas[markers_QTL_gwas$SNP%in%final_SNP,]
        #计算z值
        markers_QTL_gwas_final_SNP=transform(markers_QTL_gwas_final_SNP,z_qtl=slope/slope_se)
        head(markers_QTL_gwas_final_SNP)
        
        #创建trait相对应QTL文件夹
        trait_qtl_type_dir=paste0(trait,"/",qtl_type)
        if (! dir.exists(trait_qtl_type_dir)) {
          dir.create(trait_qtl_type_dir, recursive = T) 
        }
        
        #创建trait相对应GWAS文件夹
        trait_gwas_dir=paste0(trait,"/GWAS")
        if (! dir.exists(trait_gwas_dir)) {
          dir.create(trait_gwas_dir, recursive = T) 
        }
        
        #创建trait相对应GWAS文件夹
        trait_corr_dir=paste0(trait,"/corr")
        if (! dir.exists(trait_corr_dir)) {
          dir.create(trait_corr_dir, recursive = T) 
        }
        
        #保存全部数据
        fwrite(markers_QTL_gwas_final_SNP,
               file=paste0(trait_gwas_dir,"/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".all.txt"),
               sep = "\t")
        
        #保存GWAS中SNP的corr数据
        fwrite(corr,
               file=paste0(trait_corr_dir,"/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".corr.txt"),
               sep = "\t",
               col.names = F)
        
        #保存数据，且只保存marker，z两列并不要表头
        fwrite(markers_QTL_gwas_final_SNP[,c("SNP","z_qtl")],
               file=paste0(trait_qtl_type_dir,"/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".qtl.txt"),
               sep = "\t",
               col.names = F)
        
        #保存数据，且只保存marker，z两列并不要表头
        fwrite(markers_QTL_gwas_final_SNP[,c("SNP","z_gwas")],
               file=paste0(trait_gwas_dir,"/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".gwas.txt"),
               sep = "\t",
               col.names = F)
      }  
    }
  }
}

# ====================分割线=========================
# 感谢您使用Fast2TWAS
# 如有任何关于该包安装和使用方面的疑问请您务必联系我！(⑅˃◡˂⑅)
# 如果您有关于该包的任何意见和建议请毫不犹豫的联系我！(⑅˃◡˂⑅)
# 1.邮箱：xiechengyong@yeah.net
# 2.b站：https://space.bilibili.com/19106653
# 3.微信号：chen9477777（科研咨询）；shengxindazhaxie（开发维护）
# 4.公众号：科研牛马成长日记
# 5.淘宝商店ID：科研牛马工作室
# 6.小红书：小红薯677FE4D7
# ====================分割线==========================
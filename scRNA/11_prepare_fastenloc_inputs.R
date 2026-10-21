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
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#清除一切变量
rm(list=ls())
#释放不再使用的内存
gc()
#显式路径
getwd()
#显示R包安装位置
.libPaths()

#下载以及加载包
library(Fast2TWAS)
library(data.table)

#参数设置
inputFile = c("UC", "BD", "MDD", "SCZ")  # 本稿研究对象
#指定lead SNP的上下游范围
wind_size = 50000
#基因组参考人群
ref_pop="g1000_eur"    #c("g1000_eur","g1000_afr","g1000_amr","g1000_eas","g1000_sas")
#读入来自GTEx v8的组织列表文件
gtex_tissue_names <- fread("../03.reference/QTLEnrich/gtex_tissue_names.txt",header = F)$V1

################################################################################
#############循环构建eCAVIAR输入文件，该步骤耗时较长，注意掌握时间##############
################################################################################
write(x = "", file = "17.pre_fastenloc.log", append = F)
for (trait in inputFile) {
  #创建trait相对应QTL文件夹
  trait_qtl_type_dir=paste0(trait)
  if (! dir.exists(trait_qtl_type_dir)) {
    dir.create(trait_qtl_type_dir, recursive = T) 
  }
  
  #列出GenomicRiskLoci.txt文件位置
  GenomicRiskLoci_dir <- list.files(
    path = paste0("../13.fuma/",trait,"/"),
    pattern = "^FUMA.*",
    full.names = TRUE,
    include.dirs = TRUE
  ) %>% .[file.info(.)$isdir]
  #读入FUMA风险基因座结果
  GenomicRiskLoci=fread(paste0(GenomicRiskLoci_dir,"/GenomicRiskLoci.txt"))

  for (snp in GenomicRiskLoci$rsID) {
    #显示运行信息  
    i=which(snp == GenomicRiskLoci$rsID)
    message("\n当前性状：",trait,"，当前第",i,"个GenomicRiskLoci,LeadSNPs：",GenomicRiskLoci[i,]$rsID)
    
    #eQTL和sQTL循环
    for (qtl_type in c("eQTL","sQTL")) {
      #gtex_tissue循环
      for (tissue in gtex_tissue_names) {  
        #显示运行信息  
        message("\n当前性状：",trait,"，当前QTL类型：",qtl_type,"，当前组织：",tissue,"，当前lead SNP：",snp)
        
        #trait相对应GWAS文件夹
        trait_gwas_dir=paste0("../14.pre_eCAVIAR/",trait,"/GWAS")

        gwas_qtl_all_file <- paste0(trait_gwas_dir,"/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".all.txt")
        if (!file.exists(gwas_qtl_all_file)) {
          line=paste0("当前性状：",trait,"，当前lead SNP：",snp,"，当前QTL类型：",qtl_type,"，当前组织：",tissue,"，不存在文件")
          message(line)
          write(x = line, file = "17.pre_fastenloc.log", append = T)
          next
        }
        
        #构建fastenloc输入文件
        markers_QTL_gwas_final_SNP <- fread(gwas_qtl_all_file)
        #判断QTL和GWAS是的effect_allele和other_allele列是否相同
        for (i in 1:nrow(markers_QTL_gwas_final_SNP)) {
          #判断QTL和GWAS的效应等位基因是否相同
            if(markers_QTL_gwas_final_SNP$effect_allele[i]==markers_QTL_gwas_final_SNP$ref[i]  & 
                   markers_QTL_gwas_final_SNP$other_allele[i]==markers_QTL_gwas_final_SNP$alt[i]) {
              markers_QTL_gwas_final_SNP$slope[i]==-markers_QTL_gwas_final_SNP$slope[i]
            }
        }
        
        if (qtl_type == "sQTL") {
          # 将 variant_id 按 "_" 分割为 4 列（默认列名）
          markers_QTL_gwas_final_SNP <- markers_QTL_gwas_final_SNP %>%
            tidyr::separate(
              col = phenotype_id,
              into = c("chr", "start", "end", "cul","gene_id"),
              sep = ":",
              remove = FALSE  # 保留原始 variant_id 列
            )
        }
        
        #提取列
        markers_QTL_gwas_final_SNP <- markers_QTL_gwas_final_SNP[,c("gene_id","SNP","slope","slope_se","beta","se")]
        #保存全部数据
        fwrite(markers_QTL_gwas_final_SNP,
               file=paste0(trait_qtl_type_dir,"/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".fastenloc.txt"),
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
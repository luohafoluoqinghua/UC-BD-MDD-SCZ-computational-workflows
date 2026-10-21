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
setwd("/home/data/t230406/scRNA/15.eCAVIAR/")

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
#设置睡眠暂停时间1seconds
# sleeptime=1
#读入来自GTEx v8的组织列表文件
gtex_tissue_names <- fread("../03.reference/QTLEnrich/gtex_tissue_names.txt",header = F)$V1

################################################################################
#################循环完成eCAVIAR，该步骤耗时较长，注意掌握时间##################
################################################################################
write(x = "", file = "15.eCAVIAR.log", append = F)
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
  
  for (snp in GenomicRiskLoci$rsID) {
    #显示运行信息  
    i=which(snp == GenomicRiskLoci$rsID)
    message("\n当前性状：",trait,"，当前第",i,"个GenomicRiskLoci,LeadSNPs：",GenomicRiskLoci[i,]$rsID)
    
    #eQTL和sQTL循环
    for (qtl_type in c("eQTL","sQTL")) {
      #gtex_tissue循环
      for (tissue in gtex_tissue_names) {  
        #显示运行信息  
        message("\n当前性状：",trait,"，当前QTL类型：",qtl_type,"，当前组织：",tissue)
        
        trait_qtl_type_dir=paste0(trait,"/",qtl_type)
        trait_gwas_dir=paste0(trait,"/GWAS")
        trait_corr_dir=paste0(trait,"/corr")
        
        if (! dir.exists(paste0("result/",trait_qtl_type_dir))) {
          dir.create(paste0("result/",trait_qtl_type_dir), recursive = T) 
        }
        
        #获取QTL文件路径
        qtl_data_file <- paste0("../14.pre_eCAVIAR/",trait_qtl_type_dir,"/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".qtl.txt")
        
        #获取GWAS文件路径
        gwas_data_file <- paste0("../14.pre_eCAVIAR/",trait_gwas_dir,"/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".gwas.txt")
        
        #ld文件路径
        ld_data_file <- paste0("../14.pre_eCAVIAR/",trait_corr_dir,"/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".corr.txt")
        
        # 判断输入文件是否存在
        if (!file.exists(qtl_data_file)  || 
            !file.exists(gwas_data_file) || 
            !file.exists(ld_data_file)   || 
            length(qtl_data_file) != 1   || 
            length(gwas_data_file) != 1  ||
            length(ld_data_file) != 1) {
          line=paste0("当前性状：",trait,"，当前lead SNP：",snp,"，当前QTL类型：",qtl_type,"，当前组织：",tissue,"，无交集SNP")
          message(line)
          write(x = line, file = "15.eCAVIAR.log", append = T)
          next
        }
        
        #调用eCAVIAR进行基于贝叶斯的GWAS和eQTL共同定位
        caviar_eCAVIAR_run(gwas_z_file=gwas_data_file,
                           eqtl_z_file=qtl_data_file,
                           gwas_ld_file=ld_data_file,
                           eqtl_ld_file=ld_data_file,
                           rho_prob=0.95,
                           causal=2,
                           f=1,
                           outname=paste0(trait_qtl_type_dir,"/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".eCAVIAR"))
        
        # Sys.sleep(sleeptime)  # 暂停sleeptime秒
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
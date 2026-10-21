# ====================分割线=========================
# 感谢您使用Fast2TWAS
# 如有任何关于该包安装和使用方面的疑问请您务必联系我！(⑅˃◡˂⑅)
# 如果您有关于该包的任何意见和建议请毫不犹豫的联系我！(⑅˃◡˂⑅)
# 1.邮箱：xiechengyong@yeah.net
# 2.b站：https://space.bilibili.com/397473016?spm_id_from=..0.0
# 3.微信号：shengxindazhaxie
# 4.公众号：生信蟹道人
# 5.简书：生信蟹道人
# ====================分割线==========================


#设置当前工作路径
#setwd("D:/myproject/mylesson/6.shared_genetic_local_genetic_correlation/08.heritability_ldsc/")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#清除一切变量
rm(list=ls())
#释放不再使用的内存
gc()
#显式路径
getwd()
#显示R包安装位置
.libPaths()
#安装和加载包
library(Fast2TWAS)

#输入参数-----------------------------------------------------------------------
#读入index文件
index <- data.table::fread("../07.export_index/gwasfile_info.txt")

# 仅使用本稿四个研究对象，按UC、BD、MDD、SCZ排序。
study_traits <- c("UC", "BD", "MDD", "SCZ")
if (anyDuplicated(index$phenotype) || !all(study_traits %in% index$phenotype)) {
  stop("Regenerate 07.export_index for UC, BD, MDD and SCZ before running this step.")
}
index <- index[match(study_traits, index$phenotype), ]
if (!identical(as.character(index$name), study_traits)) {
  stop("Study input filenames do not match the four-trait configuration; rerun 07.export_index.")
}


# 生成ldsc格式数据 -------------------------------------------------------------
inputFile <-index$name
phenotype <-index$phenotype
for (i in 1:length(inputFile)) {
  message("\n当前第",i,"个gwas数据：",inputFile[i],"; phen:",phenotype[i])
  # 删除文件eaf列
  data.table::fread(paste0("../06.pre_data/",inputFile[i],".txt")) %>% subset(., select = -eaf) %>% data.table::fwrite(phenotype[i],sep = "\t")
  # ldsc格式数据
  Fast2TWAS::ldsc_format(gwasfile=phenotype[i],
                         outname = phenotype[i],
                         flag_options = "--chunksize 500000",
                         merge_alleles = "../03.reference/ldsc/w_hm3.snplist",
                         force = F) 
  
  #计算SNP-based heritability和遗传度---------------------------------------------  
  inputgwas=paste0("./result/",phenotype[i],".sumstats.gz")
  Fast2TWAS::ldsc_h2(gwasfile=inputgwas,
                     LD_ref_dir = "../03.reference/ldsc/eur_w_ld_chr/",
                     flag_options = NULL,
                     outname = paste0(phenotype[i],"_ldsc_h2"),
                     force = F)
  
  # 删除文件
  unlink(phenotype[i])
}


# ====================分割线=========================
# 感谢您使用Fast2TWAS
# 如有任何关于该包安装和使用方面的疑问请您务必联系我！(⑅˃◡˂⑅)
# 如果您有关于该包的任何意见和建议请毫不犹豫的联系我！(⑅˃◡˂⑅)
# 1.邮箱：xiechengyong@yeah.net
# 2.b站：https://space.bilibili.com/397473016?spm_id_from=..0.0
# 3.微信号：shengxindazhaxie
# 4.公众号：生信蟹道人
# 5.简书：生信蟹道人
# ====================分割线==========================
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
#setwd("D:/myproject/mylesson/6.shared_genetic_local_genetic_correlation/11.genetic_correlation_HDL/")
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


#指定性状名称
trait="UC"

#性状名
inputFile = c("BD", "MDD", "SCZ")

# 循环调用HDL包来计算rg -----------------------------------------------------------
for (phen in inputFile) {
  trait_name=index$phenotype[index$name==trait]
  other_name=index$phenotype[index$name==phen]
  message("\n当前phen1 ",trait,":",trait_name," ; ","phen2:",phen,":",other_name)
  HDL_rg_data <- HDL_rg(gwas1file=paste0("../06.pre_data/",trait,".txt"),
                        gwas2file = paste0("../06.pre_data/",phen,".txt"),
                        LD.path = "../03.reference/HDL/UKB_imputed_SVD_eigen99_extraction/",
                        Nref = 335265, #默认的UK Biobank参考样本，Nref=335265
                        outfile = paste0(trait_name,"_",other_name,".txt"))
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
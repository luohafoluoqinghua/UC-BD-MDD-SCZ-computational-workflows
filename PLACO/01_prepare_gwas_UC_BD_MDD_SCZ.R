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
#setwd("D:/myproject/mylesson/6.shared_genetic_local_genetic_correlation/06.pre_data")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#清除一切变量
rm(list=ls())
#释放不再使用的内存
gc()
#显式路径
getwd()
#显示R包安装位置
.libPaths()

#加载包
library(Fast2TWAS)

#Fast2TWAS包中建议的GWAS数据格式
show_demo_gwas()

# 输入参数：本稿四个研究对象及其OpenGWAS数据编号
# 本地原始VCF放入 ../04.rawdata/；输出简称与下游索引保持一致。
# 本次将旧疾病示例替换为本稿输入，不表示已执行或验证历史分析。
study <- data.frame(
  phenotype = c("UC", "BD", "MDD", "SCZ"),
  accession = c("ebi-a-GCST004133", "ieu-b-5110", "ieu-b-102", "ieu-b-5102"),
  stringsAsFactors = FALSE
)
inputFile <- paste0(study$accession, ".vcf.gz")
raw_files <- file.path("../04.rawdata", inputFile)
if (any(!file.exists(raw_files))) {
  stop("Missing study GWAS input: ", paste(raw_files[!file.exists(raw_files)], collapse = ", "))
}
for (i in seq_len(nrow(study))) {
  trait <- study$phenotype[i]
  data <- read_ieu_data(
      ieu_id = raw_files[i],
      download = FALSE,
      maf = 0.01,
      pvalue = NULL,
      is_fast2twas_format = TRUE)
  data.table::fwrite(data, file = paste0(trait, ".txt"), sep = "\t")
  rm(data)
  gc()
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
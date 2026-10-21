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

#加载包
library(Fast2TWAS)
library(data.table)

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
      pvalue = NULL,
      is_fast2twas_format = TRUE)
  data <- subset(data, !(CHR == 6 & BP > 25000000 & BP < 35000000))
  data.table::fwrite(data, file = paste0(trait, ".txt"), sep = "\t")
  rm(data)
  gc()
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
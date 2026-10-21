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
library(data.table)

#参数设置
inputFile = c("UC", "BD", "MDD", "SCZ")  # 本稿研究对象

# 创建文件夹 -------------------------------------------------------------------
for (trait in inputFile) {
  if (! dir.exists(trait)) {
    dir.create(trait) 
  }
}


# 放入FUMA结果之后，再运行下列代码 ---------------------------------------------
#将GenomicRiskLoci.txt数据合并起来生成stable6.txt
#创建数据框
stable6 <- data.frame()
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
  GenomicRiskLoci$Trait <- trait
  GenomicRiskLoci <- subset(GenomicRiskLoci, select=c(Trait,uniqID,rsID,pos))

  #读入GWAS数据
  gwas_data <- fread(paste0("../07.hg19Tohg38/",trait,"_hg38_qc.txt"))
  #将gwas_data数据框中CHR列的23或小写x转换为大写X
  gwas_data$CHR <- gsub("^23$|^x$", "X", as.character(gwas_data$CHR))
  #合并GenomicRiskLoci和gwas_data
  pre_stable6 <- merge(GenomicRiskLoci, gwas_data, by.x="rsID", by.y="SNP")
  #修改BP表头为hg38
  setnames(pre_stable6, "BP", "hg38")
  #修改uniqID表头为GWAS_locus
  setnames(pre_stable6, "uniqID", "GWAS_locus")
  #修改pos表头为hg19
  setnames(pre_stable6, "pos", "hg19")
  #修改pre_stable6数据中的CHR列加上chr
  pre_stable6$CHR <- paste0("chr", pre_stable6$CHR)
  #修改pre_stable6数据中的GWAS_locus列加上chr
  pre_stable6$GWAS_locus <- paste0("chr", pre_stable6$GWAS_locus)
  #将pre_stable6数据框中GWAS_locus列中以chr23开头的列转换为大写chrX开头
  pre_stable6$GWAS_locus <- gsub("^chr23", "chrX", pre_stable6$GWAS_locus)

  #合并stable6
  stable6 <- rbind(stable6, pre_stable6)
  #保存stable6
  fwrite(stable6, 
         file="stable6.txt", 
         sep="\t")
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
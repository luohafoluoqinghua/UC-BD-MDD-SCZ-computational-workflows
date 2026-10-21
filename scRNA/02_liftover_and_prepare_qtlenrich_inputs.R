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

#参数设置
inputFile = c("UC", "BD", "MDD", "SCZ")  # 本稿研究对象
maf=0.01                              #maf过滤阈值
build="hg19"                          #指定GWAS数据基因组版本"hg19"或"hg38"---记得要改！！！！！芬兰hg38---其他hg19

for (trait in inputFile) {
  # 转换基因组版本hg19TOhg38----------------------------------------------------
  # #gwas_liftover函数运行时需要下载hg38ToHg19.over.chain.gz文件，注意网速（关闭魔法/科学上网），
  #若运行失败可多试几次或增大timeout参数
  options(timeout = 600)
  if (build=="hg19") {
    data=gwas_liftover(gwasfile= paste0("../06.pre_data/",trait,".txt"),
                       build_from = "hg19",
                       build_to = "hg38",
                       # local = T,
                       savefile = paste0(trait,"_hg38.txt"))
  }else if((build=="hg38")){
    data <- fread(paste0("../06.pre_data/",trait,".txt"))
  }

  head(data)
  #删除MHC区域（IBD和它有关）
  #data<- subset(data,!(CHR==6 & BP > 25000000 & BP < 35000000))
  # MAF过滤
  if ("eaf" %in% colnames(data)) {
    data=subset(data,eaf > maf & eaf < 1-maf)
  }
  #保存数据用于FUMA
  fwrite(subset(data,select=-c(CHR,BP)),
         file = paste0(trait,"_hg38_qc_fuma.gz"),
         sep = "\t")
  
  #保存数据
  fwrite(data,
         file = paste0(trait,"_hg38_qc.txt"),
         sep = "\t")
  
  
  # 转换为QTLEnrich输入文件 ----------------------------------------------------
  #修改表头
  colnames(data)[colnames(data) == "pval"] <- "gwas_p_value"
  colnames(data)[colnames(data) == "CHR"] <- "chr"
  colnames(data)[colnames(data) == "BP"] <- "pos"
  # colnames(data)[colnames(data) == "effect_allele"] <- "effect_allele"
  colnames(data)[colnames(data) == "other_allele"] <- "non_effect_allele"
  #染色体添加chr
  data$chr <- paste0("chr",data$chr)
  head(data)
  
  #保存数据
  fwrite(data,
         file = paste0(trait,"_hg38_qc_QTLEnrich.txt"),
         sep = "\t")

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
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

#安装和加载包
if (!requireNamespace("qqman", quietly = TRUE))install.packages("qqman")
if (!requireNamespace("topr", quietly = TRUE))install.packages("topr")
library(Fast2TWAS)
library(data.table)
library(topr)

index <- fread("../07.export_index/gwasfile_info.txt")

# 仅使用本稿四个研究对象，按UC、BD、MDD、SCZ排序。
study_traits <- c("UC", "BD", "MDD", "SCZ")
if (anyDuplicated(index$phenotype) || !all(study_traits %in% index$phenotype)) {
  stop("Regenerate 07.export_index for UC, BD, MDD and SCZ before running this step.")
}
index <- index[match(study_traits, index$phenotype), ]
if (!identical(as.character(index$filename), paste0(study_traits, ".txt"))) {
  stop("Study input filenames do not match the four-trait configuration; rerun 07.export_index.")
}

width <- 14
height <- 6
sig.thres <- 5e-8

filename <- index$filename

for (trait in filename) {
  gwas_file <- paste0("../06.pre_data/", trait)
  snpManhattan <- na.omit(fread(gwas_file))
  
  # 先画图并检查返回值
  p <- topr::manhattan(snpManhattan)
  print(class(p))
  
  outname <- paste0(index$phenotype[which(filename %in% trait)], "_Manhattan.png")
  
  if (!is.null(p) && inherits(p, "ggplot")) {
    ggplot2::ggsave(outname, plot = p, width = width, height = height)
  } else {
    png(filename = outname, width = width * 300, height = height * 300, res = 300)
    topr::manhattan(snpManhattan)
    dev.off()
  }
  
  z <- snpManhattan$beta / snpManhattan$se
  lambda <- round(median(z^2) / 0.454, 3)
  print(lambda)
  
  png(filename = paste0(index$phenotype[which(filename %in% trait)], "_Q-Q.png"),
      width = height * 100, height = height * 100, units = "px")
  qqman::qq(snpManhattan$pval)
  text(2.5, 10 * ceiling(max(-log10(na.omit(snpManhattan$pval))) / 10) / 2,
       bquote(lambda == .(lambda)), cex = 1, col = "black")
  dev.off()
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
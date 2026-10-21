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
library(Fast2TWAS)
library(data.table)

#设置主性状
main_trait <- "UC"
#读入index文件
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

#筛选去掉主性状
index <- subset(index,phenotype!=main_trait)

# 循环进行pleiofdr分析 ---------------------------------------------------------
phenotype <-index$phenotype
for (trait in phenotype) {
  for (stattype in c("condfdr","conjfdr")) {
    
  #main_trait-> trait  
  path1 <- paste0("../18.condFDR_pleiofdr/result/",stattype,"_",main_trait,"_",trait)
  mat1 <- paste0(path1,"/result.mat")
  outfile1 <- paste0(stattype,"_",main_trait,"_",trait,".mat.csv")
  
  python_convert_fdrmat2csv(mat=mat1,
                            ref="../03.reference/condFDR/9545380.ref",
                            head=10,
                            is_gzip=FALSE,
                            outname=outfile1)
  #如果是conjfdr则跳过
  if (stattype=="conjfdr")next
  #trait-> main_trait 
  path2 <- paste0("../18.condFDR_pleiofdr/result/",stattype,"_",trait,"_",main_trait)
  mat2 <- paste0(path2,"/result.mat")
  outfile2 <- paste0(stattype,"_",trait,"_",main_trait,".mat.csv")
  
  python_convert_fdrmat2csv(mat=mat2,
                            ref="../03.reference/condFDR/9545380.ref",
                            head=10,
                            is_gzip=FALSE,
                            outname=outfile2)
  
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
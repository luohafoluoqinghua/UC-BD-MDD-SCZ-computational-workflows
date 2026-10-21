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
#setwd("D:/myproject/mylesson/6.shared_genetic_local_genetic_correlation/07.export_index/")
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

#输入参数-----------------------------------------------------------------------
#指定性状名称
trait="UC"

#性状名
inputFile = c("BD", "MDD", "SCZ")

# 构建input.file文件-------------------------------------------------------------
# ● name列是文件的名称，包含文件后缀
# ● phenotype是GWAS数据的性状名称
# ● cases和controls分别是GWAS数据的病例数和对照数
# ● 若GWAS数据的性状是连续性，则cases和controls设置为NA
#所有性状名
allinputFile <- c(trait,inputFile)
# 样本量依据当前稿件：UC、BD、MDD、SCZ；需与实际下载版本匹配。
#phenotype是GWAS数据的性状名称
phenotype = c("UC", "BD", "MDD", "SCZ")
#cases和controls分别是GWAS数据的病例数和对照数
#若GWAS数据的性状是连续性，则cases和controls设置为NA
cases = c(12366, 41917, 170756, 52017)
#cases和controls分别是GWAS数据的病例数和对照数
#若GWAS数据的性状是连续性，则cases和controls设置为NA
controls = c(33609, 371549, 329443, 75889)

input_info_data <- data.frame(name=allinputFile,
                              phenotype=phenotype,
                              cases=cases,
                              controls=controls)

#保存文件
data.table::fwrite(input_info_data,
                   file = "gwasfile_info.txt",
                   na = "NA",
                   quote = FALSE,
                   sep = "\t")


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
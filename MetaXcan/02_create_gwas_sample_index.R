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

################################################################################
#########构建gwasfile_info.txt 文件这个步骤决不能错，错了后面全部白瞎了#########
################################################################################
# 构建gwasfile_info.txt文件-----------------------------------------------------
# ● filename列是文件的名称，包含文件名后缀
# ● phenotype列是GWAS数据的性状名简称
# ● cases和controls分别是GWAS数据的病例数和对照数：若GWAS数据的性状是连续性，则cases和controls设置为NA
# ● samplesize列指的是GWAS数据的总样本量
################################################################################
#性状名
filename = paste0(c("UC", "BD", "MDD", "SCZ"), ".txt")

# 样本量依据当前稿件：UC、BD、MDD、SCZ；需与实际下载版本匹配。
#phenotype是GWAS数据的性状名称
phenotype = c("UC", "BD", "MDD", "SCZ")


#cases和controls分别是GWAS数据的病例数和对照数
#若GWAS数据的性状是连续性，则cases和controls设置为NA
cases = c(12366, 41917, 170756, 52017)

#cases和controls分别是GWAS数据的病例数和对照数
#若GWAS数据的性状是连续性，则cases和controls设置为NA
controls = c(33609, 371549, 329443, 75889)

#samplesize指的是GWAS数据的总样本量
samplesize = cases + controls


#构建数据框
input_info_data <- data.frame(filename=filename,
                              phenotype=phenotype,
                              cases=cases,
                              controls=controls,
                              samplesize=samplesize)

#保存文件
fwrite(input_info_data,
       file = "gwasfile_info.txt",
       na = "NA",
       quote = FALSE,
       sep = "\t")


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
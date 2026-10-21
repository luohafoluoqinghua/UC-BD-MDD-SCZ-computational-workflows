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
  traitfile1 <- paste0(main_trait,".mat")
  traitfiles <- paste0(trait,".mat")

  #conjfdr：main_trait-> trait
  pleiofdr_run(reffile="../03.reference/condFDR/ref9545380_1kgPhase3eur_LDr2p1.mat",
               traitfolder="../16.pre_data_condFDR/result/",
               traitfile1=traitfile1,
               traitname1=main_trait,
               traitfiles=traitfiles,
               traitnames=trait,
               outputdir="./result/",
               stattype="conjfdr",
               fdrthresh=0.05,
               is_randprune=TRUE,
               randprune_n=500,
               exclude_chr_pos="6 25119106 33854733; 8 7200000 12500000",
               manh_fontsize_genenames=12,
               manh_yspace=0.75,
               manh_ymargin=0.25,
               manh_colorlist="[1 0 0; 1 0.5 0 ; 0 0.75 0.75; 0 0.5 0; 0.75 0 0.75; 0 0 1; 0 1 0; 0 1 1]",
               refinfo="../03.reference/condFDR/9545380.ref",
               is_reset_pruneidx=TRUE,
               randprune_repeats="default",
               pthresh=1,
               is_perform_gc=TRUE,
               is_use_standard_gc=FALSE,
               is_randprune_gc=TRUE,
               is_exclude_from_discovery=TRUE,
               mafthresh = 0.005,
               is_exclude_ambiguous_snps=TRUE,
               is_onscreen=FALSE,
               is_dummy_zscore = FALSE,
               is_exit_matlab_upon_completion = FALSE,
               matlab_path=NULL)
  #condfdr：main_trait-> trait
  pleiofdr_run(reffile="../03.reference/condFDR/ref9545380_1kgPhase3eur_LDr2p1.mat",
               traitfolder="../16.pre_data_condFDR/result/",
               traitfile1=traitfile1,
               traitname1=main_trait,
               traitfiles=traitfiles,
               traitnames=trait,
               outputdir="./result/",
               stattype="condfdr",
               fdrthresh=0.01,
               is_randprune=TRUE,
               randprune_n=500,
               exclude_chr_pos="6 25119106 33854733; 8 7200000 12500000",
               manh_fontsize_genenames=12,
               manh_yspace=0.75,
               manh_ymargin=0.25,
               manh_colorlist="[1 0 0; 1 0.5 0 ; 0 0.75 0.75; 0 0.5 0; 0.75 0 0.75; 0 0 1; 0 1 0; 0 1 1]",
               refinfo="../03.reference/condFDR/9545380.ref",
               is_reset_pruneidx=TRUE,
               randprune_repeats="default",
               pthresh=1,
               is_perform_gc=TRUE,
               is_use_standard_gc=FALSE,
               is_randprune_gc=TRUE,
               is_exclude_from_discovery=TRUE,
               mafthresh = 0.005,
               is_exclude_ambiguous_snps=TRUE,
               is_onscreen=FALSE,
               is_dummy_zscore = FALSE,
               is_exit_matlab_upon_completion = FALSE,
               matlab_path=NULL)

  #condfdr：trait -> main_trait
  pleiofdr_run(reffile="../03.reference/condFDR/ref9545380_1kgPhase3eur_LDr2p1.mat",
               traitfolder="../16.pre_data_condFDR/result/",
               traitfile1=traitfiles,
               traitname1=trait,
               traitfiles=traitfile1,
               traitnames=main_trait,
               outputdir="./result/",
               stattype="condfdr",
               fdrthresh=0.01,
               is_randprune=TRUE,
               randprune_n=500,
               exclude_chr_pos="6 25119106 33854733; 8 7200000 12500000",
               manh_fontsize_genenames=12,
               manh_yspace=0.75,
               manh_ymargin=0.25,
               manh_colorlist="[1 0 0; 1 0.5 0 ; 0 0.75 0.75; 0 0.5 0; 0.75 0 0.75; 0 0 1; 0 1 0; 0 1 1]",
               refinfo="../03.reference/condFDR/9545380.ref",
               is_reset_pruneidx=TRUE,
               randprune_repeats="default",
               pthresh=1,
               is_perform_gc=TRUE,
               is_use_standard_gc=FALSE,
               is_randprune_gc=TRUE,
               is_exclude_from_discovery=TRUE,
               mafthresh = 0.005,
               is_exclude_ambiguous_snps=TRUE,
               is_onscreen=FALSE,
               is_dummy_zscore = FALSE,
               is_exit_matlab_upon_completion = FALSE,
               matlab_path=NULL)

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

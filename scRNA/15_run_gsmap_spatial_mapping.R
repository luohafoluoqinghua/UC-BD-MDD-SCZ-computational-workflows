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

# setwd("/home/data/t230406/scRNA/35.gsMap_quick_mode/")
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

#参数设置
inputFile = c("UC", "BD", "MDD", "SCZ")  # 本稿研究对象
#hdf5_path文件路径
hdf5_path="../34.pre_gsMap/ST/E16.5_E1S1.MOSTA.h5ad" #ST文件路径
annotation='annotation'   #根据32.pre_gsMap.R文件中的注释文件进行设置
data_layer='count'        #根据32.pre_gsMap.R文件中的数据层进行设置 

# 生成gsMap输入文件格式 --------------------------------------------------------
for (trait in inputFile) {
  gsMap_format_sumstats_run(sumstats=paste0("../07.hg19Tohg38/",trait,"_hg38_qc.txt"),
                            snp="SNP",
                            a1="effect_allele",
                            a2="other_allele",
                            info=NULL,
                            beta="beta",
                            se="se",
                            p="pval",
                            frq="eaf",
                            n="N",
                            chr="CHR",
                            pos="BP",
                            dbsnp=NULL,
                            chunksize="500000",
                            format="gsMap",
                            info_min=NULL,
                            maf_min=NULL,
                            outname=trait) 
  
}


# gsMap快速模式运行 ------------------------------------------------------------
for (trait in inputFile) {
  gsMap_quick_mode_run(sample_name=gsub(".h5ad","",basename(hdf5_path)),
                       gsMap_resource_dir='gsMap_resource',
                       hdf5_path=hdf5_path,
                       annotation=annotation,        #根据32.pre_gsMap.R文件中的注释文件进行设置
                       data_layer=data_layer,        #根据32.pre_gsMap.R文件中的数据层进行设置 
                       trait_name=trait,
                       sumstats_file=paste0('./result/',trait,".sumstats.gz"),
                       sumstats_config_file=NULL,
                       homolog_file='./gsMap_resource/homologs/mouse_human_homologs.txt', 
                       max_processes=parallel::detectCores()/2, #调用CPU一半线程数
                       latent_representation=NULL,
                       num_neighbour=NULL,
                       num_neighbour_spatial=NULL,
                       gM_slices=NULL,
                       is_pearson_residuals=FALSE,
                       argsfile=NULL)
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
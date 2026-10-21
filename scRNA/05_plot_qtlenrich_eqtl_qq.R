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
Sys.setenv(http_proxy = "http://10.10.10.215:7897")
# 设置 HTTPS 代理（通常与 HTTP 代理相同，若不同需单独设置）
Sys.setenv(https_proxy = "http://10.10.10.215:7897")
#参数设置
inputFile = c("UC", "BD", "MDD", "SCZ")  # 本稿研究对象

# 循环绘制Q-Q图
tissue_data <-fread("../03.reference/QTLEnrich/gtex_tissue_names.txt",header = F)  
for (trait in inputFile) {
  for (tissue in tissue_data$V1) {
  #显示运行信息  
  message("\n当前性状：",trait,"; 当前第",which(tissue == tissue_data$V1),"个组织：",tissue)  
  #设置关键字  
  key <- paste0(trait,"_",tissue)
  #提取输入文件路径
  qtl_file <- list.files(path = "../08.QTLEnrich_eQTL/result/output_files/", 
                         pattern = paste0("^QTLEnrich_significant.*",key),
                         full.names = T)
  pvalue_null_matrix  <- list.files(path = "../08.QTLEnrich_eQTL/result/output_files/", 
                                    pattern = paste0("^QTLEnrich_matched_null.*",key),
                                    full.names = T)
  
  # 判断输入文件是否存在
  if (!file.exists(qtl_file) || !file.exists(pvalue_null_matrix) || length(qtl_file) != 1 || length(pvalue_null_matrix) != 1) {
    stop("仔细检查QTLEnrich中的分析结果")
  }
  
  #生成Q-Q图
  QTLEnrich_qqplot_run(
    qtl_file=qtl_file,
    gwas_file=paste0("../07.hg19Tohg38/",trait,"_hg38_qc_QTLEnrich.txt"),
    pvalue_null_matrix=pvalue_null_matrix,
    is_gwas_background=T,
    is_matched_null_background=T,
    num_permutations=100,
    is_confidence_interval=T,
    upper_bound_ci=97.5,
    lower_bound_ci=2.5,
    upper_bound_clip=0,
    trait_label=trait,
    qtl_type="eQTLs",
    tissue_label=tissue,
    is_show_title=TRUE,
    qqplot_title=paste0(trait," GWAS"),
    flag_options=NULL,
    x_inches=12,
    y_inches=10,
    dpi=1000,
    exp_label=trait,
    is_help=F,
    force=F)

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
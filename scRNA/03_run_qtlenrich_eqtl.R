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
setwd("/home/data/t230406/scRNA/08.QTLEnrich_eQTL/")

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

# GTEx_Analysis_v8_eQTL --------------------------------------------------------
# 循环进行
for (trait in inputFile) {
  for (qtl_type in c("best_eqtl")) {
  QTLEnrich_run(gwas_file=paste0("../07.hg19Tohg38/",trait,"_hg38_qc_QTLEnrich.txt"),
                qtl_directory="../03.reference/QTLEnrich/GTEx_Analysis_v8_eQTL/",
                file_name=".v8.egenes.txt.gz",
                qtl_type=qtl_type,
                confounders_table="../03.reference/QTLEnrich/GTEx_v8_Confounders_Table_49Tiss_Global.txt.gz",
                null_table="../03.reference/QTLEnrich/GTEx_v8_Null_Table_49Tiss_Global.txt.gz",
                is_subset_genes=TRUE,
                is_compute_tss_distance=TRUE,
                is_keep_gene_id_suffix=TRUE,
                is_keep_null_variants=TRUE,
                is_keep_pvalue_matrix=TRUE,
                gencode_file="../03.reference/QTLEnrich/GENCODE_26_df.txt",
                null_option=NULL,
                genome_build=NULL,
                trait_name=trait,
                qtl_q_value=NULL,
                gwas_p_value=NULL,
                independent_ranking=NULL,
                lambda_factor=NULL,
                lower_bound_permutations=NULL,
                upper_bound_permutations=NULL,
                num_quantiles=NULL,
                tissue_names=NULL,
                is_subset_tissues=FALSE,
                exp_label=trait,
                is_GeneEnrich_input=TRUE,
                eGenes_directory=NULL,
                eGenes_file_name=NULL,
                flag_options=NULL,
                is_help=FALSE,
                force=FALSE)
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
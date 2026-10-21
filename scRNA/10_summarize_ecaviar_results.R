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
setwd("/home/data/t230406/scRNA/16.post_eCAVIAR/")

#清除一切变量
rm(list=ls())
#释放不再使用的内存
gc()
#显式路径
getwd()
#显示R包安装位置
.libPaths()

#下载以及加载包
library(xQTLbiolinks)
library(Fast2TWAS)
library(data.table)

#参数设置
inputFile = c("UC", "BD", "MDD", "SCZ")  # 本稿研究对象
#指定lead SNP的上下游范围
wind_size = 50000
#基因组参考人群
ref_pop="g1000_eur"    #c("g1000_eur","g1000_afr","g1000_amr","g1000_eas","g1000_sas")
#读入来自GTEx v8的组织列表文件
gtex_tissue_names <- fread("../03.reference/QTLEnrich/gtex_tissue_names.txt",header = F)$V1

################################################################################
##########################循环合并eCAVIAR分析完成之后的数据#####################
################################################################################
write(x = "", file = "16.post_eCAVIAR.log", append = F)
for (trait in inputFile) {
  #读入13.fuma文件夹中的stable6.txt文件
  stable6 <- fread(paste0("../13.fuma/stable6.txt"))
  #筛选性状
  stable6_filter <- subset(stable6, Trait == trait)
  #读入GWAS数据
  gwas_data <- fread(paste0("../07.hg19Tohg38/",trait,"_hg38_qc.txt"))
  #设置空数据框
  all_result <- data.frame()
  
  #按照stable6_filter中的rsID列进行循环
  for (snp in unique(stable6_filter$rsID)) {
    #eQTL和sQTL循环
    for (qtl_type in c("eQTL","sQTL")) {
      #gtex_tissue循环
      for (tissue in gtex_tissue_names) {
        #显示运行信息  
        message("\n当前性状：",trait,"，当前QTL类型：",qtl_type,"，当前组织：",tissue)
        eCAVIAR_col_file <- paste0("../15.eCAVIAR/result/",trait,"/",qtl_type,"/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".eCAVIAR_col")
        eCAVIAR_qtl_file <- paste0("../15.eCAVIAR/result/",trait,"/",qtl_type,"/",snp,"_",trait,"_",wind_size/1000,"kb_",ref_pop,"_markers_QTL_gwas_",tissue,"_",qtl_type,".eCAVIAR_2_set")
        #判断文件是否存在
        if (file.exists(eCAVIAR_col_file) && file.exists(eCAVIAR_qtl_file)) {
          #读入eCAVIAR_col文件
          eCAVIAR_col <- fread(eCAVIAR_col_file)
          #读入eCAVIAR_qtl文件
          eCAVIAR_qtl <- fread(eCAVIAR_qtl_file,header = F)
          #从eCAVIAR_col筛选出来
          eCAVIAR_col_filter <- subset(eCAVIAR_col, SNP_ID %in% eCAVIAR_qtl$V1)

          # 定义一个函数来执行查询并处理错误
          query_with_retry <- function(rsID_qtl, max_retries = 99) {
            retries <- 0
            while (retries < max_retries) {
              tryCatch({
                varAsso_id <- xQTLquery_varId(rsID_qtl)
                return(varAsso_id)  # 如果成功，返回结果
              }, error = function(e) {
                message("Error encountered: ", e$message)
                retries <<- retries + 1
                message("Retrying... Attempt ", retries, " of ", max_retries)
                Sys.sleep(5)  # 等待5秒后重试
              })
            }
            # return(NULL)  # 如果达到最大重试次数，返回NULL
            stop("Failed after ", max_retries, " attempts")  # 如果重试次数用完，抛出错误
          }
          
          #设置空数据框
          all_varAsso_id <- data.frame()
          for (rsID_qtl in eCAVIAR_col_filter$SNP_ID) {
            varAsso_id <-query_with_retry(rsID_qtl)
            all_varAsso_id <- rbind(all_varAsso_id, varAsso_id)
          }
          if (nrow(all_varAsso_id)==0) next
          all_varAsso_id <- subset(all_varAsso_id, select=c(variantId, snpId,ref,alt))

          #读入QTL数据
          if (qtl_type=="eQTL") {
            qtl_data_signif_variant_file <- paste0("../03.reference/QTLEnrich/GTEx_Analysis_v8_",qtl_type,"/",tissue,".v8.signif_variant_gene_pairs.txt.gz")
            qtl_data_signif_variant <- fread(qtl_data_signif_variant_file)
            
            #筛选all_varAsso数据的列
            qtl_data_signif_variant_fileter <- subset(qtl_data_signif_variant,variant_id %in% all_varAsso_id$variantId,
                                                      select=c(variant_id,gene_id,pval_nominal,slope,slope_se))
            qtl_data_signif_variant_fileter$Gene_id_intron_cluster <- ""
            }else if(qtl_type=="sQTL"){
            qtl_data_signif_variant_file <- paste0("../03.reference/QTLEnrich/GTEx_Analysis_v8_",qtl_type,"/",tissue,".v8.sqtl_signifpairs.txt.gz")
            qtl_data_signif_variant <- fread(qtl_data_signif_variant_file) 
            setnames(qtl_data_signif_variant,"phenotype_id","Gene_id_intron_cluster")
            qtl_data_signif_variant <- qtl_data_signif_variant %>% tidyr::separate(col = Gene_id_intron_cluster,
                                                        into = c("chr", "start", "end", "clu","gene_id"),
                                                        sep = ":",
                                                        remove = FALSE  # 保留原始Gene_id_intron_cluster列
                                                      )
            #筛选all_varAsso数据的列
            qtl_data_signif_variant_fileter <- subset(qtl_data_signif_variant,variant_id %in% all_varAsso_id$variantId,
                                                      select=c(variant_id,Gene_id_intron_cluster,gene_id,pval_nominal,slope,slope_se))

          }

          #合并qtl_data_signif_variant_fileter与all_varAsso_id
          all_varAsso_id_qtl_data_signif_variant_fileter <- merge(qtl_data_signif_variant_fileter, all_varAsso_id, by.x="variant_id", by.y="variantId")
          
          #将eCAVIAR_col_filter和all_varAsso_id_qtl_data_signif_variant_fileter合并
          eCAVIAR_result <- merge(eCAVIAR_col_filter, all_varAsso_id_qtl_data_signif_variant_fileter, by.x="SNP_ID", by.y="snpId")
          #eCAVIAR_result添加一列，列名为tissue
          eCAVIAR_result$Tissue <- tissue
          #合并eCAVIAR_result和GTEXv8_ensTOsymble
          GTEXv8_ensTOsymble <- GTEXv8_ensTOsymble[,c("gene_id","symble")]
          #去掉eCAVIAR_result中gene_id列中的.后面的字符，然后与GTEXv8_ensTOsymble中的数据合并
          eCAVIAR_result$gene_id <- gsub("\\..*", "", eCAVIAR_result$gene_id)
          GTEXv8_ensTOsymble_eCAVIAR_result <- merge(eCAVIAR_result, GTEXv8_ensTOsymble, by.x="gene_id", by.y="gene_id")
          #从GWAS数据中筛选SNP_ID列中的行
          gwas_data_filter <- subset(gwas_data, SNP %in% unique(eCAVIAR_result$SNP_ID),select = c(SNP,effect_allele,other_allele,pval,beta,se))
          #合并gwas_data_filter和GTEXv8_ensTOsymble_eCAVIAR_result
          GTEXv8_ensTOsymble_eCAVIAR_result_final <- merge(GTEXv8_ensTOsymble_eCAVIAR_result, gwas_data_filter, by.x="SNP_ID", by.y="SNP")
          #添加Qtl_gwas_relative_effect，判断QTL和GWAS的效应等位基因是否相同以便判断slope和beta的正负是否一致
          
          for (i in 1:nrow(GTEXv8_ensTOsymble_eCAVIAR_result_final)) {
            #判断QTL和GWAS的效应等位基因是否相同
            if (GTEXv8_ensTOsymble_eCAVIAR_result_final$effect_allele[i]==GTEXv8_ensTOsymble_eCAVIAR_result_final$alt[i] & 
                GTEXv8_ensTOsymble_eCAVIAR_result_final$other_allele[i]==GTEXv8_ensTOsymble_eCAVIAR_result_final$ref[i]) {
              #判断slope和beta的正负是否一致
              if (GTEXv8_ensTOsymble_eCAVIAR_result_final$slope[i]>0 & GTEXv8_ensTOsymble_eCAVIAR_result_final$beta[i]>0) {
                GTEXv8_ensTOsymble_eCAVIAR_result_final$Qtl_gwas_relative_effect[i] <- "SAME"
              } else if (GTEXv8_ensTOsymble_eCAVIAR_result_final$slope[i]<0 & GTEXv8_ensTOsymble_eCAVIAR_result_final$beta[i]<0) {
                GTEXv8_ensTOsymble_eCAVIAR_result_final$Qtl_gwas_relative_effect[i] <- "SAME"
              } else {
                GTEXv8_ensTOsymble_eCAVIAR_result_final$Qtl_gwas_relative_effect[i] <- "OPPOSITE"
              }
            }else if(GTEXv8_ensTOsymble_eCAVIAR_result_final$effect_allele[i]==GTEXv8_ensTOsymble_eCAVIAR_result_final$ref[i]  & 
                    GTEXv8_ensTOsymble_eCAVIAR_result_final$other_allele[i]==GTEXv8_ensTOsymble_eCAVIAR_result_final$alt[i]) {
              #判断slope和beta的正负是否一致
              if (GTEXv8_ensTOsymble_eCAVIAR_result_final$slope[i]>0 & GTEXv8_ensTOsymble_eCAVIAR_result_final$beta[i]>0) {
                GTEXv8_ensTOsymble_eCAVIAR_result_final$Qtl_gwas_relative_effect[i] <- "OPPOSITE"
              } else if (GTEXv8_ensTOsymble_eCAVIAR_result_final$slope[i]<0 & GTEXv8_ensTOsymble_eCAVIAR_result_final$beta[i]<0) {
                GTEXv8_ensTOsymble_eCAVIAR_result_final$Qtl_gwas_relative_effect[i] <- "OPPOSITE"
              } else {
                GTEXv8_ensTOsymble_eCAVIAR_result_final$Qtl_gwas_relative_effect[i] <- "SAME"
              }
            }else{
              GTEXv8_ensTOsymble_eCAVIAR_result_final$Qtl_gwas_relative_effect[i] <- "NA"
          }
          }
          
          #从stable6_filter筛选出rsID等于snp的行然后与GTEXv8_ensTOsymble_eCAVIAR_result_final合并
          stable6_filter_snp <- subset(stable6_filter, rsID==snp, select=-c(hg19,CHR,N))
          #给GTEXv8_ensTOsymble_eCAVIAR_result_final和stable6_filter_snp重新命名列
          setnames(GTEXv8_ensTOsymble_eCAVIAR_result_final, "SNP_ID", "rsID_qtl")
          setnames(GTEXv8_ensTOsymble_eCAVIAR_result_final, "variant_id", "qtl_Variant_id")
          setnames(GTEXv8_ensTOsymble_eCAVIAR_result_final, "symble", "gene_symble")
          setnames(GTEXv8_ensTOsymble_eCAVIAR_result_final, "pval_nominal", "QTL_pval_nominal")
          setnames(GTEXv8_ensTOsymble_eCAVIAR_result_final, "pval", "QTL_gwas_p_value")
          setnames(GTEXv8_ensTOsymble_eCAVIAR_result_final, "alt", "qtl_variant_ALT_(effect_allele)")
          setnames(GTEXv8_ensTOsymble_eCAVIAR_result_final, "ref", "qtl_variant_REF_(other_allele)")
          setnames(GTEXv8_ensTOsymble_eCAVIAR_result_final, "slope", "QTL_slope_relative_ALT")
          setnames(GTEXv8_ensTOsymble_eCAVIAR_result_final, "slope_se", "QTL_se")
          setnames(GTEXv8_ensTOsymble_eCAVIAR_result_final, "beta", "GWAS_effect_relative_ALT")
          setnames(GTEXv8_ensTOsymble_eCAVIAR_result_final, "se", "GWAS_se_relative_ALT")
          setnames(GTEXv8_ensTOsymble_eCAVIAR_result_final, "effect_allele", "GWAS_variant_ALT_(effect_allele)")
          setnames(GTEXv8_ensTOsymble_eCAVIAR_result_final, "other_allele", "GWAS_variant_REF_(other_allele)")
          
          setnames(stable6_filter_snp, "rsID", "rsID_GWAS")
          setnames(stable6_filter_snp, "GWAS_locus", "GWAS_variant")
          setnames(stable6_filter_snp, "beta", "GWAS_effect")
          setnames(stable6_filter_snp, "se", "GWAS_se")
          setnames(stable6_filter_snp, "pval", "GWAS_p_value")
          setnames(stable6_filter_snp, "effect_allele", "GWAS_effect_allele")
          setnames(stable6_filter_snp, "other_allele", "GWAS_other_allele")
          
          #设置stable6_filter_snp的上下游区间
          stable6_filter_snp$startbp <- stable6_filter_snp$hg38 - wind_size
          stable6_filter_snp$endbp <- stable6_filter_snp$hg38 + wind_size
          stable6_filter_snp$Number_variants <- nrow(eCAVIAR_col)
          #将GTEXv8_ensTOsymble_eCAVIAR_result_final和stable6_filter_snp拼接在一起
          stable6_filter_snp_repeated  <- stable6_filter_snp[rep(1,nrow(GTEXv8_ensTOsymble_eCAVIAR_result_final)),]
          GTEXv8_ensTOsymble_eCAVIAR_result_final_stable6 <- cbind(stable6_filter_snp_repeated,GTEXv8_ensTOsymble_eCAVIAR_result_final)
          #添加QTL Type列
          GTEXv8_ensTOsymble_eCAVIAR_result_final_stable6$QTL_Type <- qtl_type
          #
          all_result <- rbind(all_result,GTEXv8_ensTOsymble_eCAVIAR_result_final_stable6)
          
          } else {
          line=paste0("文件不存在：", eCAVIAR_col_file, "或", eCAVIAR_qtl_file)
          write(x = line, file = "16.post_eCAVIAR.log", append = T)
        }
      }
    }
  }
  #将结果保存到文件中
  fwrite(all_result, 
         file=paste0("stable_eCAVIAR_",trait,".txt"), 
         sep = "\t")
}  


for (trait in inputFile) {
  data <- fread(paste0("stable_eCAVIAR_",trait,".txt"))
  #筛选CLPP>0.01的数据
  data_filter <- subset(data,CLPP>0.01)
  #将结果保存到文件中
  fwrite(data_filter, 
         file=paste0("stable_eCAVIAR_",trait,".sig.txt"), 
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
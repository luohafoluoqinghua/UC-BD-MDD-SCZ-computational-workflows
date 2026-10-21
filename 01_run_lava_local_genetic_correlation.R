
#设置当前工作路径
setwd("~/GWAS_analysis/Lava")

#清除一切变量
rm(list=ls())
#释放不再使用的内存
gc()
#显式路径
getwd()
#BiocManager::install("snpStats")
#devtools::install_github('josefin-werme/LAVA')
library(LAVA)

### Read in summary statistics and related info
input = process.input(input.info.file="LAVA_info.txt",           # input info file
                      sample.overlap.file="sample.overlap.txt",   # sample overlap file (can be set to NULL if there is no overlap)
                      ref.prefix="./g1000_eur",                    # reference genotype data prefix
                      phenos=c("UC","BD","MDD","SCZ"))       # subset of phenotypes listed in the input info file that we want to process

loci = read.loci("./blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile")

head(loci)

for (i in 1:2495) {
  print(c(i))
  data = process.locus(loci[i,], input)
  if (class(data) == "NULL") {
    print(c("skip",i))
  }else{
    bivar = run.univ.bivar(data,target = "UC")# 修改target pheno
    save(bivar,file=paste0("./",i,".R"))
  }
}

#rbind LAVA
rm(list = ls())

library(LAVA)

setwd("./results/")

file = list.files("./")

PS = data.frame()

for (i in file) {
  load(i)
  biv = bivar$bivar
  if (class(biv)=="data.frame") {
    biv$locus = i
    PS = rbind(PS,biv)
  }
  print(i)
}

write.table(PS,"../PS.csv",row.names = F,quote = F,sep = ",")
#####################
### LAVA_FDR校正最终版
#####################
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
library(dplyr)
setwd("../")
read.csv('./PS.csv')->filterd_ALL
adjusted_p_values=p.adjust(filterd_ALL$p, method = "BH")
filterd_ALL$adjusted_p_values=adjusted_p_values
after_FDR <- dplyr::filter(filterd_ALL, adjusted_p_values < 0.05)
write.table(filterd_ALL,'./PS_FDR.xlsx', row.names=F, col.names=T, sep='\t')
up_FDR <- dplyr::filter(filterd_ALL, rho > 0)
down_FDR <- dplyr::filter(filterd_ALL, rho < 0)

# combined
filterd_ALL <-read.csv('./PS_FDR.csv')
head(filterd_ALL)
ref <- read.csv('./blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile',sep=' ')
head(ref)
final <- merge(filterd_ALL,ref)
write.csv(final,file = 'PS_LAVA_results.csv',quote = F)





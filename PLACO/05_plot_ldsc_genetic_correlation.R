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
#setwd("D:/myproject/mylesson/6.shared_genetic_local_genetic_correlation/10.plot_genetic_correlation_ldsc/")
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
library(ggplot2)
#输入参数-----------------------------------------------------------------------
#设置棒棒糖图的长和宽
width=9
height=7

#读入index文件
index <- data.table::fread("../07.export_index/gwasfile_info.txt")

# 仅使用本稿四个研究对象，按UC、BD、MDD、SCZ排序。
study_traits <- c("UC", "BD", "MDD", "SCZ")
if (anyDuplicated(index$phenotype) || !all(study_traits %in% index$phenotype)) {
  stop("Regenerate 07.export_index for UC, BD, MDD and SCZ before running this step.")
}
index <- index[match(study_traits, index$phenotype), ]
if (!identical(as.character(index$name), study_traits)) {
  stop("Study input filenames do not match the four-trait configuration; rerun 07.export_index.")
}


#指定性状名称
trait="UC"

#性状名
inputFile = c("BD", "MDD", "SCZ")

#整理ldsc_rg.log文件------------------------------------------------------------
inputFile <- c(index$phenotype[index$name==trait],index$phenotype[index$name!=trait])
# 指定文件路径
file_path <- paste0("../09.genetic_correlation_ldsc/result/",index$phenotype[index$name==trait],"_ldsc_rg.log")
# 逐行读取文件
lins <- readLines(file_path)

##编辑文件
write(x = "Trait1  Trait2  rg  se  z  p  h2_obs  h2_obs_se  h2_int  h2_int_se  gcov_int  gcov_int_se", file = "rg.txt", append = F)
for (line in lins[18:(length(lins)-3)]) {
  # 使用grep判断是否包含
  if (sum(sapply(inputFile, grepl, line)) >= 2) {
    line=gsub(" ../08.heritability_ldsc/result/|.sumstats.gz","",line)
    write(x = line, file = "rg.txt", append = T)
    print(line)
  } 
}

rg_dat=data.table::fread("rg.txt")
#添加fdr与bonferroni矫正
# p.adjust.methods
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")
rg_dat$FDR=p.adjust(p = rg_dat$p, method = "fdr")
rg_dat$bonferroni=p.adjust(p = rg_dat$p, method = "bonferroni")

#保存数据
data.table::fwrite(rg_dat,file = paste0(index$phenotype[index$name==trait],"_rg_fdr.txt"),
  sep = "\t",
  na = "NA",
  quote = FALSE)

#删除rg.txt
# unlink("rg.txt")

# 绘制棒棒糖图 -----------------------------------------------------------------
inputFile="rg.txt"
#读取输入文件
data = data.table::fread(inputFile)
# data = read.table(inputFile, header=T, sep="\t", check.names=F)

#删除NA
data=na.omit(data)

#定义圆圈颜色的函数
# p.col = c("#FFFF80","#FFFF00","#FFAA00","#FF5500","#FF0000")
p.col = c("#FDE725","#5DC863","#21908C","#3B528B","#440154")
fcolor = function(x,p.col){
  color = ifelse(x>0.8,p.col[1],ifelse(x>0.6,p.col[2],ifelse(x>0.4,p.col[3],ifelse(x>0.2,p.col[4], p.col[5])
  )))
  return(color)}

#定义设置圆圈大小的函数
p.cex = seq(2.5, 5.5, length=5)
fcex = function(x){
  x=abs(x)
  cex = ifelse(x<0.1,p.cex[1],ifelse(x<0.2,p.cex[2],ifelse(x<0.3,p.cex[3],ifelse(x<0.4,p.cex[4],p.cex[5]))))
  return(cex)}

#根据相关性检验的pvalue定义圆圈颜色
data$points.color = fcolor(x=data$p,p.col=p.col)

#根据相关系数定义圆圈大小
data$points.cex = fcex(x=data$rg)
data=data[order(data$rg),]

#生成PDF文件
pdf(file="ldsc_rg.pdf", width=width, height=height)
xlim = ceiling(max(abs(na.omit(data$rg)))*10)/10                    #定义x轴范围
layout(mat=matrix(c(1,1,1,1,1,0,2,0,3,0),nc=2),width=c(8,2.2),heights=c(1,2,1,2,1))
par(bg="white",las=1,mar=c(5,18,2,4),cex.axis=1.5,cex.lab=2)
plot(1,type="n",xlim=c(-xlim,xlim),ylim=c(0.5,nrow(data)+0.5),xlab="Genetic Correlation",ylab="",yaxt="n",yaxs="i",axes=F)
rect(par('usr')[1],par('usr')[3],par('usr')[2],par('usr')[4],col="#F5F5F5",border="#F5F5F5")
grid(ny=nrow(data),col="white",lty=1,lwd=2)
#绘制图形的线段
segments(x0=data$rg,y0=1:nrow(data),x1=0,y1=1:nrow(data),lwd=4)
#绘制图形的圆圈
points(x=data$rg,y = 1:nrow(data),col = data$points.color,pch=16,cex=data$points.cex)
#在图形中展示Trait2的名称
text(par('usr')[1],1:nrow(data),data$Trait2,adj=1,xpd=T,cex=1.5)
#展示相关性检验的p
pvalue.text=ifelse(data$p<0.001,'<0.001',sprintf("%.03f",data$p))
redcutoff_cor=0
redcutoff_pvalue=0.05
text(par('usr')[2],1:nrow(data),pvalue.text,adj=0,xpd=T,col=ifelse(abs(data$rg)>redcutoff_cor & data$p<redcutoff_pvalue,"red","black"),cex=1.5)
axis(1,tick=F)

#绘制圆圈大小的图例
par(mar=c(0,4,3,4))
plot(1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=c(0.1,0.2,0.3,0.4,0.5),col="black",pt.cex=p.cex,pch=16,bty="n",cex=2,title="abs(rg)")

#绘制圆圈颜色的图例
par(mar=c(0,6,4,6),cex.axis=1.5,cex.main=2)
barplot(rep(1,5),horiz=T,space=0,border=NA,col=p.col,xaxt="n",yaxt="n",xlab="",ylab="",main="pvalue")
axis(4,at=0:5,c(1,0.8,0.6,0.4,0.2,0),tick=F)
dev.off()



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
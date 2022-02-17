library(ggplot2)
library(ggpubr)
library(zoo)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(stringr)
library(data.table)

# setwd("C:/Users/guyha/Desktop/bisulfite sequencing")
setwd("C:/Users/owner/Desktop/bisulfite_sequencing")

RRBS_10R_1n_rep1=read.delim("1n_10R_meth.sorted_CpG.bedGraph", header = F)[-1,]
RRBS_10G_2n_rep1=read.delim("2n_10G_meth.sorted_CpG.bedGraph", header = F)[-1,]
RRBS_L_3n=read.delim("3n_L_meth.sorted_CpG.bedGraph", header = F)[-1,]

RRBS_10R_1n_rep2=read.delim("1n_10R_meth_02_10.sorted_CpG.bedGraph", header = F)[-1,]
RRBS_10R_2n_rep2=read.delim("2n_10R_meth_02_10.sorted_CpG.bedGraph", header = F)[-1,]
RRBS_K_3n=read.delim("3n_clone_K_meth_02_10.sorted_CpG.bedGraph", header = F)[-1,]
RRBS_H_3n=read.delim("3n_clone_H_meth_02_10.sorted_CpG.bedGraph", header = F)[-1,]

SwapS4=read.delim("SwapS4.bedGraph", header = F)[-1,]


colnames(RRBS_10R_1n_rep1)=c("chr","start","end","meth_percentage","num_meth","num_un_meth")
colnames(RRBS_10G_2n_rep1)=c("chr","start","end","meth_percentage","num_meth","num_un_meth")
colnames(RRBS_L_3n)=c("chr","start","end","meth_percentage","num_meth","num_un_meth")
colnames(RRBS_10R_1n_rep2)=c("chr","start","end","meth_percentage","num_meth","num_un_meth")
colnames(RRBS_10R_2n_rep2)=c("chr","start","end","meth_percentage","num_meth","num_un_meth")
colnames(RRBS_K_3n)=c("chr","start","end","meth_percentage","num_meth","num_un_meth")
colnames(RRBS_H_3n)=c("chr","start","end","meth_percentage","num_meth","num_un_meth")
colnames(SwapS4)=c("chr","start","end","meth_percentage","num_meth","num_un_meth")

RRBS_10R_1n_rep1$location=paste(RRBS_10R_1n_rep1[,1],RRBS_10R_1n_rep1[,2],sep="-")
RRBS_10G_2n_rep1$location=paste(RRBS_10G_2n_rep1[,1],RRBS_10G_2n_rep1[,2],sep="-")
RRBS_L_3n$location=paste(RRBS_L_3n[,1],RRBS_L_3n[,2],sep="-")
RRBS_10R_1n_rep2$location=paste(RRBS_10R_1n_rep2[,1],RRBS_10R_1n_rep2[,2],sep="-")
RRBS_10R_2n_rep2$location=paste(RRBS_10R_2n_rep2[,1],RRBS_10R_2n_rep2[,2],sep="-")
RRBS_K_3n$location=paste(RRBS_K_3n[,1],RRBS_K_3n[,2],sep="-")
RRBS_H_3n$location=paste(RRBS_H_3n[,1],RRBS_H_3n[,2],sep="-")
SwapS4$location=paste(SwapS4[,1],SwapS4[,2],sep="-")

RRBS_10R_1n_rep1=RRBS_10R_1n_rep1[(RRBS_10R_1n_rep1$num_meth+RRBS_10R_1n_rep1$num_un_meth)>=10,]
RRBS_10G_2n_rep1=RRBS_10G_2n_rep1[(RRBS_10G_2n_rep1$num_meth+RRBS_10G_2n_rep1$num_un_meth)>=10,]
RRBS_L_3n=RRBS_L_3n[(RRBS_L_3n$num_meth+RRBS_L_3n$num_un_meth)>=10,]
RRBS_10R_1n_rep2=RRBS_10R_1n_rep2[(RRBS_10R_1n_rep2$num_meth+RRBS_10R_1n_rep2$num_un_meth)>=10,]
RRBS_10R_2n_rep2=RRBS_10R_2n_rep2[(RRBS_10R_2n_rep2$num_meth+RRBS_10R_2n_rep2$num_un_meth)>=10,]
RRBS_K_3n=RRBS_K_3n[(RRBS_K_3n$num_meth+RRBS_K_3n$num_un_meth)>=10,]
RRBS_H_3n=RRBS_H_3n[(RRBS_H_3n$num_meth+RRBS_H_3n$num_un_meth)>=10,]
SwapS4=SwapS4[(SwapS4$num_meth+SwapS4$num_un_meth)>=10,]

gtf=read.delim("C:/Users/owner/Desktop/bisulfite sequencing/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf", skip=5, header = F)
gtf=gtf[gtf[,3]=="gene",]
gtf$gene_name=substring(gtf$V9, regexpr("gene_id", gtf$V9)+8,regexpr("; d", gtf$V9)-1)
colnames(gtf)[c(1,4:5)]=c("chr","start","end")


if (!file.exists("methylation per promoter of gene.csv")){
  
  promoter_location=gtf[,c(10,1,4,5,7)]
  colnames(promoter_location)[2:5]=c("chr","start","end", "strand")
  promoter_location$promoter_start[promoter_location$strand=="+"]=paste(promoter_location[promoter_location$strand=="+",3]-1000)
  promoter_location$promoter_end[promoter_location$strand=="+"]=paste(promoter_location[promoter_location$strand=="+",3]+1000)
  
  promoter_location$promoter_start[promoter_location$strand=="-"]=paste(promoter_location[promoter_location$strand=="-",4]-1000)
  promoter_location$promoter_end[promoter_location$strand=="-"]=paste(promoter_location[promoter_location$strand=="-",4]+1000)
  
  
  promoter_location=promoter_location[order(promoter_location$chr),]
  promoter_location=promoter_location[promoter_location$chr %in% paste("chr", c(1:22,"X"), sep=""),]
  
  
  for (i in 1:23){
    chr=paste("chr", i, sep="")
    if (i==23) chr="chrX"
    sub_chr_10R_1n_rep1=RRBS_10R_1n_rep1[RRBS_10R_1n_rep1$chr==chr, ]
    sub_chr_10R_1n_rep2=RRBS_10R_1n_rep2[RRBS_10R_1n_rep2$chr==chr, ]
    sub_chr_10G_2n_rep1=RRBS_10G_2n_rep1[RRBS_10G_2n_rep1$chr==chr, ]
    sub_chr_10R_2n_rep2=RRBS_10R_2n_rep2[RRBS_10R_2n_rep2$chr==chr, ]
    sub_chr_L_3n=RRBS_L_3n[RRBS_L_3n$chr==chr, ]    
    sub_chr_K_3n=RRBS_K_3n[RRBS_K_3n$chr==chr, ]    
    sub_chr_H_3n=RRBS_H_3n[RRBS_H_3n$chr==chr, ]    
    sub_chr_SwapS4=SwapS4[SwapS4$chr==chr, ]
    sub_promoter_location=promoter_location[promoter_location$chr==chr,]
    for (j in 1:nrow(sub_promoter_location)){
      start=sub_promoter_location$promoter_start[j]
      end=sub_promoter_location$promoter_end[j]
      promoter_location$RRBS_10R_1n_rep1[which(promoter_location$chr==chr)[1]+j-1]=mean(sub_chr_10R_1n_rep1$meth_percentage[sub_chr_10R_1n_rep1$start%in%start:end])
      promoter_location$RRBS_10R_1n_rep2[which(promoter_location$chr==chr)[1]+j-1]=mean(sub_chr_10R_1n_rep2$meth_percentage[sub_chr_10R_1n_rep2$start%in%start:end])
      promoter_location$RRBS_10G_2n_rep1[which(promoter_location$chr==chr)[1]+j-1]=mean(sub_chr_10G_2n_rep1$meth_percentage[sub_chr_10G_2n_rep1$start%in%start:end])
      promoter_location$RRBS_10R_2n_rep2[which(promoter_location$chr==chr)[1]+j-1]=mean(sub_chr_10R_2n_rep2$meth_percentage[sub_chr_10R_2n_rep2$start%in%start:end])
      promoter_location$RRBS_L_3n[which(promoter_location$chr==chr)[1]+j-1]=mean(sub_chr_L_3n$meth_percentage[sub_chr_L_3n$start%in%start:end])
      promoter_location$RRBS_K_3n[which(promoter_location$chr==chr)[1]+j-1]=mean(sub_chr_K_3n$meth_percentage[sub_chr_K_3n$start%in%start:end])
      promoter_location$RRBS_H_3n[which(promoter_location$chr==chr)[1]+j-1]=mean(sub_chr_H_3n$meth_percentage[sub_chr_H_3n$start%in%start:end])
      promoter_location$SwapS4[which(promoter_location$chr==chr)[1]+j-1]=mean(sub_chr_SwapS4$meth_percentage[sub_chr_SwapS4$start%in%start:end])
      
    }
  }
  write.csv(promoter_location, "methylation per promoter of gene.csv", row.names = F)
} else promoter_location=read.csv("methylation per promoter of gene.csv")


if (!file.exists("methylation per gene.csv")){
  
  location=gtf[,c(10,1,4,5,7)]
  colnames(location)[2:5]=c("chr","start","end", "strand")
  
  location=location[order(location$chr),]
  location=location[location$chr %in% paste("chr", c(1:22,"X"), sep=""),]
  
  
  for (i in 1:23){
    chr=paste("chr", i, sep="")
    if (i==23) chr="chrX"
    sub_chr_10R_1n_rep1=RRBS_10R_1n_rep1[RRBS_10R_1n_rep1$chr==chr, ]
    sub_chr_10R_1n_rep2=RRBS_10R_1n_rep2[RRBS_10R_1n_rep2$chr==chr, ]
    sub_chr_10G_2n_rep1=RRBS_10G_2n_rep1[RRBS_10G_2n_rep1$chr==chr, ]
    sub_chr_10R_2n_rep2=RRBS_10R_2n_rep2[RRBS_10R_2n_rep2$chr==chr, ]
    sub_chr_L_3n=RRBS_L_3n[RRBS_L_3n$chr==chr, ]    
    sub_chr_K_3n=RRBS_K_3n[RRBS_K_3n$chr==chr, ]    
    sub_chr_H_3n=RRBS_H_3n[RRBS_H_3n$chr==chr, ]    
    sub_chr_SwapS4=SwapS4[SwapS4$chr==chr, ] 
    sub_location=location[location$chr==chr,]
    for (j in 1:nrow(sub_location)){
      start=sub_location$start[j]
      end=sub_location$end[j]
      location$RRBS_10R_1n_rep1[which(location$chr==chr)[1]+j-1]=mean(sub_chr_10R_1n_rep1$meth_percentage[sub_chr_10R_1n_rep1$start%in%start:end])
      location$RRBS_10R_1n_rep2[which(location$chr==chr)[1]+j-1]=mean(sub_chr_10R_1n_rep2$meth_percentage[sub_chr_10R_1n_rep2$start%in%start:end])
      location$RRBS_10G_2n_rep1[which(location$chr==chr)[1]+j-1]=mean(sub_chr_10G_2n_rep1$meth_percentage[sub_chr_10G_2n_rep1$start%in%start:end])
      location$RRBS_10R_2n_rep2[which(location$chr==chr)[1]+j-1]=mean(sub_chr_10R_2n_rep2$meth_percentage[sub_chr_10R_2n_rep2$start%in%start:end])
      location$RRBS_L_3n[which(location$chr==chr)[1]+j-1]=mean(sub_chr_L_3n$meth_percentage[sub_chr_L_3n$start%in%start:end])
      location$RRBS_K_3n[which(location$chr==chr)[1]+j-1]=mean(sub_chr_K_3n$meth_percentage[sub_chr_K_3n$start%in%start:end])
      location$RRBS_H_3n[which(location$chr==chr)[1]+j-1]=mean(sub_chr_H_3n$meth_percentage[sub_chr_H_3n$start%in%start:end])
      location$SwapS4[which(location$chr==chr)[1]+j-1]=mean(sub_chr_SwapS4$meth_percentage[sub_chr_SwapS4$start%in%start:end])
    }
  }
  write.csv(location, "methylation per gene.csv", row.names = F)
} else ocation=read.csv("methylation per gene.csv")


promoter_meth <- read.csv("methylation per promoter of gene.csv")

promoter_meth$chr <- str_remove(promoter_meth$chr, "chr")
promoter_meth$chr[promoter_meth$chr=="X"] <- 23
promoter_meth$chr <- as.numeric(promoter_meth$chr)
promoter_meth <- na.omit(promoter_meth)
promoter_meth=promoter_meth[order(promoter_meth$chr, promoter_meth$start),]

meth_DMLs_strict_E <- promoter_meth[c(2,3,8:14,16)] #choose either meth or meth_DMLs_strict for unfiltered CpGs or filtered ones
ploidy <- c("Haploid rep1","Haploid rep2","Diploid rep1","Diploid rep2","Triploid (L)","Triploid (K)","Triploid (H)","Diploid Xi")
names(ploidy) <- c("forestgreen","green3","royalblue3","royalblue1","red1","red2","red3","blue4")

# meth_DMLs_strict_E$chr[meth_DMLs_strict_E$chr=="X"] <- 23
# meth_DMLs_strict_E$chr <- as.numeric(meth_DMLs_strict_E$chr)
# meth_DMLs_strict_E <- na.omit(meth_DMLs_strict_E)
# meth_DMLs_strict_E=meth_DMLs_strict_E[order(meth_DMLs_strict_E$chr, meth_DMLs_strict_E$start),]

# Chromosomes Size:  
chr_size=c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,156040895)

# Centromere position for each chromosomes:
centromere_pos=c(125,93.3,91,50.4,48.4,61,59.9,45.6,49,40.2,53.7,35.8,17.9,17.6,19,36.6,24,17.2,26.5,27.5,13.2,14.7,60.6)*1000000


chr_total=0
for (i in 1:23){
  chr_total=c(chr_total,sum(chr_size[1:i]))  
}

genome_size=sum(chr_size)

dev.new()
plot(1000000,1,col="white",pch=15,cex=0.4,ylim=c(-10,100), las = 1, ylab="Methylation (%)",typ="l",xlab="",xaxt = "n",xlim=c(meth_DMLs_strict_E$start[1],genome_size))


for(i in 1 :23){
  if (i>1){
    lines(c(chr_total[i],chr_total[i]),c(-10,100),col="gray48")
  }
  if (i==23){
    text(chr_total[i]+chr_size[i]/2,100,"x",cex=0.75)
  }
  else {
    text(chr_total[i]+chr_size[i]/2,100,i,cex=0.75)
  }
}

sum=0

for (i in 1:23){
  lines(c(sum+20000000,sum+chr_size[i]-20000000),c(-10,-10),lwd=10,col="black")
  lines(c(sum+20000000,sum+chr_size[i]-20000000),c(-10,-10),lwd=8,col="gray50")
  lines(c(sum+20000000,sum+chr_size[i]-20000000),c(-10,-10),lwd=7,col="gray53")
  lines(c(sum+20000000,sum+chr_size[i]-20000000),c(-10,-10),lwd=6,col="gray59")
  lines(c(sum+20000000,sum+chr_size[i]-20000000),c(-10,-10),lwd=4,col="gray75")
  lines(c(sum+20000000,sum+chr_size[i]-20000000),c(-10,-10),lwd=2,col="gray85")
  lines(c(sum+20000000,sum+chr_size[i]-20000000),c(-10,-10),lwd=1,col="gray90")
  lines(c(sum+centromere_pos[i],sum+centromere_pos[i]),c(-11,-9),col="grey13",lwd=2)
  sum=sum+chr_size[i]
}


y=meth_DMLs_strict_E$start+chr_total[meth_DMLs_strict_E$chr]

window=300 # moving average window size
for(i in c(3:10)){
  lines(frollmean(y,window),frollmean(meth_DMLs_strict_E[,i],window),col=names(ploidy)[i-2])
}
par(mai=c(21, 6, 1, 7), xpd=TRUE)
legend(x = "bottom", legend = str_replace_all(rev(ploidy),"_", " "), text.col = rev(names(ploidy)), xpd = T, ncol = 8, text.width = max(strwidth(ploidy))/10, x.intersp = c(1,1,1,1,2,2,2,2),y.intersp = 18, bty="n")

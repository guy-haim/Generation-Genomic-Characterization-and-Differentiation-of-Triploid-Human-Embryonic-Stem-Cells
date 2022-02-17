library(stringr)
library(dplyr)
library(edgeR)
library(ggplot2)
library(ggfortify)
library(pheatmap)
library(data.table)
library(zoo)
library(gplots)
library(PCAtools)
library(reshape2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(FSA)
library(ggsignif)
library(ggpubr)

#loading RNA-seq data from each sample into a data.frame and getting rid of unnecessary data (first 4 rows and second and third columns)
dsRed_4c_1 <- read.delim("dsRed-hygro4cp40_S3_counts.txt", skip = 1)[,c(1,7)]
dsRed_4c_2 <- read.delim("G1_10R_4c_p38_counts.txt", skip = 1)[,c(1,7)]
EGFP_4c_1 <- read.delim("EGFP-neo-4cp39_S10_counts.txt", skip = 1)[,c(1,7)]
EGFP_4c_2 <- read.delim("G2_10G_4c_p38_counts.txt", skip = 1)[,c(1,7)]
fusionI_3n_clone_H <- read.delim("fusionI-3n-f-4-cloneH_S8_counts.txt", skip = 1)[,c(1,7)]
fusionI_3n_clone_L_1 <- read.delim("fusion-3n-f5-cloneL_counts.txt", skip = 1)[,c(1,7)]
fusion_3n_Clone_L_2 <- read.delim("Clone-L-3n_S9_counts.txt", skip = 1)[,c(1,7)]
fusionI_3n_clone_I <- read.delim("fusionI-3n-f-4-cloneI_S6_counts.txt", skip = 1)[,c(1,7)]
fusion_3n_Clone_K <- read.delim("Clone-K-3n_S8_counts.txt", skip = 1)[,c(1,7)]
fusion_3n_Clone_J <- read.delim("Clone-J-3n_S12_counts.txt", skip = 1)[,c(1,7)]
#fusion_3n_Clone_B1 <- read.delim("fusionI-3n-f-4-cloneB_S7_counts.txt", skip = 1)[,c(1,7)]  #weired RNA-seq
fusion_3n_Clone_B <- read.delim("Clone-B-3n_S7_counts.txt", skip = 1)[,c(1,7)] #second RNA-seq of clone B
 
HuES53 <- read.delim("C:/Users/owner/Desktop/RNA-seq counts guy/Ido's samples/SRR2921778_ReadsPerGene.out.tab", skip = 1)[-c(1:4),1:2]
HuES64 <- read.delim("C:/Users/owner/Desktop/RNA-seq counts guy/Ido's samples/SRR2921779_ReadsPerGene.out.tab", skip = 1)[-c(1:4),1:2]
colnames(HuES53)[1] <- "Geneid"
colnames(HuES64)[1] <- "Geneid"
dip_pES10_rep1 <- read.delim("C:/Users/owner/Desktop/RNA-seq counts guy/Ido's samples/SRR2131926_ReadsPerGene.out.tab", skip = 1)[-c(1:4),1:2] 
dip_pES10_rep2 <- read.delim("C:/Users/owner/Desktop/RNA-seq counts guy/Ido's samples/SRR2131927_ReadsPerGene.out.tab", skip = 1)[-c(1:4),1:2] 
colnames(dip_pES10_rep1)[1] <- "Geneid"
colnames(dip_pES10_rep2)[1] <- "Geneid"

###merging all RNA-seq data into one data.frame and adding gene names according to annotation file:
gene_exp_all <- Reduce(function(x, y) merge(x, y, all=T, by = "Geneid"), 
                       list(HuES53,HuES64,dip_pES10_rep1,dip_pES10_rep2,dsRed_4c_1,dsRed_4c_2,EGFP_4c_1,EGFP_4c_2,fusionI_3n_clone_H,fusionI_3n_clone_I,fusionI_3n_clone_L_1,
                            fusion_3n_Clone_L_2,fusion_3n_Clone_K,fusion_3n_Clone_J,fusion_3n_Clone_B))
colnames(gene_exp_all) <- c("ENSEMBL","HuES53 2n","HuES64 2n","Xi 2n rep1","Xi 2n rep2","dsRed 2n rep1","dsRed 2n rep2","EGFP 2n rep1","EGFP 2n rep2","Clone H 3n","Clone I 3n","Clone L 3n rep1",
                            "Clone L 3n rep2","Clone K 3n","Clone J 3n","Clone B 3n")
gene_lengths <- read.csv("gencode.v34.transcriptLengths.csv")[,-8]

#removing the end of the ENSEMBL ID so that the same genes from different versions of the annotation file will be joined together  
gene_exp_all$ENSEMBL <- str_remove(gene_exp_all$ENSEMBL, "[.].*")
gene_lengths$ENSEMBL <- str_remove(gene_lengths$ENSEMBL, "[.].*")
gene_exp_all <- aggregate(gene_exp_all[-1], list(gene_exp_all$ENSEMBL), FUN = sum, na.rm = T)
colnames(gene_exp_all)[1] <- "ENSEMBL"

#merging all data frame to one count table with all the data needed:
gene_exp_all <- merge(gene_lengths, gene_exp_all, by = "ENSEMBL", all.y = T)
gene_exp_all <- gene_exp_all[rowSums(is.na(gene_exp_all[ ,2:22])) == 0, ] #removing empty rows
indRemoved <- which(apply(gene_exp_all[,8:22], 1, function(x) all(x == 0)) ) #saving all unexpressed genes in the tpm data.frame
gene_exp_all <- gene_exp_all[-indRemoved,]      #getting rid of rows (genes) that all of their columns are zeros
.rowNamesDF(gene_exp_all, make.names=T) <- gene_exp_all$SYMBOL

X_genes <- gene_exp_all[gene_exp_all$chr=="chrX",7]
cl=c("grey55","grey40","grey25","grey10","deepskyblue1","deepskyblue3","royalblue1","royalblue3","red1","red2","red3","red4","firebrick1","firebrick2","firebrick3")
names(cl) <- str_replace_all(colnames(gene_exp_all[,8:22]),"_", " ")


#### Differential expression (DE) analysis:

ploidy <- factor(c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2), labels = c("2n","3n")) #4 recently diploid ("normal" ESCs),and 7 triploid samples
DE_object <- DGEList(gene_exp_all[,8:22], group = ploidy, genes = gene_exp_all$gene)
keep <- filterByExpr(DE_object, min.count=2, group = ploidy)
DE_object <- DE_object[keep, , keep.lib.sizes=FALSE]
DE_TMM <- calcNormFactors(DE_object)

###PCA analysis of log2(CPM) values of each sample
cpm_normalized <- cpm(DE_TMM, log = T) #log2(cpm) after TMM normalization
cpm_normalized <- data.frame("SYMBOL" = rownames(cpm_normalized), cpm_normalized)
cpm_normalized <- merge(gene_lengths, cpm_normalized, by = "SYMBOL",all.y = T)
.rowNamesDF(cpm_normalized, make.names=T) <- cpm_normalized$SYMBOL
cpm_normalized_X <- cpm_normalized[cpm_normalized$SYMBOL %in% X_genes,]
cpm_normalized_X <- cpm_normalized_X[order(cpm_normalized_X$start, decreasing = F),]

prin_comp <- prcomp(t(cpm_normalized[,8:22]), center = T, retx = T)
gg <- cbind.data.frame(prin_comp$x, "Sample" =colnames(gene_exp_all)[8:22])
gg <- data.frame(gg, "Ploidy"=ploidy)
autoplot(prin_comp, data = gg,shape="Ploidy", colour = "Sample", label =F,size=3)+ scale_colour_manual(values=cl)+theme_bw(base_size = 18)+
  guides(shape = guide_legend(order = 1),colour = guide_legend(override.aes = list(shape=c(15))))+theme(text = element_text(size = 16))

summary(prin_comp)

pheatmap(data.matrix(cpm_normalized_X[,8:22]), scale = "row", cluster_rows=F, show_rownames = F,labels_col = str_replace_all(names(cpm_normalized_X[,8:22]),"[.]"," ")) #CPM heatmap of X-linked expressed genes, without gene clustering


################  calculate TPM and add gene names ####################

# function to calculate TPM 
tpm_calc <- function(counts, lengths) {
  rate <- counts / lengths
  rate/sum(rate)*1*10^6 
}


# calculate TPM
tpm <- gene_exp_all
for(i in names(gene_exp_all[,8:22])) {
  tpm[i] <- tpm_calc(gene_exp_all[i],gene_exp_all$length_total_exons) 
}
# write.csv(tpm,'ploidy_tpm.csv')  # write TPM table

indRemoved <- which(apply(tpm[,8:22], 1, function(x) all(x <= 1)) ) #saving all unexpressed genes in the tpm data.frame
tpm <- tpm[-indRemoved,]      #getting rid of rows (genes) that all of their columns are zeros


#############-----average_e_Karyotype-----####################################################

tpm$chr <- str_remove(tpm$chr, "chr")
tpm_E <- data.frame(tpm[,1:7], "gene_exp_2n"=(tpm[,8]+tpm[,9]+tpm[,10]+tpm[,11])/4,
                    "gene_exp_recently_2n"=(tpm[,12]+tpm[,13]+tpm[,14]+tpm[,15])/4,
                    "gene_exp_3n"=(tpm[,16]+tpm[,17]+tpm[,18]+tpm[,19]+tpm[,20]+tpm[,21]+tpm[,22])/7)

floor=1 # Expression threshold value
expressed=0.8 # Percentage of expressed samples
variable=0.1 # Percentage of most variable genes to remove
tpm_E[,-c(1:7)][tpm_E[,-c(1:7)]<1]=floor                      #?????
tpm_E=tpm_E[rowSums(tpm_E[,-c(1:7)]>floor)>=ncol(tpm_E[,-c(1:7)])* expressed,]  #?????
row_mean=rowMeans(data.matrix(tpm_E[,-c(1:7)]))
distance_mat=sweep(tpm_E[,-c(1:7)],MARGIN=1,STATS=row_mean,"-")
ssq=rowSums(distance_mat^2)
tpm_E=tpm_E[ssq<quantile(ssq,1-variable),]

med=apply(data.matrix(tpm_E[,8]),1,median) #making the median of regular diploid ESCs the baseline for expression
minus_med=sweep(tpm_E[,-c(1:7)],1,med,"/")
minus_med=log(minus_med,2)
minus_med=sweep(minus_med,2,colMeans(minus_med),"-")
CGH=cbind(tpm_E[,c(1:7)],minus_med)
CGH$chr[CGH$chr=="X"] <- 23
CGH$chr <- as.numeric(CGH$chr)
CGH <- na.omit(CGH)
CGH=CGH[order(CGH$chr, CGH$start),]


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
plot(1000000,1,col="white",pch=15,cex=0.4,ylim=c(-1,1),las=1, ylab=expression(paste("Relative  ",log[2](Expression))),typ="l",xlab="",xaxt = "n",xlim=c(CGH$start[1],genome_size))


for(i in 1 :23){
  if (i>1){
    lines(c(chr_total[i],chr_total[i]),c(-1,1),col="gray48")
  }
  if (i==23){
    text(chr_total[i]+chr_size[i]/2,1,"x",cex=0.75)
  }
  else {
    text(chr_total[i]+chr_size[i]/2,1,i,cex=0.75)
  }
}

sum=0

for (i in 1:23){
  lines(c(sum+20000000,sum+chr_size[i]-20000000),c(-1,-1),lwd=10,col="black")
  lines(c(sum+20000000,sum+chr_size[i]-20000000),c(-1,-1),lwd=8,col="gray50")
  lines(c(sum+20000000,sum+chr_size[i]-20000000),c(-1,-1),lwd=7,col="gray53")
  lines(c(sum+20000000,sum+chr_size[i]-20000000),c(-1,-1),lwd=6,col="gray59")
  lines(c(sum+20000000,sum+chr_size[i]-20000000),c(-1,-1),lwd=4,col="gray75")
  lines(c(sum+20000000,sum+chr_size[i]-20000000),c(-1,-1),lwd=2,col="gray85")
  lines(c(sum+20000000,sum+chr_size[i]-20000000),c(-1,-1),lwd=1,col="gray90")
  lines(c(sum+centromere_pos[i],sum+centromere_pos[i]),c(-1.02,-0.96),col="grey13",lwd=2)
  sum=sum+chr_size[i]
}


y=CGH$start+chr_total[CGH$chr]
j=c(0,0,0,0,0,0,0,4,1,-1)
window=300 # moving average window size
for(i in 8:10){
  lines(rollmean(y,window),rollmean(CGH[,i],window),col=cl[i-j[i]])
}
par(mai=c(21, 6, 0, 6), xpd=TRUE)
legend(x = "bottom",xjust = 0.75, legend = c("Diploid pES10","Diploid Xi pES10","Triploid pES10"), text.col = c("royalblue3","grey10","red3"), xpd = T, ncol = 3, text.width = max(strwidth(names(cl)))/130, x.intersp = c(3,6.5,5),  y.intersp = 10, bty="n")

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

setwd("C:/Users/owner/Desktop/RNA-seq counts guy/Ido's samples")
#loading RNA-seq data from each sample into a data.frame and getting rid of unnecessary data (first 4 rows and second and third columns)
pES10_h_G1_rep1 <- read.delim("SRR2131924_ReadsPerGene.out.tab", header = F, skip = 4)[,1:2]
pES10_h_G1_rep2 <- read.delim("SRR2131925_ReadsPerGene.out.tab", header = F, skip = 4)[,1:2]
# pES10_d_G1_rep1 <- read.delim("SRR2131926_ReadsPerGene.out.tab", header = F, skip = 4)[,1:2]
# pES10_d_G1_rep2 <- read.delim("SRR2131927_ReadsPerGene.out.tab", header = F, skip = 4)[,1:2]
colnames(pES10_h_G1_rep1)[1] <- "Geneid"
colnames(pES10_h_G1_rep2)[1] <- "Geneid"
# colnames(pES10_d_G1_rep1)[1] <- "Geneid"
# colnames(pES10_d_G1_rep2)[1] <- "Geneid"

setwd("C:/Users/owner/Desktop/RNA-seq counts guy")
dsRed_1c_rep1 <- read.delim("dsRed-hygro1cp40_S9_counts.txt", skip = 1)[,c(1,7)]
dsRed_4c_rep1 <- read.delim("dsRed-hygro4cp40_S3_counts.txt", skip = 1)[,c(1,7)]
dsRed_4c_rep2 <- read.delim("G1_10R_4c_p38_counts.txt", skip = 1)[,c(1,7)]
EGFP_4c_rep1 <- read.delim("EGFP-neo-4cp39_S10_counts.txt", skip = 1)[,c(1,7)]
EGFP_4c_rep2 <- read.delim("G2_10G_4c_p38_counts.txt", skip = 1)[,c(1,7)]

###merging all RNA-seq data into one data.frame and adding gene names according to annotation file:
gene_exp_all <- Reduce(function(x, y) merge(x, y, all=T, by = "Geneid"), 
                       list(pES10_h_G1_rep1, pES10_h_G1_rep2, dsRed_1c_rep1,
                            dsRed_4c_rep1, dsRed_4c_rep2, EGFP_4c_rep1, EGFP_4c_rep2))
colnames(gene_exp_all) <- c("ENSEMBL","pES10_h_G1_rep1","pES10_h_G1_rep2","dsRed_1c_rep1",
                            "dsRed_4c_rep1","dsRed_4c_rep2","EGFP_4c_rep1","EGFP_4c_rep2")

gene_lengths <- read.csv("gencode.v34.transcriptLengths.csv")[,-8]

#removing the end of the ENSEMBL ID so that the same genes from different versions of the annotation file will be joined together  
gene_exp_all$ENSEMBL <- str_remove(gene_exp_all$ENSEMBL, "[.].*")
gene_lengths$ENSEMBL <- str_remove(gene_lengths$ENSEMBL, "[.].*")
gene_exp_all <- aggregate(gene_exp_all[-1], list(gene_exp_all$ENSEMBL), FUN = sum, na.rm = T)
colnames(gene_exp_all)[1] <- "ENSEMBL"

#merging all data frame to one count table with all the data needed:
gene_exp_all <- merge(gene_lengths, gene_exp_all, by = "ENSEMBL", all.y = T)
# hs <- org.Hs.eg.db
# enterezID <- select(hs,keys = gene_lengths$SYMBOL,columns = c("ENTREZID"), keytype = "SYMBOL")
# gene_exp_all <- merge(gene_exp_all, enterezID, by = "SYMBOL", all.x = T)
gene_exp_all <- gene_exp_all[rowSums(is.na(gene_exp_all[ ,2:14])) == 0, ] #removing empty rows
indRemoved <- which(apply(gene_exp_all[,8:14], 1, function(x) all(x == 0)) ) #saving all unexpressed genes in the tpm data.frame
gene_exp_all <- gene_exp_all[-indRemoved,]      #getting rid of rows (genes) that all of their columns are zeros
.rowNamesDF(gene_exp_all, make.names=T) <- gene_exp_all$SYMBOL

X_genes <- gene_exp_all[gene_exp_all$chr=="chrX",7]
pluri_genes <- c("POU5F1","DNMT3B","LIN28A","DPPA4","TDGF1","SALL2","SALL4","SOX2","ZFP42","NANOG")
cl=c("green2","green3","green4","royalblue1","royalblue2","royalblue3","royalblue4")
names(cl) <- str_replace_all(colnames(gene_exp_all[,8:14]),"_", " ")



#### Differential expression (DE) analysis:

ploidy <- factor(c(2,2,2,1,1,1,1), labels = c("2n","1n")) 
DE_object <- DGEList(gene_exp_all[,8:14], group = ploidy, genes = gene_exp_all$gene)
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

prin_comp <- prcomp(t(cpm_normalized[,8:14]), center = T, retx = T)
gg <- cbind.data.frame(prin_comp$x, "Sample" =names(cl))
gg <- data.frame(gg, "Ploidy"=ploidy)
autoplot(prin_comp, data = gg, shape="Ploidy", colour = "Sample", label =F,size=3)+ scale_colour_manual(values=cl)+theme_bw()+
  guides(shape = guide_legend(order = 1),colour = guide_legend(override.aes = list(shape=c(15))))+theme(text = element_text(size = 16))

summary(prin_comp)

#continue DE analysis
batch <- factor(c(1,1,2,2,2,2,2), labels = c("1","2"))
design <- model.matrix(~batch+ploidy)
DE_TMM <- estimateDisp(DE_TMM, design, robust=TRUE)
fit_TMM <- glmFit(DE_TMM, design)   #fit the model
# fit_TMM <- glmQLFit(DE_TMM, design)   #fit the model (suited for low number of replicates)
# LRT_1n_Vs_2n <- glmLRT(fit_TMM)  #strict comparison of haploids and diploids
LRT_1n_Vs_2n <- glmTreat(fit_TMM, lfc=log2(2))  #strict comparison of haploids and diploids

summary(decideTests(LRT_1n_Vs_2n))

TMM_1n <- topTags(LRT_1n_Vs_2n, n=nrow(DE_TMM),adjust.method = "BH")$table #assigning fold change gene expression between diploid and haploid ESCs
TMM_1n <- data.frame("SYMBOL"=row.names(TMM_1n), TMM_1n)

colnames(TMM_1n)[2] <- "logFC_haploids"

TMM_1n <- merge(gene_lengths, TMM_1n, by = "SYMBOL", all.y = T)

sig_DE_1n <-TMM_1n[TMM_1n$FDR<=0.05,] #assigning only genes with significant fold change between haploid and diploid ESCs (by FDR) 
sig_DE_1n <- sig_DE_1n[-rowSums(is.na(sig_DE_1n[ ,2:7])) == 0, ] #removing empty rows
sig_DE_1n<- sig_DE_1n[order(sig_DE_1n$logFC_haploids, decreasing = T),]
write.csv(sig_DE_1n, file = "significant DE genes between diploid and haploid hESCs")


########### Differences in expression levels of different genes groups ########### 


## pluripotency markers ##
group <- c("1n","1n","1n","2n","2n","2n","2n")
cpm_normalized_pluri <- 2^(t(cpm_normalized[cpm_normalized$SYMBOL %in% pluri_genes,-c(1:7)]))
cpm_normalized_pluri <- data.frame(cpm_normalized_pluri, "Ploidy" = group)
melted_pluri <- melt(data = cpm_normalized_pluri, id.vars = "Ploidy", variable.name = "Gene",value.name = "Expression")
melted_pluri["log2_Expression"] <- log2(melted_pluri$Expression)

stat=compare_means(log2_Expression~Ploidy, data = melted_pluri, method="wilcox.test", p.adjust.method = "BH",group.by = "Gene")
ggbarplot(melted_pluri, x= "Gene", y= "log2_Expression",fill="Ploidy",palette = c("forestgreen","royalblue"), position = position_dodge(0.8), add="mean_se")+
  ylab(expression(log[2](CPM)))+stat_pvalue_manual(data = stat, label = "p.signif", y.position = 12, size = 3,x = "Gene",hide.ns = F)


## ion channels genes ##
ion_channels_genes <- read.delim("ion_channels_genes.txt")
ion_channels_genes <- ion_channels_genes[,c(2,13)]
colnames(ion_channels_genes) <- c("SYMBOL","Group")
ploidy <- c("1n","1n","1n","2n","2n","2n","2n")
cpm_ion_channels <- merge(cpm_normalized[,c(1,8:14)],ion_channels_genes,all.y=T, by="SYMBOL")
cpm_ion_channels <- na.omit(cpm_ion_channels)
cpm_ion_channels[,2:8] <- 2^cpm_ion_channels[,2:8]
diploid_mean <- apply(cpm_ion_channels[,5:8],1,FUN = mean)
cpm_ion_channels[,2:8] <- cpm_ion_channels[,2:8]/diploid_mean
melted_ion<-melt(cpm_ion_channels)
melted_ion$variable <- rep(ploidy, each=188)
colnames(melted_ion)[3:4] <- c("Ploidy","Expression")

stat=compare_means(Expression~Ploidy, data = melted_ion, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_ion, x= "Group", y= "Expression",fill="Ploidy",palette = c("forestgreen","royalblue"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.format", y.position = 2, size = 3,x = "Group",hide.ns = F)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

### Difference in expression (%) per gene ###
melted_ion_genes <- melted_ion
stat=compare_means(Expression~Ploidy, data = melted_ion_genes, method="wilcox.test", p.adjust.method = "BH",group.by = "SYMBOL")
# melted_ion_genes <- melted_ion_genes[stat$p.format<=0.05,]
ggbarplot(melted_ion_genes, x= "SYMBOL", y= "Expression",fill="Ploidy",palette = c("forestgreen","royalblue"),width = 0.5, position = position_dodge(0.6), add="mean_se")+
  stat_pvalue_manual(data = stat[stat$p.format<=0.05,], size=3, label = "p.format", y.position = 10, x = "SYMBOL",hide.ns = F)+xlab("")+geom_hline(yintercept = 0, linetype=2)+
  scale_y_continuous(breaks=seq(-30,10,5))+ylab("Difference in expression (%)")+theme(axis.text.x = element_text(angle = 90), text = element_text(size = 10))

# Plotting ion channels by more general subgroups (according to HGNC):
melted_ion_general <- melted_ion
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Anoctamins","Chloride channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Bestrophins","Chloride channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Chloride voltage-gated channels","Chloride channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Chloride intracellular channels","Chloride channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Chloride channels, ATP-gated CFTR","Chloride channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Volume regulated anion channel subunits","Chloride channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Sodium voltage-gated channel alpha subunits","Sodium channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Sodium voltage-gated channel beta subunits","Sodium channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Sodium channels epithelial","Sodium channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Sodium leak channels, non selective","Sodium channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Acid sensing ion channel subunits","Sodium channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Cation channels sperm associated ","Calcium channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Two pore segment channels","Calcium channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Calcium voltage-gated channel alpha1 subunits","Calcium channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Calcium voltage-gated channel auxiliary alpha2delta subunits","Calcium channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Calcium voltage-gated channel auxiliary beta subunits","Calcium channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Calcium channel auxiliary gamma subunits","Calcium channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Ryanodine receptors","Calcium channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Inositol 1,4,5-triphosphate receptors","Calcium channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Potassium voltage-gated channels","Potassium channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Potassium inwardly rectifying channel subfamily J","Potassium channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Potassium two pore domain channel subfamily K ","Potassium channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Potassium calcium-activated channels","Potassium channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Potassium sodium-activated channel subfamily T ","Potassium channels")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Aquaporins","Porins")
melted_ion_general$Group <- str_replace_all(melted_ion_general$Group,"Voltage dependent anion channels","Porins")
melted_ion_general$Expression <- (melted_ion_general$Expression-1)*100

# ### reorder dataframe by delta expression:
# melted_ion_general <- melted_ion_general[order(melted_ion_general$Group),]
# melted_ion_general$Order <- c(rep(1,407),rep(3,363),rep(5,154),rep(4,88),rep(2,682),rep(6,198))
# melted_ion_general <- melted_ion_general[order(melted_ion_general$Order),]

stat=compare_means(Expression~Ploidy, data = melted_ion_general, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_ion_general[melted_ion_general$Group!="Porins" & melted_ion_general$Group!="Gap junction proteins",], x= "Group", 
          y="Expression",fill="Ploidy",palette = c("forestgreen","royalblue"),width = 0.5, position = position_dodge(0.6),add="mean_se")+
  stat_pvalue_manual(data = stat[stat$Group!="Porins" & stat$Group!="Gap junction proteins",], size=4, label = "p.format", 
                     y.position = c(50,82,32,83), x = "Group", hide.ns = F)+xlab("")+geom_hline(yintercept = 0, linetype=2)+
  scale_y_continuous(breaks=seq(0,85,10))+
  ylab("Difference in expression (%)")+
  theme(text = element_text(size = 17),axis.text.y = element_text(size = 12))+
  scale_x_discrete(labels = c("CL-","NA+","Ca2+","K+"))

## GPCRs genes ##
GPCR_genes <- read.delim("GPCRs_genes.txt")
GPCR_genes <- GPCR_genes[,c(2,13)]
colnames(GPCR_genes) <- c("SYMBOL","Group")
ploidy <- c("1n","1n","1n","2n","2n","2n","2n")
cpm_GPCRs <- merge(cpm_normalized[,c(1,8:14)],GPCR_genes,all.y=T, by="SYMBOL")
cpm_GPCRs <- na.omit(cpm_GPCRs)
cpm_GPCRs[,2:8] <- 2^cpm_GPCRs[,2:8]
diploid_mean <- apply(cpm_GPCRs[,5:8],1,FUN = mean)
cpm_GPCRs[,2:8] <- cpm_GPCRs[,2:8]/diploid_mean
melted_GPCR<-melt(cpm_GPCRs)
melted_GPCR$variable <- rep(ploidy, each=287)
colnames(melted_GPCR)[3:4] <- c("Ploidy","Expression")

stat=compare_means(Expression~Ploidy, data = melted_GPCR, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_GPCR, x= "Group", y= "Expression",fill="Ploidy",palette = c("forestgreen","royalblue"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif", y.position = 3.5,x = "Group",hide.ns = F)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## receptors kinases genes ##
RKs_genes <- read.delim("receptors_kinases_genes.txt")
RKs_genes <- RKs_genes[,c(2,13)]
colnames(RKs_genes) <- c("SYMBOL","Group")
ploidy <- c("1n","1n","1n","2n","2n","2n","2n")
cpm_RKs <- merge(cpm_normalized[,c(1,8:14)],RKs_genes,all.y=T, by="SYMBOL")
cpm_RKs <- na.omit(cpm_RKs)
cpm_RKs[,2:8] <- 2^cpm_RKs[,2:8]
diploid_mean <- apply(cpm_RKs[,5:8],1,FUN = mean)
cpm_RKs[,2:8] <- cpm_RKs[,2:8]/diploid_mean
melted_RK<-melt(cpm_RKs)
melted_RK$variable <- rep(ploidy, each=72)
colnames(melted_RK)[3:4] <- c("Ploidy","Expression")

stat=compare_means(Expression~Ploidy, data = melted_RK, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_RK, x= "Group", y= "Expression",fill="Ploidy",palette = c("forestgreen","royalblue"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif",y.position = 7.5,x = "Group",hide.ns = F)+theme(axis.text.x = element_text(angle = 0))


## Solute carriers genes ##
SLCs_genes <- read.delim("solute_carriers_genes.txt")
SLCs_genes <- SLCs_genes[,c(2,13)]
colnames(SLCs_genes) <- c("SYMBOL","Group")
ploidy <- c("1n","1n","1n","2n","2n","2n","2n")
cpm_SLCs <- merge(cpm_normalized[,c(1,8:14)],SLCs_genes,all.y=T, by="SYMBOL")
cpm_SLCs <- na.omit(cpm_SLCs)
cpm_SLCs[,2:8] <- 2^cpm_SLCs[,2:8]
diploid_mean <- apply(cpm_SLCs[,5:8],1,FUN = mean)
cpm_SLCs[,2:8] <- cpm_SLCs[,2:8]/diploid_mean
melted_SLC<-melt(cpm_SLCs)
melted_SLC$variable <- rep(ploidy, each=364)
colnames(melted_SLC)[3:4] <- c("Ploidy","Expression")

stat=compare_means(Expression~Ploidy, data = melted_SLC, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_SLC, x= "Group", y= "Expression",fill="Ploidy",palette = c("forestgreen","royalblue"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif",y.position = 1.5,x = "Group",hide.ns = F)+theme(axis.text.x = element_text(angle = 0))


## integrins genes ##
INTs_genes <- read.delim("integrin_genes.txt")
INTs_genes <- INTs_genes[,c(2,13)]
colnames(INTs_genes) <- c("SYMBOL","Group")
ploidy <- c("1n","1n","1n","2n","2n","2n","2n")
cpm_INTs <- merge(cpm_normalized[,c(1,8:14)],INTs_genes,all.y=T, by="SYMBOL")
cpm_INTs <- na.omit(cpm_INTs)
cpm_INTs[,2:8] <- 2^cpm_INTs[,2:8]
diploid_mean <- apply(cpm_INTs[,5:8],1,FUN = mean)
cpm_INTs[,2:8] <- cpm_INTs[,2:8]/diploid_mean
melted_INT<-melt(cpm_INTs)
melted_INT$variable <- rep(ploidy, each=24)
colnames(melted_INT)[3:4] <- c("Ploidy","Expression")

stat=compare_means(Expression~Ploidy, data = melted_INT, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_INT, x= "Group", y= "Expression",fill="Ploidy",palette = c("forestgreen","royalblue"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif",y.position = 2.5,x = "Group",hide.ns = F)+theme(axis.text.x = element_text(angle = 0))


## Basic leucin zipper genes ##
BLZ_genes <- read.delim("Basic_leucin_zipper_genes.txt")
BLZ_genes <- BLZ_genes[,c(2,13)]
colnames(BLZ_genes) <- c("SYMBOL","Group")
ploidy <- c("1n","1n","1n","2n","2n","2n","2n")
cpm_BLZ <- merge(cpm_normalized[,c(1,8:14)],BLZ_genes,all.y=T, by="SYMBOL")
cpm_BLZ <- na.omit(cpm_BLZ)
cpm_BLZ[,2:8] <- 2^cpm_BLZ[,2:8]
diploid_mean <- apply(cpm_BLZ[,5:8],1,FUN = mean)
cpm_BLZ[,2:8] <- cpm_BLZ[,2:8]/diploid_mean
melted_BLZ<-melt(cpm_BLZ)
melted_BLZ$variable <- rep(ploidy, each=50)
colnames(melted_BLZ)[3:4] <- c("Ploidy","Expression")

stat=compare_means(Expression~Ploidy, data = melted_BLZ, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_BLZ, x= "Group", y= "Expression",fill="Ploidy",palette = c("forestgreen","royalblue"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif",y.position = 1.5,x = "Group",hide.ns = F)+theme(axis.text.x = element_text(angle = 0))


## Basic helix-loop-helix genes ##
BHLH_genes <- read.delim("BHLH_genes.txt")
BHLH_genes <- BHLH_genes[,c(2,13)]
colnames(BHLH_genes) <- c("SYMBOL","Group")
ploidy <- c("1n","1n","1n","2n","2n","2n","2n")
cpm_BHLH <- merge(cpm_normalized[,c(1,8:14)],BHLH_genes,all.y=T, by="SYMBOL")
cpm_BHLH <- na.omit(cpm_BHLH)
cpm_BHLH[,2:8] <- 2^cpm_BHLH[,2:8]
diploid_mean <- apply(cpm_BHLH[,5:8],1,FUN = mean)
cpm_BHLH[,2:8] <- cpm_BHLH[,2:8]/diploid_mean
melted_BHLH<-melt(cpm_BHLH)
melted_BHLH$variable <- rep(ploidy, each=89)
colnames(melted_BHLH)[3:4] <- c("Ploidy","Expression")

stat=compare_means(Expression~Ploidy, data = melted_BHLH, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_BHLH, x= "Group", y= "Expression",fill="Ploidy",palette = c("forestgreen","royalblue"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif",y.position = 5.5,x = "Group",hide.ns = F)+theme(axis.text.x = element_text(angle = 0))


## zync fingers genes ##
ZF_genes <- read.delim("zync_fingers_genes.txt")
ZF_genes <- ZF_genes[,c(2,13)]
colnames(ZF_genes) <- c("SYMBOL","Group")
ploidy <- c("1n","1n","1n","2n","2n","2n","2n")
cpm_ZF <- merge(cpm_normalized[,c(1,8:14)],ZF_genes,all.y=T, by="SYMBOL")
cpm_ZF <- na.omit(cpm_ZF)
cpm_ZF[,2:8] <- 2^cpm_ZF[,2:8]
diploid_mean <- apply(cpm_ZF[,5:8],1,FUN = mean)
cpm_ZF[,2:8] <- cpm_ZF[,2:8]/diploid_mean
melted_ZF<-melt(cpm_ZF)
melted_ZF$variable <- rep(ploidy, each=1486)
colnames(melted_ZF)[3:4] <- c("Ploidy","Expression")

stat=compare_means(Expression~Ploidy, data = melted_ZF, method="wilcox.test", methods.args=list(exact=T), p.adjust.method = "fdr",group.by = "Group")
ggbarplot(melted_ZF, x= "Group", y= "Expression",fill="Ploidy",palette = c("forestgreen","royalblue"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif",y.position = 3,x = "Group",hide.ns = T)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## general transcription factors genes ##
GTF_genes <- read.delim("general_transcription_factors_genes.txt")
GTF_genes <- GTF_genes[,c(2,13)]
colnames(GTF_genes) <- c("SYMBOL","Group")
ploidy <- c("1n","1n","1n","2n","2n","2n","2n")
cpm_GTF <- merge(cpm_normalized[,c(1,8:14)],GTF_genes,all.y=T, by="SYMBOL")
cpm_GTF <- na.omit(cpm_GTF)
cpm_GTF[,2:8] <- 2^cpm_GTF[,2:8]
diploid_mean <- apply(cpm_GTF[,5:8],1,FUN = mean)
cpm_GTF[,2:8] <- cpm_GTF[,2:8]/diploid_mean
melted_GTF<-melt(cpm_GTF)
melted_GTF$variable <- rep(ploidy, each=45)
colnames(melted_GTF)[3:4] <- c("Ploidy","Expression")

stat=compare_means(Expression~Ploidy, data = melted_GTF, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_GTF, x= "Group", y= "Expression",fill="Ploidy",palette = c("forestgreen","royalblue"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif",y.position = 1.5,x = "Group",hide.ns = F)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## TFs and DNA-binding genes ##
ploidy <- c("1n","1n","1n","2n","2n","2n","2n")
melted_TF_DB<-data.frame("SYMBOL"=c(melted_ZF$SYMBOL,melted_BLZ$SYMBOL,melted_BHLH$SYMBOL,melted_GTF$SYMBOL),
                         "Group"=c(melted_ZF$Group,melted_BLZ$Group,melted_BHLH$Group,melted_GTF$Group),
                         "Ploidy"=c(melted_ZF$Ploidy,melted_BLZ$Ploidy,melted_BHLH$Ploidy,melted_GTF$Ploidy),
                         "Expression"=c(melted_ZF$Expression,melted_BLZ$Expression,melted_BHLH$Expression,melted_GTF$Expression))

stat=compare_means(Expression~Ploidy, data = melted_TF_DB, method="wilcox.test",methods.args=list(exact=T), p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_TF_DB, x= "Group", y= "Expression",fill="Ploidy",palette = c("forestgreen","royalblue"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif",y.position = 20,x = "Group",hide.ns = F)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ scale_y_continuous(breaks=seq(0,25,5))


## general groups of genes ##
ploidy <- c("1n","1n","1n","2n","2n","2n","2n")
melted_genes<-data.frame("SYMBOL"=c(melted_ion$SYMBOL,melted_RK$SYMBOL,melted_GPCR$SYMBOL, melted_SLC$SYMBOL,melted_INT$SYMBOL, melted_ZF$SYMBOL,melted_BLZ$SYMBOL,melted_BHLH$SYMBOL,melted_GTF$SYMBOL),
                         "Group"=c(rep("Ion channels",1316),rep("RKs",504),rep("GPCRs",2009),rep("SLCs",2548),rep("Integrins",168),rep("DBPs & TFs",10402),rep("DBPs & TFs",350),rep("DBPs & TFs",623),rep("DBPs & TFs",315)),
                         "Ploidy"=c(melted_ion$Ploidy, melted_RK$Ploidy,melted_GPCR$Ploidy,melted_SLC$Ploidy, melted_INT$Ploidy, melted_ZF$Ploidy,melted_BLZ$Ploidy,melted_BHLH$Ploidy,melted_GTF$Ploidy),
                         "Expression"=c(melted_ion$Expression,melted_RK$Expression,melted_GPCR$Expression,melted_SLC$Expression,melted_INT$Expression, melted_ZF$Expression,melted_BLZ$Expression,melted_BHLH$Expression,melted_GTF$Expression))

stat=compare_means(Expression~Ploidy, data = melted_genes, method="wilcox.test",methods.args=list(exact=T), p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_genes, x= "Group", y= "Expression",fill="Ploidy",palette = c("forestgreen","royalblue"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif",y.position = 2.25,x = "Group",hide.ns = T)+ scale_y_continuous(breaks=seq(0,3,0.25))

## making a delta change plot:
melted_genes_delta <- melted_genes
melted_genes_delta$Expression <- (melted_genes_delta$Expression-1)*100

stat=compare_means(Expression~Ploidy, data = melted_genes_delta, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_genes_delta, x= "Group", y= "Expression",fill="Ploidy",palette = c("forestgreen","royalblue"), position = position_dodge(0.8), add="mean_se")+
  scale_y_continuous(breaks=seq(-10,120,25))+xlab(label = "")+
  ylab(label = "Difference in expression (%)")+
  stat_pvalue_manual(data = stat, label = "p.adj",y.position = 120,x = "Group",hide.ns = F)+
  geom_hline(yintercept = 0,linetype=2)+theme(text = element_text(size = 16),axis.text.y = element_text(size = 12))+
  scale_x_discrete(labels = c("Ion\nchannels","RKs","GPCRs","SLCs",
                              "Integrins","DBPs\nand TFs"))

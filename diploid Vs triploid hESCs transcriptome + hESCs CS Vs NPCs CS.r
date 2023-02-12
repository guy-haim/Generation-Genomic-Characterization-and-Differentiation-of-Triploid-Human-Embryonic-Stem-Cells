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

setwd("C:/Users/owner/Desktop/RNA-seq counts guy")

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
fusion_3n_Clone_B <- read.delim("Clone-B-3n_S7_counts.txt", skip = 1)[,c(1,7)] 

###merging all RNA-seq data into one data.frame and adding gene names according to annotation file:
gene_exp_all <- Reduce(function(x, y) merge(x, y, all=T, by = "Geneid"), 
                       list(dsRed_4c_1,dsRed_4c_2,EGFP_4c_1,EGFP_4c_2,fusionI_3n_clone_H,fusionI_3n_clone_I,fusionI_3n_clone_L_1,
                            fusion_3n_Clone_L_2,fusion_3n_Clone_K,fusion_3n_Clone_J,fusion_3n_Clone_B))
colnames(gene_exp_all) <- c("ENSEMBL","dsRed 2n rep1","dsRed 2n rep2","EGFP 2n rep1","EGFP 2n rep2","3n clone H","3n clone I","3n clone L rep1",
                            "3n clone L rep2","3n clone K","3n clone J","3n clone B")
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
gene_exp_all <- gene_exp_all[rowSums(is.na(gene_exp_all[ ,2:18])) == 0, ] #removing empty rows
indRemoved <- which(apply(gene_exp_all[,8:18], 1, function(x) all(x == 0)) ) #saving all unexpressed genes in the tpm data.frame
gene_exp_all <- gene_exp_all[-indRemoved,]      #getting rid of rows (genes) that all of their columns are zeros
.rowNamesDF(gene_exp_all, make.names=T) <- gene_exp_all$SYMBOL

X_genes <- gene_exp_all[gene_exp_all$chr=="chrX",7]
pluri_genes <- c("POU5F1","DNMT3B","LIN28A","DPPA4","TDGF1","SALL2","SALL4","SOX2","ZFP42","NANOG")
cl=c("deepskyblue1","deepskyblue3","royalblue1","royalblue3","red1","red2","red3","red4","firebrick1","firebrick2","firebrick3")
names(cl) <- str_replace_all(colnames(gene_exp_all[,8:18]),"_", " ")



#### Differential expression (DE) analysis:

ploidy <- factor(c(1,1,1,1,2,2,2,2,2,2,2), labels = c("2n","3n")) #4 recently diploid ("normal" ESCs),and 7 triploid samples
DE_object <- DGEList(gene_exp_all[,8:18], group = ploidy, genes = gene_exp_all$gene)
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

prin_comp <- prcomp(t(cpm_normalized[,8:18]), center = T, retx = T)
gg <- cbind.data.frame(prin_comp$x, "Sample" =colnames(gene_exp_all)[8:18])
gg <- data.frame(gg, "Ploidy"=ploidy)
autoplot(prin_comp, data = gg,shape="Ploidy", colour = "Sample", label =F,size=3)+ scale_colour_manual(values=cl)+theme_bw()+
  guides(shape = guide_legend(order = 1),colour = guide_legend(override.aes = list(shape=c(15))))+theme(text = element_text(size = 16))

summary(prin_comp)

#continue DE analysis
batch <- factor(c(2,4,2,4,2,2,1,2,2,2,3), labels = c("batch1","batch2","batch3","batch4")) #consider the batch effect when analyzing the data
design <- model.matrix(~batch + ploidy)
DE_TMM <- estimateDisp(DE_TMM, design, robust=TRUE)
fit_TMM <- glmFit(DE_TMM, design)   #fit the model
LRT_3n_Vs_2n <- glmTreat(fit_TMM, lfc=log2(2))  #strict comparison of triploids and diploids

summary(decideTests(LRT_3n_Vs_2n))

TMM_3n <- topTags(LRT_3n_Vs_2n, n=nrow(DE_TMM),adjust.method = "BH")$table #assigning fold change gene expression between diploid and triploid ESCs
TMM_3n <- data.frame("SYMBOL"=row.names(TMM_3n), TMM_3n)

colnames(TMM_3n)[2] <- "logFC_triploids"

TMM_3n <- merge(gene_lengths, TMM_3n, by = "SYMBOL", all.y = T)

sig_DE_3n <-TMM_3n[TMM_3n$FDR<=0.05,] #assigning only genes with significant fold change between triploid and diploid ESCs (by FDR) 
sig_DE_3n <- sig_DE_3n[-rowSums(is.na(sig_DE_3n[ ,2:7])) == 0, ] #removing empty rows
sig_DE_3n<- sig_DE_3n[order(sig_DE_3n$logFC_triploids, decreasing = T),]
write.csv(sig_DE_3n, file = "significant DE genes between diploid and triploid hESCs")


############# Gene Set Enrichment Analysis ############# 
library(fgsea)
library(dplyr)

ranks <- TMM_3n
ranks[, 'score'] <-ranks$logFC_triploids* (-log(ranks$FDR))
.rowNamesDF(ranks, make.names=T) <- ranks$SYMBOL
ranks <- setNames(ranks$score,ranks$SYMBOL)

go_Terms <- gmtPathways("c5.go.v7.4.symbols.gmt") 
msigdb_pathways <- go_Terms
fgseaRes <- fgseaMultilevel(msigdb_pathways, ranks, minSize=15,maxSize = 500,eps = 0)
fgseaRes <- fgseaRes[fgseaRes$padj<=0.05,]

collapsedPathways <- collapsePathways(fgseaRes[order(padj)], msigdb_pathways, ranks, nperm = 100000)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][order(-NES), pathway] 

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

fgseaResTidy <-fgseaResTidy[fgseaResTidy$padj<0.05,]
fgseaResTidy <- fgseaResTidy[fgseaResTidy$pathway %in% mainPathways,]

#Remove category from GO terms names:
fgseaResTidy$pathway <- str_remove(fgseaResTidy$pathway, "GOBP_")
fgseaResTidy$pathway <- str_remove(fgseaResTidy$pathway, "GOCC_")
fgseaResTidy$pathway <- str_remove(fgseaResTidy$pathway, "GOMF_")
fgseaResTidy$pathway <- str_replace_all(fgseaResTidy$pathway, "_"," ")

#Plotting enrichment of each GO term:
for (i in mainPathways){
  pdf(paste0(i,".pdf"),height=5,width=5)
  print(plotEnrichment(msigdb_pathways[[i]],ranks,ticksSize = 0.1) + labs(title=i))
  dev.off()
}

#Plotting GSEA of downregulated GO terms in triploids:
ggplot(data = fgseaResTidy[fgseaResTidy$NES<0 & fgseaResTidy$padj < 0.01,], aes(reorder(pathway, -log2(padj)), -log2(padj))) +
  geom_col(aes(fill=0),width = 0.6) +
  coord_flip() +
  labs(x="Pathway", y=expression(-log[2](FDR))) + 
  theme_minimal()+ 
  geom_hline(yintercept = -log2(0.05), linetype=2,color = "black", size=0.25)+
  theme(legend.position = "")


########### Differences in expression levels of different genes groups ########### 


## pluripotency markers ##
group <- c("2n","2n","2n","2n","3n","3n","3n","3n","3n","3n","3n")
cpm_normalized_pluri <- 2^(t(cpm_normalized[cpm_normalized$SYMBOL %in% pluri_genes,-c(1:7)]))
cpm_normalized_pluri <- data.frame(cpm_normalized_pluri, "Ploidy" = group)
melted_pluri <- melt(data = cpm_normalized_pluri, id.vars = "Ploidy", variable.name = "Gene",value.name = "Expression")
melted_pluri["log2_Expression"] <- log2(melted_pluri$Expression)

stat=compare_means(log2_Expression~Ploidy, data = melted_pluri, method="wilcox.test", p.adjust.method = "BH",group.by = "Gene")
ggbarplot(melted_pluri, x= "Gene", y= "log2_Expression",fill="Ploidy",palette = c("royalblue","red3"), position = position_dodge(0.8), add="mean_se")+
  ylab(expression(log[2](CPM)))


## ion channels genes ##
ion_channels_genes <- read.delim("ion_channels_genes.txt")
ion_channels_genes <- ion_channels_genes[,c(2,13)]
colnames(ion_channels_genes) <- c("SYMBOL","Group")
ploidy <- c("2n","2n","2n","2n","3n","3n","3n","3n","3n","3n","3n")
cpm_ion_channels <- merge(cpm_normalized[,c(1,8:18)],ion_channels_genes,all.y=T, by="SYMBOL")
cpm_ion_channels <- na.omit(cpm_ion_channels)
cpm_ion_channels[,2:12] <- 2^cpm_ion_channels[,2:12]
diploid_mean <- apply(cpm_ion_channels[,2:5],1,FUN = mean)
cpm_ion_channels[,2:12] <- cpm_ion_channels[,2:12]/diploid_mean
melted_ion<-melt(cpm_ion_channels)
melted_ion$variable <- rep(ploidy, each=172)
colnames(melted_ion)[3:4] <- c("Ploidy","Expression")

stat=compare_means(Expression~Ploidy, data = melted_ion, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_ion, x= "Group", y= "Expression",fill="Ploidy",palette = c("royalblue","red3"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif", y.position = 2, size = 7,x = "Group",hide.ns = T)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

### Difference in expression (%) per gene ###
melted_ion_genes <- melted_ion
stat=compare_means(Expression~Ploidy, data = melted_ion_genes, method="wilcox.test", p.adjust.method = "BH",group.by = "SYMBOL")
melted_ion_genes <- melted_ion_genes[stat$p.format<=0.05,]
ggbarplot(melted_ion_genes, x= "SYMBOL", y= "Expression",fill="Ploidy",palette = c("royalblue","red3"),width = 0.5, position = position_dodge(0.6), add="mean_se")+
  stat_pvalue_manual(data = stat[stat$p.format<=0.05,], size=4.5, label = "p.format", y.position = 3.2, x = "SYMBOL",hide.ns = F)+xlab("")+geom_hline(yintercept = 0, linetype=2)+
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

### reorder dataframe by delta expression:
melted_ion_general <- melted_ion_general[order(melted_ion_general$Group),]
melted_ion_general$Order <- c(rep(1,407),rep(3,363),rep(5,154),rep(4,88),rep(2,682),rep(6,198))
melted_ion_general <- melted_ion_general[order(melted_ion_general$Order),]

stat=compare_means(Expression~Ploidy, data = melted_ion_general, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_ion_general, x= "Group", y= "Expression",fill="Ploidy",palette = c("royalblue","red3"),width = 0.5, position = position_dodge(0.6), add="mean_se")+
  stat_pvalue_manual(data = stat, size=4.5, label = "p.adj", y.position = 10, x = "Group", hide.ns = F)+xlab("")+geom_hline(yintercept = 0, linetype=2)+
  scale_y_continuous(breaks=seq(-30,10,5))+ylab("Difference in expression (%)")+theme(text = element_text(size = 17),axis.text.y = element_text(size = 12))

### Difference in Ca channels expression (%) per gene ###
melted_calcium_channel_genes <- melted_ion_general[melted_ion_general$Group=="Calcium channels",]
stat=compare_means(Expression~Ploidy, data = melted_calcium_channel_genes, method="wilcox.test", p.adjust.method = "BH",group.by = "SYMBOL")
# melted_ion_genes <- melted_ion_genes[stat$p.format<=0.05,]
ggbarplot(melted_calcium_channel_genes, x= "SYMBOL", y= "Expression",fill="Ploidy",palette = c("royalblue","red3"),width = 0.5, position = position_dodge(0.6), add="mean_se")+
  stat_pvalue_manual(data = stat, size=4, label = "p.format", y.position = 90, x = "SYMBOL",hide.ns = F)+xlab("")+geom_hline(yintercept = 0, linetype=2)+
  scale_y_continuous(breaks=seq(-100,100,10))+ylab("Difference in expression (%)")+theme(axis.text.x = element_text(angle = 90), text = element_text(size = 12))

## GPCRs genes ##
GPCR_genes <- read.delim("GPCRs_genes.txt")
GPCR_genes <- GPCR_genes[,c(2,13)]
colnames(GPCR_genes) <- c("SYMBOL","Group")
ploidy <- c("2n","2n","2n","2n","3n","3n","3n","3n","3n","3n","3n")
cpm_GPCRs <- merge(cpm_normalized[,c(1,8:18)],GPCR_genes,all.y=T, by="SYMBOL")
cpm_GPCRs <- na.omit(cpm_GPCRs)
cpm_GPCRs[,2:12] <- 2^cpm_GPCRs[,2:12]
diploid_mean <- apply(cpm_GPCRs[,2:5],1,FUN = mean)
cpm_GPCRs[,2:12] <- cpm_GPCRs[,2:12]/diploid_mean
melted_GPCR<-melt(cpm_GPCRs)
melted_GPCR$variable <- rep(ploidy, each=264)
colnames(melted_GPCR)[3:4] <- c("Ploidy","Expression")

stat=compare_means(Expression~Ploidy, data = melted_GPCR, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_GPCR, x= "Group", y= "Expression",fill="Ploidy",palette = c("royalblue","red3"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif", y.position = 3.5,x = "Group",hide.ns = F)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## receptors kinases genes ##
RKs_genes <- read.delim("receptors_kinases_genes.txt")
RKs_genes <- RKs_genes[,c(2,13)]
colnames(RKs_genes) <- c("SYMBOL","Group")
ploidy <- c("2n","2n","2n","2n","3n","3n","3n","3n","3n","3n","3n")
cpm_RKs <- merge(cpm_normalized[,c(1,8:18)],RKs_genes,all.y=T, by="SYMBOL")
cpm_RKs <- na.omit(cpm_RKs)
cpm_RKs[,2:12] <- 2^cpm_RKs[,2:12]
diploid_mean <- apply(cpm_RKs[,2:5],1,FUN = mean)
cpm_RKs[,2:12] <- cpm_RKs[,2:12]/diploid_mean
melted_RK<-melt(cpm_RKs)
melted_RK$variable <- rep(ploidy, each=72)
colnames(melted_RK)[3:4] <- c("Ploidy","Expression")

stat=compare_means(Expression~Ploidy, data = melted_RK, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_RK, x= "Group", y= "Expression",fill="Ploidy",palette = c("royalblue","red3"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif",y.position = 1.5,x = "Group",hide.ns = F)+theme(axis.text.x = element_text(angle = 0))


## Solute carriers genes ##
SLCs_genes <- read.delim("solute_carriers_genes.txt")
SLCs_genes <- SLCs_genes[,c(2,13)]
colnames(SLCs_genes) <- c("SYMBOL","Group")
ploidy <- c("2n","2n","2n","2n","3n","3n","3n","3n","3n","3n","3n")
cpm_SLCs <- merge(cpm_normalized[,c(1,8:18)],SLCs_genes,all.y=T, by="SYMBOL")
cpm_SLCs <- na.omit(cpm_SLCs)
cpm_SLCs[,2:12] <- 2^cpm_SLCs[,2:12]
diploid_mean <- apply(cpm_SLCs[,2:5],1,FUN = mean)
cpm_SLCs[,2:12] <- cpm_SLCs[,2:12]/diploid_mean
melted_SLC<-melt(cpm_SLCs)
melted_SLC$variable <- rep(ploidy, each=346)
colnames(melted_SLC)[3:4] <- c("Ploidy","Expression")

stat=compare_means(Expression~Ploidy, data = melted_SLC, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_SLC, x= "Group", y= "Expression",fill="Ploidy",palette = c("royalblue","red3"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif",y.position = 1.5,x = "Group",hide.ns = F)+theme(axis.text.x = element_text(angle = 0))


## integrins genes ##
INTs_genes <- read.delim("integrin_genes.txt")
INTs_genes <- INTs_genes[,c(2,13)]
colnames(INTs_genes) <- c("SYMBOL","Group")
ploidy <- c("2n","2n","2n","2n","3n","3n","3n","3n","3n","3n","3n")
cpm_INTs <- merge(cpm_normalized[,c(1,8:18)],INTs_genes,all.y=T, by="SYMBOL")
cpm_INTs <- na.omit(cpm_INTs)
cpm_INTs[,2:12] <- 2^cpm_INTs[,2:12]
diploid_mean <- apply(cpm_INTs[,2:5],1,FUN = mean)
cpm_INTs[,2:12] <- cpm_INTs[,2:12]/diploid_mean
melted_INT<-melt(cpm_INTs)
melted_INT$variable <- rep(ploidy, each=24)
colnames(melted_INT)[3:4] <- c("Ploidy","Expression")

stat=compare_means(Expression~Ploidy, data = melted_INT, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_INT, x= "Group", y= "Expression",fill="Ploidy",palette = c("royalblue","red3"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif",y.position = 1.2,x = "Group",hide.ns = F)+theme(axis.text.x = element_text(angle = 0))


## Basic leucin zipper genes ##
BLZ_genes <- read.delim("Basic_leucin_zipper_genes.txt")
BLZ_genes <- BLZ_genes[,c(2,13)]
colnames(BLZ_genes) <- c("SYMBOL","Group")
ploidy <- c("2n","2n","2n","2n","3n","3n","3n","3n","3n","3n","3n")
cpm_BLZ <- merge(cpm_normalized[,c(1,8:18)],BLZ_genes,all.y=T, by="SYMBOL")
cpm_BLZ <- na.omit(cpm_BLZ)
cpm_BLZ[,2:12] <- 2^cpm_BLZ[,2:12]
diploid_mean <- apply(cpm_BLZ[,2:5],1,FUN = mean)
cpm_BLZ[,2:12] <- cpm_BLZ[,2:12]/diploid_mean
melted_BLZ<-melt(cpm_BLZ)
melted_BLZ$variable <- rep(ploidy, each=50)
colnames(melted_BLZ)[3:4] <- c("Ploidy","Expression")

stat=compare_means(Expression~Ploidy, data = melted_BLZ, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_BLZ, x= "Group", y= "Expression",fill="Ploidy",palette = c("royalblue","red3"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif",y.position = 1.5,x = "Group",hide.ns = F)+theme(axis.text.x = element_text(angle = 0))


## Basic helix-loop-helix genes ##
BHLH_genes <- read.delim("BHLH_genes.txt")
BHLH_genes <- BHLH_genes[,c(2,13)]
colnames(BHLH_genes) <- c("SYMBOL","Group")
ploidy <- c("2n","2n","2n","2n","3n","3n","3n","3n","3n","3n","3n")
cpm_BHLH <- merge(cpm_normalized[,c(1,8:18)],BHLH_genes,all.y=T, by="SYMBOL")
cpm_BHLH <- na.omit(cpm_BHLH)
cpm_BHLH[,2:12] <- 2^cpm_BHLH[,2:12]
diploid_mean <- apply(cpm_BHLH[,2:5],1,FUN = mean)
cpm_BHLH[,2:12] <- cpm_BHLH[,2:12]/diploid_mean
melted_BHLH<-melt(cpm_BHLH)
melted_BHLH$variable <- rep(ploidy, each=85)
colnames(melted_BHLH)[3:4] <- c("Ploidy","Expression")

stat=compare_means(Expression~Ploidy, data = melted_BHLH, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_BHLH, x= "Group", y= "Expression",fill="Ploidy",palette = c("royalblue","red3"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif",y.position = 1.5,x = "Group",hide.ns = F)+theme(axis.text.x = element_text(angle = 0))


## zync fingers genes ##
ZF_genes <- read.delim("zync_fingers_genes.txt")
ZF_genes <- ZF_genes[,c(2,13)]
colnames(ZF_genes) <- c("SYMBOL","Group")
ploidy <- c("2n","2n","2n","2n","3n","3n","3n","3n","3n","3n","3n")
cpm_ZF <- merge(cpm_normalized[,c(1,8:18)],ZF_genes,all.y=T, by="SYMBOL")
cpm_ZF <- na.omit(cpm_ZF)
cpm_ZF[,2:12] <- 2^cpm_ZF[,2:12]
diploid_mean <- apply(cpm_ZF[,2:5],1,FUN = mean)
cpm_ZF[,2:12] <- cpm_ZF[,2:12]/diploid_mean
melted_ZF<-melt(cpm_ZF)
melted_ZF$variable <- rep(ploidy, each=1461)
colnames(melted_ZF)[3:4] <- c("Ploidy","Expression")

stat=compare_means(Expression~Ploidy, data = melted_ZF, method="wilcox.test", methods.args=list(exact=T), p.adjust.method = "fdr",group.by = "Group")
ggbarplot(melted_ZF, x= "Group", y= "Expression",fill="Ploidy",palette = c("royalblue","red3"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif",y.position = 3,x = "Group",hide.ns = T)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## general transcription factors genes ##
GTF_genes <- read.delim("general_transcription_factors_genes.txt")
GTF_genes <- GTF_genes[,c(2,13)]
colnames(GTF_genes) <- c("SYMBOL","Group")
ploidy <- c("2n","2n","2n","2n","3n","3n","3n","3n","3n","3n","3n")
cpm_GTF <- merge(cpm_normalized[,c(1,8:18)],GTF_genes,all.y=T, by="SYMBOL")
cpm_GTF <- na.omit(cpm_GTF)
cpm_GTF[,2:12] <- 2^cpm_GTF[,2:12]
diploid_mean <- apply(cpm_GTF[,2:5],1,FUN = mean)
cpm_GTF[,2:12] <- cpm_GTF[,2:12]/diploid_mean
melted_GTF<-melt(cpm_GTF)
melted_GTF$variable <- rep(ploidy, each=45)
colnames(melted_GTF)[3:4] <- c("Ploidy","Expression")

stat=compare_means(Expression~Ploidy, data = melted_GTF, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_GTF, x= "Group", y= "Expression",fill="Ploidy",palette = c("royalblue","red3"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif",y.position = 1.5,x = "Group",hide.ns = T)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## TFs and DNA-binding genes ##
ploidy <- c("2n","2n","2n","2n","3n","3n","3n","3n","3n","3n","3n")
melted_TF_DB<-data.frame("SYMBOL"=c(melted_ZF$SYMBOL,melted_BLZ$SYMBOL,melted_BHLH$SYMBOL,melted_GTF$SYMBOL),
                         "Group"=c(melted_ZF$Group,melted_BLZ$Group,melted_BHLH$Group,melted_GTF$Group),
                         "Ploidy"=c(melted_ZF$Ploidy,melted_BLZ$Ploidy,melted_BHLH$Ploidy,melted_GTF$Ploidy),
                         "Expression"=c(melted_ZF$Expression,melted_BLZ$Expression,melted_BHLH$Expression,melted_GTF$Expression))

stat=compare_means(Expression~Ploidy, data = melted_TF_DB, method="wilcox.test",methods.args=list(exact=T), p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_TF_DB, x= "Group", y= "Expression",fill="Ploidy",palette = c("royalblue","red3"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif",y.position = 2.8,x = "Group",hide.ns = T)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ scale_y_continuous(breaks=seq(0,3,0.5))


## general groups of genes ##
ploidy <- c("2n","2n","2n","2n","3n","3n","3n","3n","3n","3n","3n")
melted_genes<-data.frame("SYMBOL"=c(melted_ion$SYMBOL,melted_RK$SYMBOL,melted_GPCR$SYMBOL, melted_SLC$SYMBOL,melted_INT$SYMBOL, melted_ZF$SYMBOL,melted_BLZ$SYMBOL,melted_BHLH$SYMBOL,melted_GTF$SYMBOL),
                         "Group"=c(rep("Ion channels",1892),rep("RKs",792),rep("GPCRs",2904),rep("SLCs",3806),rep("Integrins",264),rep("DBPs & TFs",16071),rep("DBPs & TFs",550),rep("DBPs & TFs",935),rep("DBPs & TFs",495)),
                         "Ploidy"=c(melted_ion$Ploidy, melted_RK$Ploidy,melted_GPCR$Ploidy,melted_SLC$Ploidy, melted_INT$Ploidy, melted_ZF$Ploidy,melted_BLZ$Ploidy,melted_BHLH$Ploidy,melted_GTF$Ploidy),
                         "Expression"=c(melted_ion$Expression,melted_RK$Expression,melted_GPCR$Expression,melted_SLC$Expression,melted_INT$Expression, melted_ZF$Expression,melted_BLZ$Expression,melted_BHLH$Expression,melted_GTF$Expression))

stat=compare_means(Expression~Ploidy, data = melted_genes, method="wilcox.test",methods.args=list(exact=T), p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_genes, x= "Group", y= "Expression",fill="Ploidy",palette = c("royalblue","red3"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.signif",y.position = 1.2,x = "Group",hide.ns = T)+ scale_y_continuous(breaks=seq(0,2,0.1))

## making a delta change plot:
melted_genes_delta <- melted_genes
melted_genes_delta$Expression <- (melted_genes_delta$Expression-1)*100

stat=compare_means(Expression~Ploidy, data = melted_genes_delta, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_genes_delta, x= "Group", y= "Expression",fill="Ploidy",palette = c("royalblue","red3"), position = position_dodge(0.8), add="mean_se")+
  scale_y_continuous(breaks=seq(-30,20,2.5))+xlab(label = "")+
  ylab(label = "Difference in expression (%)")+
  stat_pvalue_manual(data = stat, label = "p.adj",y.position = 11.5,x = "Group",hide.ns = F)+
  geom_hline(yintercept = 0,linetype=2)+theme(text = element_text(size = 16),axis.text.y = element_text(size = 12))+
  scale_x_discrete(labels = c("Ion\nchannels","RKs","GPCRs","SLCs",
                              "Integrins","DBPs\nand TFs"))



############# general crisprScores of membranal genes ############# 

general_groups_of_membranal_genes <- rbind(GPCR_genes,INTs_genes,ion_channels_genes,RKs_genes,SLCs_genes)
general_groups_of_membranal_genes$Group <- c(rep("GPCRs",1415),rep("Integrins",27),rep("Ion channels",220),rep("Receptor kinases",77),rep("SLCs",433))
hESCs_NPCs_crisprScores <- read.csv("1-s2.0-S1934590920302903-mmc2.csv")
hESCs_NPCs_crisprScores <- hESCs_NPCs_crisprScores[hESCs_NPCs_crisprScores$Gene.Symbol %in% unique(general_groups_of_membranal_genes$SYMBOL),]
hESCs_NPCs_crisprScores <- hESCs_NPCs_crisprScores[,-18]
hESCs_NPCs_crisprScores <- na.omit(hESCs_NPCs_crisprScores)
colnames(hESCs_NPCs_crisprScores)[2] <- "SYMBOL"

hESCs_NPCs_crisprScores <- merge(hESCs_NPCs_crisprScores[,c(2:4,9:12)],general_groups_of_membranal_genes, by="SYMBOL", all.x=T)
hESCs_NPCs_crisprScores <- hESCs_NPCs_crisprScores[,c(1,4,5,2,3,8)]
general_crisprScores_sig <- hESCs_NPCs_crisprScores[hESCs_NPCs_crisprScores$Neuroectoderm.p.value<0.05 | 
                                                      hESCs_NPCs_crisprScores$hESC.p.value <0.05,]

colnames(general_crisprScores_sig)[c(2,4)] <- c("hESCs","NPCs")
melted_general_crisprScores_sig <- melt(general_crisprScores_sig[,c(1,2,4,6)])
colnames(melted_general_crisprScores_sig)[3:4] <- c("Developmental_stage","CrisprScore")

stat=compare_means(CrisprScore~Developmental_stage, data = melted_general_crisprScores_sig, method="wilcox.test",methods.args=list(exact=T), p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_general_crisprScores_sig, x= "Group", y= "CrisprScore",fill="Developmental_stage",palette = c("darkgoldenrod1","darkturquoise","darkorchid3"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat,x= "Group",size=4.5, label = "p.format", y.position = 0.45, hide.ns = F)+
  scale_y_continuous(breaks=seq(-2,2,0.1))+geom_hline(yintercept = 0,linetype=2)+theme(text = element_text(size = 20))

stat=compare_means(CrisprScore~Developmental_stage, data = melted_general_crisprScores_sig, method="wilcox.test",methods.args=list(exact=T), p.adjust.method = "BH",group.by = "Group")
ggboxplot(melted_general_crisprScores_sig, x= "Group", y= "CrisprScore",fill="Developmental_stage",palette = c("darkgoldenrod1","darkturquoise","darkorchid3"),group.by = "Group")+  #<------ , add="mean_se"#
  stat_pvalue_manual(data = stat, label = "p.signif", y.position = 1.55, x = "Group", hide.ns = F)+
  scale_y_continuous(breaks=seq(-2,2,0.25))+geom_hline(yintercept = 0,linetype=2)


######## general crisprScores of ion channels genes (according to HGNC)  ########
general_ion_channels_genes <- ion_channels_genes
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Anoctamins","Chloride channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Bestrophins","Chloride channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Chloride voltage-gated channels","Chloride channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Chloride intracellular channels","Chloride channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Chloride channels, ATP-gated CFTR","Chloride channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Volume regulated anion channel subunits","Chloride channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Sodium voltage-gated channel alpha subunits","Sodium channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Sodium voltage-gated channel beta subunits","Sodium channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Sodium channels epithelial","Sodium channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Sodium leak channels, non selective","Sodium channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Acid sensing ion channel subunits","Sodium channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Cation channels sperm associated ","Calcium channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Two pore segment channels","Calcium channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Calcium voltage-gated channel alpha1 subunits","Calcium channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Calcium voltage-gated channel auxiliary alpha2delta subunits","Calcium channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Calcium voltage-gated channel auxiliary beta subunits","Calcium channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Calcium channel auxiliary gamma subunits","Calcium channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Ryanodine receptors","Calcium channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Inositol 1,4,5-triphosphate receptors","Calcium channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Potassium voltage-gated channels","Potassium channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Potassium inwardly rectifying channel subfamily J","Potassium channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Potassium two pore domain channel subfamily K ","Potassium channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Potassium calcium-activated channels","Potassium channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Potassium sodium-activated channel subfamily T ","Potassium channels")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Aquaporins","Porins")
general_ion_channels_genes$Group <- str_replace_all(general_ion_channels_genes$Group,"Voltage dependent anion channels","Porins")

hESCs_NPCs_crisprScores <- read.csv("1-s2.0-S1934590920302903-mmc2.csv")
ion_channels_crisprScores <- hESCs_NPCs_crisprScores[hESCs_NPCs_crisprScores$Gene.Symbol %in% unique(general_ion_channels_genes$SYMBOL),]
ion_channels_crisprScores <- ion_channels_crisprScores[,-18]
ion_channels_crisprScores <- na.omit(ion_channels_crisprScores)
colnames(ion_channels_crisprScores)[2] <- "SYMBOL"

ion_channels_crisprScores <- merge(ion_channels_crisprScores[,c(2:4,9:12)],general_ion_channels_genes, by="SYMBOL", all.x=T)
ion_channels_crisprScores <- ion_channels_crisprScores[,c(1,4,5,2,3,8)]
general_crisprScores_sig <- ion_channels_crisprScores[ion_channels_crisprScores$Neuroectoderm.p.value<0.05 | 
                                                        ion_channels_crisprScores$hESC.p.value <0.05,]

colnames(general_crisprScores_sig)[c(2,4)] <- c("hESCs","NPCs")
melted_general_crisprScores_sig <- melt(general_crisprScores_sig[,c(1,2,4,6)])
colnames(melted_general_crisprScores_sig)[3:4] <- c("Developmental_stage","CrisprScore")

### reorder dataframe by delta expression:
melted_general_crisprScores_sig <- melted_general_crisprScores_sig[order(melted_general_crisprScores_sig$Group),]
melted_general_crisprScores_sig$Order <- c(rep(2,24),rep(5,18),rep(4,6),rep(1,30),rep(3,14))
melted_general_crisprScores_sig <- melted_general_crisprScores_sig[order(melted_general_crisprScores_sig$Order),]

stat=compare_means(CrisprScore~Developmental_stage, data = melted_general_crisprScores_sig, method="wilcox.test",methods.args=list(exact=T), p.adjust.method = "BH",group.by = "Group")
ggbarplot(melted_general_crisprScores_sig, x= "Group", y= "CrisprScore",fill="Developmental_stage",palette = c("darkgoldenrod1","darkturquoise","darkorchid3"), position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat,size = 4.5, label = "p.format", y.position = 0.6, x = "Group", hide.ns = F)+
  scale_y_continuous(breaks=seq(-2,2,0.1))+geom_hline(yintercept = 0,linetype=2)+xlab("")+theme(legend.text = element_text(size = 18),text = element_text(size = 18))

stat=compare_means(CrisprScore~Developmental_stage, data = melted_general_crisprScores_sig, method="wilcox.test", p.adjust.method = "BH",group.by = "Group")
ggboxplot(melted_general_crisprScores_sig, x= "Group", y= "CrisprScore",fill="Developmental_stage",palette = c("darkgoldenrod1","darkturquoise"),group.by = "Group")+  #<------ , add="mean_se"#
  scale_y_continuous(breaks=seq(-2,2,0.25))+geom_hline(yintercept = 0,linetype=2)



################  calculate TPM and add gene names ####################

# function to calculate TPM 
tpm_calc <- function(counts, lengths) {
  rate <- counts / lengths
  rate/sum(rate)*1*10^6 
}


# calculate TPM
tpm <- gene_exp_all
for(i in names(gene_exp_all[,8:18])) {
  tpm[i] <- tpm_calc(gene_exp_all[i],gene_exp_all$length_total_exons) 
}
write.csv(tpm,'ploidy_tpm.csv')  # write TPM table


###############################-----e_Karyotype-----###############################

tpm$chr <- str_remove(tpm$chr, "chr")
tpm_E <- tpm

floor=1 # Expression threshold value
expressed=0.8 # Percentage of expressed samples
variable=0.1 # Percentage of most variable genes to remove
tpm_E[,-c(1:7)][tpm_E[,-c(1:7)]<1]=floor                      #?????
tpm_E=tpm_E[rowSums(tpm_E[,-c(1:7)]>floor)>=ncol(tpm_E[,-c(1:7)])* expressed,]  #?????
row_mean=rowMeans(data.matrix(tpm_E[,-c(1:7)]))
distance_mat=sweep(tpm_E[,-c(1:7)],MARGIN=1,STATS=row_mean,"-")
ssq=rowSums(distance_mat^2)
tpm_E=tpm_E[ssq<quantile(ssq,1-variable),]

med=apply(data.matrix(tpm_E[,8:11]),1,median) #making the median of regular diploid ESCs the baseline for expression
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
  lines(c(sum+centromere_pos[i],sum+centromere_pos[i]),c(-1.02,-0.98),col="grey13",lwd=2)
  sum=sum+chr_size[i]
}


y=CGH$start+chr_total[CGH$chr]

window=300 # moving average window size
for(i in 8:18){
  lines(rollmean(y,window),rollmean(CGH[,i],window),col=cl[i-7])
}
par(mai=c(21, 6, 0, 6), xpd=TRUE)
legend(x = "bottom",horiz = F, legend = str_replace_all(names(cl),"_", " "), text.col = cl,ncol = 6, xpd = T, text.width = max(strwidth(names(cl)))/130, x.intersp = 4, y.intersp = 4, bty="n")


roll_mean <- data.frame("genomic_loci"=rollmean(y,window),"dsRed 2n rep1"=rollmean(CGH[,8],window),"dsRed 2n rep2"=rollmean(CGH[,9],window),
                        "EGFP 2n rep1"=rollmean(CGH[,10],window),"EGFP 2n rep2"=rollmean(CGH[,11],window),"3n clone H"=rollmean(CGH[,12],window),
                        "3n clone I"=rollmean(CGH[,13],window),"3n clone L rep1"=rollmean(CGH[,14],window),"3n clone L rep2"=rollmean(CGH[,15],window),
                        "3n clone K"=rollmean(CGH[,16],window),"3n clone J"=rollmean(CGH[,17],window),"3n clone B"=rollmean(CGH[,18],window))


melted_rollmean <- melt(data = roll_mean,id.vars = "genomic_loci")
melted_rollmean <-  data.frame(melted_rollmean,"Ploidy"=c(rep("2n",52240),rep("3n",91420)))
stat <- compare_means(value~Ploidy,data = melted_rollmean,method = "wilcox.test",p.adjust.method = "BH",group.by = "genomic_loci")
sig_genomic_loci <- stat[stat$p.adj<0.05,c(1,6)]
if (nrow(sig_genomic_loci)>0){
  points(x = sig_genomic_loci$genomic_loci, y = rep(0.8,nrow(sig_genomic_loci)),col = "blue", cex = 0.2,pch=8)
}



########################### relative e-Karyotype per chr ####################
tpm$chr <- str_remove(tpm$chr, "chr")
tpm$chr[tpm$chr=="X"] <- 23
tpm$chr <- as.numeric(tpm$chr)
tpm <- na.omit(tpm)
tpm <- tpm[order(tpm$chr,tpm$start),]
for (i in 1:23){
  tpm_E <-tpm[tpm$chr==i,]
  tpm_E[,8:18] <- (tpm_E[,8:18]+1)/(((tpm_E[,8]+tpm_E[,9]+tpm_E[,10]+tpm_E[,11])/4)+1)
  tpm_E <-data.frame(tpm_E[,1:7],"Diploid"=apply(tpm_E[8:11],1,mean),"Triploid"=apply(tpm_E[12:18],1,mean))
  dev.new()
  plot(1000000,1,col="white",pch=15,cex=0.4,ylim=c(0,2), ylab="Gene Expression",typ="l",xlab="",xaxt = "n",xlim=c(0,chr_size[i]))
  xtick<-as.numeric(format(seq(0, chr_size[i], by=5000000), scientific = F))
  axis(side=1, at = xtick, labels = xtick/1000000,pos = -5.5,cex.axis=0.75)
  lines(c(0,chr_size[i]),c(0,0),lwd=10,col="black")
  lines(c(0,chr_size[i]),c(0,0),lwd=8,col="gray50")
  lines(c(0,chr_size[i]),c(0,0),lwd=7,col="gray53")
  lines(c(0,chr_size[i]),c(0,0),lwd=6,col="gray59")
  lines(c(0,chr_size[i]),c(0,0),lwd=4,col="gray75")
  lines(c(0,chr_size[i]),c(0,0),lwd=2,col="gray85")
  lines(c(0,chr_size[i]),c(0,0),lwd=1,col="gray90")
  lines(c(centromere_pos[i],centromere_pos[i]),c(-0.01,0.01),col="grey13",lwd=2)
  window=30 # moving average window size
  for(j in c(8:9)){
    lines(frollmean(tpm_E$start,window),frollmean(tpm_E[,j],window),col=c("royalblue","red3")[j-7],)
  }
  par(mai=c(21, 6, 1, 6), xpd=TRUE)
  legend(x = "bottom", legend = c("Diploid pES10","Triploid pES10"), text.col = c("royalblue","red3"), xpd = T, ncol = 4, text.width = max(strwidth(ploidy))/10, x.intersp = 3,y.intersp = 3, bty="n")
  if (i==23){
    title(main = "Chr X")
  }
  else{
    title(main = paste("Chr",i))
  }
  #statistical analysis of the differences between diploid and triploid expression levels across chromosome #i:
  tpm_i <- tpm[tpm$chr==i,]
  roll_mean <- data.frame("genomic_loci"=frollmean(tpm_i$start,window),"dsRed 2n rep1"=frollmean(tpm_i[,8],window),"dsRed 2n rep2"=frollmean(tpm_i[,9],window),
                          "EGFP 2n rep1"=frollmean(tpm_i[,10],window),"EGFP 2n rep2"=frollmean(tpm_i[,11],window),"3n clone H"=frollmean(tpm_i[,12],window),
                          "3n clone I"=frollmean(tpm_i[,13],window),"3n clone L rep1"=frollmean(tpm_i[,14],window),"3n clone L rep2"=frollmean(tpm_i[,15],window),
                          "3n clone K"=frollmean(tpm_i[,16],window),"3n clone J"=frollmean(tpm_i[,17],window),"3n clone B"=frollmean(tpm_i[,18],window))
  roll_mean <- na.omit(roll_mean)
  
  melted_rollmean <- melt(data = roll_mean,id.vars = "genomic_loci")
  melted_rollmean <-  data.frame(melted_rollmean,"Ploidy"=c(rep("2n",nrow(roll_mean)*4),rep("3n",nrow(roll_mean)*7)))
  stat <- compare_means(value~Ploidy,data = melted_rollmean,method = "wilcox.test",p.adjust.method = "BH",group.by = "genomic_loci")
  sig_genomic_loci <- stat[stat$p.adj<0.05,c(1,6)]
  points(x = sig_genomic_loci$genomic_loci, y = rep(1.9,nrow(sig_genomic_loci)),col = "blue", cex = 0.2,pch=8)
  
}


####### Cell Cycle Analysis ########
setwd("C:/Users/owner/Desktop/FACS results")
FACS_results <- read.csv("cell_cycle_analysis_2n_vs_3n (Dean Jett Fox model).csv")
# FACS_results <- FACS_results[FACS_results$Depth=="> > > ",]
# FACS_results <- data.frame(FACS_results,"Ploidy"=c(rep("3n",27),rep("2n",15)))
# FACS_results$Name <- rep(c("G1","S","G2/M"),14)

stat=compare_means(Percentage~Ploidy, data = FACS_results, method="t.test", p.adjust.method = "BH", group.by = "Cell.cycle.phase")
ggbarplot(FACS_results, x= "Cell.cycle.phase", y= "Percentage",fill="Ploidy",palette = c("royalblue","red"), position = position_dodge(0.8), add="mean_se")+
  scale_y_continuous(breaks=seq(0,70,10))+xlab(label = "")+
  ylab(label = "Cells in cell cycle (%)")+
  stat_pvalue_manual(data = stat, label = "p.format",y.position = 60,x = "Cell.cycle.phase",hide.ns = F)

proliferation_status <- data.frame("Cell.cycle.phase"=c(rep("G1",16),rep("S/G2/M",16)),"Ploidy"=rep(c(rep("3n",10),rep("2n",6)),2),
                                   "Percentage"=c(FACS_results$Percentage[FACS_results$Cell.cycle.phase=="G1"],FACS_results$Percentage[FACS_results$Cell.cycle.phase=="S"]+FACS_results$Percentage[FACS_results$Cell.cycle.phase=="G2/M"]))

stat=compare_means(Percentage~Ploidy, data = proliferation_status, method="t.test", p.adjust.method = "BH", group.by = "Cell.cycle.phase")
ggbarplot(proliferation_status, x= "Cell.cycle.phase", y= "Percentage",fill="Ploidy",palette = c("royalblue","red"), position = position_dodge(0.8), add="mean_se")+
  scale_y_continuous(breaks=seq(0,70,10))+xlab(label = "")+
  ylab(label = "Cells in cell cycle (%)")+
  stat_pvalue_manual(data = stat, label = "p.adj",y.position = 60,x = "Cell.cycle.phase",hide.ns = F)


##### cell cycle analysis w/o problematic (contact inhibited) samples #####
FACS_results <- FACS_results[FACS_results$Sample=="Fusion I_3n clone H_001.fcs" |
                               FACS_results$Sample=="Fusion I_3n clone H_002.fcs" | 
                               FACS_results$Sample=="Fusion I_3n clone L_003.fcs" |
                               FACS_results$Sample=="10R 4c p38_005.fcs" |
                               FACS_results$Sample=="exp_1_10G 4c p47_001.fcs" |
                               FACS_results$Sample=="exp_1_10R 4c p42_002.fcs"|
                               FACS_results$Sample=="exp_1_10R 4c p41_002.fcs",]

stat=compare_means(Percentage~Ploidy, data = FACS_results, method="t.test", p.adjust.method = "BH", group.by = "Cell.cycle.phase")
ggbarplot(FACS_results, x= "Cell.cycle.phase", y= "Percentage",fill="Ploidy",palette = c("royalblue","red"), position = position_dodge(0.8), add="mean_se")+
  scale_y_continuous(breaks=seq(0,70,10))+xlab(label = "")+
  ylab(label = "Cells in cell cycle (%)")+
  stat_pvalue_manual(data = stat, label = "p.format",y.position = 60,x = "Cell.cycle.phase",hide.ns = F)

proliferation_status <- data.frame("Cell.cycle.phase"=c(rep("G1",7),rep("S/G2/M",7)),"Ploidy"=rep(c(rep("3n",3),rep("2n",4)),2),
                                   "Percentage"=c(FACS_results$Percentage[FACS_results$Cell.cycle.phase=="G1"],FACS_results$Percentage[FACS_results$Cell.cycle.phase=="S"]+FACS_results$Percentage[FACS_results$Cell.cycle.phase=="G2/M"]))

#test for normality:
for (i in c("G1","S/G2/M")){
  print(shapiro.test(proliferation_status[proliferation_status$Ploidy=="2n" & proliferation_status$Cell.cycle.phase==i,3]))
}
for (i in c("G1","S/G2/M")){
  print(shapiro.test(proliferation_status[proliferation_status$Ploidy=="3n" & proliferation_status$Cell.cycle.phase==i,3]))
}

stat=compare_means(Percentage~Ploidy, data = proliferation_status, method="t.test", p.adjust.method = "BH", group.by = "Cell.cycle.phase")
ggbarplot(proliferation_status, x= "Cell.cycle.phase", y= "Percentage",fill="Ploidy",palette = c("royalblue","red"), position = position_dodge(0.8), add="mean_se")+
  scale_y_continuous(breaks=seq(0,70,10))+xlab(label = "")+
  ylab(label = "Cells in cell cycle (%)")+
  stat_pvalue_manual(data = stat, label = "p.format",y.position = 60,x = "Cell.cycle.phase",hide.ns = F)


##### growth assay #####
setwd("C:/Users/owner/Desktop/Growth assay/exp. 1")
growth <- read.csv("intensity_levels_summary_(6K_per_well).csv")
growth <- data.frame(growth,"Ploidy"=rep(c(rep("3n",12),rep("2n",12)),5))

stat=compare_means(Intensity~Ploidy, data = growth, method="t.test", p.adjust.method = "BH", group.by = "timepoint")
ggplot(growth, aes(x=timepoint, y=Intensity, group=Sample)) +
  stat_summary(fun = mean, geom='line', aes(color=Sample), lty=1, size=1.3) +
  stat_summary(fun=mean,geom='point') +
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',width=0.2) +
  theme_bw()+scale_colour_manual(values=c("deepskyblue1","deepskyblue3","royalblue1","royalblue3","red1","red2","red3","red4"))

relative_growth <- growth
day_1 <- relative_growth$Intensity[relative_growth$timepoint=="Day 1"]
Mean <- c(rep(mean(day_1[1:3]),3),rep(mean(day_1[4:6]),3),rep(mean(day_1[7:9]),3),rep(mean(day_1[10:12]),3),
          rep(mean(day_1[13:15]),3),rep(mean(day_1[16:18]),3),rep(mean(day_1[19:21]),3),rep(mean(day_1[22:24]),3))
relative_growth$Intensity <- relative_growth$Intensity/Mean

stat=compare_means(Intensity~Ploidy, data = relative_growth, method="t.test", p.adjust.method = "BH", group.by = "timepoint")
ggplot(relative_growth, aes(x=timepoint, y=Intensity, group=Sample)) +
  stat_summary(fun = mean, geom='line', aes(color=Sample), lty=1, size=1.3) +
  stat_summary(fun=mean,geom='point') +
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',width=0.2) +
  theme_bw()+scale_colour_manual(values=c("deepskyblue1","deepskyblue3","royalblue1","royalblue3","red1","red2","red3","red4"))+
  stat_compare_means(aes(group = Ploidy), label = "p.format", label.y = c(2,3,5,12,15))

###### w/o outlier Samples ######
relative_growth <- relative_growth[relative_growth$Sample!="pES10 3n clone L",]
relative_growth <- relative_growth[relative_growth$Sample!="h-pES10 10R 4c rep2",]

#test for normality:
for (i in c("Day 1","Day 2","Day 3","Day 4","Day 5")){
  print(shapiro.test(relative_growth[relative_growth$Ploidy=="2n" & relative_growth$timepoint==i,2]))
}
for (i in c("Day 1","Day 2","Day 3","Day 4","Day 5")){
  print(shapiro.test(relative_growth[relative_growth$Ploidy=="3n" & relative_growth$timepoint==i,2]))
}

stat=compare_means(Intensity~Ploidy, data = relative_growth, method="t.test", p.adjust.method = "BH", group.by = "timepoint")
ggplot(relative_growth, aes(x=timepoint, y=Intensity, group=Sample)) +
  stat_summary(fun = mean, geom='line', aes(color=Sample), lty=1, size=1.3) +
  stat_summary(fun=mean,geom='point') + ylab("Relative intensity")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',width=0.2) +
  theme_bw()+scale_colour_manual(values=c("deepskyblue1","deepskyblue3","royalblue1","red1","red2","red3"))+
  stat_compare_means(aes(group = Ploidy), label = "p.format", label.y = c(2,3,5,7,9))+
  scale_y_continuous(breaks=seq(0,10,1))

###### average, w/o outliers ######
stat=compare_means(Intensity~Ploidy, data = relative_growth, method="t.test", p.adjust.method = "BH", group.by = "timepoint")
ggplot(relative_growth, aes(x=timepoint, y=Intensity, group=Ploidy)) +
  stat_summary(fun = mean, geom='line', aes(color=Ploidy), lty=1, size=1.3) +
  stat_summary(fun=mean,geom='point') + ylab("Relative intensity")+
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',width=0.2,aes(color=Ploidy)) +
  theme_bw()+scale_colour_manual(values=c("royalblue","red"))+
  stat_compare_means(aes(group = Ploidy), label = "p.format", label.y = c(2,3,4,5.5,7.5))+
  scale_y_continuous(breaks=seq(0,10,1))


####### RNA content Analysis ########
setwd("C:/Users/owner/Desktop/AO staining 2n vs 3n")
AO_1 <- read.csv("statistics AO-1.csv")
AO_2 <- read.csv("statistics AO-2.csv")

mean_RNA_2n_1 <- mean(AO_1$mean.RNA[1:2])
relative_AO_1 <- AO_1$mean.RNA/mean_RNA_2n_1
mean_RNA_2n_2 <- mean(AO_2$mean.RNA[1:2])
relative_AO_2 <- AO_2$mean.RNA/mean_RNA_2n_2

relative_RNA <- data.frame("Sample"=rep(AO_1$Sample,2),"RNA_content"=c(relative_AO_1,relative_AO_2),"Ploidy"=c("2n","2n","3n","3n","2n","2n","3n","3n"))
#test for normality:
shapiro.test(relative_RNA[c(1,2,5,6),2])
shapiro.test(relative_RNA[c(3,4,7,8),2])

stat=compare_means(RNA_content~Ploidy, data = relative_RNA, method="t.test", p.adjust.method = "BH")
ggbarplot(relative_RNA, x= "Ploidy", y= "RNA_content",fill="Ploidy",palette = c("royalblue","red"), position = position_dodge(0.8),label = T,lab.vjust = -3,lab.nb.digits = 3, add="mean_se")+
  ylab(label = "Relative RNA content")+
  stat_pvalue_manual(data = stat, label = "p.format",y.position = 1.9,hide.ns = F)

####### DNA content Analysis ########
mean_DNA_2n_1 <- mean(AO_1$mean.DNA[1:2])
relative_AO_1 <- AO_1$mean.DNA/mean_DNA_2n_1
mean_DNA_2n_2 <- mean(AO_2$mean.DNA[1:2])
relative_AO_2 <- AO_2$mean.DNA/mean_DNA_2n_2

relative_DNA <- data.frame("Sample"=rep(AO_1$Sample,2),"DNA_content"=c(relative_AO_1,relative_AO_2),"Ploidy"=c("2n","2n","3n","3n","2n","2n","3n","3n"))

stat=compare_means(DNA_content~Ploidy, data = relative_DNA, method="t.test", p.adjust.method = "BH")
ggbarplot(relative_DNA, x= "Ploidy", y= "DNA_content",fill="Ploidy",palette = c("royalblue","red"), position = position_dodge(0.8),label = T,lab.vjust = -2.5,lab.nb.digits = 3, add="mean_se")+
  ylab(label = "Relative DNA content")+
  stat_pvalue_manual(data = stat, label = "p.format",y.position = 1.7,hide.ns = F)

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

setwd("C:/Users/owner/Desktop/RNA-seq counts guy/RNA-seq EBs and teratomas")
# setwd("C:/Users/guyha/Desktop/RNA-seq counts guy/RNA-seq EBs and teratomas")

#loading RNA-seq data from each sample into a data.frame and getting rid of unnecessary data (first 4 rows and second and third columns)
EGFP_neo_2n_T1 <- read.delim("T1-10G-4c-p41-T1_S8_counts.txt", skip = 1)[,c(1,7)]
EGFP_neo_2n_T2 <- read.delim("T2_h_pES10_EGFP_neo_p41_4c_T2_1_counts.txt", skip = 1)[,c(1,7)]
dsRed_hygro_2n_T1 <- read.delim("T1_h_pES10_dsRed_hygro_4c_p39_T1_3_counts.txt", skip = 1)[,c(1,7)]
dsRed_hygro_2n_T2 <- read.delim("T3_h_pES10_dsRed_hygro_4c_p39_T2_3_counts.txt", skip = 1)[,c(1,7)]
clone_H_3n_T <- read.delim("T5_h_pES10_fusion_I_T1_3n_f_4_clone_H_counts.txt", skip = 1)[,c(1,7)]
clone_K_3n_T <- read.delim("T4-pES10-3n-clone-K-T1_S11_counts.txt", skip = 1)[,c(1,7)]
clone_L1_3n_T <- read.delim("T5-pES10-3n-clone-L-T1_S12_counts.txt", skip = 1)[,c(1,7)] 
clone_L2_3n_T <- read.delim("T4_h_pES10_fusion_I_T1_3n_f_6_clone_L_counts.txt", skip = 1)[,c(1,7)] 
CSES4_T1 <- read.delim("teratoma_emp1ReadsPerGene.out.tab", skip = 1)[,c(1,2)]    #male diploid ESCs
CSES4_T2 <- read.delim("teratoma_emp2ReadsPerGene.out.tab", skip = 1)[,c(1,2)]    #male diploid ESCs
colnames(CSES4_T1)[1] <- "Geneid"
colnames(CSES4_T2)[1] <- "Geneid"

dip_pES10_rep1 <- read.delim("C:/Users/owner/Desktop/RNA-seq counts guy/Ido's samples/SRR2131926_ReadsPerGene.out.tab", skip = 1)[-c(1:4),1:2] 
dip_pES10_rep2 <- read.delim("C:/Users/owner/Desktop/RNA-seq counts guy/Ido's samples/SRR2131927_ReadsPerGene.out.tab", skip = 1)[-c(1:4),1:2] 
colnames(dip_pES10_rep1)[1] <- "Geneid"
colnames(dip_pES10_rep2)[1] <- "Geneid"


###merging all RNA-seq data into one data.frame and adding gene names according to annotation file:
gene_exp_all <- Reduce(function(x, y) merge(x, y, all=T, by = "Geneid"), 
                       list(dsRed_hygro_2n_T1,dsRed_hygro_2n_T2,EGFP_neo_2n_T1,EGFP_neo_2n_T2,clone_H_3n_T,clone_K_3n_T,clone_L1_3n_T,clone_L2_3n_T,
                            CSES4_T1,CSES4_T2,dip_pES10_rep1,dip_pES10_rep2))
colnames(gene_exp_all) <- c("ENSEMBL","dsRed-hygro 2n Ter1","dsRed-hygro 2n Ter2","EGFP-neo 2n Ter1","EGFP-neo 2n Ter2",
                            "clone H 3n Ter","clone K 3n Ter","clone L1 3n Ter","clone L2 3n Ter",
                            "CSES4_T1","CSES4_T2","dip_pES10_rep1","dip_pES10_rep2")
gene_lengths <- read.csv("gencode.v34.transcriptLengths.csv")[,-8]

#removing the end of the ENSEMBL ID so that the same genes from different versions of the annotation file will be joined together  
gene_exp_all$ENSEMBL <- str_remove(gene_exp_all$ENSEMBL, "[.].*")
gene_lengths$ENSEMBL <- str_remove(gene_lengths$ENSEMBL, "[.].*")
gene_exp_all <- aggregate(gene_exp_all[-1], list(gene_exp_all$ENSEMBL), FUN = sum, na.rm = T)
colnames(gene_exp_all)[1] <- "ENSEMBL"


#merging all data frame to one count table with all the data needed:
gene_exp_all <- merge(gene_lengths, gene_exp_all, by = "ENSEMBL", all.y = T)
gene_exp_all <- gene_exp_all[rowSums(is.na(gene_exp_all[ ,2:19])) == 0, ] #removing empty rows
indRemoved <- which(apply(gene_exp_all[,8:19], 1, function(x) all(x == 0)) ) #saving all unexpressed genes in the data.frame
gene_exp_all <- gene_exp_all[-indRemoved,]      #getting rid of rows (genes) that all of their columns are zeros
.rowNamesDF(gene_exp_all, make.names=T) <- gene_exp_all$SYMBOL
indRemoved_Y <- which(gene_exp_all$chr=="chrY") #saving all Y chr genes in the data.frame
gene_exp_all <- gene_exp_all[-indRemoved_Y,]      #getting rid of rows (genes) that all of their columns are zeros
write.csv(gene_exp_all,'teratomas_read_counts_table.csv')  # write read counts table
X_genes <- gene_exp_all[gene_exp_all$chr=="chrX",7]

cl=c("deepskyblue1","deepskyblue3","royalblue1","royalblue3","red1","red2","red3","red4","lightgray","darkgray","darkblue","midnightblue")
names(cl) <- str_replace_all(colnames(gene_exp_all[,8:19]),"_", " ")


#### Differential expression (DE) analysis of Teratomas (with male teratomas):

ploidy <- factor(c(1,1,1,1,2,2,2,2,1,1), labels = c("2n","3n")) #4 teratomas samples of recently diploids ("normal") and 4 teratomas samples of triploid cells
DE_object <- DGEList(gene_exp_all[,8:17], group = ploidy, genes = gene_exp_all$gene)
keep <- rowSums(cpm(DE_object)>2) >= 3
DE_object <- DE_object[keep, , keep.lib.sizes=FALSE]
DE_TMM <- calcNormFactors(DE_object)

###PCA analysis of log2(CPM) values of each sample
cpm_normalized <- cpm(DE_TMM, log = T) #log2(cpm) after TMM normalization
cpm_normalized_2 <- cpm(DE_TMM, log = F) #cpm after TMM normalization
cpm_normalized <- data.frame("SYMBOL" = rownames(cpm_normalized), cpm_normalized)
cpm_normalized_2 <- data.frame("SYMBOL" = rownames(cpm_normalized_2), cpm_normalized_2)
cpm_normalized <- merge(gene_lengths, cpm_normalized, by = "SYMBOL",all.y = T)
cpm_normalized_2 <- merge(gene_lengths, cpm_normalized_2, by = "SYMBOL",all.y = T)
.rowNamesDF(cpm_normalized, make.names=T) <- cpm_normalized$SYMBOL
.rowNamesDF(cpm_normalized_2, make.names=T) <- cpm_normalized_2$SYMBOL
cpm_normalized_X <- cpm_normalized_2[cpm_normalized_2$SYMBOL %in% X_genes,]

prin_comp <- prcomp(t(cpm_normalized[,8:17]), center = T, retx = T)
gg <- cbind.data.frame(prin_comp$x, "Sample" =colnames(gene_exp_all)[8:17])
gg <- data.frame(gg, "Ploidy"=ploidy)
autoplot(prin_comp, data = gg, colour = "Sample",shape="Ploidy", label =F,size=3)+ scale_colour_manual(values=cl[1:10])+theme_bw()+
  guides(shape = guide_legend(order = 1),colour = guide_legend(override.aes = list(shape=c(15))))

summary(prin_comp)


### X chr gene expression ###
cpm_normalized_X <- cpm_normalized_X[order(cpm_normalized_X$start, decreasing = F),]
pheatmap(data.matrix(cpm_normalized_X[,8:17]), scale = "row", cluster_rows=F, show_rownames =F,clustering_distance_cols = "maximum",color = colorRampPalette(c("darkblue", "white", "red3"))(100)) #TPM heatmap of X-linked expressed genes, without gene clustering

avg_cpm_normalized_X <- data.frame(cpm_normalized_X[,1:7], "Diploid teratomas"=(cpm_normalized_X[,8]+cpm_normalized_X[,9]+cpm_normalized_X[,10]+cpm_normalized_X[,11])/4,
                                   "Triploid teratomas"=(cpm_normalized_X[,12]+cpm_normalized_X[,13]+cpm_normalized_X[,14]+cpm_normalized_X[,15])/4,
                                   "Diploid male teratomas"=(cpm_normalized_X[,16]+cpm_normalized_X[,17])/2)
pheatmap(data.matrix(avg_cpm_normalized_X[,8:10]), scale = "row", cluster_rows=F, show_rownames =F,angle_col = 0,color = colorRampPalette(c("darkblue", "white", "red3"))(100)) #averaged TPM heatmap of X-linked expressed genes, without gene clustering
pheatmap(data.matrix(avg_cpm_normalized_X[,8:10]), scale = "row", cluster_rows=F, show_rownames =F,angle_col = 0) #averaged TPM heatmap of X-linked expressed genes, without gene clustering


#### Differential expression (DE) analysis of Teratomas (without male teratomas):

ploidy <- factor(c(1,1,1,1,2,2,2,2), labels = c("2n","3n")) #4 teratomas samples of recently diploids ("normal") and 4 teratomas samples of triploid cells
DE_object <- DGEList(gene_exp_all[,8:15], group = ploidy, genes = gene_exp_all$gene)
keep <- rowSums(cpm(DE_object)>2) >= 3
DE_object <- DE_object[keep, , keep.lib.sizes=FALSE]
DE_TMM <- calcNormFactors(DE_object)

###PCA analysis of log2(CPM) values of each sample
cpm_normalized <- cpm(DE_TMM, log = T) #log2(cpm) after TMM normalization
cpm_normalized_2 <- cpm(DE_TMM, log = F) #cpm after TMM normalization
cpm_normalized <- data.frame("SYMBOL" = rownames(cpm_normalized), cpm_normalized)
cpm_normalized_2 <- data.frame("SYMBOL" = rownames(cpm_normalized_2), cpm_normalized_2)
cpm_normalized <- merge(gene_lengths, cpm_normalized, by = "SYMBOL",all.y = T)
cpm_normalized_2 <- merge(gene_lengths, cpm_normalized_2, by = "SYMBOL",all.y = T)
.rowNamesDF(cpm_normalized, make.names=T) <- cpm_normalized$SYMBOL
.rowNamesDF(cpm_normalized_2, make.names=T) <- cpm_normalized_2$SYMBOL
cpm_normalized_X <- cpm_normalized_2[cpm_normalized_2$SYMBOL %in% X_genes,]

prin_comp <- prcomp(t(cpm_normalized[,8:15]), center = T, retx = T)
gg <- cbind.data.frame(prin_comp$x, "Sample" =colnames(gene_exp_all)[8:15])
gg <- data.frame(gg, "Ploidy"=ploidy)
autoplot(prin_comp, data = gg, colour = "Sample",shape="Ploidy", label =F,size=3)+ scale_colour_manual(values=cl[1:8])+theme_bw()+
  guides(shape = guide_legend(order = 1),colour = guide_legend(override.aes = list(shape=c(15))))


### continue DE analysis ###
batch <- factor(c(1,1,2,1,1,2,2,1), labels = c("batch1","batch2"))
design <- model.matrix(~batch + ploidy)
DE_TMM <- estimateDisp(DE_TMM, design, robust=TRUE)
fit_TMM <- glmFit(DE_TMM, design)   #fit the model
LRT_3n_Vs_2n <- glmLRT(fit_TMM)  #comparing triploids to diploids

summary(decideTests(LRT_3n_Vs_2n))  #decide by FDR
TMM_3n <- topTags(LRT_3n_Vs_2n, n=nrow(DE_TMM))$table #assigning fold change gene expression between diploid and triploid ESCs
TMM_3n <- data.frame("SYMBOL"=row.names(TMM_3n), TMM_3n)
colnames(TMM_3n)[2] <- "logFC_triploids"
TMM_3n <- merge(gene_lengths, TMM_3n, by = "SYMBOL", all.y = T)
sig_DE_3n <-TMM_3n[TMM_3n$FDR<=0.05,] #assigning only genes with significant fold change between triploid and diploid ESCs (by FDR) 
sig_DE_3n <- sig_DE_3n[-rowSums(is.na(sig_DE_3n[ ,2:7])) == 0, ] #removing empty rows
sig_DE_3n<- sig_DE_3n[order(sig_DE_3n$logFC_triploids, decreasing = T),]
sig_DE_3n<- data.frame(sig_DE_3n[,1:8],"Fold.Change"=2^sig_DE_3n$logFC_triploids,sig_DE_3n[,9:12])


####Perform GO analysis via GO Enrichment Analysis in THE GENE ONTOLOGY RESOURCE, STRING Database and GSEA-MSigDB website:
library(fgsea)
library(dplyr)

# ranks <- TMM_3n[TMM_3n$FDR<=0.05 & abs(TMM_3n$logFC_triploids)>1,]
ranks <- TMM_3n
ranks[, 'score'] <-ranks$logFC_triploids* (-log(ranks$FDR))

# ranks[ranks$logFC_triploids < 0, 'score'] <-
#   ranks[ranks$logFC_triploids < 0,'logFC_triploids'] + log(ranks[ranks$logFC_triploids < 0,'FDR'])

.rowNamesDF(ranks, make.names=T) <- ranks$SYMBOL

ranks <- setNames(ranks$score,ranks$SYMBOL)
cell_types <- gmtPathways("c8.all.v7.4.symbols.gmt") 
fgseaRes <- fgseaMultilevel(cell_types, ranks, minSize=20, maxSize=500,eps = 0)
fgseaRes <- na.omit(fgseaRes) ##removing pathways which do not have an estimate of the error deviation of the calculated p-values (where log2err=NA)
fgseaRes <- fgseaRes[fgseaRes$padj < 0.05,]

collapsedPathways <- collapsePathways(fgseaRes[order(padj)][fgseaRes$padj < 0.05,], cell_types, ranks, nperm = 100000,pval.threshol = 0.01)
mainPathways_fdr <- fgseaRes[pathway %in% collapsedPathways$mainPathways][order(padj), pathway]
dev.new()
plotGseaTable(cell_types[mainPathways_fdr], ranks, fgseaRes, gseaParam = 0.5)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

fgseaResTidy <-fgseaResTidy[fgseaResTidy$padj<0.01,]
fgseaResTidy <- fgseaResTidy[fgseaResTidy$pathway %in% mainPathways_fdr,]

ggplot(data = fgseaResTidy[fgseaResTidy$NES<0,], aes(reorder(pathway, -padj), -log2(padj))) +
  geom_col(aes(fill=NES)) +
  coord_flip() +
  labs(x="Cell Type", y=expression(-log[2](FDR))) + 
  theme_minimal()+ 
  geom_hline(yintercept = -log2(0.05), linetype=2,color = "red", size=1)

fgseaResTidy <- fgseaResTidy[order(fgseaResTidy$padj),]

#select the top 3 pathways(could be different numbers based on the permutation test results and number of overlapping cell types):
ggplot(data = fgseaResTidy[c(1,9,10),], aes(reorder(pathway, -padj), -log2(padj))) +
  geom_col(aes(fill=0),width = 0.5) + theme(axis.text.y=element_text(size= 8,face = "bold",hjust = 0.5))+
  coord_flip() +
  labs(x="Cell Type", y=expression(-log[2](FDR))) + 
  geom_hline(yintercept = -log2(0.05), linetype=2,color = "black", size=0.25)+theme(legend.position = "")+
  scale_x_discrete(labels=rev(c("Midbrain\nNeurotypes","PFC\nMicroglia","Visceral\nNeurons")))


############## visualization of GSEA ###########################
library(enrichplot)
library(DOSE)
library(clusterProfiler)
library(ggnewscale)

de <-sort(ranks,decreasing = T)

gene.df <- bitr(names(de), fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
de2 <- data_frame("SYMBOL"=names(de),"score"=de)
de2 <- merge(gene.df,de2, by="SYMBOL", all.x = T)
de2 <- de2[,3:4]
de2 <- sort(setNames(de2$score,de2$ENTREZID),decreasing = T)

edo <- gseDGN(de2,pvalueCutoff = 0.01)
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
edox <- pairwise_termsim(edox)
cnetplot(edox, categorySize="qvalue",showCategory = 7, foldChange=ranks, colorEdge = TRUE, cex_label_gene = 0.5)+
  scale_color_gradientn(name = expression(paste(Log[2],"(FC)")),
                        colours = c("blue4","lightskyblue","red"),
                        limits= c(-90, 10),values = c(0,0.9,1) )


ego <- gseGO(geneList = de,
                 OrgDb = org.Hs.eg.db,
                 keyType = "SYMBOL",
                 ont = "CC", eps = 0,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.01)

ego <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
ego <- pairwise_termsim(ego)
cnetplot(ego, categorySize="qvalue",showCategory = 5, foldChange=ranks, colorEdge = TRUE, cex_label_gene = 0.5)+
  scale_color_gradientn(name = expression(paste(Log[2],"(FC)")),
                        colours = c("blue4","lightskyblue","red"),
                        limits= c(-90, 10),values = c(0,0.9,1) )


ego2 <- gseGO(geneList = de,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP", eps = 0,
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01)
ego2 <- simplify(ego2, cutoff=0.7, by="p.adjust", select_fun=min)
ego2 <- pairwise_termsim(ego2)
cnetplot(ego2, foldChange=ranks,showCategory = 5,categorySize="qvalue", cex_label_gene = 0.5, colorEdge = TRUE)+
  scale_color_gradientn(name = expression(paste(Log[2],"(FC)")),
                        colours = c("blue4","lightskyblue","red"),
                        limits= c(-90, 10),values = c(0,0.9,1) )


ego3 <- gseGO(geneList = de,
                 OrgDb = org.Hs.eg.db,
                 keyType = "SYMBOL",
                 ont = "MF", eps = 0,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05)
ego3 <- pairwise_termsim(ego3)

### Ranks, pvalues and padjust are too similar so the simplification is done manually:

#taking the MFs that were the most depleted:
top5_dep_MFs <- ego3@result[["Description"]][which(ego3@result[["rank"]] %in% sort(unique(ego3@result[["rank"]]),decreasing = T)[1:5])]
ego3@result <- ego3@result[ego3@result[["Description"]] %in% top5_dep_MFs,]

#getting rid of redundant MFs:
dependent_MF <- c("guanyl nucleotide binding","guanyl ribonucleotide binding","nucleoside binding","purine ribonucleoside binding",
                  "purine nucleoside binding","ribonucleoside binding","microtubule binding","DNA-binding transcription repressor activity, RNA polymerase II-specific")
ego3@result <- ego3@result[!(ego3@result[["Description"]] %in% dependent_MF),]

cnetplot(ego3, foldChange=ranks, showCategory = 6,categorySize="qvalue", cex_label_gene = 0.5, colorEdge = TRUE,layout = "dh")+
  scale_color_gradientn(name = expression(paste(Log[2],"(FC)")),
                        colours = c("blue4","lightskyblue","red"),
                        limits= c(-70, 10),values = c(0,0.875,1) )

#Layout of the map, e.g. 'star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 'randomly', 'fr', 'kk', 'drl' or 'lgl'


# function to calculate TPM 
tpm_calc <- function(counts, lengths) {
  rate <- counts / lengths
  rate/sum(rate)*1*10^6 
}

#######TeratoScore#######
gct <- read.delim(file= "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", skip=2)
#removing the end of the ENSEMBL ID so that the same genes from different versions of the annotation file will be joined together  

gct <- data.frame(gct[,-c(3,4,13:15,21,28,29,34:37,47:51,56)],"Adipose.Tissue"=(gct$Adipose...Subcutaneous+gct$Adipose...Visceral..Omentum.)/2,
                  "Nervous.System"=(gct$Brain...Cerebellar.Hemisphere+gct$Brain...Cerebellum+gct$Brain...Cortex+gct$Brain...Spinal.cord..cervical.c.1.)/4,
                  "Gut"=(gct$Colon...Sigmoid+gct$Colon...Transverse+gct$Small.Intestine...Terminal.Ileum+gct$Stomach)/4,
                  "Heart"=(gct$Heart...Atrial.Appendage+gct$Heart...Left.Ventricle)/2,
                  "Kidney"=(gct$Kidney...Cortex+gct$Kidney...Medulla)/2,
                  "Skin"=(gct$Skin...Not.Sun.Exposed..Suprapubic.+gct$Skin...Sun.Exposed..Lower.leg.)/2,
                  "Blood"=(gct$Spleen+gct$Whole.Blood)/2)

gct$Name <- str_remove(gct$Name, "[.].*")
colnames(gct)[1] <- "ENSEMBL"

ES_Placenta <- read.delim("Placenta_ES_sorted.txt")
ES_Placenta <- ES_Placenta[,c(1,6:15)]
ES_Placenta <- data.frame(ES_Placenta[,-c(2:11)],"Placenta"=apply(ES_Placenta[,7:11],1,median), "Pluripotent"=apply(ES_Placenta[,2:6],1,median))
colnames(ES_Placenta)[1] <- "ENSEMBL"

# calculate TPM for all samples and tissues:
teratomas_gene_exp <- gene_exp_all[,c(1,8:15)]
gtex_tissues_exp <- gct[,-2]
unified_exp <- merge(ES_Placenta,teratomas_gene_exp, by="ENSEMBL",all.y=T)
unified_exp <- merge(unified_exp,gtex_tissues_exp, by="ENSEMBL",all.x=T)
unified_exp <- merge(gene_lengths[,c(1,6,7)], unified_exp, by="ENSEMBL",all.y=T)
unified_exp <- na.omit(unified_exp)

for(i in names(unified_exp[,4:56])) {
  unified_exp[i] <- tpm_calc(unified_exp[i],unified_exp$length_total_exons) 
}
unified_exp <- unified_exp[-duplicated(unified_exp$SYMBOL),]

#divide to samples' TPM and GTEX tissues +ES+Placenta:
tissues_exp_profile <- unified_exp[,c(1:5,14:56)]
samples_TPM <- merge(gene_exp_all[,1:7],unified_exp[,c(1,6:13)], by="ENSEMBL",all.y=T)
write.csv(samples_TPM,'samples_TPM.csv')  # write TPM table

#removing all non protein coding genes:
genes <- tissues_exp_profile$ENSEMBL
library(biomaRt)
mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "ensembl.org")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart = mart,
  attributes = c(
    "hgnc_symbol",
    "entrezgene_id",
    "ensembl_gene_id",
    "gene_biotype"),
  filter = "ensembl_gene_id",
  values = genes,
  uniqueRows=TRUE)
annotLookup_2 <- annotLookup[annotLookup$gene_biotype=="protein_coding",]
tissues_exp_profile <- tissues_exp_profile[tissues_exp_profile$ENSEMBL %in% annotLookup_2$ensembl_gene_id,]
###################

colnames(tissues_exp_profile)[32] <- "Skeletal.Muscle"
tissues_exp_profile <- tissues_exp_profile[,c(1,3,6:48,4:5)]
tissues_exp_profile <- tissues_exp_profile[!duplicated(tissues_exp_profile$ENSEMBL),]
tissues_exp_profile <- tissues_exp_profile[!duplicated(tissues_exp_profile$SYMBOL),]
row.names(tissues_exp_profile) <- tissues_exp_profile$SYMBOL


## filtering the list ##
keep <- c()
for (i in row.names(tissues_exp_profile)){
  keep <- c(keep,(max(as.numeric(tissues_exp_profile[i,3:47]))/Rfast::nth(as.numeric(tissues_exp_profile[i,3:47]), 2, descending = T))>10)
}#keeping only genes that are expressed 10-fold or higher in their specific tissue than in all other tissues
tissues_exp_profile <- tissues_exp_profile[keep,]
tissues_exp_profile <- na.omit(tissues_exp_profile)

library(reshape2)
a <- melt(tissues_exp_profile[, -1], id.vars = 'SYMBOL')

## assign tissue specificity by assigning each gene the tissue at which it has the maximum expression levles:
tissues_exp_profile <- data.frame(tissues_exp_profile,"Tissue"=rep("",nrow(tissues_exp_profile)))
tissues <- c()
for (i in row.names(tissues_exp_profile)){
  tissues <- c(tissues,colnames(tissues_exp_profile)[which(as.numeric(tissues_exp_profile[i,3:47])==max(as.numeric(tissues_exp_profile[i,3:47])))+2] )
}
tissues_exp_profile$Tissue <- tissues


tissues_exp_profile <- tissues_exp_profile[,c(1:2,48,3:47)]
tissues_exp_profile <- tissues_exp_profile[order(tissues_exp_profile$Tissue),]
tissues_exp_profile <- tissues_exp_profile[apply(tissues_exp_profile[,-c(1:3)],1,max)>=5,] #filtering out genes which their highest expression level in their specific tissue is less than 5
all_tissue_specific_genes <- tissues_exp_profile[,1:3]
write.csv(all_tissue_specific_genes, "all_tissue_specific_genes.csv")
table(tissues_exp_profile$Tissue)

### taking only genes that are specific to our tissues of interest:
tissue_specific_genes <- all_tissue_specific_genes[all_tissue_specific_genes$Tissue=="Nervous.System" |
                                                     all_tissue_specific_genes$Tissue=="Skin" |
                                                     all_tissue_specific_genes$Tissue=="Blood" |
                                                     all_tissue_specific_genes$Tissue=="Kidney" |
                                                     all_tissue_specific_genes$Tissue=="Heart" |
                                                     all_tissue_specific_genes$Tissue=="Skeletal.Muscle" |
                                                     all_tissue_specific_genes$Tissue=="Gut" |
                                                     all_tissue_specific_genes$Tissue=="Liver" |
                                                     all_tissue_specific_genes$Tissue=="Lung" |
                                                     all_tissue_specific_genes$Tissue=="Pancreas" |
                                                     all_tissue_specific_genes$Tissue=="Placenta" |
                                                     all_tissue_specific_genes$Tissue=="Pluripotent",]

tissues_exp_profile_2 <- tissues_exp_profile[tissues_exp_profile$ENSEMBL %in% tissue_specific_genes$ENSEMBL,]
tissues <- tissue_specific_genes[!duplicated(tissue_specific_genes$Tissue),3]
tissues <- tissues[c(6,11,2,3,10,12,1,4,5,7,8,9)]
avg_lineage_gene_exp <- data.frame("dsRed-hygro 2n Ter1"=rep(0,12),"dsRed-hygro 2n Ter2"=rep(0,12),"EGFP-neo 2n Ter1"=rep(0,12),"EGFP-neo 2n Ter2"=rep(0,12),
                                   "clone H 3n Ter"=rep(0,12),"clone K 3n Ter"=rep(0,12),"clone L1 3n Ter"=rep(0,12),"clone L2 3n Ter"=rep(0,12))
row.names(avg_lineage_gene_exp)<- tissues


samples_TPM <- merge(samples_TPM[,c(1,7:15)], tissue_specific_genes[,-2], by = "ENSEMBL", all.y = T)
samples_TPM <- samples_TPM[order(samples_TPM$Tissue),]
samples_TPM <- samples_TPM[,c(1,2,11,3:10)]
samples_TPM <- na.omit(samples_TPM)
samples_TPM <- samples_TPM[!duplicated(samples_TPM$SYMBOL),]
row.names(samples_TPM) <- samples_TPM$SYMBOL

#keeping only genes that are expressed in diploid teratomas:
keep <- rowSums(samples_TPM[,4:7]>1) >= 1   #saving all genes expressed in at least 1 diploid teratoma
samples_TPM <- samples_TPM[keep,]

#dividing sapmles' TPM values by the average expression values of each gene in its specific tissue
samples_TPM_div <- samples_TPM
for (i in tissues){
  a<- row.names(samples_TPM[samples_TPM$Tissue==i,])
  samples_TPM_div[a,4:11] <- samples_TPM[a,4:11]/tissues_exp_profile_2[a,i]
}
colnames(samples_TPM_div)[4:11] <- colnames(avg_lineage_gene_exp)
for (i in colnames(avg_lineage_gene_exp)){
  for (j in tissues){
    avg_lineage_gene_exp[j,i] <- mean(samples_TPM_div[samples_TPM_div$Tissue==j,i])*100
  }
}

## average gene expression by ploidy and by lineage ##
group <- c("2n","2n","2n","2n","3n","3n","3n","3n")
samples_TPM_div_2 <- t(samples_TPM_div)
lineage <- samples_TPM_div_2[3,]
lineage <- rep(lineage,each=8)
samples_TPM_div_2 <- samples_TPM_div_2[-c(1:3),]
samples_TPM_div_2 <- data.frame(samples_TPM_div_2, "Ploidy" = group)
melted_gene_exp <- melt(data = samples_TPM_div_2, id.vars = "Ploidy", variable.name = "gene",value.name = "gene_expression",)
melted_gene_exp$gene_expression <- as.numeric(melted_gene_exp$gene_expression)*100 #make the values the percentage expression of each gene in its related tissue
melted_gene_exp <- data.frame(melted_gene_exp,"Lineage"=str_replace_all(lineage,"[.]"," "))
melted_gene_exp <- melted_gene_exp[order(melted_gene_exp$Lineage),]
melted_gene_exp$Order <- c(rep(8,56),rep(6,152),rep(11,64),rep(3,16),rep(10,296),rep(5,56),
                           rep(2,72),rep(9,120),rep(4,232),rep(1,248),rep(12,224),rep(7,184))
melted_gene_exp <- melted_gene_exp[order(melted_gene_exp$Order),]

ggboxplot(melted_gene_exp, x= "Lineage", y= "gene_expression",fill="Ploidy",
          outlier.shape = NA, palette = c("royalblue","red3"),fun="mean_sd1")+
  stat_compare_means(aes(group = Ploidy), label = "p.signif", method = "kruskal.test",label.y = 175, hide.ns = T,paired = F)+
  coord_cartesian(ylim=c(0,175))

stat=compare_means(gene_expression~Ploidy, data = melted_gene_exp, method="wilcox.test", p.adjust.method = "BH",group.by = "Lineage")
ggbarplot(melted_gene_exp, x= "Lineage", y= "gene_expression",fill="Ploidy",palette = c("royalblue","red3"),width = 0.7, position = position_dodge(0.8), add="mean_se")+
  stat_pvalue_manual(data = stat, x = "Lineage",size = 4.5, label = "p.format", y.position = 65,hide.ns = F)+
  theme_classic()+ylab("Relative Lineage-Specific Gene Expression (%)")+xlab(NULL)+
  theme(text = element_text(size = 13.5,face = "bold"), legend.position = "top")+
  scale_x_discrete(labels = c("Pluripotent","Nervous\nSystem","Kidney","Placenta","Lung","Gut","Skin","Blood","Pancreas","Liver","Heart","Skeletal\nMuscle"))


################ calculate TPM and add gene names ####################

# function to calculate TPM 
tpm_calc <- function(counts, lengths) {
  rate <- counts / lengths
  rate/sum(rate)*1*10^6 
}


# calculate TPM
tpm <- gene_exp_all
for(i in names(gene_exp_all[,8:19])) {
  tpm[i] <- tpm_calc(gene_exp_all[i],gene_exp_all$length_total_exons) 
}

#plot XIST expression levels in the teratomas
xist_exp <- data.frame("XIST"=t(tpm[tpm$SYMBOL=="XIST",8:19]),"Ploidy"=c("2n", "2n", "2n", "2n", "3n", "3n", "3n", "3n", "2n male", "2n male", "2n Xi ESCs", "2n Xi ESCs"))
stat=compare_means(XIST~Ploidy, data = xist_exp, method="wilcox.test", p.adjust.method = "BH")
ggbarplot(xist_exp, x= "Ploidy", y= "XIST",fill="Ploidy",palette = c("royalblue","red3","grey","darkblue"), position = position_dodge(0.8), add="mean_se")+
  ylab("XIST expressiom levels (TPM)")+geom_hline(yintercept = 1,linetype=2)

#plot XIST expression levels in the teratomas
xist_exp <- data.frame("XIST"=t(tpm[tpm$SYMBOL=="XIST",8:19]),"Ploidy"=c("2n", "2n", "2n", "2n", "3n", "3n", "3n", "3n", "2n male", "2n male", "2n Xi ESCs", "2n Xi ESCs"))
xist_exp$XIST[xist_exp$XIST==0] <- 0.01
xist_exp$XIST <- log2(xist_exp$XIST)
stat=compare_means(XIST~Ploidy, data = xist_exp, method="wilcox.test", p.adjust.method = "BH")
ggbarplot(xist_exp, x= "Ploidy", y= "XIST",fill="Ploidy",palette = c("royalblue","red3","grey","darkblue"), position = position_dodge(0.8), add="mean_se")+
  ylab(expression(Log[2](TPM)))+geom_hline(yintercept = 0,linetype=2)



################ nuclear area analysis 2n vs 3n teratomas ####################
setwd("C:/Users/owner/Desktop/nuclear size teratomas/area measurments")
pacinian_curpsule <- read.csv("pacinian_curpsule_nuclear_measurements.csv")
colnames(pacinian_curpsule)[1] <- "Nuclear_area"
#general nuclear volume in teratomas:
Nuclear_volume <- pacinian_curpsule$Nuclear_area/pi
Nuclear_volume <- Nuclear_volume^0.5
Nuclear_volume <- Nuclear_volume^3*pi*4/3          #in µm^3
nuclear_measurments <- data.frame("Measurments"=c(pacinian_curpsule$Nuclear_area,Nuclear_volume),"type"=c(rep("Nuclear_area",100),rep("Nuclear_volume",100)), "Ploidy"=rep(pacinian_curpsule$ploidy,2))

relative_nuclear_measurments <- nuclear_measurments
relative_nuclear_measurments$Measurments[relative_nuclear_measurments$type=="Nuclear_area"] <- relative_nuclear_measurments$Measurments[relative_nuclear_measurments$type=="Nuclear_area"]/mean(relative_nuclear_measurments$Measurments[relative_nuclear_measurments$type=="Nuclear_area" & relative_nuclear_measurments$Ploidy=="2n"])
relative_nuclear_measurments$Measurments[relative_nuclear_measurments$type=="Nuclear_volume"] <- relative_nuclear_measurments$Measurments[relative_nuclear_measurments$type=="Nuclear_volume"]/mean(relative_nuclear_measurments$Measurments[relative_nuclear_measurments$type=="Nuclear_volume" & relative_nuclear_measurments$Ploidy=="2n"])

stat=compare_means(Measurments~Ploidy, data = relative_nuclear_measurments, method="t.test", p.adjust.method = "BH",group.by = "type")
ggbarplot(relative_nuclear_measurments, x= "type", y= "Measurments",fill="Ploidy",palette = c("royalblue","red3"), position = position_dodge(0.8),lab.vjust = c(-1.6,-1.5,-2.5,-1.6),label = T,lab.nb.digits = 2, add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.format",x= "type", y.position = 1.7,tip.length = 0.01,hide.ns = F)+ylab("Relative size")

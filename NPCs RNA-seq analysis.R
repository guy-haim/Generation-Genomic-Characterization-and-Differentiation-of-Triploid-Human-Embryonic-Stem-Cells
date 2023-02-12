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

setwd("C:/Users/owner/Desktop/RNA-seq counts guy/NPCs")
# setwd("C:/Users/guyha/Desktop/RNA-seq counts guy/RNA-seq EBs and teratomas")

#loading RNA-seq data from each sample into a data.frame and getting rid of unnecessary data (first 4 rows and second and third columns)
NPCs_2n_rep1 <- read.delim("10G_p46_4c_NPCs_day_7_ReadsPerGene.out.tab", header = F, skip = 4)[,c(1,4)]
NPCs_2n_rep2 <- read.delim("10R_4c_p44_NPCs_day_7_ReadsPerGene.out.tab", header = F, skip = 4)[,c(1,4)]
colnames(NPCs_2n_rep1)[1] <- "Geneid"
colnames(NPCs_2n_rep2)[1] <- "Geneid"

NPCs_3n_clone_H <- read.delim("3n_f10_clone_H_NPCs_day_7_ReadsPerGene.out.tab", header = F, skip = 4)[,c(1,4)]
NPCs_3n_clone_B <- read.delim("3n_f8_clone_B_NPCs_day_7_ReadsPerGene.out.tab", header = F, skip = 4)[,c(1,4)]
colnames(NPCs_3n_clone_H)[1] <- "Geneid"
colnames(NPCs_3n_clone_B)[1] <- "Geneid"



###merging all RNA-seq data into one data.frame and adding gene names according to annotation file:
gene_exp_all <- Reduce(function(x, y) merge(x, y, all=T, by = "Geneid"), 
                       list(NPCs_2n_rep1,NPCs_2n_rep2,NPCs_3n_clone_H,NPCs_3n_clone_B))
colnames(gene_exp_all) <- c("ENSEMBL","NPCs_2n_rep1","NPCs_2n_rep2",
                            "NPCs_3n_clone_H","NPCs_3n_clone_B")
gene_lengths <- read.csv("gene_lengths.genecode.human.V38.csv")[,-1]

#removing the end of the ENSEMBL ID so that the same genes from different versions of the annotation file will be joined together  
gene_exp_all$ENSEMBL <- str_remove(gene_exp_all$ENSEMBL, "[.].*")
gene_lengths$ENSEMBL <- str_remove(gene_lengths$ENSEMBL, "[.].*")
gene_exp_all <- aggregate(gene_exp_all[-1], list(gene_exp_all$ENSEMBL), FUN = sum, na.rm = T)
colnames(gene_exp_all)[1] <- "ENSEMBL"


#merging all data frame to one count table with all the data needed:
gene_exp_all <- merge(gene_lengths, gene_exp_all, by = "ENSEMBL", all.y = T)
gene_exp_all <- gene_exp_all[rowSums(is.na(gene_exp_all[ ,2:12])) == 0, ] #removing empty rows
indRemoved <- which(apply(gene_exp_all[,8:11], 1, function(x) all(x == 0)) ) #saving all unexpressed genes in the data.frame
gene_exp_all <- gene_exp_all[-indRemoved,]      #getting rid of rows (genes) that all of their columns are zeros
.rowNamesDF(gene_exp_all, make.names=T) <- gene_exp_all$SYMBOL
indRemoved_Y <- which(gene_exp_all$chr=="chrY") #saving all Y chr genes in the data.frame
gene_exp_all <- gene_exp_all[-indRemoved_Y,]      #getting rid of rows (genes) that all of their columns are zeros
write.csv(gene_exp_all,'NPCs_read_counts_table.csv')  # write read counts table
X_genes <- gene_exp_all[gene_exp_all$chr=="chrX",7]

cl=c("royalblue1","royalblue3","red1","red3")
names(cl) <- str_replace_all(colnames(gene_exp_all[,8:11]),"_", " ")


#### Differential expression (DE) analysis of Teratomas (with male teratomas):

ploidy <- factor(c(1,1,2,2), labels = c("2n","3n")) #4 teratomas samples of recently diploids ("normal") and 4 teratomas samples of triploid cells
DE_object <- DGEList(gene_exp_all[,8:11], group = ploidy, genes = gene_exp_all$gene)
keep <- rowSums(cpm(DE_object)>2) >= 2
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

prin_comp <- prcomp(t(cpm_normalized[,8:11]), center = T, retx = T)
gg <- cbind.data.frame(prin_comp$x, "Sample" =colnames(gene_exp_all)[8:11])
gg <- data.frame(gg, "Ploidy"=ploidy)
autoplot(prin_comp, data = gg, colour = "Sample", shape="Ploidy", label =F,size=3)+ scale_colour_manual(values=cl)+theme_bw()+
  guides(shape = guide_legend(order = 1),colour = guide_legend(override.aes = list(shape=c(15))))

summary(prin_comp)

pheatmap(data.matrix(cpm_normalized[,8:11]), scale = "row", cluster_rows=F, show_rownames =F,color = colorRampPalette(c("darkblue", "white", "red3"))(100)) #TPM heatmap of X-linked expressed genes, without gene clustering

### X chr gene expression ###
cpm_normalized_X <- cpm_normalized_X[order(cpm_normalized_X$start, decreasing = F),]
pheatmap(data.matrix(cpm_normalized_X[,8:11]), scale = "row", cluster_rows=F, show_rownames =F,color = colorRampPalette(c("darkblue", "white", "red3"))(100)) #TPM heatmap of X-linked expressed genes, without gene clustering



### continue DE analysis ###
design <- model.matrix(~ploidy)
DE_TMM <- estimateDisp(DE_TMM, design, robust=TRUE)
# fit_TMM <- glmFit(DE_TMM, design)   #fit the model
fit_TMM <- glmQLFit(DE_TMM, design)   #fit the model (suited for low number of replicates)
# LRT_3n_Vs_2n <- glmLRT(fit_TMM)  #comparing triploids to diploids
LRT_3n_Vs_2n <- glmTreat(fit_TMM, lfc=log2(2))  #strict comparison of triploids and diploids
summary(decideTests(LRT_3n_Vs_2n))  #decide by FDR

TMM_3n <- topTags(LRT_3n_Vs_2n, n=nrow(DE_TMM))$table #assigning fold change gene expression between diploid and triploid ESCs
TMM_3n <- data.frame("SYMBOL"=row.names(TMM_3n), TMM_3n)
colnames(TMM_3n)[2] <- "logFC_triploids"
TMM_3n <- merge(gene_lengths, TMM_3n, by = "SYMBOL", all.y = T)
sig_DE_3n <-TMM_3n[TMM_3n$FDR<=0.05,] #assigning only genes with significant fold change between triploid and diploid ESCs (by FDR) 
sig_DE_3n <- sig_DE_3n[-rowSums(is.na(sig_DE_3n[ ,2:7])) == 0, ] #removing empty rows
sig_DE_3n<- sig_DE_3n[order(sig_DE_3n$logFC_triploids, decreasing = T),]
sig_DE_3n<- data.frame(sig_DE_3n[,1:8],"Fold.Change"=2^sig_DE_3n$logFC_triploids,sig_DE_3n[,9:12])

TMM_3n<- TMM_3n[order(TMM_3n$logFC_triploids, decreasing = T),]
TMM_3n<- data.frame(TMM_3n[,1:8],"Fold.Change"=2^TMM_3n$logFC_triploids,TMM_3n[,9:12])

writeLines(sig_DE_3n[1:300,1])
writeLines(sig_DE_3n[1825:2124,1])


####Perform GO analysis via GO Enrichment Analysis in THE GENE ONTOLOGY RESOURCE, STRING Database and GSEA-MSigDB website:
library(fgsea)
library(dplyr)

ranks <- TMM_3n
ranks[, 'score'] <-ranks$logFC_triploids* (-log(ranks$FDR))
ranks <- ranks[!duplicated(ranks$SYMBOL),]
.rowNamesDF(ranks, make.names=T) <- ranks$SYMBOL

ranks <- setNames(ranks$score,ranks$SYMBOL)
cell_types <- gmtPathways("c8.all.v2022.1.Hs.symbols.gmt") 
fgseaRes_c8 <- fgseaMultilevel(cell_types, ranks, minSize=20, maxSize=500,eps = 0)
fgseaRes_c8 <- na.omit(fgseaRes_c8) ##removing pathways which do not have an estimate of the error deviation of the calculated p-values (where log2err=NA)
fgseaRes_c8 <- fgseaRes_c8[fgseaRes_c8$padj < 0.05,]

collapsedPathways_c8 <- collapsePathways(fgseaRes_c8[order(padj)][fgseaRes_c8$padj < 0.05,], cell_types, ranks, nperm = 10000,pval.threshol = 0.01)
mainPathways_fdr_c8 <- fgseaRes_c8[pathway %in% collapsedPathways_c8$mainPathways][order(padj), pathway]
dev.new()
plotGseaTable(cell_types[mainPathways_fdr_c8], ranks, fgseaRes_c8, gseaParam = 0.5)

fgseaResTidy_c8 <- fgseaRes_c8 %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy_c8 %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

fgseaResTidy_c8 <-fgseaResTidy_c8[fgseaResTidy_c8$padj<0.05,]
fgseaResTidy_c8 <- fgseaResTidy_c8[fgseaResTidy_c8$pathway %in% mainPathways_fdr_c8,]
fgseaResTidy_c8 <- fgseaResTidy_c8[order(fgseaResTidy_c8$padj),]

ggplot(data = fgseaResTidy_c8[fgseaResTidy_c8$NES<0,], aes(reorder(pathway, -padj), -log2(padj))) +
  geom_col(aes(fill=NES)) +
  coord_flip() +
  labs(x="Cell Type", y=expression(-log[2](FDR))) + 
  theme_minimal()+ 
  geom_hline(yintercept = -log2(0.05), linetype=2,color = "red", size=1)


ggplot(data = fgseaResTidy_c8[fgseaResTidy_c8$NES>0,], aes(reorder(pathway, -padj), -log2(padj))) +
  geom_col(aes(fill=NES)) +
  coord_flip() +
  labs(x="Cell Type", y=expression(-log[2](FDR))) + 
  theme_minimal()+ 
  geom_hline(yintercept = -log2(0.05), linetype=2,color = "red", size=1)


##### Hallmark gene sets GSEA #####
hallmark_gene_sets <- gmtPathways("h.all.v2022.1.Hs.symbols.gmt")
fgseaRes_h <- fgseaMultilevel(hallmark_gene_sets, ranks, minSize=20, maxSize=500,eps = 0)
fgseaRes_h <- na.omit(fgseaRes_h) ##removing pathways which do not have an estimate of the error deviation of the calculated p-values (where log2err=NA)
fgseaRes_h <- fgseaRes_h[fgseaRes_h$padj < 0.05,]

collapsedPathways_h <- collapsePathways(fgseaRes_h[order(padj)][fgseaRes_h$padj < 0.05,], hallmark_gene_sets, ranks, nperm = 10000,pval.threshol = 0.01)
mainPathways_fdr_h <- fgseaRes_h[pathway %in% collapsedPathways_h$mainPathways][order(padj), pathway]
dev.new()
plotGseaTable(hallmark_gene_sets[mainPathways_fdr_h], ranks, fgseaRes_h, gseaParam = 0.5)

fgseaResTidy_h <- fgseaRes_h %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy_h %>% 
  dplyr::select(-leadingEdge, -ES) %>% 
  arrange(padj) %>% 
  DT::datatable()

fgseaResTidy_h <-fgseaResTidy_h[fgseaResTidy_h$padj<0.05,]
fgseaResTidy_h <- fgseaResTidy_h[fgseaResTidy_h$pathway %in% mainPathways_fdr_h,]
fgseaResTidy_h <- fgseaResTidy_h[order(fgseaResTidy_h$padj),]

ggplot(data = fgseaResTidy_h[fgseaResTidy_h$NES<0,], aes(reorder(pathway, -padj), -log2(padj))) +
  geom_col(aes(fill=NES)) +
  coord_flip() +
  labs(x="Pathway", y=expression(-log[2](FDR))) + 
  theme_minimal()+ 
  geom_hline(yintercept = -log2(0.05), linetype=2,color = "red", size=1)


ggplot(data = fgseaResTidy_h[fgseaResTidy_h$NES>0,], aes(reorder(pathway, -padj), -log2(padj))) +
  geom_col(aes(fill=NES)) +
  coord_flip() +
  labs(x="Pathway", y=expression(-log[2](FDR))) + 
  theme_minimal()+ 
  geom_hline(yintercept = -log2(0.05), linetype=2,color = "red", size=1)


plotGseaTable(c(hallmark_gene_sets[mainPathways_fdr_h[1:2]],cell_types[mainPathways_fdr_c8[5]]), ranks, rbind(fgseaRes_c8,fgseaRes_h), gseaParam = 0.5)





################  calculate TPM and add gene names ####################

# function to calculate TPM 
tpm_calc <- function(counts, lengths) {
  rate <- counts / lengths
  rate/sum(rate)*1*10^6 
}


# calculate TPM
tpm <- gene_exp_all
tpm <- na.omit(tpm)
for(i in names(tpm[,8:11])) {
  tpm[i] <- tpm_calc(tpm[i],tpm$length_total_exons) 
}

mean_2n_TPM <- apply(X = tpm[,8:9], MARGIN = 1, FUN = mean)
mean_3n_TPM <- apply(X = tpm[,10:11], MARGIN = 1, FUN = mean)
mean_TPM <- data.frame(tpm[1:7],mean_2n_TPM,mean_3n_TPM)
relative_mean_exp <- data.frame(mean_TPM[1:7],mean_TPM[,8:9]/mean_2n_TPM)
colnames(relative_mean_exp)[8:9] <- c("2n_relative_exp","3n_relative_exp")
relative_exp <- data.frame(tpm[1:7],tpm[,8:11]/mean_2n_TPM)
write.csv(relative_mean_exp,"relative_mean_exp.csv")
write.csv(mean_TPM,"mean_TPM.csv")
NPCs_markers <- c("NES","SOX1","PAX6","POU3F2","NEUROD4","NEUROG2","ASCL1","MEIS1","LMX1A")
npcs_markers_relative_mean_exp <- relative_mean_exp[relative_mean_exp$SYMBOL %in% NPCs_markers,]

npcs_markers_tpms <- tpm[tpm$SYMBOL %in% NPCs_markers,]
npcs_markers_tpms <- npcs_markers_tpms[,7:11]
npcs_markers_tpms[,2:5] <- apply(format(x = npcs_markers_tpms[,2:5], scientific = F,digits = 2),2,FUN = as.numeric)
keep <- rowSums(npcs_markers_tpms[2:3]) >= 2 | rowSums(npcs_markers_tpms[4:5]) >= 2 #keeping genes that are expressed [mean(tpm)>=1] in either diploid or triploid  
npcs_markers_tpms <- npcs_markers_tpms[keep,]

melted_NPCs_tpm <- melt(data = npcs_markers_tpms, id.vars = "SYMBOL")
melted_NPCs_tpm <- melted_NPCs_tpm[order(melted_NPCs_tpm$SYMBOL),]
melted_NPCs_tpm <- data.frame(melted_NPCs_tpm,"Ploidy"=c("2n","2n","3n","3n"))
melted_NPCs_tpm[,3] <- as.numeric(format(x = melted_NPCs_tpm[,3], scientific = F))

ggbarplot(melted_NPCs_tpm[melted_NPCs_tpm$SYMBOL %in% c("NES","POU3F2","SOX1"),], x= "SYMBOL", y= "value",fill="Ploidy",palette = c("royalblue","red3"),width = 0.7, position = position_dodge(0.8), add="mean")+
  theme_classic()+ylab("Neural markers Expression (TPM)")+xlab(NULL)+
  theme(text = element_text(size = 13.5,face = "bold"), legend.position = "top")

ggbarplot(melted_NPCs_tpm[melted_NPCs_tpm$SYMBOL %in% c("ASCL1","PAX6"),], x= "SYMBOL", y= "value",fill="Ploidy",palette = c("royalblue","red3"),width = 0.7, position = position_dodge(0.8), add="mean")+
  theme_classic()+ylab("Neural markers Expression (TPM)")+xlab(NULL)+
  theme(text = element_text(size = 13.5,face = "bold"), legend.position = "top")


melted_NPCs_tpm[melted_NPCs_tpm$value==0,3] <-0.1

for (i in names(mean_2n_TPM[names(mean_2n_TPM) %in% melted_NPCs_tpm$SYMBOL])){
  melted_NPCs_tpm[melted_NPCs_tpm$SYMBOL==i,3] <- melted_NPCs_tpm[melted_NPCs_tpm$SYMBOL==i,3]/mean_2n_TPM[names(mean_2n_TPM) %in% melted_NPCs_tpm$SYMBOL][i]
}
ggbarplot(melted_NPCs_tpm, x= "SYMBOL", y= "value",fill="Ploidy",palette = c("royalblue","red3"),width = 0.7, position = position_dodge(0.8), add="mean")+
  theme_classic()+ylab("Neural markers Expression (TPM)")+xlab(NULL)+
  theme(text = element_text(size = 13.5,face = "bold"), legend.position = "top")

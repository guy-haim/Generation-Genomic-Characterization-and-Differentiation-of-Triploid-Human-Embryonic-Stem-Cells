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
growth <- data.frame(growth,"Ploidy"=rep(c(rep("3n",9),rep("2n",9)),5))

stat=compare_means(Intensity~Ploidy, data = growth, method="t.test", p.adjust.method = "BH", group.by = "timepoint")
ggplot(growth, aes(x=timepoint, y=Intensity, group=Sample)) +
  stat_summary(fun = mean, geom='line', aes(color=Sample), lty=1, size=1.3) +
  stat_summary(fun=mean,geom='point') +
  stat_summary(fun.data=mean_cl_boot,geom='errorbar',width=0.2) +
  theme_bw()+scale_colour_manual(values=c("deepskyblue1","royalblue1","royalblue3","red1","red3","red4"))

relative_growth <- growth
day_1 <- relative_growth$Intensity[relative_growth$timepoint=="Day 1"]
Mean <- c(rep(mean(day_1[1:3]),3),rep(mean(day_1[4:6]),3),rep(mean(day_1[7:9]),3),rep(mean(day_1[10:12]),3),
          rep(mean(day_1[13:15]),3),rep(mean(day_1[16:18]),3))
relative_growth$Intensity <- relative_growth$Intensity/Mean

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

###### average ######
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

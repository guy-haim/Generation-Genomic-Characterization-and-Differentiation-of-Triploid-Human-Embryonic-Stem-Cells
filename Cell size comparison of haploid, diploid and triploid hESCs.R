### assigning diameter values (in µm) to cell line:
diameter_3n_H <- c(13.5, 12.5, 12.6, 12.7)                      
diameter_3n_B <- c(12.4, 11.7, 12.8, 12.0)                      
diameter_3n_I <- c(13.0, 12.6, 13.2, 13.4)                      

diameter_2n_EGFP <- c(11.0, 11.3, 11.5, 11.6)                   
diameter_2n_dsRed <- c(11.7, 11.3, 11.8, 11.6,
                       11.4, 11.2, 11.8, 11.9)                  
diameter_1n_dsRed <- c(10.6, 10.9, 10.1, 10.0, 9.8, 10.0,
                       10.7, 10.6, 10.1, 10.5, 10.1, 10.0)      

### combining diameter values according to ploidy:
diameter_3n <- c(diameter_3n_B, diameter_3n_H, diameter_3n_I)
diameter_3n_AVG <- mean(diameter_3n)
diameter_2n <- c(diameter_2n_EGFP, diameter_2n_dsRed)
diameter_2n_AVG <- mean(diameter_2n)
diameter_1n <- diameter_1n_dsRed
diameter_1n_AVG <- mean(diameter_1n)

#calculating surface area:
surface_Area_3n <- (diameter_3n/2)^2*pi*4      #in µm^2
surface_Area_2n <- (diameter_2n/2)^2*pi*4      #in µm^2
surface_Area_1n <- (diameter_1n/2)^2*pi*4      #in µm^2

#calculating volume:
volume_3n <- (diameter_3n/2)^3*pi*4/3          #in µm^3
volume_2n <- (diameter_2n/2)^3*pi*4/3          #in µm^3
volume_1n <- (diameter_1n/2)^3*pi*4/3          #in µm^3


#test normality of the data
shapiro.test(diameter_3n)
shapiro.test(diameter_2n)
shapiro.test(diameter_1n)


#organize the result in a data.frame
cell_size <- data.frame("diameter_in_µm" = c(diameter_1n, diameter_2n, diameter_3n), 
                        "surface_area_in_µm2" = c(surface_Area_1n, surface_Area_2n, surface_Area_3n), 
                        "volume_in_µm3" = c(volume_1n, volume_2n, volume_3n), 
                        "ploidy" = rep(c('1n','2n','3n'), c(12,12,12)))

#anova test of diameter, surface area and volume between ploidy states
aov_diameter <- aov(cell_size$diameter_in_µm~cell_size$ploidy)
summary(aov_diameter)[[1]][1,5] #get the anova p.val
aov_surface_area <- aov(cell_size$surface_area_in_µm2~cell_size$ploidy)
summary(aov_surface_area)[[1]][1,5] #get the anova p.val
aov_volume <- aov(cell_size$volume_in_µm3~cell_size$ploidy)
summary(aov_volume)[[1]][1,5] #get the anova p.val

#post hoc test
TukeyHSD(aov_diameter)
TukeyHSD(aov_surface_area)
TukeyHSD(aov_volume)
#multiple comparison and correction by FDR
library(ggstatsplot)
diameter_comparison <- pairwiseComparisons::pairwise_comparisons(cell_size, x=ploidy, y=diameter_in_µm, p.adjust.method = "fdr",)
surface_area_comparison <- pairwiseComparisons::pairwise_comparisons(cell_size, x=ploidy, y=surface_area_in_µm2, p.adjust.method = "fdr")
vol_comparison <- pairwiseComparisons::pairwise_comparisons(cell_size, x=ploidy, y=volume_in_µm3, p.adjust.method = "fdr")

#visualize cell size differences between ploidy states
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(ggsignif)
library(ggpubr)
stat=compare_means(diameter_in_µm~ploidy, data = cell_size, method="t.test", p.adjust.method = "BH",ref.group = "2n")
g1=ggbarplot(cell_size, x= "ploidy", y= "diameter_in_µm",fill="ploidy",palette = c("green4","royalblue","red4"),width = 0.5, position = position_dodge(0.5), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.adj", y.position = c(12,13.1),hide.ns = F)+
  ylab("Diameter (µm)")+theme_classic()+scale_y_continuous(breaks=seq(0,15,2))+xlab("Ploidy")+theme(legend.position = "none",text = element_text(size = 15))

stat=compare_means(surface_area_in_µm2~ploidy, data = cell_size, method="t.test", p.adjust.method = "BH",ref.group = "2n")
g2=ggbarplot(cell_size, x= "ploidy", y= "surface_area_in_µm2",fill="ploidy",palette = c("green4","royalblue","red4"), position = position_dodge(0.4), add="mean_se")+
  stat_pvalue_manual(data = stat, size=4.5, label = "p.adj", y.position = c(450,540),hide.ns = F)+
  ylab(expression(paste("Surface Area ",(µm^2))))+theme_classic()+scale_y_continuous(breaks=seq(0,1000,100))+xlab("Ploidy")+theme(legend.position = "none",text = element_text(size = 16))

stat=compare_means(volume_in_µm3~ploidy, data = cell_size, method="t.test", p.adjust.method = "BH",ref.group = "2n")
g3=ggbarplot(cell_size, x= "ploidy", y= "volume_in_µm3",fill="ploidy",palette = c("green4","royalblue","red4"), position = position_dodge(0.4), add="mean_se")+
  stat_pvalue_manual(data = stat, size=4.5, label = "p.adj", y.position = c(900,1160),hide.ns = F)+
  ylab(expression(paste("Volume ",(µm^3))))+theme_classic()+scale_y_continuous(breaks=seq(0,2000,250))+xlab("Ploidy")+theme(legend.position = "none",text = element_text(size = 16))

grid.arrange(g2, g3, ncol=2)


library(plyr)
library(stringr)

cell_size_2 <- cell_size
cell_size_2$ploidy <- as.numeric(str_replace_all(cell_size_2$ploidy,"n",""))

r1=ggplot(data = cell_size_2,mapping = aes(x= ploidy, y= volume_in_µm3,group=ploidy))
r1=r1+geom_dotplot(method = "histodot",binaxis = "y",stackdir = "center",dotsize = 0.5,fill=c(rep("forestgreen",12),rep("royalblue",12),rep("red3",12)))+
  geom_smooth(method = "lm", se = T,group="ploidy",color="black")+
  theme_classic()+
  xlab("Ploidy")+
  ylab(expression(paste("Volume ",(µm^3))))+
  scale_x_continuous(name=expression(Ploidy), breaks=c(1,2,3), labels=c("1n","2n","3n"), limits=c(0.9,3.1))
summary(lm(cell_size_2$volume_in_µm3 ~ cell_size_2$ploidy))
r1=r1+ annotate(geom="text", x=1.5, y=1000, label=expression("Volume = 310 + 253*n\n"),size=4.5)+
  annotate(geom="text", x=1.5, y=980, label="R = 0.92",size=4.5)+theme(legend.position = "none",text = element_text(size = 15))


cell_size <- data.frame(cell_size[,1:3],"SA_Vol"=cell_size$surface_area_in_µm2/cell_size$volume_in_µm3,"ploidy"=cell_size$ploidy)

ggplot(data = cell_size, aes(x= ploidy, y=SA_Vol/mean_SA_vol_2n))+stat_summary(fun = "mean",geom = "bar", fill=c("green4","royalblue","red4"))+
  geom_signif(test = "t.test",test.args = "less",textsize = 6, comparisons = list(c("2n", "1n"), c("3n","2n"),c("3n","1n")), map_signif_level = T,y_position = c(1.2,1.1,1.3))+
  geom_errorbar(stat='summary',fun.data=mean_se,width=.3)+ylab("Relative Surface Area / Volume")+theme_classic()



cell_size$relative_SV_ratio <- cell_size$SA_Vol/mean_SA_vol_2n
## plotting with corrected p-values by FDR (only significant) ##
stat=compare_means(relative_SV_ratio~ploidy, data = cell_size, method="t.test",  p.adjust.method = "BH",ref.group = "2n")
SV1=ggbarplot(cell_size, x= "ploidy", y= "relative_SV_ratio",fill="ploidy",palette = c("green4","royalblue","red4"),width = 0.5, position = position_dodge(0.6), add="mean_se")+
  stat_pvalue_manual(data = stat, label = "p.adj", y.position = c(1.17,1.05),hide.ns = F)+
  ylab("Relative Surface Area / Volume")+theme_classic()+scale_y_continuous(breaks=seq(0,2,0.1))+xlab("Ploidy")+theme(legend.position = "none",text = element_text(size = 15))
SV1

ggarrange(g1,r1,SV1,ncol = 3,nrow = 1, widths = c(1,2,1))

relative_cell_size <- data.frame("relative_diameter" = cell_size$diameter_in_µm/diameter_2n, 
                        "relative_surface_area" = cell_size$surface_area_in_µm2/surface_Area_2n, 
                        "relative_volume" = cell_size$volume_in_µm3/volume_2n, 
                        "ploidy" = rep(c('1n','2n','3n'), c(12,12,12)))

#visualize relative cell size differences between ploidy states:
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(ggsignif)

fun_mean <- function(x){
  return(round(x=(data.frame(y=mean(x),label=mean(x,na.rm=T))), digits=2))
  }
stat=compare_means(relative_diameter~ploidy, data = relative_cell_size, method="t.test", p.adjust.method = "BH",ref.group = "2n")
d=ggbarplot(relative_cell_size, x= "ploidy", y= "relative_diameter",fill="ploidy",palette = c("green4","royalblue","red4"),width = 0.6, position = position_dodge(0.4), add="mean_se")+
  stat_pvalue_manual(data = stat, size=4.5, label = "p.adj", y.position = c(1.12,1.22),hide.ns = F)+
  ylab("Relative Diameter")+theme_classic()+xlab("Ploidy")+theme(legend.position = "none",text = element_text(size = 16))+
  stat_summary(fun.data = fun_mean, geom="text", vjust=c(-1.05,-0.4,-0.9))

stat=compare_means(relative_surface_area~ploidy, data = relative_cell_size, method="t.test", p.adjust.method = "BH",ref.group = "2n")
s=ggbarplot(relative_cell_size, x= "ploidy", y= "relative_surface_area",fill="ploidy",palette = c("green4","royalblue","red4"),width = 0.6, position = position_dodge(0.4), add="mean_se")+
  stat_pvalue_manual(data = stat, size=4.5, label = "p.adj", y.position = c(1.15,1.38),hide.ns = F)+
  ylab("Relative Surface Area")+theme_classic()+xlab("Ploidy")+theme(legend.position = "none",text = element_text(size = 16))+
  stat_summary(fun.data = fun_mean, geom="text", vjust=c(-1.1,-0.5,-1))

stat=compare_means(relative_volume~ploidy, data = relative_cell_size, method="t.test", p.adjust.method = "BH",ref.group = "2n")
v=ggbarplot(relative_cell_size, x= "ploidy", y= "relative_volume",fill="ploidy",palette = c("green4","royalblue","red4"),width = 0.6, position = position_dodge(0.4), add="mean_se")+
  stat_pvalue_manual(data = stat, size=4.5, label = "p.adj", y.position = c(1.15,1.55),hide.ns = F)+
  ylab("Relative Volume")+theme_classic()+xlab("Ploidy")+theme(legend.position = "none",text = element_text(size = 16))+
  stat_summary(fun.data = fun_mean, geom="text", vjust=c(-1.05,-0.5,-1.3))

ggarrange(d,s,v, ncol=3,labels = "AUTO")
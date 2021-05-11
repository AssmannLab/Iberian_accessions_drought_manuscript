list.of.packages <- c("ggplot2", "readr","cowplot","magick","png","pdftools","grid","ggrepel","tidyverse","ggpubr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)library(tidyverse)

library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(Hmisc)
library(ggpmisc)
library(qvalue)
library(ggpubr)
library(cowplot)
library(magick)
library(png)
library(pdftools)
library(grid)
library(tidyr)

##Before running code, unzip data/annotation_GWAS.csv, data/GWAS/GWAS_Leaf_area_potential_plasticity_index_drought_LM.csv and data/GWAS/GWAS_SAUR26_eGWAS_LM.csv
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
##Fig. 2A (GO network): Gene Ontology (GO) biological process analysis was performed using the open-source software ShinyGO v0.61: http://bioinformatics.sdstate.edu/go60/ (Ge et al., 2020) with the following settings: the search species was “Arabidopsis thaliana,”  the P-value cutoff (FDR) was 0.05, and the number of most significant terms to show was 10. 
##To replicate these results, the lists of genes that were input into these analyses based on significance thresholds can be selected from the file "Annotated_GWAS_Leaf_area_potential_plasticity_index_drought_LM.csv" within the "Tables" folder. See below for the code and filtering steps applied to the raw output of GWAS analysis.
##GWAS analysis can be replicated using the online tool GWAPP (http://gwas.gmi.oeaw.ac.at/) was employed using a linear regression model (Seren et al., 2012). The phenotype that we input in the analysis (leaf area potential plasticity, labeled as "LAP_plasticity_index") can be retrieved from the "data" folder (phenotypes.csv) 
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
####Fig. 2A: Manhattan plot. GWAS on leaf area plasticity
dat1 = read.csv("data/GWAS/GWAS_Leaf_area_potential_plasticity_index_drought_LM.csv")

head(dat1)

p <- dat1$score
p=2*pnorm(-abs(p))
nullRatio <- pi0est(p)
nullRatioS <- pi0est(p, lambda=seq(0.40, 0.95, 0.05), smooth.log.pi0="TRUE")
nullRatioM <- pi0est(p, pi0.method="bootstrap")
qobj = qvalue(p, lambda=seq(0.05, 0.95, 0.1), smooth.log.pi0="TRUE")
qvalues <- qobj$qvalues
lfdr <- qobj$lfdr

#qobj_fdrlevel
dat1 <-cbind(dat1,qvalues,lfdr)


dat2 = read_csv("data/annotation_GWAS.csv")
head(dat2)
total <- merge(dat1,dat2,by=c("chr", "pos"))
head(total)
total$GVE.x <- NULL
total$score.y <- NULL
total$maf.y <- NULL
total$mac.y <- NULL
total$maf.y <- NULL
total$GVE.y <- NULL
total$QUAL <- NULL
total$FILTER <- NULL
total$ID <- NULL
total <-total[!(total$INFO == "downstream_gene_variant"), ]
colnames(total)[colnames(total)=="score.x"] <- "score"

colnames(total)[colnames(total)=="maf.x"] <- "maf"
colnames(total)[colnames(total)=="mac.x"] <- "mac"

total <- total[order(-total$score) ,  ]
total <- total[total$maf > 0.1, ]
head(total)
duplicated_position_column <-total$pos
total <-cbind(total,duplicated_position_column)

total<-total %>% unite(label, symbol, duplicated_position_column, sep = " pos. ")
total$X1 <- NULL
total$GVE <- NULL
colnames(total)[colnames(total)=="chr"] <- "CHR"
colnames(total)[colnames(total)=="pos"] <- "BP"

head(total)
###save annotated file for GWAS analysis
write.csv(total, "Tables/Annotated_GWAS_Leaf_area_potential_plasticity_index_drought_LM.csv", row.names = FALSE)
head(don)

don <- total %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(total, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

axisdf = don %>% group_by(CHR) %>% summarise(center=( max(BPcum) + min(BPcum) ) / 2 )

p2<-don %>% 
  
  ggplot( aes(x=BPcum, y=score)) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.4, size=3) +
  scale_color_manual(values = rep(c("#F8B195", "#F67280", "#C06C84","#6C5B7B","#355C7D"), 5 )) +

  #geom_hline(yintercept = 4, color = "#A9A9A9", linetype = "longdash")+
  # custom X axis:
  scale_x_continuous(labels=c("1" = "Chr. 1", "2" = "Chr. 2","3" = "Chr. 3","4" = "Chr. 4","5" = "Chr. 5"), breaks= axisdf$center )+
  #  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  #scale_x_continuous( label = "axisdf$CHR", breaks= axisdf$center ) +
  #scale_y_continuous(limits = c(0, 7.5))+
  #xlab(expression(bold("Cyclic nucleotide-gated channel 16 (FPKM)")))+ 
  theme(legend.position = "none")+
  # remove space between plot area and x axis
  
  # Custom the theme:
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm")) +
  theme(panel.background = element_rect(fill = 'white')) +
  ylab(expression(bold(atop(paste("Leaf area plasticity index:"), paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P")))))))+ 
  theme(axis.text.y = element_text(family = "Arial", color="black", size=25, face="bold")) +
  theme(axis.text.x = element_text(family = "Arial", color="black", size=25, face="bold")) +
  theme(axis.title.y = element_text(family = "Arial", color="black", size=25, face="bold")) +
  theme(axis.title.x=element_blank())+
theme(plot.title = element_text(family = "Arial", color="#696969",face="bold",  size=35)) +  
 ggtitle("GWAS")
p2
png("panels/p2.png", width = 14, height = 7, units = 'in', res = 350)
p2
dev.off() 

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#####Fig. 2B. Created on keynote. files and figure available in the "SAUR26_gene_model" folder
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

####Fig. 2B: Manhattan plot eGWAS for SAUR26

dat1 = read.csv("data/GWAS/GWAS_SAUR26_eGWAS_LM.csv")

head(dat1)

p <- dat1$score
p=2*pnorm(-abs(p))
nullRatio <- pi0est(p)
nullRatioS <- pi0est(p, lambda=seq(0.40, 0.95, 0.05), smooth.log.pi0="TRUE")
nullRatioM <- pi0est(p, pi0.method="bootstrap")
qobj = qvalue(p, lambda=seq(0.05, 0.95, 0.1), smooth.log.pi0="TRUE")
qvalues <- qobj$qvalues
lfdr <- qobj$lfdr

#qobj_fdrlevel
dat1 <-cbind(dat1,qvalues,lfdr)


dat2 = read_csv("data/annotation_GWAS.csv")
head(dat2)
total <- merge(dat1,dat2,by=c("chr", "pos"))
head(total)
total$GVE.x <- NULL
total$score.y <- NULL
total$maf.y <- NULL
total$mac.y <- NULL
total$maf.y <- NULL
total$GVE.y <- NULL
total$QUAL <- NULL
total$FILTER <- NULL
total$ID <- NULL
total <-total[!(total$INFO == "downstream_gene_variant"), ]
colnames(total)[colnames(total)=="score.x"] <- "score"

colnames(total)[colnames(total)=="maf.x"] <- "maf"
colnames(total)[colnames(total)=="mac.x"] <- "mac"

total <- total[order(-total$score) ,  ]
total <- total[total$maf > 0.1, ]
head(total)
duplicated_position_column <-total$pos
total <-cbind(total,duplicated_position_column)

total<-total %>% unite(label, symbol, duplicated_position_column, sep = " pos. ")
total$X1 <- NULL
total$GVE <- NULL
colnames(total)[colnames(total)=="chr"] <- "CHR"
colnames(total)[colnames(total)=="pos"] <- "BP"

head(total)

write.csv(total, "Tables/annotation_GWAS_SAUR26_eGWAS_LM.csv", row.names = FALSE)

don1 <- total %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(total, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

axisdf = don1 %>% group_by(CHR) %>% summarise(center=( max(BPcum) + min(BPcum) ) / 2 )

p3<-don1 %>% 
  
  ggplot( aes(x=BPcum, y=score)) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.4, size=3) +
  scale_color_manual(values = rep(c("#F8B195", "#F67280", "#C06C84","#6C5B7B","#355C7D"), 5 )) +
  
  #geom_hline(yintercept = 4, color = "#A9A9A9", linetype = "longdash")+
  # custom X axis:
  scale_x_continuous(labels=c("1" = "Chr. 1", "2" = "Chr. 2","3" = "Chr. 3","4" = "Chr. 4","5" = "Chr. 5"), breaks= axisdf$center )+
  #  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  #scale_x_continuous( label = "axisdf$CHR", breaks= axisdf$center ) +
 # scale_y_continuous(limits = c(0, 14))+
  #xlab(expression(bold("Cyclic nucleotide-gated channel 16 (FPKM)")))+ 
  theme(legend.position = "none")+
  # remove space between plot area and x axis
  
  # Custom the theme:
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm")) +
  theme(panel.background = element_rect(fill = 'white')) +
  ylab(expression(bold(atop(paste("SAUR26 expression:"), paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P")))))))+ 
  theme(axis.text.y = element_text(family = "Arial", color="black", size=25, face="bold")) +
  theme(axis.text.x = element_text(family = "Arial", color="black", size=25, face="bold")) +
  theme(axis.title.y = element_text(family = "Arial", color="black", size=25, face="bold")) +
  theme(axis.title.x=element_blank())+
theme(plot.title = element_text(family = "Arial", color="#696969",face="bold",  size=35)) +  
  ggtitle("eGWAS")
p3
png("panels/p3.png", width = 14, height = 7, units = 'in', res = 350)
p3
dev.off() 


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
####Fig. 2D

phenotypes = read.csv("data/phenotypes.csv")
head(phenotypes)

Phenotypesselected<- select(phenotypes, accession_id, SAUR26_haplotypes,LAP_plasticity_index)
head(Phenotypesselected)

my_comparisons <- c("Major", "Minor")
my_font_size <- 7
p4 <- ggplot(Phenotypesselected, aes(x=SAUR26_haplotypes, y=LAP_plasticity_index, fill=SAUR26_haplotypes)) + geom_boxplot(width=0.3,alpha=1) +  geom_violin(trim=FALSE,alpha=0.2) + 
  scale_fill_manual(values=c("#909ea4", "#D3D3D3")) + scale_x_discrete(limits=c("Major","Minor"), labels = c('Major haplotype','Minor haplotype'), expand = c(0.1, 0.15)) +
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.25,"cm")) +
  theme(axis.text.x = element_text(family = "Arial", color="black", face="bold", size=25)) +
  theme(axis.text.y = element_text(family = "Arial", color="black", face="bold", size=15)) +
  theme(axis.title.x = element_text(family = "Arial", color="black", face="bold", size=30)) +
  theme(axis.title.y = element_text(family = "Arial", color="black", face="bold", size=30)) +
  theme(plot.title = element_text(hjust = 0.5, vjust =  3, family = "Arial", color="#43464B", face="bold", size=25)) +
  theme(axis.title.x=element_blank())+ 
  theme(plot.margin=unit(c(1,1,1,1),"cm"))+
 # ggtitle("Our sample") +
  guides(fill = FALSE) + 
  ylab("Leaf area plasticity index")+  stat_compare_means(size = my_font_size,label.x.npc = "left",label.y.npc= "top",vjust=-3.2, hjust=0.5) 
p4$layers[[2]]$aes_params$textsize <- my_font_size
p4

p4
png("panels/p4.png", width = 14, height = 7, units = 'in', res = 350)
p4
dev.off() 


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
####Fig. 2Da. Transcriptome data is retrieved from the GEO dataset with accession number GSE80744 and SRA study SRP074107 (Kawakatsu et al., 2016)
phenotypes = read.csv("data/phenotypes.csv")
head(phenotypes)

Phenotypesselected<- select(phenotypes, accession_id, SAUR26_haplotypes)
head(Phenotypesselected)
transposedfinal = read_csv("data/transposedfinal.csv")
Transcriptselected<- select(transposedfinal, accession_id,AT3G03850)
head(Transcriptselected)
final<-merge(Phenotypesselected, Transcriptselected, by = "accession_id", sort = TRUE)
head(final)

my_comparisons <- c("Major", "Minor")
my_font_size <- 7
p5 <- ggplot(final, aes(x=SAUR26_haplotypes, y=AT3G03850, fill=SAUR26_haplotypes))  + geom_boxplot(width=0.20,alpha=1) +  geom_violin(trim=FALSE,alpha=0.2) + 
  scale_fill_manual(values=c("#909ea4", "#D3D3D3")) + scale_x_discrete(limits=c("Major","Minor"), labels = c('Major haplotype','Minor haplotype'), expand = c(0.1, 0.1)) +
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
    theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.25,"cm")) +
  theme(axis.text.x = element_text(family = "Arial", color="black", face="bold", size=25)) +
  theme(axis.text.y = element_text(family = "Arial", color="black", face="bold", size=15)) +
  theme(axis.title.x = element_text(family = "Arial", color="black", face="bold", size=30)) +
  theme(axis.title.y = element_text(family = "Arial", color="black", face="bold", size=30)) +
  theme(plot.title = element_text(hjust = 0.5, vjust =  3, family = "Arial", color="#43464B", face="bold", size=25)) +
  theme(axis.title.x=element_blank())+ 
  theme(plot.margin=unit(c(1,1,1,1),"cm"))+
 ggtitle("Our sample") +
  guides(fill = FALSE) + 
  ylab("SAUR26 expression (RPKM)")+  stat_compare_means(size = my_font_size,label.x.npc = "left",label.y.npc= "top",vjust=-8.5, hjust=0.5) 
p5$layers[[2]]$aes_params$textsize <- my_font_size
p5

p5
png("panels/p5.png", width = 8, height = 10, units = 'in', res = 350)
p5
dev.off() 

######################################################################################################################################
#######################################################################################################################################
####Fig. 2Eb. Transcriptome data is retrieved from the GEO dataset with accession number GSE80744 and SRA study SRP074107 (Kawakatsu et al., 2016)

phenotypes = read.csv("data/Iberian_SAUR26.csv")
head(phenotypes)

Phenotypesselected<- select(phenotypes, accession_id, SAUR26_haplotypes)
head(Phenotypesselected)
transposedfinal = read_csv("data/transposedfinal.csv")
Transcriptselected<- select(transposedfinal, accession_id,AT3G03850)
head(Transcriptselected)
final<-merge(Phenotypesselected, Transcriptselected, by = "accession_id", sort = TRUE)
head(final)

my_comparisons <- c("Major", "Minor")
my_font_size <- 7
p6 <- ggplot(final, aes(x=SAUR26_haplotypes, y=AT3G03850, fill=SAUR26_haplotypes)) + geom_boxplot(width=0.3,alpha=1) +  geom_violin(trim=FALSE,alpha=0.2) + 
  scale_fill_manual(values=c("#909ea4", "#D3D3D3"))  + scale_x_discrete(limits=c("Major","Minor"), labels = c('Major haplotype','Minor haplotype'), expand = c(0.1, 0.1)) +
  geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.25,"cm")) +
  theme(axis.text.x = element_text(family = "Arial", color="black", face="bold", size=25)) +
  theme(axis.text.y = element_text(family = "Arial", color="black", face="bold", size=15)) +
  theme(axis.title.x = element_text(family = "Arial", color="black", face="bold", size=30)) +
  theme(axis.title.y = element_text(family = "Arial", color="black", face="bold", size=30)) +
  theme(plot.title = element_text(hjust = 0.5, vjust =  3, family = "Arial", color="#43464B", face="bold", size=25)) +
  theme(axis.title.x=element_blank())+ 
  theme(plot.margin=unit(c(1,1,1,1),"cm"))+
  ggtitle("Iberian population") +
  guides(fill = FALSE) + 
  ylab("SAUR expression (RPKM)")+  stat_compare_means(size = my_font_size,label.x.npc = "left",label.y.npc= "top",vjust=-5, hjust=0.5) 
p6$layers[[2]]$aes_params$textsize <- my_font_size
p6

p6
png("panels/p6.png", width = 8, height = 10, units = 'in', res = 350)
p6
dev.off() 
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
####Fig. 2E. Transcriptome data is retrieved from the GEO dataset with accession number GSE80744 and SRA study SRP074107 (Kawakatsu et al., 2016)

df7=read_csv("data/phenotypes.csv")
head(df7)

haplotypes=read_csv("data/Iberian_SAUR26.csv")
head(haplotypes)

Phenotypesselected<- select(df7, accession_id, LAP_plasticity_index, name)
head(Phenotypesselected)
transposedfinal = read_csv("data/transposedfinal.csv")
Transcriptselected<- select(transposedfinal, accession_id,AT3G03850)
head(Transcriptselected)
final<-merge(Phenotypesselected, Transcriptselected, by = "accession_id", sort = TRUE)
final2<-merge(final, haplotypes, by = "accession_id", sort = TRUE)

head(final2)

###Plot1= altitude vs leaf area time potential, short days well-watered
formula <- y ~ x
p7 <- ggplot(final2, aes(x=AT3G03850 , y= LAP_plasticity_index)) +
  geom_point(aes(color=Haplotype),alpha = 0.9,  size=16)+   geom_text_repel(aes(AT3G03850, LAP_plasticity_index, label = name), size=6) +
  geom_smooth(method = "lm", se = T, fill="black", colour="black", size=0.5, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                       label.x.npc = "left", label.y.npc = "bottom",
                                                                                                       formula = formula, parse = TRUE, size = 7, colour="black") + 
  scale_color_manual(values=c('#909ea4','#D3D3D3'))+ theme_bw()# +

p7<- p7 + theme(axis.line = element_line(size=1, colour = "black"),
                panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm")) +  
  theme(legend.position = c(0.9, 0.17),legend.background = element_rect(fill = "white", color = "black")) + theme(legend.text = element_text(colour="black", size=14,face="bold"))  + theme(legend.title = element_text(colour="black", size=18,face="bold"))     
p7<- p7 + 
  theme(axis.text = element_text(family = "Arial", color="black", face="bold",size=20)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=30)) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(5, 5, 5,5)
  )+ ylim(-1,0.6)+
  ylab("Leaf area plasticity index") + 
  xlab("SAUR26 expression (FPKM)")
p7
png("panels/p7.png", width = 10, height = 7, units = 'in', res = 350)
p7
dev.off() 

#########################################################################################################
#########################################################################################################
###Fig. 2G cartoon was created in keynote and the configuration of that figure can be found in the "SAUR26 cartoon" folder
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################



M5 <- image_read("panels/p5.png")
M6 <- image_read("panels/p6.png")
plot_grid(rasterGrob(M5),rasterGrob(M6), greedy = TRUE, ncol = 2, nrow = 1, align = 'v', hjust=0, label_size=40)
png("panels/right.png", width = 10, height = 7, units = 'in', res = 350)
plot_grid(rasterGrob(M5),rasterGrob(M6), greedy = TRUE, ncol = 2, nrow = 1, align = 'v', hjust=0, label_size=40)
dev.off()

M2 <- image_read("panels/p2.png")
M3 <- image_read("panels/p3.png")
M4 <- image_read("panels/p4.png")
right <- image_read("panels/right.png")
M7 <- image_read("panels/p7.png")
M8 <- image_read_pdf("panels/p8.pdf")

plot_grid(rasterGrob(M2),rasterGrob(M3),rasterGrob(M4),rasterGrob(right),rasterGrob(M7),rasterGrob(M8), labels=c("", "C", "D","E","F","G"), greedy = TRUE, ncol = 2, nrow = 3, align = 'v', hjust=0, label_size=20)
png("panels/main.png", width = 11, height = 10, units = 'in', res = 350)
plot_grid(rasterGrob(M2),rasterGrob(M3),rasterGrob(M4),rasterGrob(right),rasterGrob(M7),rasterGrob(M8), labels=c("", "C", "D","E","F","G"), greedy = TRUE, ncol = 2, nrow = 3, align = 'v', hjust=0, label_size=20)
dev.off()

GO <- image_read("panels/GOfigure.png")
genemodel <- image_read_pdf("panels/p1.pdf")

plot_grid(rasterGrob(GO),rasterGrob(genemodel), labels=c("A", "B"), greedy = TRUE, ncol = 2, nrow = 1, align = 'v', hjust=0, label_size=105,rel_heights = c(2, 6))
png("panels/top.png", width = 60, height = 18, units = 'in', res = 350)
plot_grid(rasterGrob(GO),rasterGrob(genemodel), labels=c("A", "B"), greedy = TRUE, ncol = 2, nrow = 1, align = 'v', hjust=0, label_size=105,rel_heights = c(2, 6))
dev.off()


top <- image_read("panels/top.png")

main <- image_read("panels/main.png")

plot_grid(rasterGrob(top),rasterGrob(main), labels=c("", ""), greedy = TRUE, ncol = 1, nrow = 2, align = 'v', hjust=0, label_size=22,rel_heights = c(2, 6))
png("Figure 2.png", width = 12, height = 13, units = 'in', res = 350)
plot_grid(rasterGrob(top),rasterGrob(main), labels=c("", ""), greedy = TRUE, ncol = 1, nrow = 2, align = 'v', hjust=0, label_size=22,rel_heights = c(2, 6))
dev.off()


#########################################################################################################
#########################################################################################################
#########################################################################################################

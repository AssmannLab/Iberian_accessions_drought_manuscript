list.of.packages <- c("ggplot2", "readr","cowplot","magick","png","pdftools","grid","ggrepel","tidyverse","ggpubr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

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
#####Go to https://bioconductor.org/packages/release/bioc/html/qvalue.html for instructions on how to install the qvalue package.


dat2b <-dat2[!(dat2$INFO == "downstream_gene_variant"), ]
head(dat2b)
dat2b$maf <- NULL
dat2b$mac <- NULL
dat2b$GVE <- NULL
dat2b$ID <- NULL
dat2b$QUAL <- NULL
dat2b$FILTER <- NULL
nrow(dat2b)

dat3 = read_csv("data/gene_description.csv")
head(dat3)
total <- merge(dat2b,dat3,by="gene")
total$score <- NULL
total$description <- NULL

head(total)
final <-total[,c("chr","pos","REF","ALT","gene","annotation","symbol","INFO")]
colnames(final)[colnames(final)=="INFO"] <- "effect"

head(final)

write.csv(final, "data/annotation_GWAS.csv", row.names = FALSE)


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
df=read_csv("data/phenotypesb.csv")
head(df)
###Plot1= Nppspring vs leaf area plasticity in response to drought
formula <- y ~ x
p1 <- ggplot(df, aes(x=NPP , y= PI)) +
  geom_point(color='grey',alpha = 0.8,  size=12)+   geom_text_repel(aes(NPP, PI, label = accession_name, fontface="bold")) +
  geom_smooth(method = "lm", se = T, fill="black", colour="black", size=0.5, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                         label.x.npc = "right", label.y.npc = "top",
                                                                                                         formula = formula, parse = TRUE, size = 7, colour="black") +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.25,"cm")) +
  theme(axis.text.x = element_text(family = "Arial", color="black", face="bold", size=30)) +
  theme(axis.text.y = element_text(family = "Arial", color="black", face="bold", size=30)) +
  theme(axis.title.x = element_text(family = "Arial", color="black", face="bold", size=35)) +
  theme(axis.title.y = element_text(family = "Arial", color="black", face="bold", size=40)) +
  theme(plot.title = element_text(hjust = 0.5, vjust =  3, family = "Arial", color="#43464B", face="bold", size=25)) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))+
  theme(plot.title = element_text(family = "Arial", color="#696969",face="bold",  size=35)) +  
  theme(legend.position="none") + 
  ylab(expression(bold("Leaf area plasticity index (PI)")))+ 
  xlab(expression(bold(atop("MOD17A2 Net primary productivity spring", paste((gC/m^"2"/day))))))+
  ggtitle("Phenotype x environment ")+
  scale_y_continuous(limits = c(-0.75, 0.75), breaks = seq(-0.75, 0.75, by = 0.25))

p1
coef
linearMod <- lm(NPP ~ PI, data=df)  # build linear regression model on full data
print(linearMod)
summary(linearMod)
summary(linearMod)$r.squared
p1
png("panels/p1.png", width = 14, height = 10, units = 'in', res = 350)
p1
dev.off()
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
##GWAS analysis can be replicated using the online tool GWAPP (http://gwas.gmi.oeaw.ac.at/) was employed using a linear regression model (Seren et al., 2012). The phenotype that we input in the analysis (leaf area potential plasticity, labeled as "LAP_plasticity_index") can be retrieved from the "data" folder (phenotypes.csv) 
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
####Fig. 2A: Manhattan plot. GWAS on leaf area plasticity
dat1 = read.csv("data/GWAS/PI_LM.csv")

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
total$GVE <- NULL

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
write.csv(total, "Tables/Annotated_GWAS_PI_LM.csv", row.names = FALSE)
total<-total %>% unite(label, symbol, duplicated_position_column, sep = " pos. ")
head(total)

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
  mutate(quadrant = case_when(score >= 5.49   ~ "Q1")) %>% 
  mutate(label = case_when(score >= 5.49  ~  paste(label))) %>% 
  
  ggplot( aes(x=BPcum, y=score)) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.4, size=3) +
  scale_color_manual(values = rep(c("#F8B195", "#F67280", "#C06C84","#6C5B7B","#355C7D"), 5 )) +
  
  #  geom_text_repel(aes(BPcum, score, label = label), fontface = 'bold.italic', size=2, force=0.1,segment.size = 0.1,arrow = arrow(length = unit(0.01, 'npc')),family = 'Arial',
  #                 box.padding = unit(0.15, "lines"),
  #                point.padding = unit(0.45, "lines"),
  #               segment.color = 'grey50',max.overlaps = Inf ) +  
  
  geom_label_repel(aes(BPcum, score, label = label,fill = factor(quadrant)),
                   fontface = 'bold.italic', color = 'white', size =3,
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.5, "lines"),
                   segment.color = 'grey50',max.overlaps = Inf 
  ) +
  
  
  #geom_hline(yintercept = 4, color = "#A9A9A9", linetype = "longdash")+
  # custom X axis:
  scale_x_continuous(labels=c("1" = "Chr. 1", "2" = "Chr. 2","3" = "Chr. 3","4" = "Chr. 4","5" = "Chr. 5"), breaks= axisdf$center )+
  #  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  #scale_x_continuous( label = "axisdf$CHR", breaks= axisdf$center ) +
  # scale_y_continuous(limits = c(0, 7.5))+
  #xlab(expression(bold("Cyclic nucleotide-gated channel 16 (FPKM)")))+ 
  theme(legend.position = "none")+
  # remove space between plot area and x 
  
  # Custom the theme:
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm")) +
  theme(panel.background = element_rect(fill = 'white')) +
  ylab(expression(bold(atop(paste("Leaf area plasticity index (PI):"), paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P")))))))+ 
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
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

####Fig. 2B: Manhattan plot eGWAS for SAUR26

dat1 = read.csv("data/GWAS/eGWAS_SAUR26_LM.csv")
dat1<-dat1 %>% unite(label, symbol, duplicated_position_column, sep = " pos. ")
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
head(don1)
write.csv(total, "Tables/a.csv", row.names = FALSE)

p3<-don1 %>% 
  mutate(quadrant = case_when(score <= 39 & score >= 35  ~ "Q1")) %>% 
  mutate(label = case_when(score <= 39 & score >= 35  ~  paste(label))) %>%   
  
  ggplot( aes(x=BPcum, y=score)) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.4, size=3) +
  scale_color_manual(values = rep(c("#F8B195", "#F67280", "#C06C84","#6C5B7B","#355C7D"), 5 )) +
  
  #  geom_text_repel(aes(BPcum, score, label = label), fontface = 'bold.italic', size=2, force=0.1,segment.size = 0.1,arrow = arrow(length = unit(0.01, 'npc')),family = 'Arial',
  #                 box.padding = unit(0.15, "lines"),
  #                point.padding = unit(0.45, "lines"),
  #               segment.color = 'grey50',max.overlaps = Inf ) +  
  
  geom_label_repel(aes(BPcum, score, label = label,fill = factor(quadrant)),
                   fontface = 'bold.italic', color = 'white', size =2,
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.5, "lines"),
                   segment.color = 'grey50',max.overlaps = Inf 
  ) +
  
  
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

phenotypes = read.csv("data/phenotypesb.csv")
head(phenotypes)

Phenotypesselected<- select(phenotypes, accession_id, SAUR26_haplotypes,PI)
Phenotypesselected[complete.cases(Phenotypesselected), ]
head(Phenotypesselected)

my_comparisons <- c("Major", "Minor")
my_font_size <- 7
p4 <- ggplot(Phenotypesselected, aes(x=SAUR26_haplotypes, y=PI, fill=SAUR26_haplotypes)) + geom_boxplot(width=0.3,alpha=1) +  geom_violin(trim=FALSE,alpha=0.2) + 
  scale_fill_manual(values=c("#909ea4", "#D3D3D3")) + scale_x_discrete(limits=c("Major","Minor"), labels = c('Major haplotype','Minor haplotype'), expand = c(0.1, 0.15)) +
  geom_point( shape = 21,size=4, position = position_jitterdodge(), color="black",alpha=0.8)+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.25,"cm")) +
  theme(axis.text.x = element_text(family = "Arial", color="black", face="bold", size=30)) +
  theme(axis.text.y = element_text(family = "Arial", color="black", face="bold", size=30)) +
  theme(axis.title.x = element_text(family = "Arial", color="black", face="bold", size=35)) +
  theme(axis.title.y = element_text(family = "Arial", color="black", face="bold", size=40)) +
  theme(plot.title = element_text(hjust = 0.5, vjust =  3, family = "Arial", color="#43464B", face="bold", size=25)) +
  theme(axis.title.x=element_blank())+ 
  theme(plot.margin=unit(c(1,1,1,1),"cm"))+
 # ggtitle("Our sample") +
  guides(fill = "none") + 
  ylab("Leaf area plasticity index (PI)")+  stat_compare_means(size = my_font_size,label.x.npc = "left",label.y.npc= "top",vjust=-5.8, hjust=0.5) 
p4$layers[[2]]$aes_params$textsize <- my_font_size
p4

p4
png("panels/p4.png", width = 14, height = 10, units = 'in', res = 350)
p4
dev.off() 


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
####Fig. 2Da. Transcriptome data is retrieved from the GEO dataset with accession number GSE80744 and SRA study SRP074107 (Kawakatsu et al., 2016)
phenotypes = read.csv("data/phenotypesb.csv")
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
  geom_point( shape = 21,size=4, position = position_jitterdodge(), color="black",alpha=0.8)+
    theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.25,"cm")) +
  theme(axis.text.x = element_text(family = "Arial", color="black", face="bold", size=30)) +
  theme(axis.text.y = element_text(family = "Arial", color="black", face="bold", size=30)) +
  theme(axis.title.x = element_text(family = "Arial", color="black", face="bold", size=35)) +
  theme(axis.title.y = element_text(family = "Arial", color="black", face="bold", size=40)) +
  theme(plot.title = element_text(hjust = 0.5, vjust =  3, family = "Arial", color="#43464B", face="bold", size=25)) +
  theme(axis.title.x=element_blank())+ 
  theme(plot.margin=unit(c(1,1,1,1),"cm"))+
 ggtitle("Our sample") +
  guides(fill = "none") + 
  ylab("SAUR26 expression (RPKM)")+  stat_compare_means(size = my_font_size,label.x.npc = "left",label.y.npc= "top",vjust=-8.0, hjust=0.5) 
p5$layers[[2]]$aes_params$textsize <- my_font_size
p5

p5
png("panels/p5.png", width = 10, height = 11, units = 'in', res = 350)
p5
dev.off() 

######################################################################################################################################
#######################################################################################################################################
####Fig. 2Eb. Transcriptome data is retrieved from the GEO dataset with accession number GSE80744 and SRA study SRP074107 (Kawakatsu et al., 2016)

Phenotypesselected = read.csv("data/Iberian_SAUR26.csv")
head(phenotypes)

head(Phenotypesselected)
transposedfinal = read_csv("data/transposedfinal.csv")
Transcriptselected<- select(transposedfinal, accession_id,AT3G03850)
head(Transcriptselected)
final<-merge(Phenotypesselected, Transcriptselected, by = "accession_id", sort = TRUE)
head(final)

my_comparisons <- c("Major", "Minor")
my_font_size <- 7
p6 <- ggplot(final, aes(x=Haplotype, y=AT3G03850, fill=Haplotype)) + geom_boxplot(width=0.3,alpha=1) +  geom_violin(trim=FALSE,alpha=0.2) + 
  scale_fill_manual(values=c("#909ea4", "#D3D3D3"))  + scale_x_discrete(limits=c("Major","Minor"), labels = c('Major haplotype','Minor haplotype'), expand = c(0.1, 0.1)) +
  geom_point( shape = 21,size=4, position = position_jitterdodge(), color="black",alpha=0.8)+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.25,"cm")) +
  theme(axis.text.x = element_text(family = "Arial", color="black", face="bold", size=30)) +
  theme(axis.text.y = element_text(family = "Arial", color="black", face="bold", size=30)) +
  theme(axis.title.x = element_text(family = "Arial", color="black", face="bold", size=35)) +
  theme(axis.title.y = element_text(family = "Arial", color="black", face="bold", size=40)) +
  theme(plot.title = element_text(hjust = 0.5, vjust =  3, family = "Arial", color="#43464B", face="bold", size=25)) +
  theme(axis.title.x=element_blank())+ 
  theme(plot.margin=unit(c(1,1,1,1),"cm"))+
  ggtitle("Iberian population") +
  guides(fill = "none") + 
  ylab("SAUR expression (RPKM)")+  stat_compare_means(size = my_font_size,label.x.npc = "left",label.y.npc= "top",vjust=-6, hjust=0.5) 
p6$layers[[2]]$aes_params$textsize <- my_font_size
p6

p6
png("panels/p6.png", width = 10, height = 11, units = 'in', res = 350)
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

df7=read_csv("data/phenotypesb.csv")
head(df7)

head(df7)
transposedfinal = read_csv("data/transposedfinal.csv")
Transcriptselected<- select(transposedfinal, accession_id,AT3G03850)
head(Transcriptselected)
final<-merge(df7, Transcriptselected, by = "accession_id", sort = TRUE)
head(final)

###Plot1= altitude vs leaf area time potential, short days well-watered
formula <- y ~ x
p7 <- ggplot(final, aes(x=AT3G03850 , y= PI)) +
  geom_point(aes(color=SAUR26_haplotypes),alpha = 0.9,  size=16)+   geom_text_repel(aes(AT3G03850, PI, label = accession_name), size=6) +
  geom_smooth(method = "lm", se = T, fill="black", colour="black", size=0.5, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                       label.x.npc = "left", label.y.npc = "bottom",
                                                                                                       formula = formula, parse = TRUE, size = 10, colour="black") + 
  scale_color_manual(values=c('#909ea4','#D3D3D3'), name="Haplotypes")+ theme_bw() + 
  theme(legend.position = c(0.9, 0.17),legend.background = element_rect(fill = "white", color = "black")) + 
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.25,"cm")) +  theme(axis.text.x = element_text(family = "Arial", color="black", face="bold", size=25)) +
  theme(axis.text.x = element_text(family = "Arial", color="black", face="bold", size=30)) +
  theme(axis.text.y = element_text(family = "Arial", color="black", face="bold", size=30)) +
  theme(axis.title.x = element_text(family = "Arial", color="black", face="bold", size=35)) +
  theme(axis.title.y = element_text(family = "Arial", color="black", face="bold", size=40)) +
  theme(plot.title = element_text(hjust = 0.5, vjust =  3, family = "Arial", color="#43464B", face="bold", size=25)) +
  theme(plot.margin=unit(c(1,1,1,1),"cm"))+ ylim(-1,0.6)+
  ylab("Leaf area plasticity index (PI)") + xlab("SAUR26 expression (FPKM)")
p7
png("panels/p7.png", width = 14, height = 10, units = 'in', res = 350)
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
png("panels/right.png", width = 18, height = 10, units = 'in', res = 350)
plot_grid(rasterGrob(M5),rasterGrob(M6), greedy = TRUE, ncol = 2, nrow = 1, align = 'v', hjust=0, label_size=40)
dev.off()
M1 <- image_read("panels/p1.png")
genemodel <- image_read_pdf("panels/p1.pdf")
M2 <- image_read("panels/p2.png")
M3 <- image_read("panels/p3.png")
M4 <- image_read("panels/p4.png")
right <- image_read("panels/right.png")
M7 <- image_read("panels/p7.png")
M8 <- image_read_pdf("panels/p8.pdf")

plot_grid(rasterGrob(M1),rasterGrob(M2),rasterGrob(genemodel),rasterGrob(M3),rasterGrob(M4),rasterGrob(right),rasterGrob(M7),rasterGrob(M8), labels=c("a", "b", "c","d","e","f","g","h"), greedy = TRUE, ncol = 2, nrow = 4, align = 'v', hjust=0, label_size=20)
png("Figure 2.png", width = 9, height = 10, units = 'in', res = 350)
plot_grid(rasterGrob(M1),rasterGrob(M2),rasterGrob(genemodel),rasterGrob(M3),rasterGrob(M4),rasterGrob(right),rasterGrob(M7),rasterGrob(M8), labels=c("a", "b", "c","d","e","f","g","h"), greedy = TRUE, ncol = 2, nrow = 4, align = 'v', hjust=0, label_size=20)
dev.off()







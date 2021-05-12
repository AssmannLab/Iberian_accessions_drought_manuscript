#READ TWO FILES
library(tidyverse)
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(Hmisc)
library(ggpmisc)
library(qvalue)
library(ggpubr)
library(cowplot)
library(patchwork)
library(magick)
library(png)
library(grid)

###################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################



####Create an annotation file compatible with the files we are creating.
annotation = read_csv("data/annotation.csv")
annotation <-annotation[!(annotation$INFO == "downstream_gene_variant"), ]
annotation$gene <- gsub("(.*)..*", "\\1", annotation$gene)
annotation$gene <-gsub("\\.","",annotation$gene)
nrow(annotation)
head(annotation)
gene_description = read_csv("data/gene_description.csv")
head(gene_description)
annotation <- merge(annotation,gene_description,by= "gene")
head(annotation)



df71a=read_csv("data/GWAS/GWAS_d13C_Spain_LM.csv")
df71a <- df71a[df71a$maf > 0.1, ]
df71a <- merge(df71a,annotation,by=c("chr", "pos"))
head(df71a)

df1b=read_csv("data/GWAS/GWAS_d13C_Our sample_Spain_LM.csv")
df1b <- df1b[df1b$maf > 0.1, ]
df1b <- merge(df1b,annotation,by=c("chr", "pos"))
head(df1b)

head(df1)
df1 <- merge(df71a,df1b,by=c("chr", "pos"))
head(df1)
nrow(df1)

df1<-df1 %>% unite(annotation, symbol.y, pos, sep = " pos. ")
df1 <- df1[order(df1$score.x) ,  ]
head(df1)


df1$GVE.x <- NULL
df1$gene.x <- NULL
df1$INFO.x <- NULL
df1$annotation.x <- NULL
df1$Gene_id.x <- NULL
df1$GVE.y <- NULL
head(df1)

write.csv(df1, "/Tables/GWAS_d13C.csv", row.names = FALSE)
df1 <- df1[order(-df1$score.x) ,  ]
head(df1)
df1 <- df1 %>% 
  rename(
    Iberian_score.d13C = score.x,
    Our_sample_score.d13C = score.y)
head(df1)
formula <- y ~ x
p1<-df1 %>% 
  mutate(quadrant = case_when(Iberian_score.d13C > 4 & Our_sample_score.d13C > 4   ~ "Q1",
                              Iberian_score.d13C <= 4 & Our_sample_score.d13C > 4  ~ "Q2",
                              Iberian_score.d13C <= 4 & Our_sample_score.d13C <= 4 ~ "Q3",
                              TRUE                                         ~ "Q4")) %>% 
  mutate(label = case_when(Iberian_score.d13C >= 4 & Our_sample_score.d13C >= 4 ~  paste(annotation))) %>% 
  ggplot(aes(x = Iberian_score.d13C, y = Our_sample_score.d13C)) + 
  geom_point(aes(color=quadrant),alpha = 0.5,  size=3) + 
  geom_smooth(method = "lm", se = T, fill="black", colour="#696969", size=0.7, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                           label.x.npc = "left", label.y.npc = "top",
                                                                                                           formula = formula, parse = TRUE, size = 5, colour="#3792cb") +
  scale_colour_manual(values = c( "#47a569", "#6689c3", "#c5c5c5","#e9af41"))+
  geom_text_repel(aes(Iberian_score.d13C, Our_sample_score.d13C, label = label), fontface = 'bold', size=1, force=0.1,segment.size = 0.1,arrow = arrow(length = unit(0.01, 'npc')),family = 'Arial',
                  fontface = 'bold',
                  box.padding = unit(0.15, "lines"),
                  point.padding = unit(0.45, "lines"),
                  segment.color = 'grey50') + 
  geom_vline(xintercept = 4, color = "black", linetype = "dashed")+
  geom_hline(yintercept = 4, color = "black", linetype = "dashed")+
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=17)) +  
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm"))  +
  theme(axis.text = element_text(family = "Arial", color="black",  face="bold",size=16)) +
  theme(plot.title = element_text(family = "Arial", color="#696969", face="bold",size=28)) +  theme(panel.background = element_rect(fill = 'white')) +
  ggtitle("GWAS")+
  theme(legend.position="none") + 
  ylab(expression(bold(atop(paste("Our sample: Iberian intrinsic water use efficiency"," (\u03B4"^"13", "C)"), paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P")))))))+ 
  xlab(expression(bold(atop(paste("Iberian intrinsic water use efficiency"," (\u03B4"^"13", "C)"), paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P")))))))+ 
  scale_x_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2))+
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2))

p1
png("p1.png", width = 7, height = 7, units = 'in', res = 350)
p1
dev.off() 

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

###Open dataframe with phenotypes
####Create an annotation file compatible with the files we are creating.
annotation = read_csv("data/gene_description.csv")
head(annotation)

###Open dataframe with phenotypes



##res<-rcorr(as.matrix(phenotypes))
##AA<-signif(res$r, 2)
##write.csv(AA, "/Users/angel_admin/Box/Assmann_Lab/2_Researcher_Data/1_Current/Angel_Ferrero-Serrano/Spanish accessions/Fig28_transcriptome_potentialvsstability/Iberian_PxExP_corr.csv")


#####Run autocorrelation matrix Phenotype A
phenotypes = read_csv("data/Iberian_d13C.csv")
head(phenotypes)
phenotypesselected<- select(phenotypes, Dittberner_moleco18_Delta_13C)
accession_id1<-phenotypes$accession_id
selected <-cbind(accession_id1, phenotypesselected)
head(selected)
transposedfinal = read_csv("data/transposedfinal.csv")
head(transposedfinal)
colnames(transposedfinal)[colnames(transposedfinal)=="accession_id"] <- "accession_id1"
final<-merge(selected, transposedfinal, by = "accession_id1", sort = TRUE)
head(final)
accession_idsorted<-final$accession_id1
final$accession_id1 <- NULL
res<-rcorr(as.matrix(final), type = c("pearson","spearman"))
AA<-signif(res$r, 2)
BB<-signif(res$P,2)
outputB1 <-cbind(accession_idsorted,AA)
outputB1 <- outputB1[-1, 1:2]
outputB1 <- as.data.frame(outputB1)
outputB1 <- outputB1 %>% 
  rownames_to_column(var = "gene")
outputB2 <-cbind(accession_idsorted,BB)
outputB2 <- outputB2[-1, 1:2]
outputB2 <- as.data.frame(outputB2)
outputB2 <- outputB2 %>% 
  rownames_to_column(var = "gene")
outputB <- merge(outputB1,outputB2,by="gene")
outputB <- outputB[order(-outputB$Dittberner_moleco18_Delta_13C.x) ,  ]
head(outputB)
outputB <- merge(outputB,annotation,by="gene")
outputB <- outputB %>% 
  rename(
    Dittberner_moleco18_Delta_13C_rs = Dittberner_moleco18_Delta_13C.x,
    Dittberner_moleco18_Delta_13C_P_val = Dittberner_moleco18_Delta_13C.y)
outputB$accession_idsorted.x <- NULL
outputB$accession_idsorted.y <- NULL
outputB <- outputB[order(-outputB$Dittberner_moleco18_Delta_13C_rs) ,  ]
head(outputB, n=100)
tail(outputB, n=50)
head(outputB)
write.csv(outputB, "Tables/TWAS_Iberian_Dittberner_moleco18_Delta_13C.csv", row.names = FALSE)

head(outputB)


#####Run autocorrelation matrix condition B

###Open dataframe with phenotypes
Our_phenotypes = read.csv("data/Our_sample_phenotypes.csv")
head(Our_phenotypes)
head(Our_phenotypes)
phenotypesselected<- select(Our_phenotypes, Dittberner_moleco18_Delta_13C)
accession_id1<-Our_phenotypes$accession_id
selected <-cbind(accession_id1, phenotypesselected)
head(selected)
transposedfinal = read_csv("data/transposedfinal.csv")
head(transposedfinal)
colnames(transposedfinal)[colnames(transposedfinal)=="accession_id"] <- "accession_id1"
final<-merge(selected, transposedfinal, by = "accession_id1", sort = TRUE)
head(final)
accession_idsorted<-final$accession_id1
final$accession_id1 <- NULL
res<-rcorr(as.matrix(final), type = c("pearson","spearman"))
AA<-signif(res$r, 2)
BB<-signif(res$P,2)
outputC1 <-cbind(accession_idsorted,AA)
outputC1 <- outputC1[-1, 1:2]
outputC1 <- as.data.frame(outputC1)
outputC1 <- outputC1 %>% 
  rownames_to_column(var = "gene")
outputC2 <-cbind(accession_idsorted,BB)
outputC2 <- outputC2[-1, 1:2]
outputC2 <- as.data.frame(outputC2)
outputC2 <- outputC2 %>% 
  rownames_to_column(var = "gene")
outputC <- merge(outputC1,outputC2,by="gene")
outputC <- outputC[order(-outputC$Dittberner_moleco18_Delta_13C.x) ,  ]
head(outputC)
outputC <- merge(outputC,annotation,by="gene")
outputC <- outputC %>% 
  rename(
    Our_sample_Dittberner_moleco18_Delta_13C_rs = Dittberner_moleco18_Delta_13C.x,
    Our_sample_Dittberner_moleco18_Delta_13C_P_val = Dittberner_moleco18_Delta_13C.y)
outputC$accession_idsorted.x <- NULL
outputC$accession_idsorted.y <- NULL
outputC <- outputC[order(-outputC$Our_sample_Dittberner_moleco18_Delta_13C_rs) ,  ]
head(outputC, n=100)
tail(outputC, n=50)
head(outputC)
write.csv(outputC, "Tables/TWAS_OUR_Iberian_Dittberner_moleco18_Delta_13C.csv", row.names = FALSE)

head(outputC)


#####Merge df3

df2 <- merge(outputB,outputC,by="gene")
nrow(df2)
head(df2)
df8[df2 == "NAN"] <- NA
df8[complete.cases(df2),]

df2$accession_idsorted.y <- NULL
df2$annotation.x <- NULL
df2$accession_idsorted.x <- NULL
df2$Gene_id.x <- NULL
df2$symbol.x <- NULL

head(df2)
df2 <- df2 %>% 
  rename(
    symbol = symbol.y,
    annotation = annotation.y)
head(df2)

formula <- y ~ x
p2<-df2 %>% 
  mutate(quadrant = case_when(Dittberner_moleco18_Delta_13C_rs <= -0.4 & Our_sample_Dittberner_moleco18_Delta_13C_rs >= 0.4 ~ "Q1",
                              Dittberner_moleco18_Delta_13C_rs <= -0.4 & Our_sample_Dittberner_moleco18_Delta_13C_rs <= -0.4 ~ "Q2",
                              Dittberner_moleco18_Delta_13C_rs >= 0.4 & Our_sample_Dittberner_moleco18_Delta_13C_rs >= 0.4 ~ "Q3",
                              Dittberner_moleco18_Delta_13C_rs >= 0.4 & Our_sample_Dittberner_moleco18_Delta_13C_rs <= -0.4 ~ "Q4",
                              TRUE                                         ~ "Q5")) %>% 
  mutate(label = case_when(Dittberner_moleco18_Delta_13C_rs <= -0.4 & Our_sample_Dittberner_moleco18_Delta_13C_rs >= 0.4 ~  paste(symbol),
                           Dittberner_moleco18_Delta_13C_rs <= -0.4 & Our_sample_Dittberner_moleco18_Delta_13C_rs <= -0.4 ~  paste(symbol),
                           Dittberner_moleco18_Delta_13C_rs >= 0.4 & Our_sample_Dittberner_moleco18_Delta_13C_rs >= 0.4 ~  paste(symbol),
                           Dittberner_moleco18_Delta_13C_rs >= 0.4 & Our_sample_Dittberner_moleco18_Delta_13C_rs <= -0.4 ~  paste(symbol))) %>% 
  
  mutate(shade = case_when(Dittberner_moleco18_Delta_13C_rs <= -0.4 & Our_sample_Dittberner_moleco18_Delta_13C_rs >= 0.4 ~  "S1",
                           Dittberner_moleco18_Delta_13C_rs <= -0.4 & Our_sample_Dittberner_moleco18_Delta_13C_rs <= -0.4 ~  "S2",
                           Dittberner_moleco18_Delta_13C_rs >= 0.4 & Our_sample_Dittberner_moleco18_Delta_13C_rs >= 0.4 ~  "S3",
                           Dittberner_moleco18_Delta_13C_rs >= 0.4 & Our_sample_Dittberner_moleco18_Delta_13C_rs <= -0.4 ~  "S4")) %>% 
  
  ggplot(aes(x = Dittberner_moleco18_Delta_13C_rs, y = Our_sample_Dittberner_moleco18_Delta_13C_rs)) + 
  geom_point(aes(color=quadrant,alpha = shade) , size=12, shape = "-") +   scale_alpha_discrete(0.3,0.9,0.3,0.3,0.3)+
  
  geom_smooth(method = "lm", se = T, fill="black", colour="#696969", size=0.7, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                             label.x.npc = "left", label.y.npc = "top",
                                                                                                             formula = formula, parse = TRUE, size = 5, colour="#3792cb") +
  scale_colour_manual(values = c("#ce494a","#47a569", "#c5c5c5"))+  
  geom_text_repel(aes(Dittberner_moleco18_Delta_13C_rs, Our_sample_Dittberner_moleco18_Delta_13C_rs, label = label),  size=2, force=0.4,segment.size = 0.5,arrow = arrow(length = unit(0.01, 'npc')),family = 'Arial',
                  fontface = 'bold',
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.5, "lines"),
                  segment.color = 'grey50',max.overlaps = Inf) + 
  geom_vline(xintercept = 0.4, color = "black", linetype = "dashed")+
  geom_vline(xintercept = -0.4, color = "black", linetype = "dashed")+
  geom_hline(yintercept = 0.4, color = "black", linetype = "dashed")+
  geom_hline(yintercept = -0.4, color = "black", linetype = "dashed")+
  
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm"))  +
  theme(axis.text = element_text(family = "Arial", color="black", size=16)) +
  theme(axis.title = element_text(family = "Arial", color="black", size=17)) +
  theme(legend.position="none") + 
  theme(plot.title = element_text(family = "Arial", color="#696969", face="bold",size=28)) +  theme(panel.background = element_rect(fill = 'white')) +
  ylab(expression(bold(atop(paste("Our sample: Iberian intrinsicwater use efficiency (\u03B4"^"13", "C)"), paste("transcriptome-wide association ", (r[s]))))))+ 
  xlab(expression(bold(atop(paste("Iberian intrinsic water use efficiency (\u03B4"^"13", "C)"), paste("transcriptome-wide association ", (r[s]))))))+ 
   ggtitle("TWAS")+
  scale_y_continuous(limits = c(-0.8, 0.8), breaks = seq(-0.8, 0.8, by = 0.4))+  scale_x_continuous(limits = c(-0.8, 0.8), breaks = seq(-0.8, 0.8, by = 0.4))

p2 
png("p2.png", width = 7, height = 7, units = 'in', res = 350)
p2
dev.off() 




####################################################################################################################################################
####################


M1 <- image_read("p1.png")
M2 <- image_read("p2.png")


plot_grid(rasterGrob(M1),rasterGrob(M2),labels=c('A', 'B'), greedy = TRUE, ncol =2, nrow = 1, align = 'v', hjust=0, label_size=45)
png("Fig.S2345.png", width = 14, height = 6, units = 'in', res = 350)
plot_grid(rasterGrob(M1),rasterGrob(M2),labels=c('A', 'B'), greedy = TRUE, ncol =2, nrow = 1, align = 'v', hjust=0, label_size=45)
dev.off()

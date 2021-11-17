list.of.packages <- c("ggplot2", "viridis", "ggrepel","readr", "tidyverse","tidyr","cowplot","qvalue","dplyr","viridis","ggpubr","ggthemes", "Hmisc", "cowplot","magick","grid", "ggpmisc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(readr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(qvalue)
library(Hmisc)
library(cowplot)
library(dplyr)
library(ggpmisc)
library(magick)
library(grid)

####Plot1
####################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################

df1=read_csv("data/phenotypesb.csv")
head(df1)

####Plot1
meanx <- mean(df1$iWUE,na.rm = TRUE)
print(meanx)
meany <- mean(df1$PI,na.rm = TRUE)
print(meany)
formula <- y ~ x
p1<-df1 %>% 
  mutate(quadrant = case_when(iWUE > -35.54932923 & PI > 0   ~ "Q1",
                              iWUE <= -35.54932923 & PI > 0  ~ "Q2",
                              iWUE <= -35.54932923 & PI <= 0 ~ "Q3",
                              TRUE                                         ~ "Q4")) %>% 
  ggplot(aes(x = iWUE, y = PI)) + 
  geom_point(aes(color=quadrant),alpha = 0.7,  size=9) +  
  geom_smooth(method = "lm", se = T, fill="black", colour="black", size=0.2, alpha = 0.1) +  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")),
                                                                                                          label.x.npc = "right", label.y.npc = "top",
                                                                                                          formula = formula, parse = TRUE, size = 5, colour="black")+
  scale_colour_manual(values = c("#47a569", "#6689c3", "#ce494a", "#e9af41"))+
  geom_text_repel(aes(iWUE, PI, label = accession_name),  size=4, force=0.1,segment.size = 0.5,arrow = arrow(length = unit(0.01, 'npc')),family = 'Arial',
                  fontface = 'bold',colour = "#696969",
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.5, "lines"),
                  segment.color = 'grey50',max.overlaps = Inf) + 
  geom_vline(xintercept = -35.57922, color = "black", linetype = "dashed")+
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm"))  +
  theme(plot.title = element_text(family = "Arial", color="#696969",face="bold",  size=35))  +  theme(panel.background = element_rect(fill = 'white')) +
   ggtitle("Phenotypes")+
  xlab(expression(bold(paste("WUE"[i]," (\u03B4"^"13", "C (\u2030))"))))+  
  ylab("Leaf area plasticity index")+
  theme(axis.text = element_text(family = "Arial", color="black", size=16)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=18)) +  
  theme(legend.position="none") 
p1

png("panels/p1.png", width = 7, height = 7, units = 'in', res = 350)
p1
dev.off() 
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###Plot4= leaf area plasticity index vs d13C GWAS
###################################################################################################################################################
###################################################################################################################################################

####Create an annotation file compatible with the files we are creating.
annotation = read_csv("data/annotation_GWAS.csv")
df2a=read_csv("data/GWAS/iWUE_LM.csv")
df2a <- df2a[df2a$maf > 0.1, ]
df2a <- merge(df2a,annotation,by=c("chr", "pos"))
head(df2a)

df2b=read_csv("data/GWAS/PI_LM.csv")
df2b <- df2b[df2b$maf > 0.1, ]
df2b <- merge(df2b,annotation,by=c("chr", "pos"))
head(df2b)

head(df2)
df2 <- merge(df2a,df2b,by=c("chr", "pos"))
head(df2)
nrow(df2)

p <- df2$score.x
p=2*pnorm(-abs(p))
nullRatio <- pi0est(p)
nullRatioS <- pi0est(p, lambda=seq(0.40, 0.95, 0.05), smooth.log.pi0="TRUE")
nullRatioM <- pi0est(p, pi0.method="bootstrap")
qobj = qvalue(p, lambda=seq(0.05, 0.95, 0.1), smooth.log.pi0="TRUE")
qvalues.x <- qobj$qvalues
lfdr <- qobj$lfdr
df2 <-cbind(df2,qvalues.x)

p <- df2$score.y
p=2*pnorm(-abs(p))
nullRatio <- pi0est(p)
nullRatioS <- pi0est(p, lambda=seq(0.40, 0.95, 0.05), smooth.log.pi0="TRUE")
nullRatioM <- pi0est(p, pi0.method="bootstrap")
qobj = qvalue(p, lambda=seq(0.05, 0.95, 0.1), smooth.log.pi0="TRUE")
qvalues.y <- qobj$qvalues
lfdr <- qobj$lfdr
df2 <-cbind(df2,qvalues.y)

df2<-df2 %>% unite(annotationb, gene.y, pos, sep = " pos. ")
df2 <- df2[order(df2$qvalues.x) ,  ]
head(df2)


df2$GVE.x <- NULL
df2$gene.x <- NULL
df2$INFO.x <- NULL
df2$annotation.x <- NULL
df2$Gene_id.x <- NULL
df2$GVE.y <- NULL

write_csv(df2, "/Tables/GWAS_d13C vs leaf area stability index drought.csv")
df2 <- df2[order(-df2$score.x) ,  ]
head(df2)
df2 <- df2 %>% 
  rename(
    score.d13C = score.x,
    score.LA_plasticity = score.y,
    annotation = annotation.y)
head(df2)

p2<-df2 %>% 
  mutate(quadrant = case_when(score.d13C > 4 & score.LA_plasticity > 4   ~ "Q1",
                              score.d13C <= 4 & score.LA_plasticity > 4  ~ "Q2",
                              score.d13C <= 4 & score.LA_plasticity <= 4 ~ "Q3",
                              TRUE                                         ~ "Q4")) %>% 
  mutate(label = case_when(score.d13C >= 4 & score.LA_plasticity >= 4 ~  paste(annotationb))) %>% 
  ggplot(aes(x = score.d13C, y = score.LA_plasticity)) + 
  geom_point(aes(color=quadrant),alpha = 0.5,  size=3) + 
  scale_colour_manual(values = c("#6689c3", "#c5c5c5", "#e9af41"))+
  geom_text_repel(aes(score.d13C, score.LA_plasticity, label = label), fontface = 'bold', size=1, force=0.1,segment.size = 0.1,arrow = arrow(length = unit(0.01, 'npc')),family = 'Arial',
                  fontface = 'bold',
                  box.padding = unit(0.15, "lines"),
                  point.padding = unit(0.45, "lines"),
                  segment.color = 'grey50') + 
  geom_vline(xintercept = 4, color = "black", linetype = "dashed")+
  geom_hline(yintercept = 4, color = "black", linetype = "dashed")+
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=18)) +  
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm"))  +
  theme(axis.text = element_text(family = "Arial", color="black",  face="bold",size=16)) +
  theme(plot.title = element_text(family = "Arial", color="#696969",face="bold",  size=35))  +  theme(panel.background = element_rect(fill = 'white')) +
  ggtitle("GWAS")+
  theme(legend.position="none") + 
  ylab(expression(bold(atop("Leaf area plasticity index", paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P")))))))+
  xlab(expression(bold(atop(paste("WUE"[i]," (\u03B4"^"13", "C)"), paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P")))))))+ 
  scale_x_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2))+
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2))

p2
png("panels/p2.png", width = 7, height = 7, units = 'in', res = 350)
p2
dev.off() 
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
####Create an annotation file compatible with the files we are creating.
annotation = read_csv("data/gene_description.csv")
head(annotation)
#annotation<-annotation %>% 
# rename(
#  Gene_id = gene
#)
head(annotation)

###Open dataframe with phenotypes
phenotypes=read_csv("data/phenotypesb.csv")
spec(phenotypes)


head(phenotypes)
#####Run autocorrelation matrix condition B
head(phenotypes)
phenotypesselected<- select(phenotypes, PI)
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
outputB <-cbind(accession_idsorted,AA)
outputB <- outputB[-1, 1:2]
outputB <- as.data.frame(outputB)
outputB <- outputB %>% 
  rownames_to_column(var = "gene")
outputB <-cbind(accession_idsorted,BB)
outputB <- outputB[order(-outputB$PI) ,  ]
head(outputB)
outputB <- merge(outputB,annotation,by="gene")
outputB$accession_idsorted.x <- NULL
outputB$accession_idsorted.y <- PI
outputB <- outputB[order(-outputB$PI) ,  ]
head(outputB, n=100)
tail(outputB, n=50)
head(outputB)
write.csv(outputB, "Tables/TWAS_PI.csv", row.names = FALSE)

#####Run autocorrelation matrix condition C
head(phenotypes)
phenotypesselected<- select(phenotypes, iWUE)
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
outputC <-cbind(accession_idsorted,AA)
outputC <- outputC[-1, 1:2]
outputC <- as.data.frame(outputC)
outputC <- outputC %>% 
  rownames_to_column(var = "gene")
outputC <-cbind(accession_idsorted,BB)
outputC <- outputC[order(-outputC$iWUE) ,  ]
head(outputC)
outputC <- merge(outputC,annotation,by="gene")
outputC$accession_idsorted.x <- NULL
outputC$accession_idsorted.y <- NULL
outputC <- outputC[order(-outputC$iWUE) ,  ]
head(outputC, n=100)
tail(outputC, n=50)
head(outputC)
write.csv(outputC, "Tables/TWAS_iWUE.csv", row.names = FALSE)

#####Merge df4 
df3 <- merge(outputC,outputB,by="gene")
nrow(df3)
df3[df3 == "NAN"] <- NA
df3[complete.cases(df3),]

df3$accession_idsorted.y <- NULL
df3$annotation.x <- NULL
df3$accession_idsorted.x <- NULL
df3$Gene_id.x <- NULL
df3$symbol.x <- NULL
df3$annotation.x <- NULL
head(df3)

df3 <- df3 %>% 
  rename(
    symbol = symbol.y,
    annotation = annotation.y)
head(df3)


p3<-df3 %>% 
  mutate(quadrant = case_when(iWUE <= -0.4 & PI >= 0.4 ~ "Q1",
                              iWUE <= -0.4 & PI <= -0.4 ~ "Q2",
                              iWUE >= 0.4 & PI >= 0.4 ~ "Q3",
                              iWUE >= 0.4 & PI <= -0.4 ~ "Q4",
                              TRUE                                         ~ "Q5")) %>% 
  mutate(label = case_when(iWUE <= -0.4 & PI >= 0.4 ~  paste(symbol),
                           iWUE <= -0.4 & PI <= -0.4 ~  paste(symbol),
                           iWUE >= 0.4 & PI >= 0.4 ~  paste(symbol),
                           iWUE >= 0.4 & PI <= -0.4 ~  paste(symbol))) %>% 
  
  ggplot(aes(x = iWUE, y = PI)) + 
  geom_point(aes(color=quadrant),alpha = 0.8,  size=12, shape = "-") +  
  scale_colour_manual(values = c( "#6689c3","#ce494a", "#47a569", "#e9af41", "#c5c5c5"))+
  geom_text_repel(aes(iWUE, PI, label = label),  size=2, force=0.1,segment.size = 0.3,arrow = arrow(length = unit(0.01, 'npc')),family = 'Arial',
                  fontface = 'bold',
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.35, "lines"),
                  segment.color = 'grey50',max.overlaps = Inf) + 
  geom_vline(xintercept = 0.4, color = "black", linetype = "dashed")+
  geom_vline(xintercept = -0.4, color = "black", linetype = "dashed")+
  geom_hline(yintercept = 0.4, color = "black", linetype = "dashed")+
  geom_hline(yintercept = -0.4, color = "black", linetype = "dashed")+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm"))  +
  theme(axis.text = element_text(family = "Arial", color="black", size=18)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=16)) +  
  theme(legend.position="none") + 
  theme(plot.title = element_text(family = "Arial", color="#696969",face="bold",  size=35)) +  theme(panel.background = element_rect(fill = 'white')) +
  ggtitle("TWAS")+
  xlab(expression(bold(atop(paste("WUE"[i]," (\u03B4"^"13", "C)"), paste("transcriptome-wide association ", (r[s]))))))+ 
  ylab(expression(bold(atop("Leaf area plasticity index", paste("transcriptome-wide association ", (r[s]))))))+ 
  scale_y_continuous(limits = c(-0.8, 0.8), breaks = seq(-0.8, 0.8, by = 0.4))


p3

png("panels/p3.png", width = 7, height = 7, units = 'in', res = 350)
p3
dev.off() 

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################

phenotypes=read_csv("data/phenotypesb.csv")
head(phenotypes)
phenotypesselected<- select(phenotypes, accession_name,PI)
accession_id1<-phenotypes$accession_id
selected <-cbind(accession_id1, phenotypesselected)
selected <-selected %>% 
  rename( accession_id = accession_id1)
head(selected)
transposedfinal = read_csv("data/transposedfinal.csv")
accession_id2<-transposedfinal$accession_id

transposedfinal<- select(transposedfinal, AT3G48010)
transposedfinal <-cbind(accession_id2, transposedfinal)
transposedfinal <-transposedfinal %>% 
  rename( accession_id = accession_id2)
head(transposedfinal)

df4<-merge(selected, transposedfinal, by = "accession_id", sort = TRUE)
head(df1)
formula <- y ~ x

head(df1)
###Plot2= d13C vs leaf area plasticity in response to drought
p4 <- ggplot(df4, aes(x=AT3G48010 , y= PI)) +
  geom_point(color="#47a569",alpha = 0.7,  size=12) +  
    geom_smooth(method = "lm", se = T, fill="#696969", colour="#696969", size=0.7, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                             label.x.npc = "right", label.y.npc = "top",
                                                                                                             formula = formula, parse = TRUE, size = 7, colour="black") +
  geom_text_repel(aes(AT3G48010, PI, label = accession_name), fontface = 'bold', size=5, force=0.3,segment.size = 1,family = 'Arial',
                  fontface = 'bold',
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.35, "lines"),
                  segment.color = 'grey50') + 
  
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=35)) +  
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm"))  +
  theme(axis.text = element_text(family = "Arial", color="black",  face="bold",size=20)) +
  theme(plot.title = element_text(family = "Arial", color="black", face="bold",size=28)) +  theme(panel.background = element_rect(fill = 'white')) +
  theme(legend.position="none") + 
  ylab(expression(bold("Leaf area plasticity index")))+ 
  xlab(expression(bold(paste(bolditalic("CNGC16")~"(FPKM)"))))
#ggtitle("Phenotype x phenotype")
# scale_y_continuous(limits = c(-0.75, 0.75), breaks = seq(-0.75, 0.75, by = 0.25))+
#scale_x_continuous(limits = c(-37, -33), breaks = seq(-37, -33, by = 1))

p4
png("panels/p4.png", width = 8, height = 7, units = 'in', res = 350)
p4
dev.off() 



###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
phenotypesselected<- select(phenotypes, accession_name,iWUE)
accession_id1<-phenotypes$accession_id
selected <-cbind(accession_id1, phenotypesselected)
selected <-selected %>% 
  rename( accession_id = accession_id1)
head(selected)
transposedfinal = read_csv("data/transposedfinal.csv")
accession_id2<-transposedfinal$accession_id

transposedfinal<- select(transposedfinal, AT3G48010)
transposedfinal <-cbind(accession_id2, transposedfinal)
transposedfinal <-transposedfinal %>% 
  rename( accession_id = accession_id2)
head(transposedfinal)

df5<-merge(selected, transposedfinal, by = "accession_id", sort = TRUE)
head(df5)
formula <- y ~ x
p5 <- ggplot(df5, aes(x=AT3G48010 , y= iWUE)) +
  geom_point(color="#47a569",alpha = 0.7,  size=12) +  
  geom_smooth(method = "lm", se = T, fill="#696969", colour="#696969", size=0.7, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                             label.x.npc = "left", label.y.npc = "top",
                                                                                                             formula = formula, parse = TRUE, size = 7, colour="black") +
  geom_text_repel(aes(AT3G48010, iWUE, label = accession_name), fontface = 'bold', size=5, force=0.3,segment.size = 1,family = 'Arial',
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.35, "lines"),
                  segment.color = 'grey50') + 
  
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=35)) +  
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm"))  +
  theme(axis.text = element_text(family = "Arial", color="black",  face="bold",size=20)) +
  theme(plot.title = element_text(family = "Arial", color="black", face="bold",size=28)) +  theme(panel.background = element_rect(fill = 'white')) +
  theme(legend.position="none") + 
  ylab(expression(bold(paste("WUE"[i]," (\u03B4"^"13", "C (\u2030))"))))+   
  xlab(expression(bold(paste(bolditalic("CNGC16")~"(FPKM)"))))

p5
png("panels/p5.png", width = 8, height = 7, units = 'in', res = 350)
p5
dev.off() 

##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################

phenotypes=read_csv("data/phenotypesb.csv")
head(phenotypes)
phenotypesselected<- select(phenotypes, accession_name,PI)
accession_id1<-phenotypes$accession_id
selected <-cbind(accession_id1, phenotypesselected)
selected <-selected %>% 
  rename( accession_id = accession_id1)
head(selected)
transposedfinal = read_csv("data/transposedfinal.csv")
accession_id2<-transposedfinal$accession_id

transposedfinal<- select(transposedfinal, AT5G08180)
transposedfinal <-cbind(accession_id2, transposedfinal)
transposedfinal <-transposedfinal %>% 
  rename( accession_id = accession_id2)
head(transposedfinal)

df6<-merge(selected, transposedfinal, by = "accession_id", sort = TRUE)
head(df6)
formula <- y ~ x

head(df6)
###Plot2= d13C vs leaf area plasticity in response to drought
p6 <- ggplot(df6, aes(x=AT5G08180 , y= PI)) +
  geom_point(color="#ce494a",alpha = 0.7,  size=12) +  
  
  geom_smooth(method = "lm", se = T, fill="#696969", colour="#696969", size=0.7, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                             label.x.npc = "left", label.y.npc = "top",
                                                                                                             formula = formula, parse = TRUE, size = 7, colour="black") +
  geom_text_repel(aes(AT5G08180, PI, label = accession_name), fontface = 'bold', size=5, force=0.3,segment.size = 1,family = 'Arial',
                  fontface = 'bold',
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.35, "lines"),
                  segment.color = 'grey50') + 
  
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=35)) +  
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm"))  +
  theme(axis.text = element_text(family = "Arial", color="black",  face="bold",size=20)) +
  theme(plot.title = element_text(family = "Arial", color="black", face="bold",size=28)) +  theme(panel.background = element_rect(fill = 'white')) +
  theme(legend.position="none") + 
  ylab(expression(bold("Leaf area plasticity index")))+ 
  xlab(expression(bold(paste(bolditalic("CBF5")~"(FPKM)"))))#+
#scale_x_continuous(limits = c(0, 250), breaks = seq(0, 250, by = 40))

p6
png("panels/p6.png", width = 8, height = 7, units = 'in', res = 350)
p6
dev.off() 
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
phenotypesselected<- select(phenotypes, accession_name,iWUE)
accession_id1<-phenotypes$accession_id
selected <-cbind(accession_id1, phenotypesselected)
selected <-selected %>% 
  rename( accession_id = accession_id1)
head(selected)
transposedfinal = read_csv("data/transposedfinal.csv")
accession_id2<-transposedfinal$accession_id

transposedfinal<- select(transposedfinal, AT3G57150)
transposedfinal <-cbind(accession_id2, transposedfinal)
transposedfinal <-transposedfinal %>% 
  rename( accession_id = accession_id2)
head(transposedfinal)

df7<-merge(selected, transposedfinal, by = "accession_id", sort = TRUE)
head(df7)
formula <- y ~ x
p7 <- ggplot(df7, aes(x=AT3G57150 , y= iWUE)) +
  geom_point(color="#ce494a",alpha = 0.7,  size=12) +  
  
  geom_smooth(method = "lm", se = T, fill="#696969", colour="#696969", size=0.7, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                             label.x.npc = "right", label.y.npc = "top",
                                                                                                             formula = formula, parse = TRUE, size = 7, colour="black") +
  geom_text_repel(aes(AT3G57150, iWUE, label = accession_name), fontface = 'bold', size=5, force=0.3,segment.size = 1,family = 'Arial',
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.35, "lines"),
                  segment.color = 'grey50') + 
  
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=35)) +  
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm"))  +
  theme(axis.text = element_text(family = "Arial", color="black",  face="bold",size=20)) +
  theme(plot.title = element_text(family = "Arial", color="black", face="bold",size=28)) +  theme(panel.background = element_rect(fill = 'white')) +
  theme(legend.position="none") + 
  xlab(expression(bold(paste(bolditalic("CBF5")~"(FPKM)"))))+
  ylab(expression(bold(paste("WUE"[i]," (\u03B4"^"13", "C (\u2030))"))))


p7
png("panels/p7.png", width = 8, height = 7, units = 'in', res = 350)
p7
dev.off() 

M4 <- image_read("panels/p4.png")
M5 <- image_read("panels/p5.png")
M6 <- image_read("panels/p6.png")
M7 <- image_read("panels/p7.png")



plot_grid(rasterGrob(M6),rasterGrob(M4),rasterGrob(M7),rasterGrob(M5), labels=c('d', 'f',"e","g"), greedy = TRUE, ncol =2, nrow = 2, align = 'v', hjust=0.05, label_size=30)
png("panels/examples.png", width = 9, height = 7, units = 'in', res = 350)
plot_grid(rasterGrob(M6),rasterGrob(M4),rasterGrob(M7),rasterGrob(M5), labels=c('d', 'f',"e","g"), greedy = TRUE, ncol =2, nrow = 2, align = 'v', hjust=0.05, label_size=30)
dev.off()




M1 <- image_read("panels/p1.png")
M2 <- image_read("panels/p2.png")
M3 <- image_read("panels/p3.png")
examples <- image_read("panels/examples.png")



plot_grid(rasterGrob(M1), rasterGrob(M3),rasterGrob(M2),rasterGrob(examples),labels=c('a', 'c',"b",""), greedy = TRUE, ncol =2, nrow = 2, align = 'v', hjust=0, label_size=20)
png("Figure 4.png", width = 8, height = 7, units = 'in', res = 350)
plot_grid(rasterGrob(M1), rasterGrob(M3),rasterGrob(M2),rasterGrob(examples),labels=c('a', 'c',"b",""), greedy = TRUE, ncol =2, nrow = 2, align = 'v', hjust=0, label_size=20)
dev.off()


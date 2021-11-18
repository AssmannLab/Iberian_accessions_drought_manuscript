list.of.packages <- c("ggplot2", "ggrepel","readr", "tidyverse","tidyr","cowplot","qvalue","dplyr","viridis","ggpubr","ggthemes")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(Hmisc)
library(readr)
library(tidyverse)
library(dplyr)
library(tibble)
library(scico)
library(tidyr)
library(cowplot)
library(grid)
library(png)
library(magick)
library(qvalue)
library(dplyr)
library(ggplot2)
library(readr)
library(viridis)
library(tidyr)
library(ggrepel)
library(ggpubr)
library(tidyr)
library(ggthemes)
library(viridis)
library(ggpmisc)
library(ggpubr)
library(cowplot)

############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################

####Create an annotation file compatible with the files we are creating.
annotation = read_csv("data/gene_description.csv")
head(annotation)
annotation<-annotation %>% 
  rename(
    Gene_id = gene
  )
head(annotation)

###Open dataframe with phenotypes
phenotypes = read.csv("data/phenotypesb.csv", header = TRUE, sep = ",", quote = "\"",
                      dec = ".", fill = TRUE,comment.char = "")
head(phenotypes)
#####Run autocorrelation matrix: Leaf area potential Well-watered
phenotypesselected<- select(phenotypes, LAP_well_watered)
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
outputA <-cbind(accession_idsorted,AA)
outputA <- outputA[-1, 1:2]
outputA <- as.data.frame(outputA)
outputA <- outputA[order(-outputA$LAP_well_watered) ,  ]
head(outputA)
outputA <- outputA %>% 
rownames_to_column(var = "Gene_id")
outputA <- merge(outputA,annotation,by="Gene_id")
head(outputA, n=50)
tail(outputA, n=50)
outputA$accession_idsorted <- NULL
head(outputA)
# Basic histogram plot with threshold line 
p1<- gghistogram(outputA, x = "LAP_well_watered", bins = 40, 
            fill = "#008764", color = "#008764",
            rug = FALSE)+geom_vline(xintercept = c(0.4,-0.4), color = "black", linetype = "dashed")+
  font("xlab", size = 30, face = "bold")+
  font("ylab", size = 30, face = "bold")+
  font("xy.text", size = 20, face = "bold")  +
             xlim(-1, 1)+ xlab(expression(bold(atop("Relative leaf area (rLA, well-watered):", paste("transcriptome-wide association ", (r[s])))))) +
  ylab(expression(bold(Count)))

p1
png("panels/p1.png", width = 15, height = 7, units = 'in', res = 350)
p1
dev.off() 
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
#####Run autocorrelation matrix: Leaf area potential drought
phenotypesselected<- select(phenotypes, LAP_drought)
accession_id1<-phenotypes$accession_id
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
outputB <- outputB[order(-outputB$LAP_drought) ,  ]
head(outputB)
outputB <- outputB %>% 
  rownames_to_column(var = "Gene_id")
outputB <- merge(outputB,annotation,by="Gene_id")
head(outputB, n=50)
tail(outputB, n=50)
# Basic histogram plot with threshold line 
p2<- gghistogram(outputB, x = "LAP_drought", bins = 40, 
                 fill = "#9f4c49", color = "#9f4c49",
                 rug = FALSE)+geom_vline(xintercept = c(0.4,-0.4), color = "black", linetype = "dashed")+
  font("xlab", size = 30, face = "bold")+
  font("ylab", size = 30, face = "bold")+
  font("xy.text", size = 20, face = "bold")  +
  xlim(-1, 1)+ xlab(expression(bold(atop("Relative leaf area (rLA, drought):", paste("transcriptome-wide association ", (r[s])))))) +
  ylab(expression(bold(Count)))

p2
png("panels/p2.png", width = 15, height = 7, units = 'in', res = 350)
p2
dev.off() 


#####Run autocorrelation matrix: Leaf area plasticity index (drought)

phenotypesselected<- select(phenotypes, PI)
accession_id1<-phenotypes$accession_id
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
outputC <- outputC[order(-outputC$PI) ,  ]
head(outputC)
outputC <- outputC %>% 
  rownames_to_column(var = "Gene_id")
outputC <- merge(outputC,annotation,by="Gene_id")
head(outputC, n=50)
tail(outputC, n=50)
outputC$accession_idsorted <- NULL

# Basic histogram plot with threshold line 
p3<- gghistogram(outputC, x = "PI", bins = 40, 
                 fill = "#ffce00", color = "#ffce00",
                 rug = FALSE)+geom_vline(xintercept = c(0.4,-0.4), color = "black", linetype = "dashed")+
  font("xlab", size = 30, face = "bold")+
  font("ylab", size = 30, face = "bold")+
  font("xy.text", size = 20, face = "bold")  +
  xlim(-1, 1)+ xlab(expression(bold(atop("Leaf area plasticity index (PI):", paste("transcriptome-wide association ", (r[s])))))) +
  ylab(expression(bold(Count)))

p3
png("panels/p3.png", width = 15, height = 7, units = 'in', res = 350)
p3
dev.off() 

#####Run autocorrelation matrix: Yield potential short days Well-watered

phenotypesselected<- select(phenotypes, iWUE)
accession_id1<-phenotypes$accession_id
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
outputD <-cbind(accession_idsorted,AA)
outputD <- outputD[-1, 1:2]
outputD <- as.data.frame(outputD)
outputD <- outputD[order(-outputD$iWUE) ,  ]
head(outputD)
outputD <- outputD %>% 
  rownames_to_column(var = "Gene_id")
outputD <- merge(outputD,annotation,by="Gene_id")
head(outputD, n=50)
tail(outputD, n=50)
outputD$accession_idsorted <- NULL
# Basic histogram plot with threshold line 
p4<- gghistogram(outputD, x = "iWUE", bins = 40, 
                 fill = "#00a7ff", color = "#00a7ff",
                 rug = FALSE)+geom_vline(xintercept = c(0.4,-0.4), color = "black", linetype = "dashed")+
  font("xlab", size = 30, face = "bold")+
  font("ylab", size = 30, face = "bold")+
  font("xy.text", size = 20, face = "bold")  +
  xlim(-1, 1)+ 
  xlab(expression(bold(atop(paste("Our sample WUE"[i]," (\u03B4"^"13", "C)"), paste("transcriptome-wide association ", (r[s]))))))+ 
    ylab(expression(bold(Count)))

p4
png("panels/p4.png", width = 15, height = 7, units = 'in', res = 350)
p4
dev.off() 


#####Run autocorrelation matrix: Yield potential short days drought
phenotypes2 = read.csv("data/Iberian_d13C.csv", header = TRUE, sep = ",", quote = "\"",
                      dec = ".", fill = TRUE,comment.char = "")
head(phenotypes2)
phenotypesselected<- select(phenotypes2, Iberian_d13C)
accession_id1<-phenotypes2$accession_id
accession_id1<-phenotypes2$accession_id
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
outputE <-cbind(accession_idsorted,AA)
outputE <- outputE[-1, 1:2]
outputE <- as.data.frame(outputE)
outputE <- outputE[order(-outputE$Iberian_d13C) ,  ]
head(outputE)
outputE <- outputE %>% 
  rownames_to_column(var = "Gene_id")
outputE <- merge(outputE,annotation,by="Gene_id")
head(outputE, n=50)
tail(outputE, n=50)
outputE$accession_idsorted <- NULL

p5<- gghistogram(outputE, x = "Iberian_d13C", bins = 40, 
                 fill = "#074589", color = "#074589",
                 rug = FALSE)+geom_vline(xintercept = c(0.4,-0.4), color = "black", linetype = "dashed")+
  font("xlab", size = 30, face = "bold")+
  font("ylab", size = 30, face = "bold")+
  font("xy.text", size = 20, face = "bold")  +
  xlim(-1, 1)+ 
  xlab(expression(bold(atop(paste("1001 Genomes. Iberian population WUE"[i]," (\u03B4"^"13", "C)"), paste("transcriptome-wide association ", (r[s]))))))+ 
  ylab(expression(bold(Count)))

p5
png("panels/p5.png", width = 15, height = 7, units = 'in', res = 350)
p5
dev.off() 

##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
annotation = read_csv("data/annotation_GWAS.csv")

dat1 = read.csv("data/GWAS/LAP_WW_LM.csv")

head(dat1)
dat1 <- merge(dat1,annotation,by=c("chr", "pos"))
head(dat1)
dat1$X1 <- NULL
dat1$GVE<- NULL

dat1 <- dat1[order(-dat1$score) ,  ]
dat1 <- dat1[dat1$maf > 0.1, ]
# Basic histogram plot with threshold line 
p6<- gghistogram(dat1, x = "score", bins = 40, 
                 fill = "#008764", color = "#008764",
                 rug = FALSE)+geom_vline(xintercept = 4, color = "black", linetype = "dashed")+
  font("xlab", size = 30, face = "bold")+
  font("ylab", size = 30, face = "bold")+
  font("xy.text", size = 20, face = "bold")  +
  xlim(0, 10)+ xlab(expression(bold(atop(paste("Relative leaf area (rLA, well-watered):"), paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P")))))))+ 
  ylab(expression(bold(Count)))

p6
png("panels/p6.png", width = 15, height = 7, units = 'in', res = 350)
p6
dev.off() 

##########################################################################################################################################################
##########################################################################################################################################################
annotation = read_csv("data/annotation_GWAS.csv")

dat2 = read.csv("data/GWAS/LAP_Drought_LM.csv")

head(dat2)
dat2 <- merge(dat2,annotation,by=c("chr", "pos"))
head(dat1)
dat2$X1 <- NULL
dat2$GVE<- NULL

dat2 <-dat2[!(dat2$INFO == "downstream_gene_variant"), ]

dat2 <- dat2[order(-dat2$score) ,  ]
dat2 <- dat2[dat2$maf > 0.1, ]
# Basic histogram plot with threshold line 
p7<- gghistogram(dat2, x = "score", bins = 40, 
                 fill = "#9f4c49", color = "#9f4c49",
                 rug = FALSE)+geom_vline(xintercept = 4, color = "black", linetype = "dashed")+
  font("xlab", size = 30, face = "bold")+
  font("ylab", size = 30, face = "bold")+
  font("xy.text", size = 20, face = "bold")  +
  xlim(0, 10)+  xlab(expression(bold(atop(paste("Relative leaf area (rLA, drought)"), paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P")))))))+ 
  
  ylab(expression(bold(Count)))

p7
png("panels/p7.png", width = 15, height = 7, units = 'in', res = 350)
p7
dev.off() 
##########################################################################################################################################################
##########################################################################################################################################################
annotation = read_csv("data/annotation_GWAS.csv")

dat3 = read.csv("data/GWAS/PI.csv")

dat3 <- dat3[order(-dat3$score) ,  ]
dat3 <- dat3[dat3$maf > 0.1, ]
# Basic histogram plot with threshold line 
p8<- gghistogram(dat3, x = "score", bins = 40, 
                 fill = "#ffce00", color = "#ffce00",
                 rug = FALSE)+geom_vline(xintercept = 4, color = "black", linetype = "dashed")+
  font("xlab", size = 30, face = "bold")+
  font("ylab", size = 30, face = "bold")+
  font("xy.text", size = 20, face = "bold")  +
  xlim(0, 10)+  xlab(expression(bold(atop(paste("Leaf area plasticity index"), paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P")))))))+ 
  
  ylab(expression(bold(Count)))

p8
png("panels/p8.png", width = 15, height = 7, units = 'in', res = 350)
p8
dev.off() 

##########################################################################################################################################################
##########################################################################################################################################################
annotation = read_csv("data/annotation_GWAS.csv")

dat4 = read.csv("data/GWAS/Our_sample_d13C_LM.csv")

head(dat4)
dat4 <- merge(dat4,annotation,by=c("chr", "pos"))
head(dat4)
dat4$X1 <- NULL
dat4$GVE<- NULL

dat4 <-dat4[!(dat4$INFO == "downstream_gene_variant"), ]

dat4 <- dat4[order(-dat4$score) ,  ]
dat4 <- dat4[dat4$maf > 0.1, ]
# Basic histogram plot with threshold line 
p9<- gghistogram(dat4, x = "score", bins = 40, 
                 fill = "#00a7ff", color = "#00a7ff",
                 rug = FALSE)+geom_vline(xintercept = 4, color = "black", linetype = "dashed")+
  font("xlab", size = 30, face = "bold")+
  font("ylab", size = 30, face = "bold")+
  font("xy.text", size = 20, face = "bold")  +
  xlim(0, 10)+ 
  xlab(expression(bold(atop(paste("Our sample WUE"[i]," (\u03B4"^"13", "C)"), paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P")))))))+ 
  ylab(expression(bold(Count)))

p9
png("panels/p9.png", width = 15, height = 7, units = 'in', res = 350)
p9
dev.off() 
##########################################################################################################################################################
##########################################################################################################################################################
annotation = read_csv("data/annotation_GWAS.csv")

dat5 = read.csv("data/GWAS/Iberian_d13C_LM.csv")

head(dat5)
dat5 <- merge(dat5,annotation,by=c("chr", "pos"))
head(dat5)
dat5$GVE<- NULL

head (dat5)
dat5 <- dat5[order(-dat5$score) ,  ]
dat5 <- dat5[dat5$maf > 0.1, ]
# Basic histogram plot with threshold line 
p10<- gghistogram(dat5, x = "score", bins = 40, 
                 fill = "#074589", color = "#074589",
                 rug = FALSE)+geom_vline(xintercept = 4, color = "black", linetype = "dashed")+
  font("xlab", size = 30, face = "bold")+
  font("ylab", size = 30, face = "bold")+
  font("xy.text", size = 20, face = "bold")  +
  xlim(0, 10)+ 
  xlab(expression(bold(atop(paste("1001 Genomes. Iberian population WUE"[i]," (\u03B4"^"13", "C)"), paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P")))))))+ 
  ylab(expression(bold(Count)))

p10
png("panels/p10.png", width = 15, height = 7, units = 'in', res = 350)
p10
dev.off() 
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################

M6 <- image_read("panels/p1.png")
M7 <- image_read("panels/p2.png")
M8 <- image_read("panels/p3.png")
M9 <- image_read("panels/p4.png")
M10 <- image_read("panels/p5.png")
M1 <- image_read("panels/p6.png")
M2 <- image_read("panels/p7.png")
M3 <- image_read("panels/p8.png")
M4 <- image_read("panels/p9.png")
M5 <- image_read("panels/p10.png")


plot_grid(rasterGrob(M1),rasterGrob(M6),rasterGrob(M2),rasterGrob(M7),rasterGrob(M3),rasterGrob(M8),rasterGrob(M4),rasterGrob(M9),rasterGrob(M5),rasterGrob(M10),labels=c('a', 'f','b', 'g','c', 'h','d', 'i','e', 'j'), greedy = TRUE, ncol =2, nrow = 5, align = 'v', hjust=0, label_size=12)
png("Figure S3.png", width = 7, height = 7, units = 'in', res = 350)
plot_grid(rasterGrob(M1),rasterGrob(M6),rasterGrob(M2),rasterGrob(M7),rasterGrob(M3),rasterGrob(M8),rasterGrob(M4),rasterGrob(M9),rasterGrob(M5),rasterGrob(M10),labels=c('a', 'f','b', 'g','c', 'h','d', 'i','e', 'j'), greedy = TRUE, ncol =2, nrow = 5, align = 'v', hjust=0, label_size=12)
dev.off()


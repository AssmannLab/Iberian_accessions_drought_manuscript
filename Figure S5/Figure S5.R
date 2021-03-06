list.of.packages <- c("ggplot2", "ggrepel","readr", "tidyverse","tidyr","cowplot","qvalue","dplyr","grid","ggpubr","ggthemes","png")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(readr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(qvalue)
library(viridis)
library(Hmisc)
library(cowplot)
library(patchwork)
library(magick)
library(png)
library(grid)
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################Fig. S1A 
####Create an annotation file compatible with the files we are creating.
annotation = read_csv("data/annotation_GWAS.csv")
nrow(annotation)
head(annotation)
gene_description = read_csv("data/gene_description.csv")
head(gene_description)
annotation <- merge(annotation,gene_description,by= "gene")
head(annotation)
#########Step 1: annotate and calculate q-values for both conditions (well-watered (WW) & drought(LW))
df2a=read_csv("data/GWAS/LAP_well_watered_LM.csv")
df2a <- df2a[df2a$maf > 0.1, ]
df2b=read_csv("data/GWAS/LAP_drought_LM.csv")
df2b <- df2b[df2b$maf > 0.1, ]

####Step 2: Merge both dataframes into one dataset based on common genomic location (chromosome & location) & annotate based on the TAIR11 genome release
df2 <- merge(df2a,df2b,by=c("chr", "pos"))
df2 <- merge(df2,annotation,by=c("chr", "pos"))
df2[complete.cases(df2), ]
head(df2)

####Step 3: Calculate phenotypic potentials and stability index. Add then to the dataframe.

df2$plasticity_index <- (df2$score.y) - (df2$score.x)

####Step 4: Create a unique ID for each snp. This is done merging the chromosome & position columns.
df2<- df2 %>%  
  unite(ID, c("chr", "pos"), remove = FALSE)

####Step 5: Clean-up the dataframe.
df2$maf.x <- NULL
df2$mac.x <- NULL
df2$GVE.x <- NULL
df2$mac.y <- NULL
df2$GVE.y <- NULL
df2$symbol.x <- NULL
df2$gene.y <- NULL
df2$INFO.y <- NULL
df2$annotation.y <- NULL
df2<-df2 %>% 
  rename(
    Score_WW =score.x,
    Score_LW=score.y,
    Locus_id = gene,
    MAF = maf.y,
    Description = annotation.x,
    Symbol = symbol.y
  )
df2 <- df2[order(-df2$plasticity_index) ,  ]
head(df2)

#df2 <-df2 %>%
 # filter(Score_WW > 3.99999 | Score_LW>3.99999)##### render plot with high genomic plasticity index
####Step 6: Save the dataframe in the /Tables directory.
write_csv(x = df2, "Tables/GWAS_LAP_plasticity_drought.csv")

####Step 7: Plot the data.
df2<- df2 %>%
  gather(TRT, Score, Score_WW:Score_LW)
head(df2)
tail(df2)
df2 <- df2[order(-df2$plasticity_index) ,  ]
pal <- c( "#008764","#9f4c49")

p2 <- ggplot(df2, aes(factor(TRT), Score, fill=TRT)) +
  scale_x_discrete(limits=c("Score_WW","Score_LW"), labels = c('Well-watered','Drought'), expand = c(0.1, 0.1)) +
  geom_point(aes(color=plasticity_index),size=2,
             alpha = 0.75,
             position=position_jitter(width=0.005)) +
  geom_line(aes(group=ID, color=plasticity_index),size = 0.3, 
            alpha = 0.75)+ scale_colour_gradientn(name="Plasticity index",limits = c(-8, 8), colours = pal) +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.25,"cm")) +
  theme(axis.text.x = element_text(family = "Arial", color="black", face="bold", size=18)) +
  theme(axis.text.y = element_text(family = "Arial", color="black", face="bold", size=10)) +
  theme(axis.title.x = element_text(family = "Arial", color="black", face="bold", size=25)) +
  theme(axis.title.y = element_text(family = "Arial", color="black", face="bold", size=12)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = - 5, family = "Arial", color="#696969", face="bold", size=25)) +
  theme(axis.title.x=element_blank())+ 
  theme(plot.margin=unit(c(1,1,1,1),"cm"))+
  guides(fill = "none") + 
  ylim(0, 6)+
  ylab("Relative leaf area (rLA)\n(-log(p-value))") + 
  xlab("Treatment") + 
  theme(legend.position = "none")+
  geom_hline(yintercept = 4, linetype = "dashed", colour = "grey20", size=0.3)+  ggtitle("GWAS: reaction norms")



p2
png("panels/p2.png", width = 10, height = 7, units = 'in', res = 350)
p2
dev.off() 

############################################################################################################################################################################Fig. S1B
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################

############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################
###########################################################################################################################################################################
############################################################################################################################################################################

####Create an annotation file compatible with the files we are creating.
annotation = read_csv("data/gene_description.csv")
head(annotation)


###Open dataframe with phenotypes
phenotypes = read.csv("data/phenotypesb.csv", header = TRUE, sep = ",", quote = "\"",
                      dec = ".", fill = TRUE,comment.char = "")
head(phenotypes)
#####Run autocorrelation matrix Phenotype A
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
AA<-signif(res$r, 3)
BB<-signif(res$P,3)
outputA1 <-cbind(accession_idsorted,AA)
outputA1 <- outputA1[-1, 1:2]
outputA1 <- as.data.frame(outputA1)
outputA1 <- outputA1 %>% 
  rownames_to_column(var = "gene")

outputA2 <-cbind(accession_idsorted,BB)
outputA2 <- outputA2[-1, 1:2]
outputA2 <- as.data.frame(outputA2)
outputA2 <- outputA2 %>% 
  rownames_to_column(var = "gene")
outputA <- merge(outputA1,outputA2,by="gene")
outputA <- outputA %>% 
  rename(
    LA_WW_rs = LAP_well_watered.x,
    LA_WW_P_val = LAP_well_watered.y)
outputA$accession_idsorted.x <- NULL
outputA$accession_idsorted.y <- NULL

head(outputA)

outputA <- outputA[order(-outputA$LA_WW_rs) ,  ]
head(outputA)

#####Run autocorrelation matrix condition B
phenotypesselected<- select(phenotypes, LAP_drought)
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
AA<-signif(res$r, 3)
BB<-signif(res$P,3)
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
outputB <- outputB %>% 
  rename(
    LA_LW_rs = LAP_drought.x,
    LA_LW_P_val = LAP_drought.y)
outputB$accession_idsorted.x <- NULL
outputB$accession_idsorted.y <- NULL

head(outputB)

outputB <- outputB[order(-outputB$LA_LW_rs) ,  ]
head(outputB)
#####Merge df1

df3 <- merge(outputA,outputB,by="gene")
df3 <- merge(df3,annotation,by="gene")
head(df3)

df3[complete.cases(df3), ]
df3$plasticity_index <- (df3$LA_LW_rs) - (df3$LA_WW_rs)


df3<-df3 %>% 
  rename(
    description = annotation
  )

nrow(df3)
head(df3)

######
#########
#df4 <-filter(df4, ((LAP_SD_WW > 0.399 & LAP_SD_LW < 0.399)| (LAP_SD_WW < -0.399 & LAP_SD_LW > -0.399)| (LAP_SD_LW > 0.399 & LAP_SD_WW < 0.399)| (LAP_SD_LW < -0.399 & LAP_SD_WW > -0.399)))

#df3 <-filter(df3, ((LA_WW_rs > 0.399 & LA_LW_rs < 0)| (LA_WW_rs < -0.399 & LA_LW_rs > 0)| (LA_LW_rs > 0.399 & LA_WW_rs < 0)| (LA_LW_rs < -0.399 & LA_WW_rs > 0)))
df3<-df3[complete.cases(df3), ]
write_csv(x = df3, "Tables/TWAS_Leaf_area_plasticity_drought.csv")
head(df3)
tail(df3)

######Plot data
df3<- df3 %>%
  gather(TRT, rs, LA_LW_rs:LA_WW_rs)
head(df3)
df3 <- df3[order(-df3$plasticity_index) ,  ]
df3 <- na.omit(df3, cols="Phenotypic_potential")

pal <- c("#9f4c49", "#008764")

p3 <- ggplot(df3, aes(factor(TRT), rs, fill=TRT)) +
  scale_x_discrete(limits=c("LA_LW_rs","LA_WW_rs"), labels = c('Well-watered','Drought'), expand = c(0.1, 0.1)) +
  geom_point(aes(color=plasticity_index),size=2,
             alpha = 0.75) +
  geom_line(aes(group=gene, color=plasticity_index),size = 0.3,
            alpha = 0.75)+ scale_colour_gradientn(name="Plasticity index",limits = c(-1, 1), colours = pal) +
  theme(axis.line = element_line(size=0.5, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.25,"cm")) +
  theme(axis.text.x = element_text(family = "Arial", color="black", face="bold", size=18)) +
  theme(axis.text.y = element_text(family = "Arial", color="black", face="bold", size=10)) +
  theme(axis.title.x = element_text(family = "Arial", color="black", face="bold", size=25)) +
  theme(axis.title.y = element_text(family = "Arial", color="black", face="bold", size=12)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = - 5, family = "Arial", color="#696969", face="bold", size=25)) +
  theme(axis.title.x=element_blank())+ 
  guides(fill = FALSE) + 
  theme(legend.position = "none")+
  
    ylim(-0.7, 0.7)+
  ylab(expression(bold(atop("Relative leaf area (rLA)", paste("transcriptome-wide association ", (r[s]))))))+ 
  geom_hline(yintercept = 0, linetype = "solid", colour = "grey20", size=0.3)+
geom_hline(yintercept = 0.4, linetype = "dashed", colour = "grey20", size=0.3)+
  geom_hline(yintercept = -0.4, linetype = "dashed", colour = "grey20", size=0.3)+ggtitle("TWAS: reaction norms")
p3

png("panels/p3.png", width = 10, height = 7, units = 'in', res = 350)
p3
dev.off() 

M4 <- image_read("panels/p2.png")
M5 <- image_read("panels/p3.png")


plot_grid(rasterGrob(M4),rasterGrob(M5),labels=c('a', 'b'), greedy = TRUE, ncol =1, nrow = 2, align = 'v', hjust=0, label_size=45)
png("Figure S5.png", width = 7, height = 12, units = 'in', res = 350)
plot_grid(rasterGrob(M4),rasterGrob(M5),labels=c('a', 'b'), greedy = TRUE, ncol =1, nrow = 2, align = 'v', hjust=0, label_size=45)
dev.off()

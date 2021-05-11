list.of.packages <- c("ggplot2", "ggrepel","readr", "tidyverse","tidyr","cowplot","qvalue","dplyr","viridis","ggpubr","ggthemes","Hmisc","magick","png","grid","tidyr","ggpmisc")
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
library(ggpubr)
library(magick)
library(png)
library(grid)
library(tidyr)
library(ggpmisc)
##############################################################################################################################################################
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

####Create an annotation file compatible with the files we are creating.
annotation = read_csv("data/gene_description.csv")
head(annotation)
annotation<-annotation %>% 
  rename(
    Gene_id = gene
  )
head(annotation)

###Open dataframe with phenotypes
phenotypes = read.csv("data/phenotypes.csv", header = TRUE, sep = ",", quote = "\"",
                      dec = ".", fill = TRUE,comment.char = "")
head(phenotypes)
#####Run autocorrelation matrix Phenotype A
phenotypesselected<- select(phenotypes, LAP_SD_WW)
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
outputA <- outputA[order(-outputA$LAP_SD_WW) ,  ]
head(outputA)
outputA <- outputA %>% 
rownames_to_column(var = "Gene_id")
outputA <- merge(outputA,annotation,by="Gene_id")
head(outputA, n=50)
tail(outputA, n=50)
outputA$accession_idsorted <- NULL
outputA <- outputA %>%    ##Impose r = +-0.39 threshold
filter_at(vars(contains("LAP_SD_WW")), .vars_predicate =  any_vars(abs(.) > 0.39))
outputA <- outputA[order(-outputA$LAP_SD_WW) ,  ]
head(outputA)
write_csv(x = outputA, "Tables/TWAS_LAP_SD_WW.csv")

#####Run autocorrelation matrix condition B
phenotypesselected<- select(phenotypes, LAP_SD_LW)
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
outputB <- outputB[order(-outputB$LAP_SD_LW) ,  ]
head(outputB)
outputB <- outputB %>% 
  rownames_to_column(var = "Gene_id")
outputB <- merge(outputB,annotation,by="Gene_id")
head(outputB, n=50)
tail(outputB, n=50)
outputA$accession_idsorted <- NULL
outputB <- outputB %>%    ##Impose r = +-0.39 threshold
filter_at(vars(contains("LAP_SD_LW")), .vars_predicate =  any_vars(abs(.) > 0.39))
outputB <- outputB[order(-outputB$LAP_SD_LW) ,  ]
tail(outputB)

write_csv(x = outputB, "Tables/TWAS_LAP_SD_LW.csv")


#####Merge df1

df1 <- merge(outputA,outputB,by="Gene_id")
head(df1)
df1[complete.cases(df1), ]
df1$stability_index <- (df1$LAP_SD_LW) - (df1$LAP_SD_WW)

df1$accession_idsorted.y <- NULL
df1$annotation.x <- NULL
df1$accession_idsorted.x <- NULL
df1$symbol.x <- NULL
df1$accession_idsorted <- NULL

df1<-df1 %>% 
  rename(
    description = annotation.y
  )
df1<-df1 %>% 
  rename(
    symbol =  symbol.y

  )
nrow(df1)
head(df1)
######Plot data
df1<- df1 %>%
  gather(TRT, Phenotypic_potential, LAP_SD_LW:LAP_SD_WW)
head(df1)
df1 <- df1[order(-df1$stability_index) ,  ]
df1 <- na.omit(df1, cols="Phenotypic_potential")
write_csv(x = df1, "Tables/TWAS_Leaf_area_homeostasis_drought.csv")

pal <- c("#d53e4f", "#fc8d59", "#fee08b", "#e6f598", "#99d594", "#3288bd")

df1a <-df1[which(df1$stability_index>=-0.4), ]###Here we filter out reaction norms that high high slope and not indicative of homeostasis.
tail(df1a)

#######Plot
data_ends <- df1a %>% filter(TRT == "LAP_SD_LW")

p1a <- ggplot(df1a, aes(factor(TRT), Phenotypic_potential, fill=TRT)) +
  scale_x_discrete(limits=c("LAP_SD_LW","LAP_SD_WW"), labels = c('Well-watered','Drought'), expand = c(0.1, 0.1)) +
  geom_point(color="black",size=8,
             alpha = 1,
             ) +
  geom_text_repel(
    aes(label = symbol), data = data_ends,
    fontface ="bold.italic", color = "#696969", size = 4,box.padding = 1,force = 0.5,
    nudge_x      = -0.15,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2,segment.color = '#696969')+
  geom_line(aes(group=Gene_id, color=stability_index),size = 2,
            alpha = 1, colour = "black") +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.25,"cm")) +
  theme(axis.text.x = element_text(family = "Arial", color="black", face="bold", size=18)) +
  theme(axis.text.y = element_text(family = "Arial", color="black", face="bold", size=7)) +
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_text(family = "Arial", color="black", face="bold", size=12)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = - 5, family = "Arial", color="#696969", face="bold", size=30)) +
  guides(fill = FALSE) +  theme(legend.position = "none")+
  theme(plot.margin=unit(c(1,1,0.1,1),"cm"))+
    ylim(0.4, 0.6)+
  ylab(expression(bold(atop("Leaf area potential:", paste("transcriptome-wide association ", (r[s])))))) +
  xlab("Treatment") +  ggtitle("TWAS: reaction norms (homeostasis)")


p1a


p1b <- ggplot(df1a, aes(factor(TRT), Phenotypic_potential, fill=TRT)) +
  scale_x_discrete(limits=c("LAP_SD_LW","LAP_SD_WW"), labels = c('Well-watered','Drought'), position = "top",expand = c(0.1, 0.1)) +
  geom_point(color = "black",size=8,
             alpha = 1,
  ) +
  geom_text_repel(
    aes(label = symbol), data = data_ends,
    fontface ="bold.italic", color = "#696969", size = 4,box.padding = 1,force = 0.5,
    nudge_x      = -0.15,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2,segment.color = '#696969')+
  geom_line(aes(group=Gene_id, color=stability_index),size = 2,
            alpha = 1, colour = "black") +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.25,"cm")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  theme(axis.text.y = element_text(family = "Arial", color="black", face="bold", size=7)) +
  theme(axis.title.y = element_text(family = "Arial", color="black", face="bold", size=12)) +
  theme(legend.title = element_blank())+ 
  guides(fill = FALSE) +theme(legend.justification=c(1,0), legend.position=c(1,0))+
  theme(plot.margin=unit(c(0.2,1,1,1),"cm"))+
    theme(legend.key.size = unit(0.25, 'cm'))+
  ylim(-0.43, -0.40)+   ylab(expression(bold(atop("Leaf area potential:", paste("transcriptome-wide association ", (r[s])))))) 

p1b

p1 <- plot_grid(
                 p1a, p1b,
                 labels = "", ncol = 1)

p1
png("panels/p1.png", width = 14, height = 7, units = 'in', res = 350)
p1
dev.off() 
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


phenotypes=read_csv("data/phenotypes.csv")
head(phenotypes)
phenotypesselected<- select(phenotypes, name,LAP_SD_WW,LAP_SD_LW)
accession_id1<-phenotypes$accession_id
selected <-cbind(accession_id1, phenotypesselected)
selected <-selected %>% 
  rename( accession_id = accession_id1)
head(selected)
transposedfinal = read_csv("data/transposedfinal.csv")
accession_id2<-transposedfinal$accession_id

transposedfinal<- select(transposedfinal, AT2G43060)
transposedfinal <-cbind(accession_id2, transposedfinal)
transposedfinal <-transposedfinal %>% 
  rename( accession_id = accession_id2)
head(transposedfinal)

df2<-merge(selected, transposedfinal, by = "accession_id", sort = TRUE)
head(df2)

df2<- df2 %>%
  gather(Treatment, LA, LAP_SD_WW:LAP_SD_LW)

head (df2)

df2$Treatment <- as.factor(df2$Treatment)
head(df2)


df2$Treatment <- factor(df2$Treatment, levels=c("LAP_SD_WW", "LAP_SD_LW"))

library(mdthemes)
library(ggplot2)
library(ggpubr)
levels(df2$Treatment) <- c("Well-watered", "Drought")
F = as.formula("y~x")
p2 <-ggplot(df2, aes(AT2G43060, LA, color = Treatment,shape=Treatment)) +
  geom_point(alpha = 0.7,  size=10) + 
  scale_colour_manual(name = "Treatment",
                      labels = c("Well-watered", "Drought"),
                      values = c('#000000','#808080')) +   
  scale_shape_manual(name = "Treatment",
                     labels = c("Well-watered", "Drought"),
                     values = c(16, 15))+
  stat_smooth(aes(fill = Treatment, color = Treatment), method = "lm", formula = F,fullrange = TRUE, se = FALSE) +
  stat_regline_equation(
    label.x = c(200,200), label.y = c(1.75,1.65),
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula=F,size=6,show.legend = FALSE )+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.5,"cm")) +
  theme(axis.text = element_text(family = "Arial", color="black", face="bold", size=12)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=24)) +
  theme(plot.title = element_text(family = "Arial", color="black", face="bold", size=45)) +  theme(panel.background = element_rect(fill = 'white')) +
  ylab("Leaf area potential") + 
  xlab(expression(bold(paste(bolditalic("IBH1")~"(FPKM)"))))+
  theme(legend.position = c(0.80, 0.92))+
  theme(legend.key.size = unit(1, "cm"))+ xlim(200,800)+ ylim(0.1,1.75) + theme(legend.position = c(0.9, 0.15),
                                                                          legend.background = element_rect(fill = "white", color = "black"))       

p2
png("panels/p2.png", width = 8, height = 7, units = 'in', res = 350)
p2
dev.off() 

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
library("tidyr")
library(ggpmisc)

phenotypes=read_csv("data/phenotypes.csv")
head(phenotypes)
phenotypesselected<- select(phenotypes, name,LAP_SD_WW,LAP_SD_LW)
accession_id1<-phenotypes$accession_id
selected <-cbind(accession_id1, phenotypesselected)
selected <-selected %>% 
  rename( accession_id = accession_id1)
head(selected)
transposedfinal = read_csv("data/transposedfinal.csv")
accession_id2<-transposedfinal$accession_id

transposedfinal<- select(transposedfinal, AT4G14220)
transposedfinal <-cbind(accession_id2, transposedfinal)
transposedfinal <-transposedfinal %>% 
  rename( accession_id = accession_id2)
head(transposedfinal)

df3<-merge(selected, transposedfinal, by = "accession_id", sort = TRUE)
head(df3)

df3<- df3 %>%
  gather(Treatment, LA, LAP_SD_WW:LAP_SD_LW)

head (df3)

df3$Treatment <- as.factor(df3$Treatment)
head(df3)


df3$Treatment <- factor(df3$Treatment, levels=c("LAP_SD_WW", "LAP_SD_LW"))

library(mdthemes)
library(ggplot2)
library(ggpubr)
levels(df3$Treatment) <- c("Well-watered", "Drought")
F = as.formula("y~x")
p3 <-ggplot(df3, aes(AT4G14220, LA, color = Treatment,shape=Treatment)) +
  geom_point(alpha = 0.7,  size=10) + 
  scale_colour_manual(name = "Treatment",
                      labels = c("Well-watered", "Drought"),
                      values = c('#000000','#808080')) +   
  scale_shape_manual(name = "Treatment",
                     labels = c("Well-watered", "Drought"),
                     values = c(16, 15))+
  stat_smooth(aes(fill = Treatment, color = Treatment), method = "lm", formula = F,fullrange = TRUE, se = FALSE) +
  stat_regline_equation(
    label.x = c(400,400), label.y = c(1.75,1.65),
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula=F,size=6,show.legend = FALSE)+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.5,"cm")) +
  theme(axis.text = element_text(family = "Arial", color="black", face="bold", size=12)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=24)) +
  theme(plot.title = element_text(family = "Arial", color="black", face="bold", size=45)) +  theme(panel.background = element_rect(fill = 'white')) +
  ylab("Leaf area potential") + 
  xlab(expression(bold(paste(bolditalic("RHF1A")~"(FPKM)"))))+
  theme(legend.position = c(0.80, 0.92))+xlim(400,1600)+ ylim(0.1,1.75) +
  theme(legend.key.size = unit(1, "cm"))+  theme(legend.position = c(0.9, 0.9),legend.background = element_rect(fill = "white", color = "black"))       

p3
png("panels/p3.png", width = 8, height = 7, units = 'in', res = 350)
p3
dev.off() 

############################################################################################################################################################################arrange panels
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

M1 <- image_read("panels/p1.png")
M2 <- image_read("panels/p2.png")
M3 <- image_read("panels/p3.png")
plot_grid(rasterGrob(M2),rasterGrob(M3),labels=c( 'B',"C"), greedy = TRUE, ncol =2, nrow = 1, align = 'v', hjust=0, label_size=47)
png("panels/bottom.png", width = 18, height = 7, units = 'in', res = 350)
plot_grid(rasterGrob(M2),rasterGrob(M3),labels=c( 'B',"C"), greedy = TRUE, ncol =2, nrow = 1, align = 'v', hjust=0, label_size=47)
dev.off()

M4 <- image_read("panels/bottom.png")
plot_grid(rasterGrob(M1),rasterGrob(M4),labels=c( 'A',""), greedy = TRUE, ncol =1, nrow = 2, align = 'v', hjust=0, label_size=17)
png("Figure 3.png", width = 7, height = 7, units = 'in', res = 350)
plot_grid(rasterGrob(M1),rasterGrob(M4),labels=c( 'A',""), greedy = TRUE, ncol =1, nrow = 2, align = 'v', hjust=0, label_size=17)
dev.off()



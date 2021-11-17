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
outputA <- outputA %>%    ##Impose r = +-0.39 threshold
filter_at(vars(contains("LAP_well_watered")), .vars_predicate =  any_vars(abs(.) > 0.39))
outputA <- outputA[order(-outputA$LAP_well_watered) ,  ]
head(outputA)
write_csv(x = outputA, "Tables/TWAS_LAP_well_watered.csv")

#####Run autocorrelation matrix condition B
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
outputA$accession_idsorted <- NULL
outputB <- outputB %>%    ##Impose r = +-0.39 threshold
filter_at(vars(contains("LAP_drought")), .vars_predicate =  any_vars(abs(.) > 0.39))
outputB <- outputB[order(-outputB$LAP_drought) ,  ]
tail(outputB)

write_csv(x = outputB, "Tables/TWAS_LAP_SD_LW.csv")


#####Merge df1

df1 <- merge(outputA,outputB,by="Gene_id")
head(df1)
df1[complete.cases(df1), ]
df1$stability_index <- (df1$LAP_drought) - (df1$LAP_well_watered)

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
  gather(TRT, Phenotypic_potential, LAP_drought:LAP_well_watered)
head(df1)
df1 <- df1[order(-df1$stability_index) ,  ]
df1 <- na.omit(df1, cols="Phenotypic_potential")
write_csv(x = df1, "Tables/TWAS_Leaf_area_homeostasis_drought.csv")

pal <- c("#d53e4f", "#fc8d59", "#fee08b", "#e6f598", "#99d594", "#3288bd")

df1a <-df1[which(df1$stability_index>=-0.4), ]###Here we filter out reaction norms that high high slope and not indicative of homeostasis.
tail(df1a)

#######Plot
data_ends <- df1a %>% filter(TRT == "LAP_drought")

p1a <- ggplot(df1a, aes(factor(TRT), Phenotypic_potential, fill=TRT)) +
  scale_x_discrete(limits=c("LAP_drought","LAP_well_watered"), labels = c('Well-watered','Drought'), expand = c(0.1, 0.1)) +
  geom_point(color="black",size=8,
             alpha = 1,
             ) +
  geom_text_repel(
    aes(label = symbol), data = data_ends,
    fontface ="bold.italic", color = "#696969", size = 5,box.padding = 1,force = 0.5,
    nudge_x      = -0.15,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2,segment.color = '#696969')+
  geom_line(aes(group=Gene_id, color=stability_index),size = 2,
            alpha = 1, colour = "black") +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.25,"cm")) +
  theme(axis.text.x = element_text(family = "Arial", color="black", face="bold", size=18)) +
  theme(axis.text.y = element_text(family = "Arial", color="black", face="bold", size=12)) +
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_text(family = "Arial", color="black", face="bold", size=15)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = - 5, family = "Arial", color="#696969", face="bold", size=20)) +
  guides(fill = "none") +  theme(legend.position = "none")+
  theme(plot.margin=unit(c(1,1,0.1,1),"cm"))+
    ylim(0.4, 0.6)+
  ylab(expression(bold(atop("Relative leaf area (rLA):", paste("transcriptome-wide association ", (r[s])))))) +
  xlab("Treatment") +  ggtitle("TWAS: reaction norms (homeostasis)")


p1a


p1b <- ggplot(df1a, aes(factor(TRT), Phenotypic_potential, fill=TRT)) +
  scale_x_discrete(limits=c("LAP_drought","LAP_well_watered"), labels = c('Well-watered','Drought'), position = "top",expand = c(0.1, 0.1)) +
  geom_point(color = "#7f7f7f",size=8,
             alpha = 1,
  ) +
  geom_text_repel(
    aes(label = symbol), data = data_ends,
    fontface ="bold.italic", color = "#696969", size = 5,box.padding = 1,force = 0.5,
    nudge_x      = -0.15,
    direction    = "y",
    hjust        = 0,
    segment.size = 0.2,segment.color = '#7f7f7f')+
  geom_line(aes(group=Gene_id, color=stability_index),size = 2,
            alpha = 1, colour = "#7f7f7f") +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.25,"cm")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank())+
  theme(axis.text.y = element_text(family = "Arial", color="black", face="bold", size=12)) +
  theme(axis.title.y = element_text(family = "Arial", color="black", face="bold", size=15)) +
  theme(legend.title = element_blank())+ 
  guides(fill = FALSE) +theme(legend.justification=c(1,0), legend.position=c(1,0))+
  theme(plot.margin=unit(c(0.2,1,1,1),"cm"))+
    theme(legend.key.size = unit(0.25, 'cm'))+
  ylim(-0.43, -0.40)+   ylab(expression(bold(atop("Relative leaf area (rLA):", paste("transcriptome-wide association ", (r[s])))))) 

p1b

p1 <- plot_grid(
                 p1a, p1b,
                 labels = "", ncol = 1)

p1
png("panels/p1.png", width = 12, height = 12, units = 'in', res = 350)
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


phenotypes=read_csv("data/phenotypesb.csv")
head(phenotypes)
phenotypesselected<- select(phenotypes, accession_name,LAP_well_watered,LAP_drought)
accession_id1<-phenotypes$accession_id
selected <-cbind(accession_id1, phenotypesselected)
selected <-selected %>% 
  rename( accession_id = accession_id1)
head(selected)
transposedfinal = read_csv("data/transposedfinal.csv")
accession_id2<-transposedfinal$accession_id

transposedfinal<- select(transposedfinal, AT5G39471)
transposedfinal <-cbind(accession_id2, transposedfinal)
transposedfinal <-transposedfinal %>% 
  rename( accession_id = accession_id2)
head(transposedfinal)

df2<-merge(selected, transposedfinal, by = "accession_id", sort = TRUE)
head(df2)

df2<- df2 %>%
  gather(Treatment, LA, LAP_well_watered:LAP_drought)

head (df2)

df2$Treatment <- as.factor(df2$Treatment)
head(df2)


df2$Treatment <- factor(df2$Treatment, levels=c("LAP_well_watered", "LAP_drought"))

library(mdthemes)
library(ggplot2)
library(ggpubr)
levels(df2$Treatment) <- c("Well-watered", "Drought")
F = as.formula("y~x")
p2 <-ggplot(df2, aes(AT5G39471, LA, color = Treatment,shape=Treatment)) +
  geom_point(alpha = 0.8,  size=10) + 
  scale_colour_manual(name = "Treatment",
                      labels = c("Well-watered", "Drought"),
                      values = c('#000000','#808080')) +   
  scale_shape_manual(values=c(19, 1))+

  stat_smooth(aes(fill = Treatment, color = Treatment,linetype=Treatment), method = "lm", formula = F,fullrange = TRUE, se = FALSE ) +
  scale_color_manual(values = c("Well-watered" = "black",
                                "Drought"="black")) +
  stat_regline_equation(
    label.x = c(8,8), label.y = c(0.23,0.11),
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula=F,size=6,show.legend = FALSE )+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.5,"cm")) +
  theme(axis.text = element_text(family = "Arial", color="black", face="bold", size=14)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=24)) +
  theme(plot.title = element_text(family = "Arial", color="black", face="bold", size=45)) +  theme(panel.background = element_rect(fill = 'white')) +
  ylab("Relative leaf area (rLA)") + 
  xlab(expression(bold(paste(bolditalic("AT5G39471")~"(FPKM)"))))+
  theme(legend.position = c(0.80, 0.92))+ ylim(0.1,1.75) +xlim(0,45)+
  theme(legend.key.size = unit(1, "cm"))+  theme(legend.position = c(0.1, 0.1))   

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
############################################################################################################################################################################


phenotypes=read_csv("data/phenotypesb.csv")
head(phenotypes)
phenotypesselected<- select(phenotypes, accession_name,LAP_well_watered,LAP_drought)
accession_id1<-phenotypes$accession_id
selected <-cbind(accession_id1, phenotypesselected)
selected <-selected %>% 
  rename( accession_id = accession_id1)
head(selected)
transposedfinal = read_csv("data/transposedfinal.csv")
accession_id2<-transposedfinal$accession_id

transposedfinal<- select(transposedfinal, AT1G68870)
transposedfinal <-cbind(accession_id2, transposedfinal)
transposedfinal <-transposedfinal %>% 
  rename( accession_id = accession_id2)
head(transposedfinal)

df3<-merge(selected, transposedfinal, by = "accession_id", sort = TRUE)
head(df3)

df3<- df3 %>%
  gather(Treatment, LA, LAP_well_watered:LAP_drought)

head (df3)

df3$Treatment <- as.factor(df3$Treatment)
head(df3)


df3$Treatment <- factor(df3$Treatment, levels=c("LAP_well_watered", "LAP_drought"))

library(mdthemes)
library(ggplot2)
library(ggpubr)
levels(df3$Treatment) <- c("Well-watered", "Drought")
F = as.formula("y~x")
p3 <-ggplot(df3, aes(AT1G68870, LA, color = Treatment,shape=Treatment)) +
  geom_point(alpha = 0.8,  size=10) + 
  scale_colour_manual(name = "Treatment",
                      labels = c("Well-watered", "Drought"),
                      values = c('#000000','#808080')) +   
  scale_shape_manual(values=c(19, 1))+
  
  stat_smooth(aes(fill = Treatment, color = Treatment,linetype=Treatment), method = "lm", formula = F,fullrange = TRUE, se = FALSE ) +
  scale_color_manual(values = c("Well-watered" = "black",
                                "Drought"="black")) +
  stat_regline_equation(
    label.x = c(120,120), label.y = c(0.23,0.11),
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula=F,size=6,show.legend = FALSE )+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.5,"cm")) +
  theme(axis.text = element_text(family = "Arial", color="black", face="bold", size=14)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=24)) +
  theme(plot.title = element_text(family = "Arial", color="black", face="bold", size=45)) +  theme(panel.background = element_rect(fill = 'white')) +
  ylab("Relative leaf area (rLA)") + 
  xlab(expression(bold(paste(bolditalic("SOFL2")~"(FPKM)"))))+
  theme(legend.position = c(0.80, 0.92))+ ylim(0.1,1.75) +xlim(0,700)+
  theme(legend.key.size = unit(1, "cm"))+  theme(legend.position = c(0.1, 0.1))   

p3
png("panels/p3.png", width = 8, height = 7, units = 'in', res = 350)
p3
dev.off() 

############################################################################################################################################################################
############################################################################################################################################################################
#########
###########################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################


phenotypes=read_csv("data/phenotypesb.csv")
head(phenotypes)
phenotypesselected<- select(phenotypes, accession_name,LAP_well_watered,LAP_drought)
accession_id1<-phenotypes$accession_id
selected <-cbind(accession_id1, phenotypesselected)
selected <-selected %>% 
  rename( accession_id = accession_id1)
head(selected)
transposedfinal = read_csv("data/transposedfinal.csv")
accession_id2<-transposedfinal$accession_id

transposedfinal<- select(transposedfinal, AT1G55490)
transposedfinal <-cbind(accession_id2, transposedfinal)
transposedfinal <-transposedfinal %>% 
  rename( accession_id = accession_id2)
head(transposedfinal)

df4<-merge(selected, transposedfinal, by = "accession_id", sort = TRUE)
head(df4)

df4<- df4 %>%
  gather(Treatment, LA, LAP_well_watered:LAP_drought)

head (df4)

df4$Treatment <- as.factor(df4$Treatment)
head(df4)


df4$Treatment <- factor(df4$Treatment, levels=c("LAP_well_watered", "LAP_drought"))

levels(df4$Treatment) <- c("Well-watered", "Drought")
F = as.formula("y~x")
p4 <-ggplot(df4, aes(AT1G55490, LA, color = Treatment,shape=Treatment)) +
  geom_point(alpha = 0.8,  size=10) + 
  scale_colour_manual(name = "Treatment",
                      labels = c("LAP_well_watered", "LAP_drought"),
                      values = c('#000000','#808080')) +   
  scale_shape_manual(values=c(19, 1))+
  
  stat_smooth(aes(fill = Treatment, color = Treatment,linetype=Treatment), method = "lm", formula = F,fullrange = TRUE, se = FALSE ) +
  scale_color_manual(values = c("Well-watered" = "black",
                                "Drought"="black")) +
  stat_regline_equation(
    label.x = c(11500,11500), label.y = c(0.23,0.11),
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula=F,size=6,show.legend = FALSE )+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.5,"cm")) +
  theme(axis.text = element_text(family = "Arial", color="black", face="bold", size=14)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=24)) +
  theme(plot.title = element_text(family = "Arial", color="black", face="bold", size=45)) +  theme(panel.background = element_rect(fill = 'white')) +
  ylab("Relative leaf area (rLA)") + 
  xlab(expression(bold(paste(bolditalic("CPN60B")~"(FPKM)"))))+
  theme(legend.position = c(0.80, 0.92))+ ylim(0.1,1.75) +#xlim(0,700)+
  theme(legend.key.size = unit(1, "cm"))+  theme(legend.position = c(0.1, 0.1))   

p4
png("panels/p4.png", width = 8, height = 7, units = 'in', res = 350)
p4
dev.off() 


############################################################################################################################################################################arrange panels
############################################################################################################################################################################
############################################################################################################################################################################
###########################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################


phenotypes=read_csv("data/phenotypesb.csv")
head(phenotypes)
phenotypesselected<- select(phenotypes, accession_name,LAP_well_watered,LAP_drought)
accession_id1<-phenotypes$accession_id
selected <-cbind(accession_id1, phenotypesselected)
selected <-selected %>% 
  rename( accession_id = accession_id1)
head(selected)
transposedfinal = read_csv("data/transposedfinal.csv")
accession_id2<-transposedfinal$accession_id

transposedfinal<- select(transposedfinal, AT5G43320)
transposedfinal <-cbind(accession_id2, transposedfinal)
transposedfinal <-transposedfinal %>% 
  rename( accession_id = accession_id2)
head(transposedfinal)

df5<-merge(selected, transposedfinal, by = "accession_id", sort = TRUE)
head(df5)

df5<- df5 %>%
  gather(Treatment, LA, LAP_well_watered:LAP_drought)

head (df5)

df5$Treatment <- as.factor(df5$Treatment)
head(df5)


df5$Treatment <- factor(df5$Treatment, levels=c("LAP_well_watered", "LAP_drought"))

levels(df5$Treatment) <- c("Well-watered", "Drought")
F = as.formula("y~x")
p5 <-ggplot(df5, aes(AT5G43320, LA, color = Treatment,shape=Treatment)) +
  geom_point(alpha = 0.8,  size=10) + 
  scale_colour_manual(name = "Treatment",
                      labels = c("LAP_well_watered", "LAP_drought"),
                      values = c('#000000','#808080')) +   
  scale_shape_manual(values=c(19, 1))+
  
  stat_smooth(aes(fill = Treatment, color = Treatment,linetype=Treatment), method = "lm", formula = F,fullrange = TRUE, se = FALSE ) +
  scale_color_manual(values = c("Well-watered" = "black",
                                "Drought"="black")) +
  stat_regline_equation(
    label.x = c(2400,2400), label.y = c(0.23,0.11),
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula=F,size=6,show.legend = FALSE )+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.5,"cm")) +
  theme(axis.text = element_text(family = "Arial", color="black", face="bold", size=14)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=24)) +
  theme(plot.title = element_text(family = "Arial", color="black", face="bold", size=45)) +  theme(panel.background = element_rect(fill = 'white')) +
  ylab("Relative leaf area (rLA)") + 
  xlab(expression(bold(paste(bolditalic("CKL8")~"(FPKM)"))))+
  theme(legend.position = c(0.80, 0.92))+ ylim(0.1,1.75) +xlim(2000,4500)+
  theme(legend.key.size = unit(1, "cm"))+  theme(legend.position = c(0.1, 0.1))   

p5
png("panels/p5.png", width = 8, height = 7, units = 'in', res = 350)
p5
dev.off() 



############################################################################################################################################################################
####
M2 <- image_read("panels/p2.png")
M3 <- image_read("panels/p3.png")
M4 <- image_read("panels/p4.png")
M5 <- image_read("panels/p5.png")

#####Figure Final version

plot_grid(rasterGrob(M2),rasterGrob(M3), rasterGrob(M4), rasterGrob(M5),labels=c('b',"c","d","e"), greedy = TRUE, ncol =2, nrow = 2, align = 'v', hjust=0, label_size=30)
png("panels/topright.png", width = 8, height = 7, units = 'in', res = 350)
plot_grid(rasterGrob(M2),rasterGrob(M3), rasterGrob(M4), rasterGrob(M5),labels=c('b',"c","d","e"), greedy = TRUE, ncol =2, nrow = 2, align = 'v', hjust=0, label_size=30)
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


phenotypes=read_csv("data/phenotypesb.csv")
head(phenotypes)
phenotypesselected<- select(phenotypes, accession_name,LAP_well_watered,LAP_drought)
accession_id1<-phenotypes$accession_id
selected <-cbind(accession_id1, phenotypesselected)
selected <-selected %>% 
  rename( accession_id = accession_id1)
head(selected)
transposedfinal = read_csv("data/transposedfinal.csv")
accession_id2<-transposedfinal$accession_id

transposedfinal<- select(transposedfinal, AT5G60160)
transposedfinal <-cbind(accession_id2, transposedfinal)
transposedfinal <-transposedfinal %>% 
  rename( accession_id = accession_id2)
head(transposedfinal)

df6<-merge(selected, transposedfinal, by = "accession_id", sort = TRUE)
head(df6)

df6<- df6 %>%
  gather(Treatment, LA, LAP_well_watered:LAP_drought)

head (df6)

df6$Treatment <- as.factor(df6$Treatment)
head(df6)


df6$Treatment <- factor(df6$Treatment, levels=c("LAP_well_watered", "LAP_drought"))
library(mdthemes)
library(ggplot2)
library(ggpubr)
levels(df6$Treatment) <- c("Well-watered", "Drought")
F = as.formula("y~x")
p6 <-ggplot(df6, aes(AT5G60160, LA, color = Treatment,shape=Treatment)) +
  geom_point(alpha = 0.4,  size=10) + 
  scale_colour_manual(name = "Treatment",
                      labels = c("Well-watered", "Drought"),
                      values = c('#7f7f7f','#7f7f7f')) +   
  scale_shape_manual(values=c(19, 1))+
  
  stat_smooth(aes(fill = Treatment, color = Treatment,linetype=Treatment), method = "lm", formula = F,fullrange = TRUE, se = FALSE ) +
  scale_color_manual(values = c("Well-watered" = "black",
                                "Drought"="black")) +
  stat_regline_equation(
    label.x = c(1100,1100), label.y = c(0.23,0.11),
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula=F,size=6,show.legend = FALSE )+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.5,"cm")) +
  theme(axis.text = element_text(family = "Arial", color="black", face="bold", size=14)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=24)) +
  theme(plot.title = element_text(family = "Arial", color="black", face="bold", size=45)) +  theme(panel.background = element_rect(fill = 'white')) +
  ylab("Relative leaf area (rLA)") + 
  xlab(expression(bold(paste(bolditalic("AT5G60160")~"(FPKM)"))))+
  theme(legend.position = c(0.80, 0.92))+ ylim(0.1,1.75) +xlim(1000,1600)+
  theme(legend.key.size = unit(1, "cm"))+  theme(legend.position = c(0.1, 0.1))   

p6
png("panels/p6.png", width = 8, height = 7, units = 'in', res = 350)
p6
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
############################################################################################################################################################################

phenotypes=read_csv("data/phenotypesb.csv")
head(phenotypes)
phenotypesselected<- select(phenotypes, accession_name,LAP_well_watered,LAP_drought)
accession_id1<-phenotypes$accession_id
selected <-cbind(accession_id1, phenotypesselected)
selected <-selected %>% 
  rename( accession_id = accession_id1)
head(selected)
transposedfinal = read_csv("data/transposedfinal.csv")
accession_id2<-transposedfinal$accession_id

transposedfinal<- select(transposedfinal, AT4G15260)
transposedfinal <-cbind(accession_id2, transposedfinal)
transposedfinal <-transposedfinal %>% 
  rename( accession_id = accession_id2)
head(transposedfinal)

df7<-merge(selected, transposedfinal, by = "accession_id", sort = TRUE)
head(df7)

df7<- df7 %>%
  gather(Treatment, LA, LAP_well_watered:LAP_drought)

head (df7)

df7$Treatment <- as.factor(df7$Treatment)
head(df7)


df7$Treatment <- factor(df7$Treatment, levels=c("LAP_well_watered", "LAP_drought"))

library(mdthemes)
library(ggplot2)
library(ggpubr)
levels(df7$Treatment) <- c("Well-watered", "Drought")
F = as.formula("y~x")
p7 <-ggplot(df7, aes(AT4G15260, LA, color = Treatment,shape=Treatment)) +
  geom_point(alpha = 0.4,  size=10) + 
  scale_colour_manual(name = "Treatment",
                      labels = c("Well-watered", "Drought"),
                      values = c('#7f7f7f','#7f7f7f')) +   
  scale_shape_manual(values=c(19, 1))+
  
  stat_smooth(aes(fill = Treatment, color = Treatment,linetype=Treatment), method = "lm", formula = F,fullrange = TRUE, se = FALSE ) +
  scale_color_manual(values = c("Well-watered" = "black",
                                "Drought"="black")) +
  stat_regline_equation(
    label.x = c(600,600), label.y = c(0.23,0.11),
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula=F,size=6,show.legend = FALSE )+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.5,"cm")) +
  theme(axis.text = element_text(family = "Arial", color="black", face="bold", size=14)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=24)) +
  theme(plot.title = element_text(family = "Arial", color="black", face="bold", size=45)) +  theme(panel.background = element_rect(fill = 'white')) +
  ylab("Relative leaf area (rLA)") + 
  xlab(expression(bold(paste(bolditalic("AT4G15260")~"(FPKM)"))))+
  theme(legend.position = c(0.80, 0.92))+ ylim(0.1,1.75)+# +xlim(0,700)+

  theme(legend.key.size = unit(1, "cm"))+  theme(legend.position = c(0.1, 0.1))   

p7
png("panels/p7.png", width = 8, height = 7, units = 'in', res = 350)
p7
dev.off() 

############################################################################################################################################################################
############################################################################################################################################################################
#########
###########################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################


phenotypes=read_csv("data/phenotypesb.csv")
head(phenotypes)
phenotypesselected<- select(phenotypes, accession_name,LAP_well_watered,LAP_drought)
accession_id1<-phenotypes$accession_id
selected <-cbind(accession_id1, phenotypesselected)
selected <-selected %>% 
  rename( accession_id = accession_id1)
head(selected)
transposedfinal = read_csv("data/transposedfinal.csv")
accession_id2<-transposedfinal$accession_id

transposedfinal<- select(transposedfinal, AT5G16800)
transposedfinal <-cbind(accession_id2, transposedfinal)
transposedfinal <-transposedfinal %>% 
  rename( accession_id = accession_id2)
head(transposedfinal)

df8<-merge(selected, transposedfinal, by = "accession_id", sort = TRUE)
head(df8)

df8<- df8 %>%
  gather(Treatment, LA, LAP_well_watered:LAP_drought)

head (df8)

df8$Treatment <- as.factor(df8$Treatment)
head(df8)


df8$Treatment <- factor(df8$Treatment, levels=c("LAP_well_watered", "LAP_drought"))

levels(df8$Treatment) <- c("Well-watered", "Drought")
F = as.formula("y~x")
p8 <-ggplot(df8, aes(AT5G16800, LA, color = Treatment,shape=Treatment)) +
  geom_point(alpha = 0.4,  size=10) + 
  scale_colour_manual(name = "Treatment",
                      labels = c("LAP_well_watered", "LAP_drought"),
                      values = c('#7f7f7f','#7f7f7f')) +   
  scale_shape_manual(values=c(19, 1))+
  
  stat_smooth(aes(fill = Treatment, color = Treatment,linetype=Treatment), method = "lm", formula = F,fullrange = TRUE, se = FALSE ) +
  scale_color_manual(values = c("Well-watered" = "black",
                                "Drought"="black")) +
  stat_regline_equation(
    label.x = c(675,675), label.y = c(0.23,0.11),
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula=F,size=6,show.legend = FALSE )+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.5,"cm")) +
  theme(axis.text = element_text(family = "Arial", color="black", face="bold", size=14)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=24)) +
  theme(plot.title = element_text(family = "Arial", color="black", face="bold", size=45)) +  theme(panel.background = element_rect(fill = 'white')) +
  ylab("Relative leaf area (rLA)") + 
  xlab(expression(bold(paste(bolditalic("AT5G16800")~"(FPKM)"))))+
  theme(legend.position = c(0.80, 0.92))+ ylim(0.1,1.75) +xlim(500,1500)+
  theme(legend.key.size = unit(1, "cm"))+  theme(legend.position = c(0.1, 0.1))   

p8
png("panels/p8.png", width = 8, height = 7, units = 'in', res = 350)
p8
dev.off() 


############################################################################################################################################################################arrange panels
############################################################################################################################################################################
############################################################################################################################################################################
###########################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################

phenotypes=read_csv("data/phenotypesb.csv")
head(phenotypes)
phenotypesselected<- select(phenotypes, accession_name,LAP_well_watered,LAP_drought)
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

df9<-merge(selected, transposedfinal, by = "accession_id", sort = TRUE)
head(df9)

df9<- df9 %>%
  gather(Treatment, LA, LAP_well_watered:LAP_drought)

head (df9)

df9$Treatment <- as.factor(df9$Treatment)
head(df9)


df9$Treatment <- factor(df9$Treatment, levels=c("LAP_well_watered", "LAP_drought"))

levels(df9$Treatment) <- c("Well-watered", "Drought")
F = as.formula("y~x")
p9 <-ggplot(df9, aes(AT4G14220, LA, color = Treatment,shape=Treatment)) +
  geom_point(alpha = 0.4,  size=10) + 
  scale_colour_manual(name = "Treatment",
                      labels = c("LAP_well_watered", "LAP_drought"),
                      values = c('#7f7f7f','#7f7f7f')) +   
  scale_shape_manual(values=c(19, 1))+
  
  stat_smooth(aes(fill = Treatment, color = Treatment,linetype=Treatment), method = "lm", formula = F,fullrange = TRUE, se = FALSE ) +
  scale_color_manual(values = c("Well-watered" = "black",
                                "Drought"="black")) +
  stat_regline_equation(
    label.x = c(410,410), label.y = c(0.23,0.11),
    aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~")),
    formula=F,size=6,show.legend = FALSE )+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.5,"cm")) +
  theme(axis.text = element_text(family = "Arial", color="black", face="bold", size=14)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=24)) +
  theme(plot.title = element_text(family = "Arial", color="black", face="bold", size=45)) +  theme(panel.background = element_rect(fill = 'white')) +
  ylab("Relative leaf area (rLA)") + 
  xlab(expression(bold(paste(bolditalic("RHF1A")~"(FPKM)"))))+
  theme(legend.position = c(0.80, 0.92))+ ylim(0.1,1.75)+ xlim(250,1200)+
  theme(legend.key.size = unit(1, "cm"))+  theme(legend.position = c(0.1, 0.1))   

p9
png("panels/p9.png", width = 8, height = 7, units = 'in', res = 350)
p9
dev.off() 



############################################################################################################################################################################
####
M6 <- image_read("panels/p6.png")
M7 <- image_read("panels/p7.png")
M8 <- image_read("panels/p8.png")
M9 <- image_read("panels/p9.png")

#####Figure Final version

plot_grid(rasterGrob(M6),rasterGrob(M7), rasterGrob(M8), rasterGrob(M9),labels=c('g',"h","i","j"), greedy = TRUE, ncol =2, nrow = 2, align = 'v', hjust=0, label_size=30)
png("panels/bottomright.png", width = 8, height = 7, units = 'in', res = 350)
plot_grid(rasterGrob(M6),rasterGrob(M7), rasterGrob(M8), rasterGrob(M9),labels=c('g',"h","i","j"), greedy = TRUE, ncol =2, nrow = 2, align = 'v', hjust=0, label_size=30)
dev.off()


F1 <- image_read("panels/p1.png")
F2 <- image_read("panels/topright.png")
F3 <- image_read("panels/bottomright.png")

plot_grid(rasterGrob(F2),rasterGrob(F3),labels=c('',""), greedy = TRUE, ncol =1, nrow = 2, align = 'v', hjust=0, label_size=30)
png("panels/right.png", width = 8, height = 7, units = 'in', res = 350)
plot_grid(rasterGrob(F2),rasterGrob(F3),labels=c('',""), greedy = TRUE, ncol =1, nrow = 2, align = 'v', hjust=0, label_size=30)
dev.off()

Right <- image_read("panels/right.png")


plot_grid(rasterGrob(F1),rasterGrob(Right),labels=c('a',""), greedy = TRUE, ncol =2, nrow = 1, align = 'v', hjust=0, label_size=20)
png("Figure 3.png", width = 15, height = 7, units = 'in', res = 350)
plot_grid(rasterGrob(F1),rasterGrob(Right),labels=c('a',""), greedy = TRUE, ncol =2, nrow = 1, align = 'v', hjust=0, label_size=20)
dev.off()
list.of.packages <- c("ggplot2", "viridis", "ggrepel","readr", "tidyverse","tidyr","cowplot","qvalue","dplyr","ggpmisc","ggpubr","ggthemes", "Hmisc", "cowplot","magick","grid")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(readr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(ggpubr)
library(Hmisc)
library(cowplot)
library(dplyr)
library(ggpmisc)
library(magick)
library(grid)

phenotypes=read_csv("data/phenotypes.csv")
head(phenotypes)
phenotypesselected<- select(phenotypes, name,LAPdSI_drought)
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

df1<-merge(selected, transposedfinal, by = "accession_id", sort = TRUE)
head(df1)
formula <- y ~ x

head(df1)
###Plot2= d13C vs leaf area plasticity in response to drought
p1 <- ggplot(df1, aes(x=AT3G48010 , y= LAPdSI_drought)) +
  geom_point(color='black',alpha = 0.5,  size=12)+   
  
  geom_smooth(method = "lm", se = T, fill="#696969", colour="#696969", size=0.7, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                       label.x.npc = "left", label.y.npc = "top",
                                                                                                       formula = formula, parse = TRUE, size = 5, colour="#696969") +
  geom_text_repel(aes(AT3G48010, LAPdSI_drought, label = name), fontface = 'bold', size=5, force=0.3,segment.size = 1,family = 'Arial',
                  fontface = 'bold',
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.35, "lines"),
                  segment.color = 'grey50') + 
  
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=35)) +  
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm"))  +
  theme(axis.text = element_text(family = "Arial", color="black",  face="bold",size=16)) +
  theme(plot.title = element_text(family = "Arial", color="black", face="bold",size=28)) +  theme(panel.background = element_rect(fill = 'white')) +
  theme(legend.position="none") + 
  ylab(expression(bold("Leaf area plasticity index")))+ 
  xlab(expression(bold(paste(bolditalic("CNGC16")~"(FPKM)"))))
#ggtitle("Phenotype x phenotype")
# scale_y_continuous(limits = c(-0.75, 0.75), breaks = seq(-0.75, 0.75, by = 0.25))+
#scale_x_continuous(limits = c(-37, -33), breaks = seq(-37, -33, by = 1))

p1
png("panels/p1.png", width = 10, height = 7, units = 'in', res = 350)
p1
dev.off() 


###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
phenotypesselected<- select(phenotypes, name,Dittberner_moleco18_Delta_13C)
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

df2<-merge(selected, transposedfinal, by = "accession_id", sort = TRUE)
head(df1)
formula <- y ~ x
p2 <- ggplot(df2, aes(x=AT3G48010 , y= Dittberner_moleco18_Delta_13C)) +
  geom_point(color='black',alpha = 0.5,  size=12)+   
  
  geom_smooth(method = "lm", se = T, fill="#696969", colour="#696969", size=0.7, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                       label.x.npc = "left", label.y.npc = "top",
                                                                                                       formula = formula, parse = TRUE, size = 5, colour="#696969") +
  geom_text_repel(aes(AT3G48010, Dittberner_moleco18_Delta_13C, label = name), fontface = 'bold', size=5, force=0.3,segment.size = 1,family = 'Arial',
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.35, "lines"),
                  segment.color = 'grey50') + 
  
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=35)) +  
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm"))  +
  theme(axis.text = element_text(family = "Arial", color="black",  face="bold",size=16)) +
  theme(plot.title = element_text(family = "Arial", color="black", face="bold",size=28)) +  theme(panel.background = element_rect(fill = 'white')) +
  theme(legend.position="none") + 
  ylab(expression(bold(paste("WUE"[i]," (\u03B4"^"13", "C (\u2030))"))))+   
  xlab(expression(bold(paste(bolditalic("CNGC16")~"(FPKM)"))))

  #ggtitle("Phenotype x phenotype")
# scale_y_continuous(limits = c(-0.75, 0.75), breaks = seq(-0.75, 0.75, by = 0.25))+
#scale_x_continuous(limits = c(-37, -33), breaks = seq(-37, -33, by = 1))

p2
png("panels/p2.png", width = 10, height = 7, units = 'in', res = 350)
p2
dev.off() 

##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################

phenotypes=read_csv("data/phenotypes.csv")
head(phenotypes)
phenotypesselected<- select(phenotypes, name,LAPdSI_drought)
accession_id1<-phenotypes$accession_id
selected <-cbind(accession_id1, phenotypesselected)
selected <-selected %>% 
  rename( accession_id = accession_id1)
head(selected)
transposedfinal = read_csv("data/transposedfinal.csv")
accession_id2<-transposedfinal$accession_id

transposedfinal<- select(transposedfinal, AT3G47950)
transposedfinal <-cbind(accession_id2, transposedfinal)
transposedfinal <-transposedfinal %>% 
  rename( accession_id = accession_id2)
head(transposedfinal)

df3<-merge(selected, transposedfinal, by = "accession_id", sort = TRUE)
head(df3)
formula <- y ~ x

head(df3)
###Plot2= d13C vs leaf area plasticity in response to drought
p3 <- ggplot(df3, aes(x=AT3G47950 , y= LAPdSI_drought)) +
  geom_point(color='black',alpha = 0.5,  size=12)+   
  
  geom_smooth(method = "lm", se = T, fill="#696969", colour="#696969", size=0.7, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                       label.x.npc = "right", label.y.npc = "top",
                                                                                                       formula = formula, parse = TRUE, size = 5, colour="#696969") +
  geom_text_repel(aes(AT3G47950, LAPdSI_drought, label = name), fontface = 'bold', size=5, force=0.3,segment.size = 1,family = 'Arial',
                  fontface = 'bold',
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.35, "lines"),
                  segment.color = 'grey50') + 
  
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=35)) +  
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm"))  +
  theme(axis.text = element_text(family = "Arial", color="black",  face="bold",size=16)) +
  theme(plot.title = element_text(family = "Arial", color="black", face="bold",size=28)) +  theme(panel.background = element_rect(fill = 'white')) +
  theme(legend.position="none") + 
  ylab(expression(bold("Leaf area plasticity index")))+ 
  xlab(expression(bold(paste(bolditalic("AHA4")~"(FPKM)"))))
#ggtitle("Phenotype x phenotype")
# scale_y_continuous(limits = c(-0.75, 0.75), breaks = seq(-0.75, 0.75, by = 0.25))+
#scale_x_continuous(limits = c(-37, -33), breaks = seq(-37, -33, by = 1))

p3
png("panels/p3.png", width = 10, height = 7, units = 'in', res = 350)
p3
dev.off() 


###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
phenotypesselected<- select(phenotypes, name,Dittberner_moleco18_Delta_13C)
accession_id1<-phenotypes$accession_id
selected <-cbind(accession_id1, phenotypesselected)
selected <-selected %>% 
  rename( accession_id = accession_id1)
head(selected)
transposedfinal = read_csv("data/transposedfinal.csv")
accession_id2<-transposedfinal$accession_id

transposedfinal<- select(transposedfinal, AT3G47950)
transposedfinal <-cbind(accession_id2, transposedfinal)
transposedfinal <-transposedfinal %>% 
  rename( accession_id = accession_id2)
head(transposedfinal)

df4<-merge(selected, transposedfinal, by = "accession_id", sort = TRUE)
head(df4)
formula <- y ~ x
p4 <- ggplot(df4, aes(x=AT3G47950 , y= Dittberner_moleco18_Delta_13C)) +
  geom_point(color='black',alpha = 0.5,  size=12)+   
  
  geom_smooth(method = "lm", se = T, fill="#696969", colour="#696969", size=0.7, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                       label.x.npc = "right", label.y.npc = "top",
                                                                                                       formula = formula, parse = TRUE, size = 5, colour="#696969") +
  geom_text_repel(aes(AT3G47950, Dittberner_moleco18_Delta_13C, label = name), fontface = 'bold', size=5, force=0.3,segment.size = 1,family = 'Arial',
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.35, "lines"),
                  segment.color = 'grey50') + 
  
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=35)) +  
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm"))  +
  theme(axis.text = element_text(family = "Arial", color="black",  face="bold",size=16)) +
  theme(plot.title = element_text(family = "Arial", color="black", face="bold",size=28)) +  theme(panel.background = element_rect(fill = 'white')) +
  theme(legend.position="none") + 
  xlab(expression(bold(paste(bolditalic("AHA4")~"(FPKM)"))))+
  ylab(expression(bold(paste("WUE"[i]," (\u03B4"^"13", "C (\u2030))"))))+   
  xlim(limits = c(20, 140))

#ggtitle("Phenotype x phenotype")
# scale_y_continuous(limits = c(-0.75, 0.75), breaks = seq(-0.75, 0.75, by = 0.25))+

p4
png("panels/p4.png", width = 10, height = 7, units = 'in', res = 350)
p4
dev.off() 


M1 <- image_read("panels/p1.png")
M2 <- image_read("panels/p2.png")
M3 <- image_read("panels/p3.png")
M4 <- image_read("panels/p4.png")


plot_grid(rasterGrob(M1),rasterGrob(M3),rasterGrob(M2), rasterGrob(M4),labels=c('A', 'C','B', 'D'), greedy = TRUE, ncol =2, nrow = 2, label_size=35)
png("Figure 4.png", width = 19, height = 12, units = 'in', res = 350)
plot_grid(rasterGrob(M1),rasterGrob(M3),rasterGrob(M2), rasterGrob(M4),labels=c('A', 'C','B', 'D'), greedy = TRUE, ncol =2, nrow = 2,   label_size=35)
dev.off()


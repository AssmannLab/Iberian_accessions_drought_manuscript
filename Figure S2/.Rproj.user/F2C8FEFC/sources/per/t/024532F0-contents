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
library(readr)
library(dplyr)
library(tidyr)
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
####Create an annotation file compatible with the files we are creating.
annotation = read_csv("data/annotation_GWAS.csv")
nrow(annotation)
head(annotation)
gene_description = read_csv("data/gene_description.csv")
head(gene_description)
annotation <- merge(annotation,gene_description,by= "gene")
head(annotation)
df1a=read_csv("data/Linear model/PI_LM.csv")
head(df1a)
df1b=read_csv("data/AMM Mixed model/PI_AMM.csv")
df1b <- df1b[order(df1b$score) ,  ]
head(df1b)
df1 <- merge(df1a,df1b,by=c("chr", "pos"))
df1 <- merge(df1,annotation,by=c("chr", "pos"))
head(df1)

df1 <- df1 %>% 
  rename(
    Leaf_Area_Plasticity_Index_LM_score = score.x,
    Leaf_Area_Plasticity_Index_LM_MAF = maf.x,
    Leaf_Area_Plasticity_Index_LM_MAC = mac.x,
    Leaf_Area_Plasticity_Index_AMM_score = score.y,
    Leaf_Area_Plasticity_Index_AMM_MAF = maf.y,
    Leaf_Area_Plasticity_Index_AMM_MAC = mac.y, annotation = annotation.y, symbol= symbol.y )
df1$GVE.x <- NULL
df1$GVE.y <- NULL
df1$symbol.x <- NULL
df1$annotation.x <- NULL

head(df1)

df1<-df1 %>% unite(label, symbol, pos, sep = " pos. ")
df1 <- df1[order(df1$Leaf_Area_Plasticity_Index_LM_score) ,  ]
head(df1)
df1 <- df1[df1$Leaf_Area_Plasticity_Index_LM_MAF > 0.1, ]

write.csv(df1, "Tables/GWAS_Leaf_area_plasticity_Linear_vs_AMM_models.csv", row.names = TRUE)


p1<-df1 %>% 
  mutate(quadrant = case_when(Leaf_Area_Plasticity_Index_LM_score > 4 & Leaf_Area_Plasticity_Index_AMM_score > 4   ~ "Q1",
                              Leaf_Area_Plasticity_Index_LM_score <= 4 & Leaf_Area_Plasticity_Index_AMM_score > 4  ~ "Q2",
                              Leaf_Area_Plasticity_Index_LM_score <= 4 & Leaf_Area_Plasticity_Index_AMM_score <= 4 ~ "Q3",
                              TRUE                                         ~ "Q4")) %>% 
  mutate(label = case_when(Leaf_Area_Plasticity_Index_LM_score >= 4 & Leaf_Area_Plasticity_Index_AMM_score >= 4 ~  paste(label))) %>% 
  ggplot(aes(x = Leaf_Area_Plasticity_Index_LM_score, y = Leaf_Area_Plasticity_Index_AMM_score)) + 
  geom_point(aes(color=quadrant),alpha = 0.5,  size=3) + 
  geom_text_repel(aes(Leaf_Area_Plasticity_Index_LM_score, Leaf_Area_Plasticity_Index_AMM_score, label = label),  size=2, force=0.1,segment.size = 0.5,arrow = arrow(length = unit(0.01, 'npc')),family = 'Arial',
                  fontface = 'bold.italic',
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.5, "lines"),
                  segment.color = 'grey50',max.overlaps = Inf) + 
  scale_colour_manual(values = c( "#47a569", "#6689c3", "#c5c5c5","#e9af41"))+
  geom_vline(xintercept = 4, color = "black", linetype = "dashed")+
  geom_hline(yintercept = 4, color = "black", linetype = "dashed")+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm"))  +
  theme(axis.text = element_text(family = "Arial", color="black", size=16)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=13)) +  
  theme(plot.title = element_text(family = "Arial", color="#696969",face="bold",  size=20))  +  theme(panel.background = element_rect(fill = 'white')) +
  theme(legend.position="none") +
  scale_x_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2))+
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2))+
  ggtitle("Linear vs. Accelerated Mixed Model")+
  xlab(expression(bold(atop(paste("Leaf Area Plasticity Index: Linear Model"), paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P")))))))+ 
  ylab(expression(bold(atop(paste("Leaf Area Plasticity Index: Accelerated Mixed Model"), paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P")))))))
p1
png("panels/p1.png", width = 7, height = 7, units = 'in', res = 350)
p1
dev.off() 


##################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

####Create an annotation file compatible with the files we are creating.

df2a=read_csv("data/Linear model/Iberian_d13C_LM.csv")
head(df2a)
head(df2a)
df2b=read_csv("data/AMM Mixed model/Iberian_d13C_AMM.csv")
df2b <- df2b[order(df2b$score) ,  ]
head(df2b)
df2 <- merge(df2a,df2b,by=c("chr", "pos"))
df2 <- merge(df2,annotation,by=c("chr", "pos"))
head(df2)

df2 <- df2 %>% 
  rename(
    Iberian_d13C_LM_score = score.x,
    Iberian_d13C_LM_MAF = maf.x,
    Iberian_d13C_LM_MAC = mac.x,
    Iberian_d13C_AMM_score = score.y,
    Iberian_d13C_AMM_MAF = maf.y,
    Leaf_Area_Plasticity_Index_AMM_MAC = mac.y)
df2$GVE.x <- NULL
df2$GVE.y <- NULL
df2$annotation.y <- NULL


head(df2)


df2<-df2 %>% unite(label, symbol, pos, sep = " pos. ")
df2 <- df2[order(df2$Iberian_d13C_LM_score) ,  ]
head(df2)
df2 <- df2[df2$Iberian_d13C_LM_MAF > 0.1, ]

write.csv(df2, "Tables/GWAS_Iberian_d13C_Linear_vs_AMM_models.csv", row.names = TRUE)


p2<-df2 %>% 
  mutate(quadrant = case_when(Iberian_d13C_LM_score > 4 & Iberian_d13C_AMM_score > 4   ~ "Q1",
                              Iberian_d13C_LM_score <= 4 & Iberian_d13C_AMM_score > 4  ~ "Q2",
                              Iberian_d13C_LM_score <= 4 & Iberian_d13C_AMM_score <= 4 ~ "Q3",
                              TRUE                                         ~ "Q4")) %>% 
  mutate(label = case_when(Iberian_d13C_LM_score >= 4 & Iberian_d13C_AMM_score >= 4 ~  paste(label))) %>% 
  ggplot(aes(x = Iberian_d13C_LM_score, y = Iberian_d13C_AMM_score)) + 
  geom_point(aes(color=quadrant),alpha = 0.5,  size=3) + 
  geom_text_repel(aes(Iberian_d13C_LM_score, Iberian_d13C_AMM_score, label = label),  size=1.5, force=1,segment.size = 0.5,arrow = arrow(length = unit(0.01, 'npc')),family = 'Arial',
                  fontface = 'bold.italic',
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.5, "lines"),
                  segment.color = 'grey50',max.overlaps = Inf) + 
  scale_colour_manual(values = c( "#47a569",  "#c5c5c5","#e9af41"))+
  geom_vline(xintercept = 4, color = "black", linetype = "dashed")+
  geom_hline(yintercept = 4, color = "black", linetype = "dashed")+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm"))  +
  theme(axis.text = element_text(family = "Arial", color="black", size=16)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=13)) +  
  theme(plot.title = element_text(family = "Arial", color="#696969",face="bold",  size=20))  +  theme(panel.background = element_rect(fill = 'white')) +
  theme(legend.position="none") +
  scale_x_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2))+
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2))+
  ggtitle("Linear vs. Accelerated Mixed Model")+
  xlab(expression(bold(atop(paste("Iberian iWUE"," (\u03B4"^"13", "C): Linear Model"), paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P")))))))+ 
  ylab(expression(bold(atop(paste("Iberian iWUE"," (\u03B4"^"13", "C): Accelerated Mixed Model"), paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P"))))))) 

p2
png("panels/p2.png", width = 7, height = 7, units = 'in', res = 350)
p2
dev.off() 


##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################


df3a=read_csv("data/Linear model/PI_LM.csv")
head(df3a)
head(df3a)
df3b=read_csv("data/GWAS_flow/GWAS_Flow_PI.csv")
head(df3b)
df3b$score <- (-log10(df3b$pval))
df3b <- df3b[order(-df3b$score) ,  ]
head(df3b)
df3 <- merge(df3a,df3b,by=c("chr", "pos"), all = TRUE)
df3 <- merge(df3,annotation,by=c("chr", "pos"))
head(df3)
df3 <- df3 %>% 
  rename(
    Leaf_Area_Plasticity_Index_LM_score = score.x,
    Leaf_Area_Plasticity_Index_LM_MAF = maf,
    Leaf_Area_Plasticity_Index_LM_MAC = mac.x,
    Leaf_Area_Plasticity_Index_GWAS_Flow_score = score.y,
    pval_GWAS_Flow = pval,
    Leaf_Area_Plasticity_Index_GWAS_Flow_MAC = mac.y)
df3$GVE <- NULL
df3$symbol.x <- NULL
df3$annotation.x <- NULL
df3 <- df3[df3$Leaf_Area_Plasticity_Index_LM_MAF > 0.1, ]

head(df3)


df3<-df3 %>% unite(label, symbol, pos, sep = " pos. ")
df3 <- df3[order(-df3$Leaf_Area_Plasticity_Index_LM_score) ,  ]
head(df3)

write.csv(df3, "Tables/GWAS_lAP_plasticity_GWAS_Flow_vs_Linear_models.csv", row.names = TRUE)


p3<-df3 %>% 
  mutate(quadrant = case_when(Leaf_Area_Plasticity_Index_LM_score > 4 & Leaf_Area_Plasticity_Index_GWAS_Flow_score > 4   ~ "Q1",
                              Leaf_Area_Plasticity_Index_LM_score <= 4 & Leaf_Area_Plasticity_Index_GWAS_Flow_score > 4  ~ "Q2",
                              Leaf_Area_Plasticity_Index_LM_score <= 4 & Leaf_Area_Plasticity_Index_GWAS_Flow_score <= 4 ~ "Q3",
                              TRUE                                         ~ "Q4")) %>% 
  mutate(label = case_when(Leaf_Area_Plasticity_Index_LM_score >= 4 & Leaf_Area_Plasticity_Index_GWAS_Flow_score >= 4 ~  paste(label))) %>% 
  ggplot(aes(x = Leaf_Area_Plasticity_Index_LM_score, y = Leaf_Area_Plasticity_Index_GWAS_Flow_score)) + 
  geom_point(aes(color=quadrant),alpha = 0.5,  size=3) + 
  geom_text_repel(aes(Leaf_Area_Plasticity_Index_LM_score, Leaf_Area_Plasticity_Index_GWAS_Flow_score, label = label),  size=1.5, force=1,segment.size = 0.5,arrow = arrow(length = unit(0.01, 'npc')),family = 'Arial',
                  fontface = 'bold.italic',
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.5, "lines"),
                  segment.color = 'grey50',max.overlaps = Inf) + 
  scale_colour_manual(values = c( "#47a569",  "#c5c5c5","#e9af41"))+
  geom_vline(xintercept = 4, color = "black", linetype = "dashed")+
  geom_hline(yintercept = 4, color = "black", linetype = "dashed")+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm"))  +
  theme(axis.text = element_text(family = "Arial", color="black", size=16)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=13)) +  
  theme(plot.title = element_text(family = "Arial", color="#696969",face="bold",  size=20))  +  theme(panel.background = element_rect(fill = 'white')) +
  theme(legend.position="none") +
  scale_x_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2))+
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2))+
  ggtitle("Linear vs.GWAS-Flow")+
  ylab(expression(bold(atop(paste("Leaf area plasticity: GWAS-Flow"), paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P")))))))+ 
  xlab(expression(bold(atop(paste("Leaf Area Plasticity Index: Linear Model"), paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P")))))))
p3
png("panels/p3.png", width = 7, height = 7, units = 'in', res = 350)
p3
dev.off() 


##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################


df4a=read_csv("data/Linear model/Iberian_d13C_LM.csv")
head(df4a)
head(df4a)
df4b=read_csv("data/GWAS_flow/1001pop_Delta_13C_GWAS_Flow.csv")
head(df4b)
df4b$score <- (-log10(df4b$pval))
df4b <- df4b[order(-df4b$score) ,  ]
head(df4b)
df4 <- merge(df4a,df4b,by=c("chr", "pos"), all = TRUE)
df4 <- merge(df4,annotation,by=c("chr", "pos"))
head(df4)

df4 <- df4 %>% 
  rename(
    Iberian_d13C_LM_score = score.x,
    Iberian_d13C_LM_MAF = maf,
    Iberian_d13C_LM_MAC = mac.x,
    Iberian_d13C_GWAS_Flow_score = score.y,
    pval_GWAS_Flow = pval,
    Iberian_d13C_GWAS_Flow_MAC = mac.y)
df4$GVE <- NULL
df4$symbol.x <- NULL
df4$annotation.x <- NULL
head(df4)


df4<-df4 %>% unite(label, symbol, pos, sep = " pos. ")
df4 <- df4[order(-df4$Iberian_d13C_LM_score) ,  ]
df4 <- df4[df4$Iberian_d13C_LM_MAF > 0.1, ]

head(df4)

write.csv(df4, "Tables/GWAS_1001G_d13C_GWAS_Flow_vs_AMM_models.csv", row.names = TRUE)


p4<-df4 %>% 
  mutate(quadrant = case_when(Iberian_d13C_LM_score > 4 & Iberian_d13C_GWAS_Flow_score > 4   ~ "Q1",
                              Iberian_d13C_LM_score <= 4 & Iberian_d13C_GWAS_Flow_score > 4  ~ "Q2",
                              Iberian_d13C_LM_score <= 4 & Iberian_d13C_GWAS_Flow_score <= 4 ~ "Q3",
                              TRUE                                         ~ "Q4")) %>% 
  mutate(label = case_when(Iberian_d13C_LM_score >= 4 & Iberian_d13C_GWAS_Flow_score >= 4 ~  paste(label))) %>% 
  ggplot(aes(x = Iberian_d13C_LM_score, y = Iberian_d13C_GWAS_Flow_score)) + 
  geom_point(aes(color=quadrant),alpha = 0.5,  size=3) + 
  geom_text_repel(aes(Iberian_d13C_LM_score, Iberian_d13C_GWAS_Flow_score, label = label),  size=1.5, force=1,segment.size = 0.5,arrow = arrow(length = unit(0.01, 'npc')),family = 'Arial',
                  fontface = 'bold.italic',
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.5, "lines"),
                  segment.color = 'grey50',max.overlaps = Inf) + 
  scale_colour_manual(values = c( "#47a569",  "#c5c5c5","#e9af41"))+
  geom_vline(xintercept = 4, color = "black", linetype = "dashed")+
  geom_hline(yintercept = 4, color = "black", linetype = "dashed")+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm"))  +
  theme(axis.text = element_text(family = "Arial", color="black", size=16)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=13)) +  
  theme(plot.title = element_text(family = "Arial", color="#696969",face="bold",  size=20))  +  theme(panel.background = element_rect(fill = 'white')) +
  theme(legend.position="none") +
  scale_x_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2))+
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2))+
  ggtitle("Linear vs.GWAS-Flow")+
ylab(expression(bold(atop(paste("Iberian iWUE"," (\u03B4"^"13", "C): GWAS-Flow"), paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P")))))))+ 
xlab(expression(bold(atop(paste("Iberian iWUE"," (\u03B4"^"13", "C): Linear Model"), paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P"))))))) 
  

p4
png("panels/p4.png", width = 7, height = 7, units = 'in', res = 350)
p4
dev.off() 






M1 <- image_read("panels/p1.png")
M2 <- image_read("panels/p2.png")
M3 <- image_read("panels/p3.png")
M4 <- image_read("panels/p4.png")



plot_grid(rasterGrob(M1),rasterGrob(M2),rasterGrob(M3),rasterGrob(M4),labels=c('a', 'b','c', 'd'), greedy = TRUE, ncol =2, nrow = 2, align = 'v', hjust=0, label_size=30)
png("Figure.S2.png", width = 12, height = 12, units = 'in', res = 350)
plot_grid(rasterGrob(M1),rasterGrob(M2),rasterGrob(M3),rasterGrob(M4),labels=c('a', 'b','c', 'd'), greedy = TRUE, ncol =2, nrow = 2, align = 'v', hjust=0, label_size=30)
dev.off()



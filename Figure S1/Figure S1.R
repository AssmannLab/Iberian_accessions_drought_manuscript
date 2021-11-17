list.of.packages <- c("ggplot2", "readr","cowplot","magick","png","pdftools","grid","ggrepel","tidyverse","ggpubr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ggplot2)
library(readr)
library(cowplot)
library(magick)
library(png)
library(pdftools)
library(grid)
library(ggrepel)
library(tidyverse)
library(ggpubr)

#### Fig. 1A can be found on the "panels" folder ("p1.jpg")

######################################################################################################## #Figures 1B & 1C
#########################################################################################################
#########################################################################################################
#########################################################################################################
## Load leaf area phenotypes
## LA_SD_WW = Leaf area under well-watered conditions (cm2)
## LA_SD_LW = Leaf area under drought (cm2)

data<-read.csv("data/phenotypes.csv")
head(data)


######Plot barplots for well-watered conditions (Fig. 1B)
p3 <- ggplot(subset(data, !is.na(LA_SD_WW))) +
  geom_bar( aes(x=accession_name, y=LA_SD_WW), stat="identity", fill="#008764", alpha=0.7) +
  geom_errorbar( aes(x=accession_name, ymin=LA_SD_WW-LA_SD_WW_SE, ymax=LA_SD_WW+LA_SD_WW_SE), width=0.4, colour="black", alpha=0.9, size=1.3) +
  # scale_y_continuous(breaks = c(5,10,15,20), limits = c(0,300)) +
  scale_y_continuous( limits = c(0,200)) +
  theme_bw() +
  theme(axis.line = element_line( colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.25,"cm")) +
  theme(panel.background=element_rect(fill='white')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text = element_text(family = "sans", color="#666666", face="bold", size=18)) +
  theme(axis.title = element_text(family = "sans", color="#666666", face="bold", size=30)) +
    theme(plot.title = element_text(family = "Arial", color="#008764", face="bold", size=55)) +  theme(panel.background = element_rect(fill = 'white')) +
  ggtitle("Well-watered")+
  ylab(expression(bold("Leaf area"~cm^2)))+ 
  xlab(expression(bold("Accession")))
p3

######Plot barplots for well-watered conditions (Fig. 1B)
p4 <- ggplot(subset(data, !is.na(LA_SD_LW))) +
  geom_bar( aes(x=accession_name, y=LA_SD_LW), stat="identity", fill="#9f4c49", alpha=0.7) +
  geom_errorbar( aes(x=accession_name, ymin=LA_SD_LW-LA_SD_LW_SE, ymax=LA_SD_LW+LA_SD_LW_SE), width=0.4, colour="black", alpha=0.9, size=1.3) +
  # scale_y_continuous(breaks = c(5,10,15,20), limits = c(0,300)) +
  scale_y_continuous( limits = c(0,200)) +
  theme_bw() +
  theme(axis.line = element_line( colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.25,"cm")) +
  theme(panel.background=element_rect(fill='white')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text = element_text(family = "sans", color="#666666", face="bold", size=18)) +
  theme(axis.title = element_text(family = "sans", color="#666666", face="bold", size=30)) +
  theme(plot.title = element_text(family = "Arial", color="#9f4c49", face="bold", size=55)) +  theme(panel.background = element_rect(fill = 'white')) +
  ggtitle("Drought")+
  ylab(expression(bold("Leaf area"~cm^2)))+ 
  xlab(expression(bold("Accession")))
p4

########################################################################################################
#########################################################################################################
#########################################################################################################
#######################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
plot_grid(p3,p4, labels=c("a", "b"), greedy = TRUE, ncol = 1, nrow =2, align = 'v', hjust=0, label_size=65)
png("Figure S1.png", width = 25, height = 24, units = 'in', res = 350)
plot_grid(p3,p4, labels=c("a", "b"), greedy = TRUE, ncol = 1, nrow =2, align = 'v', hjust=0, label_size=65)
dev.off()
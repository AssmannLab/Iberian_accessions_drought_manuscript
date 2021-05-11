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

############################################################################################################# Fig. 1D can be found on the "panels" folder ("p2b.pdf.jpg"). This image has been modified from the original that can be seen in a keynote file ("cartoon rosettes.key")

######################################################################################################### Reaction norm figures (Fig. 1E & 1F)
#########################################################################################################
#########################################################################################################
#########################################################################################################

## Load leaf area  phenotypes. This file can be found the the "data" folder. 
## LA_SD_WW = Leaf area under well-watered conditions (cm2)
## LA_SD_LW = Leaf area under drought (cm2)

###Open dataframe with phenotypes
data = read_csv("data/phenotypes.csv")
df5a <- data[c("accession_id", "accession_name", "LA_SD_WW", "LA_SD_LW")]
head(df5a)


######Plot reaction norms for leaf area (well-watered vs. drought) (Fig. 1E)

head(df5a)

df5a<-df5a[!is.na(df5a$LA_SD_LW), ]
df5a<- df5a %>%
  gather(TRT, phenotype, LA_SD_LW:LA_SD_WW)
df5a<-df5a[complete.cases(df5a), ]

head(df5a)

my_comparisons <- list( c("LA_SD_WW", "LA_SD_LW") )

data_ends <- df5a %>% filter(TRT == "LA_SD_WW")

p5a <- ggplot(df5a, aes(factor(TRT), phenotype, fill=TRT)) + geom_violin(alpha=0.5,color = NA) +
  scale_fill_manual(values=c("#9f4c49", "#008764")) + scale_x_discrete(limits=c("LA_SD_WW","LA_SD_LW"), labels = c('Well-watered','Drought'), expand = c(0.1, 0.1)) +
  geom_point(color= "black",size=5,
             alpha =0.6,
             position=position_jitter(width=0.005)) +    
  geom_text_repel(aes(label = accession_name), data = data_ends,
                  fontface ="bold", color = "black", size = 7, box.padding = 1, force = 0.5, segment.colour="#575d5d",
                  nudge_x      = -0.30,
                  direction    = "y",
                  hjust        = 0,
                  segment.size = 0.2,max.overlaps = Inf
  ) +
  geom_line(aes(group=accession_id), color="black",size = 1.2,
            alpha = 0.7)+ 
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.25,"cm")) +
  theme(axis.text.x = element_text(family = "Arial", color="black", face="bold", size=35)) +
  theme(axis.text.y = element_text(family = "Arial", color="black", face="bold", size=15)) +
  theme(axis.title.x = element_text(family = "Arial", color="black", face="bold", size=40)) +
  theme(axis.title.y = element_text(family = "Arial", color="black", face="bold", size=40)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = - 5, family = "Arial", color="black", face="bold", size=50)) +
  theme(axis.title.x=element_blank())+ 
  theme(plot.margin=unit(c(1,1,1,1),"cm"))+
  guides(fill = FALSE) + 
  ylab(expression(bold(paste("Leaf area (", cm^2,")"))))+
  xlab("Treatment")+  ggtitle("Leaf area: reaction norms")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") 
p5a


######Plot eaction norms for leaf area potentials (well-watered vs. drought) (Fig. 1E). In this plot and further analysis, we remove the values of non-Iberian Col-0, which we used as an internal reference, as well as the the values for accession IP-Coc1 (accession id = 9535), which were consistently unrealistically influential across different phenotypes based on their extreme Cook’s distance. 

######Leaf area reaction norm
## Load leaf area  phenotypes. This file can be found the the "data" folder. 
## LA_SD_WW = Leaf area under well-watered conditions (cm2)
## LA_SD_LW = Leaf area under drought (cm2)

###Open dataframe with phenotypes
df5b = read_csv("data/phenotypic_potential.csv")
df5b <- df5b[c("accession_id", "name", "LAP_SD_WW", "LAP_SD_LW")]
head(df5b)

df5b<-df5b[!is.na(df5b$LAP_SD_LW), ]
df5b<- df5b %>%
  gather(TRT, phenotype, LAP_SD_LW:LAP_SD_WW)
df5b<-df5b[complete.cases(df5), ]

head(df5b)

my_comparisons <- list( c("LAP_SD_WW", "LAP_SD_LW") )

data_ends <- df5b %>% filter(TRT == "LAP_SD_WW")
library(ggpubr)

p5b <- ggplot(df5b, aes(factor(TRT), phenotype, fill=TRT)) + geom_violin(alpha=0.5,color = NA) +
  scale_fill_manual(values=c("#9f4c49", "#008764")) + scale_x_discrete(limits=c("LAP_SD_WW","LAP_SD_LW"), labels = c('Well-watered','Drought'), expand = c(0.1, 0.1)) +
  geom_point(color= "black",size=5,
             alpha =0.6,
             position=position_jitter(width=0.005)) +    
  geom_text_repel(aes(label = name), data = data_ends,
                  fontface ="bold", color = "black", size = 7, box.padding = 1, force = 0.5, segment.colour="#575d5d",
                  nudge_x      = -0.30,
                  direction    = "y",
                  hjust        = 0,
                  segment.size = 0.2,max.overlaps = Inf
  ) +
  geom_line(aes(group=accession_id), color="black",size = 1.2,
            alpha = 0.7)+ 
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.25,"cm")) +
  theme(axis.text.x = element_text(family = "Arial", color="black", face="bold", size=35)) +
  theme(axis.text.y = element_text(family = "Arial", color="black", face="bold", size=15)) +
  theme(axis.title.x = element_text(family = "Arial", color="black", face="bold", size=40)) +
  theme(axis.title.y = element_text(family = "Arial", color="black", face="bold", size=40)) +
  theme(plot.title = element_text(hjust = 0.5, vjust = - 5, family = "Arial", color="black", face="bold", size=50)) +
  theme(axis.title.x=element_blank())+ 
  theme(plot.margin=unit(c(1,1,1,1),"cm"))+
  guides(fill = FALSE) + 
  ylab(expression(bold(paste("Leaf area potential"))))+
  xlab("Treatment")+ 
   ggtitle("Leaf area potential: reaction norms")

p5b


#########################################################################################################Arrange panels and save figure
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
plot_grid(p5a,p5b, labels=c("E", "F"), greedy = TRUE, ncol = 1, nrow =2, align = 'v', hjust=0, label_size=75)
png("panels/BottomRight.png", width = 22, height = 24, units = 'in', res = 350)
plot_grid(p5a,p5b, labels=c("E", "F"), greedy = TRUE, ncol = 1, nrow =2, align = 'v', hjust=0, label_size=75)
dev.off()
BottomRight <- image_read("panels/BottomRight.png")

plot_grid(p3,p4, labels=c("B", "C"), greedy = TRUE, ncol = 1, nrow = 2, align = 'v', hjust=0, label_size=75)
png("panels/TopRight.png", width = 22, height = 24, units = 'in', res = 350)
plot_grid(p3,p4, labels=c("B", "C"), greedy = TRUE, ncol = 1, nrow = 2, align = 'v', hjust=0, label_size=75)
dev.off()
TopRight <- image_read("panels/TopRight.png")

M1 <- image_read("panels/p1.jpg")
M2 <- image_read_pdf("panels/p2b.pdf")

plot_grid(rasterGrob(M1),rasterGrob(TopRight),labels=c("A", ""),  greedy = TRUE, ncol = 2, nrow = 1, align = 'v', hjust=0.03, label_size=40)
png("panels/top.png", width = 22, height = 11.5, units = 'in', res = 350)
plot_grid(rasterGrob(M1),rasterGrob(TopRight),labels=c("A", ""),  greedy = TRUE, ncol = 2, nrow = 1, align = 'v', hjust=0.03, label_size=40)
dev.off()

TOP <- image_read("panels/top.png")

plot_grid(rasterGrob(M2),rasterGrob(BottomRight), labels=c("D", ""), greedy = TRUE, ncol = 2, nrow = 1, align = 'v', hjust=0, label_size=35)
png("panels/Bottom.png", width = 22, height = 11.5, units = 'in', res = 350)
plot_grid(rasterGrob(M2),rasterGrob(BottomRight), labels=c("D", ""), greedy = TRUE, ncol = 2, nrow = 1, align = 'v', hjust=0, label_size=35)
dev.off()

Bottom <- image_read("panels/Bottom.png")

plot_grid(rasterGrob(TOP),rasterGrob(Bottom), labels=c("", ""), greedy = TRUE, ncol = 1, nrow = 2, align = 'v', hjust=0, label_size=40)
png("panels/Figure 1.png", width = 20, height = 22, units = 'in', res = 350)
plot_grid(rasterGrob(TOP),rasterGrob(Bottom), labels=c("", ""), greedy = TRUE, ncol = 1, nrow = 2, align = 'v', hjust=0, label_size=40)
dev.off()




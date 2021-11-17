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
library(ggpmisc)

######################################################################################################### Reaction norm figures (Fig. 1E & 1F)
#########################################################################################################
#########################################################################################################
#########################################################################################################

## Load leaf area  phenotypes. This file can be found the the "data" folder. 
## LA_SD_WW = Leaf area under well-watered conditions (cm2)
## LA_SD_LW = Leaf area under drought (cm2)

###Open dataframe with phenotypes
data = read_csv("data/phenotypesb.csv")
df1a <- data[c("accession_id", "accession_name", "Leaf_area_well_watered", "Leaf_area_drought")]
head(df1a)


######Plot reaction norms for leaf area (well-watered vs. drought) (Fig. 1E)

head(df1a)

df1a<-df1a[!is.na(df1a$Leaf_area_well_watered), ]
df1a<- df1a %>%
  gather(TRT, phenotype, Leaf_area_well_watered:Leaf_area_drought)
df1a<-df1a[complete.cases(df5a), ]

head(df1a)

my_comparisons <- list( c("Leaf_area_well_watered", "Leaf_area_drought") )

data_ends <- df1a %>% filter(TRT == "Leaf_area_well_watered")

p1a <- ggplot(df1a, aes(factor(TRT), phenotype, fill=TRT)) + geom_violin(alpha=0.5,color = NA) +
  scale_fill_manual(values=c("#9f4c49", "#008764")) + scale_x_discrete(limits=c("Leaf_area_well_watered","Leaf_area_drought"), labels = c('Well-watered','Drought'), expand = c(0.1, 0.1)) +
  geom_point(color= "black",size=5,
             alpha =0.6,
             position=position_jitter(width=0.005)) +    
  geom_text_repel(aes(label = accession_name), data = data_ends,
                  fontface ="bold", color = "black", size = 5, box.padding = 1, force = 0.5, segment.colour="#575d5d",
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
  theme(plot.title = element_text(hjust = 0.5, vjust = - 5, family = "Arial", color="black", face="bold", size=30)) +
  theme(axis.title.x=element_blank())+ 
  theme(plot.margin=unit(c(1,1,1,1),"cm"))+
  guides(fill = "none") + 
  ylab(expression(bold(paste("Leaf area (", cm^2,")"))))+
  xlab("Treatment")+  ggtitle("Leaf area: reaction norms")+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") 
p1a


######Plot eaction norms for leaf area potentials (well-watered vs. drought) (Fig. 1E). In this plot and further analysis, we remove the values of non-Iberian Col-0, which we used as an internal reference, as well as the the values for accession IP-Coc1 (accession id = 9535), which were consistently unrealistically influential across different phenotypes based on their extreme Cook’s distance. 

######Leaf area reaction norm
## Load leaf area  phenotypes. This file can be found the the "data" folder. 
## LA_SD_WW = Leaf area under well-watered conditions (cm2)
## LA_SD_LW = Leaf area under drought (cm2)

###Open dataframe with phenotypes
data = read_csv("data/phenotypesb.csv")
df1b <- data[c("accession_id", "accession_name", "LAP_well_watered", "LAP_drought")]
head(data)
head(df1b)
df1b <- df1b[c("accession_id", "accession_name", "LAP_well_watered", "LAP_drought")]
head(df1b)

df1b<-df1b[!is.na(df1b$LAP_drought), ]
df1b<- df1b %>%
  gather(TRT, phenotype, LAP_drought:LAP_well_watered)
df1b<-df1b[complete.cases(df1b), ]

head(df1b)

my_comparisons <- list( c("LAP_well_watered", "LAP_SD_LW") )

data_ends <- df1b %>% filter(TRT == "LAP_well_watered")

p1b <- ggplot(df1b, aes(factor(TRT), phenotype, fill=TRT)) + geom_violin(alpha=0.5,color = NA) +
  scale_fill_manual(values=c("#9f4c49", "#008764")) + scale_x_discrete(limits=c("LAP_well_watered","LAP_drought"), labels = c('Well-watered','Drought'), expand = c(0.1, 0.1)) +
  geom_point(color= "black",size=5,
             alpha =0.6,
             position=position_jitter(width=0.005)) +    
  geom_text_repel(aes(label = accession_name), data = data_ends,
                  fontface ="bold", color = "black", size = 5, box.padding = 1, force = 0.5, segment.colour="#575d5d",
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
  theme(plot.title = element_text(hjust = 0.5, vjust = - 5, family = "Arial", color="black", face="bold", size=30)) +
  theme(axis.title.x=element_blank())+ 
  theme(plot.margin=unit(c(1,1,1,1),"cm"))+
  guides(fill = "none") + 
  ylab(expression(bold(paste("Relative leaf area (rLA)"))))+
  xlab("Treatment")+ 
   ggtitle("Relative leaf area (rLA): reaction norms")

p1b


#########################################################################################################Arrange panels and save figure
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
df2 = read_csv("data/phenotypesb.csv")
head(df2)


head(df3)

####Plot6

formula <- y ~ x
p2<-df2 %>% 
  mutate(quadrant = case_when(LAP_well_watered > 1 & PI > 0   ~ "Q1",
                              LAP_well_watered <= 1 & PI > 0  ~ "Q2",
                              LAP_well_watered <= 1 & PI <= 0 ~ "Q3",
                              TRUE                                         ~ "Q4")) %>% 
  ggplot(aes(x = LAP_well_watered, y = PI)) + 
  geom_point(aes(color=quadrant),alpha = 0.7,  size=11) +  
  geom_smooth(method = "lm", se = T, fill="black", colour="black", size=0.2, alpha = 0.1) +  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")),
                                                                                                          label.x.npc = "right", label.y.npc = "top",
                                                                                                          formula = formula, parse = TRUE, size = 6, colour="black")+
  scale_colour_manual(values = c("#47a569", "#6689c3", "#ce494a", "#e9af41"))+
  geom_text_repel(aes(LAP_well_watered, PI, label = accession_name),  size=4, force=0.1,segment.size = 0.5,arrow = arrow(length = unit(0.01, 'npc')),family = 'Arial',
                  fontface = 'bold',colour = "#696969",
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.5, "lines"),
                  segment.color = 'grey50',max.overlaps = Inf) + 
  geom_vline(xintercept = 1, color = "black", linetype = "dashed")+
  geom_hline(yintercept = 0, color = "black", linetype = "dashed")+
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm"))  +
  theme(plot.title = element_text(family = "Arial", color="#696969",face="bold",  size=35))  +  theme(panel.background = element_rect(fill = 'white')) +
  xlab("Relative leaf area (rLA) under well-watered conditions ") + 
  ylab("Leaf area plasticity index (PI)")+
  theme(axis.text = element_text(family = "Arial", color="black", size=16)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=16)) +  
  theme(legend.position="none") 
p2
png("panels/p2.png", width = 7, height = 7, units = 'in', res = 350)
p2
dev.off() 

#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################

M2 <- image_read_pdf("panels/I2.pdf")

plot_grid(p1a,p1b,rasterGrob(M1), labels=c("b", "c","e"), greedy = TRUE, ncol = 1, nrow =3, align = 'v',  label_size=60)
png("panels/Right.png", width = 12, height = 33, units = 'in', res = 350)
plot_grid(p1a,p1b,rasterGrob(M2), labels=c("b", "c","e"), greedy = TRUE, ncol = 1, nrow =3, align = 'v',  label_size=60)
dev.off()

M1 <- image_read("panels/I1.jpg")
p2 <- image_read("panels/p2.png")

plot_grid(rasterGrob(M1),rasterGrob(p2), labels=c("a", "d"), greedy = TRUE, ncol = 1, nrow =2, align = 'v',  label_size=40)
png("panels/left.png", width = 9, height = 21, units = 'in', res = 350)
plot_grid(rasterGrob(M1),rasterGrob(p2), labels=c("a", "d"), greedy = TRUE, ncol = 1, nrow =2, align = 'v',  label_size=40)
dev.off()

left <- image_read("panels/left.png")
right <- image_read("panels/right.png")

plot_grid(rasterGrob(left),rasterGrob(right), labels=c("", ""), greedy = TRUE, ncol = 2, nrow =1, align = 'v',  label_size=50)
png("Figure1.png", width = 21, height = 24, units = 'in', res = 350)
plot_grid(rasterGrob(left),rasterGrob(right), labels=c("", ""), greedy = TRUE, ncol = 2, nrow =1, align = 'v',  label_size=50)
dev.off()


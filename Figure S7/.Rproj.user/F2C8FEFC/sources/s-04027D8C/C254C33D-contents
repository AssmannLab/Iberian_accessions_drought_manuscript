########Scatter plots: our phenotypes against altitude
#scale_fill_manual(values=c("#82d6b6", "#008764"))  +  scale_x_discrete(limits=c("SD","LD")) +

# scale_fill_manual(values=c("#9f4c49", "#008764")) +  scale_x_discrete(limits=c("WW","LW")) +

# color #008764 for well watered treatment under short days
# color #82d6b6 for drought treatment under short days
# color #82d6b6 for long day treatment
#  scale_fill_manual(values=c("#32adb2", "#ff9a00")) +  scale_x_discrete(limits=c("drought","photoperiod")) +
library(dplyr) 
library(ggplot2)
library(ggthemes)
library(viridis)
library(ggpmisc)
library(ggrepel)
library(readr)
library(cowplot)
library(Hmisc)
library(readr)
library(tidyverse)


df=read_csv("data/data.csv")
head(df)
###Plot1= Nppspring vs leaf area plasticity in response to drought
formula <- y ~ x
p1 <- ggplot(df, aes(x=PI , y= Percentage_reduction)) +
  geom_point(color='grey',alpha = 0.8,  size=12)+   geom_text_repel(aes(PI, Percentage_reduction, label = accession_name, fontface="bold")) +
  geom_smooth(method = "lm", se = T, fill="black", colour="black", size=0.5, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                         label.x.npc = "right", label.y.npc = "bottom",
                                                                                                         formula = formula, parse = TRUE, size = 7, colour="black") +
  
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=25)) +  
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm"))  +
  theme(axis.text = element_text(family = "Arial", color="black",  face="bold",size=10)) +
  theme(plot.title = element_text(family = "Arial", color="black", face="bold",size=28)) +  theme(panel.background = element_rect(fill = 'white')) +
  theme(legend.position="none") + 
  xlab(expression(bold("Leaf area plasticity index")))+ 
  ylab(expression(bold("Leaf area change (%)")))

p1
coef
linearMod <- lm(Leaf_area_plasticity_index ~ Percentage_reduction, data=df)  # build linear regression model on full data
print(linearMod)
summary(linearMod)
summary(linearMod)$r.squared

png("Figure S4.png", width = 7, height = 7, units = 'in', res = 250)
p1
dev.off()
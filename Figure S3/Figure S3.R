#READ TWO FILES
library(tidyverse)
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(Hmisc)
library(ggpmisc)
library(qvalue)
library(ggpubr)
library(cowplot)
library(patchwork)
library(magick)
library(png)
library(grid)

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

df1 = read_csv("data/Iberian_d13C.csv")
head(phenotypes)


Phenotypesselected<- select(df1, accession_id, Dittberner_moleco18_Delta_13C)
head(Phenotypesselected)
transposedfinal = read_csv("data/transposedfinal.csv")
Transcriptselected<- select(transposedfinal, accession_id,AT3G55800)
head(Transcriptselected)
df1<-merge(Phenotypesselected, Transcriptselected, by = "accession_id", sort = TRUE)
head(df1)

###Plot1= altitude vs leaf area time potential, short days well-watered
formula <- y ~ x
p1 <- ggplot(df1, aes(x=AT3G55800 , y= Dittberner_moleco18_Delta_13C)) +
  geom_point(color='black',alpha = 0.4,  size=5)+  
  geom_smooth(method = "lm", se = T, fill="black", colour="black", size=0.5, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                         label.x.npc = "right", label.y.npc = "bottom",
                                                                                                         formula = formula, parse = TRUE, size = 4, colour="black") +  theme_bw()# +

p1<- p1 + theme(axis.line = element_line(size=1, colour = "black"),
                panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm")) 
p1<- p1 + 
  theme(axis.text = element_text(family = "Arial", color="black", face="bold",size=10)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=15)) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(5, 5, 5,5)
  )+ 
  #ylim(-1,1)+
  ylab(expression(bold(paste("WUE"[i]," (\u03B4"^"13", "C (\u2030))"))))+   
  xlab(expression(bold(paste(bolditalic("SBPASE")~"(FPKM)"))))
p1
png("panels/p1.png", width = 7, height = 7, units = 'in', res = 350)
p1
dev.off() 

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

df2 = read_csv("data/Iberian_d13C.csv")
head(phenotypes)


Phenotypesselected<- select(df2, accession_id, Dittberner_moleco18_Delta_13C)
head(Phenotypesselected)
transposedfinal = read_csv("data/transposedfinal.csv")
Transcriptselected<- select(transposedfinal, accession_id,AT4G34830)
head(Transcriptselected)
df2<-merge(Phenotypesselected, Transcriptselected, by = "accession_id", sort = TRUE)
head(df2)

###Plot1= altitude vs leaf area time potential, short days well-watered
formula <- y ~ x
p2 <- ggplot(df2, aes(x=AT4G34830 , y= Dittberner_moleco18_Delta_13C)) +
  geom_point(color='black',alpha = 0.4,  size=5)+  
  geom_smooth(method = "lm", se = T, fill="black", colour="black", size=0.5, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                         label.x.npc = "right", label.y.npc = "bottom",
                                                                                                         formula = formula, parse = TRUE, size = 4, colour="black") +  theme_bw()# +

p2<- p2 + theme(axis.line = element_line(size=1, colour = "black"),
                panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm")) 
p2<- p2 + 
  theme(axis.text = element_text(family = "Arial", color="black", face="bold",size=10)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=15)) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(5, 5, 5,5)
  )+ 
  #ylim(-1,1)+
  ylab(expression(bold(paste("WUE"[i]," (\u03B4"^"13", "C)"))))+  
  xlab(expression(bold(paste(bolditalic("MRL1")~"(FPKM)"))))
  
p2
png("panels/p2.png", width = 7, height = 7, units = 'in', res = 350)
p2
dev.off() 

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

df3 = read_csv("data/Iberian_d13C.csv")
head(phenotypes)


Phenotypesselected<- select(df3, accession_id, Dittberner_moleco18_Delta_13C)
head(Phenotypesselected)
transposedfinal = read_csv("data/transposedfinal.csv")
Transcriptselected<- select(transposedfinal, accession_id,AT4G37925)
head(Transcriptselected)
df3<-merge(Phenotypesselected, Transcriptselected, by = "accession_id", sort = TRUE)
head(df3)

###Plot1= altitude vs leaf area time potential, short days well-watered
formula <- y ~ x
p3 <- ggplot(df3, aes(x=AT4G37925 , y= Dittberner_moleco18_Delta_13C)) +
  geom_point(color='black',alpha = 0.4,  size=5)+  
  geom_smooth(method = "lm", se = T, fill="black", colour="black", size=0.5, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                         label.x.npc = "right", label.y.npc = "bottom",
                                                                                                         formula = formula, parse = TRUE, size = 4, colour="black") +  theme_bw()# +

p3<- p3 + theme(axis.line = element_line(size=1, colour = "black"),
                panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm")) 
p3<- p3 + 
  theme(axis.text = element_text(family = "Arial", color="black", face="bold",size=10)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=15)) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(5, 5, 5,5)
  )+ 
  #ylim(-1,1)+
  ylab(expression(bold(paste("WUE"[i]," (\u03B4"^"13", "C)"))))+  
  xlab(expression(bold(paste(bolditalic("NdhM")~"(FPKM)"))))

p3
png("panels/p3.png", width = 7, height = 7, units = 'in', res = 350)
p3
dev.off() 

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

df4 = read_csv("data/Iberian_d13C.csv")
head(phenotypes)


Phenotypesselected<- select(df4, accession_id, Dittberner_moleco18_Delta_13C)
head(Phenotypesselected)
transposedfinal = read_csv("data/transposedfinal.csv")
Transcriptselected<- select(transposedfinal, accession_id,AT3G54050)
head(Transcriptselected)
df4<-merge(Phenotypesselected, Transcriptselected, by = "accession_id", sort = TRUE)
head(df4)

###Plot1= altitude vs leaf area time potential, short days well-watered
formula <- y ~ x
p4 <- ggplot(df4, aes(x=AT3G54050, y= Dittberner_moleco18_Delta_13C)) +
  geom_point(color='black',alpha = 0.4,  size=5)+  
  geom_smooth(method = "lm", se = T, fill="black", colour="black", size=0.5, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                         label.x.npc = "right", label.y.npc = "bottom",
                                                                                                         formula = formula, parse = TRUE, size = 4, colour="black") +  theme_bw()# +

p4<- p4 + theme(axis.line = element_line(size=1, colour = "black"),
                panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm")) 
p4<- p4 + 
  theme(axis.text = element_text(family = "Arial", color="black", face="bold",size=10)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=15)) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(5, 5, 5,5)
  )+ 
  #ylim(-1,1)+
  ylab(expression(bold(paste("WUE"[i]," (\u03B4"^"13", "C)"))))+  
  xlab(expression(bold(paste(bolditalic("HCEF1")~"(FPKM)"))))

p4
png("panels/p4.png", width = 7, height = 7, units = 'in', res = 350)
p4
dev.off() 

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

df5 = read_csv("data/Iberian_d13C.csv")
head(phenotypes)


Phenotypesselected<- select(df5, accession_id, Dittberner_moleco18_Delta_13C)
head(Phenotypesselected)
transposedfinal = read_csv("data/transposedfinal.csv")
Transcriptselected<- select(transposedfinal, accession_id,AT4G09650)
head(Transcriptselected)
df5<-merge(Phenotypesselected, Transcriptselected, by = "accession_id", sort = TRUE)
head(df5)

###Plot1= altitude vs leaf area time potential, short days well-watered
formula <- y ~ x
p5 <- ggplot(df5, aes(x=AT4G09650, y= Dittberner_moleco18_Delta_13C)) +
  geom_point(color='black',alpha = 0.4,  size=5)+  
  geom_smooth(method = "lm", se = T, fill="black", colour="black", size=0.5, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                         label.x.npc = "right", label.y.npc = "bottom",
                                                                                                         formula = formula, parse = TRUE, size = 4, colour="black") +  theme_bw()# +

p5<- p5 + theme(axis.line = element_line(size=1, colour = "black"),
                panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm")) 
p5<- p5 + 
  theme(axis.text = element_text(family = "Arial", color="black", face="bold",size=10)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=15)) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(5, 5, 5,5)
  )+ 
  #ylim(-1,1)+
  ylab(expression(bold(paste("WUE"[i]," (\u03B4"^"13", "C)"))))+  
  xlab(expression(bold(paste(bolditalic("ATPD")~"(FPKM)"))))

p5
png("panels/p5.png", width = 7, height = 7, units = 'in', res = 350)
p5
dev.off() 

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################


df8 = read_csv("data/Iberian_d13C.csv")
head(phenotypes)


Phenotypesselected<- select(df8, accession_id, Dittberner_moleco18_Delta_13C)
head(Phenotypesselected)
transposedfinal = read_csv("data/transposedfinal.csv")
Transcriptselected<- select(transposedfinal, accession_id,AT5G64040)
head(Transcriptselected)
df8<-merge(Phenotypesselected, Transcriptselected, by = "accession_id", sort = TRUE)
head(df8)

###Plot1= altitude vs leaf area time potential, short days well-watered
formula <- y ~ x
p8 <- ggplot(df8, aes(x=AT5G64040, y= Dittberner_moleco18_Delta_13C)) +
  geom_point(color='black',alpha = 0.4,  size=5)+  
  geom_smooth(method = "lm", se = T, fill="black", colour="black", size=0.5, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                         label.x.npc = "right", label.y.npc = "bottom",
                                                                                                         formula = formula, parse = TRUE, size = 4, colour="black") +  theme_bw()# +

p8<- p8 + theme(axis.line = element_line(size=1, colour = "black"),
                panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm")) 
p8<- p8 + 
  theme(axis.text = element_text(family = "Arial", color="black", face="bold",size=10)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=15)) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(5, 5, 5,5)
  )+ 
  #ylim(-1,1)+
  ylab(expression(bold(paste("WUE"[i]," (\u03B4"^"13", "C)"))))+  
  xlab(expression(bold(paste(bolditalic("PSAN")~"(FPKM)"))))

p8
png("panels/p8.png", width = 7, height = 7, units = 'in', res = 350)
p8
dev.off() 

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

df9 = read_csv("data/Iberian_d13C.csv")
head(phenotypes)


Phenotypesselected<- select(df9, accession_id, Dittberner_moleco18_Delta_13C)
head(Phenotypesselected)
transposedfinal = read_csv("data/transposedfinal.csv")
Transcriptselected<- select(transposedfinal, accession_id,AT3G50820)
head(Transcriptselected)
df9<-merge(Phenotypesselected, Transcriptselected, by = "accession_id", sort = TRUE)
head(df9)

###Plot1= altitude vs leaf area time potential, short days well-watered
formula <- y ~ x
p9 <- ggplot(df9, aes(x=AT3G50820, y= Dittberner_moleco18_Delta_13C)) +
  geom_point(color='black',alpha = 0.4,  size=5)+  
  geom_smooth(method = "lm", se = T, fill="black", colour="black", size=0.5, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                         label.x.npc = "right", label.y.npc = "bottom",
                                                                                                         formula = formula, parse = TRUE, size = 4, colour="black") +  theme_bw()# +

p9<- p9 + theme(axis.line = element_line(size=1, colour = "black"),
                panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm")) 
p9<- p9 + 
  theme(axis.text = element_text(family = "Arial", color="black", face="bold",size=10)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=15)) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(5, 5, 5,5)
  )+ 
  #ylim(-1,1)+
  ylab(expression(bold(paste("WUE"[i]," (\u03B4"^"13", "C)"))))+  
  xlab(expression(bold(paste(bolditalic("PSBO2")~"(FPKM)"))))

p9
png("panels/p9.png", width = 7, height = 7, units = 'in', res = 350)
p9
dev.off() 
#######################################################################################################################################
#######################################################################################################################################
########################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

df10 = read_csv("data/Iberian_d13C.csv")
head(phenotypes)


Phenotypesselected<- select(df10, accession_id, Dittberner_moleco18_Delta_13C)
head(Phenotypesselected)
transposedfinal = read_csv("data/transposedfinal.csv")
Transcriptselected<- select(transposedfinal, accession_id,AT1G48370)
head(Transcriptselected)
df10<-merge(Phenotypesselected, Transcriptselected, by = "accession_id", sort = TRUE)
head(df10)

###Plot1= altitude vs leaf area time potential, short days well-watered
formula <- y ~ x
p10 <- ggplot(df10, aes(x=AT1G48370, y= Dittberner_moleco18_Delta_13C)) +
  geom_point(color='black',alpha = 0.4,  size=5)+  
  geom_smooth(method = "lm", se = T, fill="black", colour="black", size=0.5, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                         label.x.npc = "right", label.y.npc = "bottom",
                                                                                                         formula = formula, parse = TRUE, size = 4, colour="black") +  theme_bw()# +

p10<- p10 + theme(axis.line = element_line(size=1, colour = "black"),
                panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm")) 
p10<- p10 + 
  theme(axis.text = element_text(family = "Arial", color="black", face="bold",size=10)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=15)) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(5, 5, 5,5)
  )+ 
  #ylim(-1,1)+
  ylab(expression(bold(paste("WUE"[i]," (\u03B4"^"13", "C)"))))+  
  xlab(expression(bold(paste(bolditalic("YSL8")~"(FPKM)"))))

p10
png("panels/p10.png", width = 7, height = 7, units = 'in', res = 350)
p10
dev.off() 
###################################################################################################################################################
###################################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

df11 = read_csv("data/Iberian_d13C.csv")
head(phenotypes)


Phenotypesselected<- select(df11, accession_id, Dittberner_moleco18_Delta_13C)
head(Phenotypesselected)
transposedfinal = read_csv("data/transposedfinal.csv")
Transcriptselected<- select(transposedfinal, accession_id,AT2G26510)
head(Transcriptselected)
df11<-merge(Phenotypesselected, Transcriptselected, by = "accession_id", sort = TRUE)
head(df11)

###Plot1= altitude vs leaf area time potential, short days well-watered
formula <- y ~ x
p11 <- ggplot(df11, aes(x=AT2G26510, y= Dittberner_moleco18_Delta_13C)) +
  geom_point(color='black',alpha = 0.4,  size=5)+  
  geom_smooth(method = "lm", se = T, fill="black", colour="black", size=0.5, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                         label.x.npc = "right", label.y.npc = "bottom",
                                                                                                         formula = formula, parse = TRUE, size = 4, colour="black") +  theme_bw()# +

p11<- p11 + theme(axis.line = element_line(size=1, colour = "black"),
                  panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
                  panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm")) 
p11<- p11 + 
  theme(axis.text = element_text(family = "Arial", color="black", face="bold",size=10)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=15)) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(5, 5, 5,5)
  )+ 
  #ylim(-1,1)+
  ylab(expression(bold(paste("WUE"[i]," (\u03B4"^"13", "C)"))))+  
  xlab(expression(bold(paste(bolditalic("PDE135")~"(FPKM)"))))

p11
png("panels/p11.png", width = 7, height = 7, units = 'in', res = 350)
p11
dev.off() 
####################################################################################################################################################
###################################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################



plot_grid(p1,p2,p3,p4,p5,p8,p9,p10,p11,labels=c('A', 'B',"C","D","E","F","G","H","I"), greedy = TRUE, ncol =3, nrow = 3, align = 'v', hjust=0, label_size=35)
png("Figure S3.png", width = 15, height = 15, units = 'in', res = 350)
plot_grid(p1,p2,p3,p4,p5,p8,p9,p10,p11,labels=c('A', 'B',"C","D","E","F","G","H","I"), greedy = TRUE, ncol =3, nrow = 3, align = 'v', hjust=0, label_size=35)
dev.off()

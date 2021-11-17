list.of.packages <- c("ggplot2", "viridis", "ggrepel","readr", "tidyverse","tidyr","cowplot","qvalue","dplyr","ggpmisc","ggpubr","ggthemes", "Hmisc", "cowplot","magick","grid", "ggpmisc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


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
library(magick)
library(png)
library(grid)

#######################################################################################################################################Fig. 6A d13C GWAS
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################Figure 1A: before running the script, pleas unzip files from the /data folder and /data/GWAS f
dat1 = read.csv("data/GWAS/iWUE_Iberian_LM.csv")

head(dat1)
nrow(dat1)
dat2 = read_csv("data/annotation_GWAS.csv")
head(dat2)
total <- merge(dat1,dat2,by=c("chr", "pos"))
head(total)
total$X1 <- NULL
total$GVE<- NULL

total <- total[order(-total$score) ,  ]
total <- total[total$maf > 0.1, ]
head(total)
duplicated_position_column <-total$pos
total <-cbind(total,duplicated_position_column)

total<-total %>% unite(label, symbol, duplicated_position_column, sep = " pos. ")
total$X1 <- NULL
total$GVE <- NULL
colnames(total)[colnames(total)=="chr"] <- "CHR"
colnames(total)[colnames(total)=="pos"] <- "BP"

head(total)

write.csv(total, "Tables/Annotated_iWUE_Iberian_LM.csv", row.names = FALSE)

don <- total %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(total, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)
head(don)

axisdf = don %>% group_by(CHR) %>% summarise(center=( max(BPcum) + min(BPcum) ) / 2 )

p1<-don %>% 
  mutate(quadrant = case_when(score >= 6.5   ~ "Q1")) %>% 
  mutate(label = case_when(score >= 6.5  ~  paste(label))) %>% 
  
  ggplot( aes(x=BPcum, y=score)) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.4, size=3) +
  scale_color_manual(values = rep(c("#F8B195", "#F67280", "#C06C84","#6C5B7B","#355C7D"), 5 )) +
  
#  geom_text_repel(aes(BPcum, score, label = label), fontface = 'bold.italic', size=2, force=0.1,segment.size = 0.1,arrow = arrow(length = unit(0.01, 'npc')),family = 'Arial',
 #                 box.padding = unit(0.15, "lines"),
  #                point.padding = unit(0.45, "lines"),
   #               segment.color = 'grey50',max.overlaps = Inf ) +  
  
  geom_label_repel(aes(BPcum, score, label = label,fill = factor(quadrant)),
    fontface = 'bold.italic', color = 'white', size =2,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines"),
    segment.color = 'grey50',max.overlaps = Inf 
  ) +
  
  
  #geom_hline(yintercept = 4, color = "#A9A9A9", linetype = "longdash")+
  # custom X axis:
  scale_x_continuous(labels=c("1" = "Chr. 1", "2" = "Chr. 2","3" = "Chr. 3","4" = "Chr. 4","5" = "Chr. 5"), breaks= axisdf$center )+
  #  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  #scale_x_continuous( label = "axisdf$CHR", breaks= axisdf$center ) +
 # scale_y_continuous(limits = c(0, 7.5))+
  #xlab(expression(bold("Cyclic nucleotide-gated channel 16 (FPKM)")))+ 
  theme(legend.position = "none")+
  # remove space between plot area and x axis
  
  # Custom the theme:
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm")) +
  theme(panel.background = element_rect(fill = 'white')) +
  ylab(expression(bold(atop(paste("WUE"[i]," (\u03B4"^"13", "C)"), paste("Genome-wide" ~ association ~ (-log["10"] ~ bolditalic("P")))))))+ 
  theme(axis.text.y = element_text(family = "Arial", color="black", size=25, face="bold")) +
  theme(axis.text.x = element_text(family = "Arial", color="black", size=25, face="bold")) +
  theme(axis.title.y = element_text(family = "Arial", color="black", size=25, face="bold")) +
  theme(axis.title.x=element_blank()) +
  theme(plot.title = element_text(family = "Arial", color="#696969",face="bold",  size=35)) +  
  ggtitle("GWAS")+ylim(0,8)
p1


#######################################################################################################################################Fig. 6A (inset)
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
phenotypes = read.csv("data/ClpX3_haplotypes.csv")
head(phenotypes)

my_comparisons <- c("Major", "Minor")
my_font_size <- 2.5
inset <- ggplot(phenotypes, aes(x=ClpX3_haplotypes, y=Dittberner_moleco18_Delta_13C, fill=ClpX3_haplotypes)) + geom_boxplot(width=0.3,alpha=1) +  geom_violin(trim=FALSE,alpha=0.2) + 
  scale_fill_manual(values=c("#fa7b49", "#fdd5c5")) + scale_x_discrete(limits=c("Major","Minor"), labels = c('Major haplotype','Minor haplotype'), expand = c(0.1, 0.15)) +
  geom_point( shape = 21,size=0.75, position = position_jitterdodge(), color="black",alpha=1)+
  theme(axis.line = element_line(size=0.75, colour = "black"),
        panel.grid.major = element_line(size=0.05,colour = "#d3d3d3"), panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.25,"cm")) +
  theme(axis.text.x = element_text(family = "Arial", color="black", face="bold", size=8)) +
  theme(axis.text.y = element_text(family = "Arial", color="black", face="bold", size=7)) +
  theme(axis.title.y = element_text(family = "Arial", color="black", face="bold", size=8)) +
  theme(axis.title.x=element_blank())+ 
  # ggtitle("Our sample") +
  guides(fill = FALSE) + 
  ylab(expression(bold(paste("WUE"[i]," (\u03B4"^"13", "C (\u2030))"))))+   
  stat_compare_means(size = my_font_size,label.x.npc = "left",label.y.npc= "top",vjust=-1.8, hjust=0.5) +
theme(plot.title = element_text(family = "Arial", color="#696969",face="bold.italic",  size=13)) +  
  ggtitle("ClpX3")
inset$layers[[2]]$aes_params$textsize <- my_font_size
inset

final <- p1 + inset_element(inset, left = 0.28, bottom = 0.7, right = 0.45, top = 0.99)

png("panels/p1.png", width = 17, height = 7, units = 'in', res = 350)
final
dev.off() 

#######################################################################################################################################Fig. 6B d13C TWAS
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################


annotation = read_csv("data/gene_description.csv")
head(annotation)
phenotypes = read_csv("data/Iberian_d13C.csv")
head(phenotypes)
#####Run autocorrelation matrix
head(phenotypes)
phenotypesselected<- select(phenotypes, Dittberner_moleco18_Delta_13C)
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
BB<-signif(res$P,2)
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
outputA <- outputA[order(-outputA$Dittberner_moleco18_Delta_13C.x) ,  ]
head(outputA)
outputA <- merge(outputA,annotation,by="gene")
outputA <- outputA %>% 
  rename(
    Dittberner_moleco18_Delta_13C_rs = Dittberner_moleco18_Delta_13C.x,
    Dittberner_moleco18_Delta_13C_P_val = Dittberner_moleco18_Delta_13C.y)
outputA$accession_idsorted.x <- NULL
outputA$accession_idsorted.y <- NULL
outputA <- outputA[order(-outputA$Dittberner_moleco18_Delta_13C_rs) ,  ]
head(outputA, n=100)
tail(outputA, n=50)
head(outputA)
write.csv(outputA, "Tables/TWAS_Iberian_Dittberner_moleco18_Delta_13C.csv", row.names = FALSE)



#data = read_tsv("/Users/angel_admin/Documents/stability paper october 18/ANNOTATED FULL/WUE_fixed_appended.tsv")
datb = read.csv("data/mapping.csv")
head(datb)


dat2 <- merge(outputA,datb,by= "gene")
head(dat2)
nrow(dat2)
#dat1 <- dat1 %>% 
# mutate(CHR = ifelse(grepl("2$", Trtmt), 2, 1))

labelcolor <-  c("#003f5c", "#58508d", "#bc5090","#ffa600")


don2 <- dat2 %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(dat2, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

axisdf = don2 %>% group_by(CHR) %>% summarise(center=( max(BPcum) + min(BPcum) ) / 2 )

p2<-don2 %>% 
  mutate(quadrant = case_when(Dittberner_moleco18_Delta_13C_rs <= -0.4   ~ "Q1",
                              Dittberner_moleco18_Delta_13C_rs >= 0.4  ~ "Q4")) %>% 
  mutate(label = case_when(Dittberner_moleco18_Delta_13C_rs <= -0.4  ~  paste(symbol),
                           Dittberner_moleco18_Delta_13C_rs >= 0.4  ~  paste(symbol))) %>% 
  ggplot( aes(x=BPcum, y=Dittberner_moleco18_Delta_13C_rs)) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.6, size=12, shape = "-") +
  scale_color_manual(values = rep(c("#F8B195", "#F67280", "#C06C84","#6C5B7B","#355C7D","#31d314","#d32e14" ), 5 )) +
  
  geom_label_repel(
    aes(BPcum, Dittberner_moleco18_Delta_13C_rs, label = label,fill = factor(quadrant)),
    fontface = 'bold.italic', color = 'white', size =3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines"),
    segment.color = 'grey50',max.overlaps = Inf 
  ) +
  
  geom_hline(yintercept = 0.4, color = "black", linetype = "longdash")+
  geom_hline(yintercept = 0, color = "black", linetype = "solid")+
  
  geom_hline(yintercept = -0.4, color = "black", linetype = "longdash")+
  
  # custom X axis:
  scale_x_continuous(labels=c("1" = "Chr. 1", "2" = "Chr. 2","3" = "Chr. 3","4" = "Chr. 4","5" = "Chr. 5", "C"="C", "M" = "M"), breaks= axisdf$center )+
  
  #scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  #scale_x_continuous( label = "axisdf$CHR", breaks= axisdf$center ) +
  scale_y_continuous(limits = c(-0.6, 0.6))+
  #xlab(expression(bold("Cyclic nucleotide-gated channel 16 (FPKM)")))+ 
  theme(legend.position = "none")+
  # remove space between plot area and x axis
  
  # Custom the theme:
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm")) +
  theme(panel.background = element_rect(fill = 'white')) +
  ylab(expression(bold(atop(paste("WUE"[i]," (\u03B4"^"13", "C)"), paste("transcriptome-wide association ", (r[s]))))))+ 
  theme(axis.text.y = element_text(family = "Arial", color="black", size=20, face="bold")) +
  theme(axis.text.x = element_text(family = "Arial", color="black", size=25, face="bold")) +
  theme(axis.title.y = element_text(family = "Arial", color="black", size=25, face="bold")) +
  theme(axis.title.x=element_blank())+
  theme(plot.title = element_text(family = "Arial", color="#696969",face="bold",  size=35)) +  
  ggtitle("TWAS")
p2
png("panels/p2.png", width = 17, height = 7, units = 'in', res = 350)
p2
dev.off() 

#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################Fig. 6C


df4 = read_csv("data/Iberian_d13C.csv")
head(phenotypes)


Phenotypesselected<- select(df4, accession_id, Dittberner_moleco18_Delta_13C)
head(Phenotypesselected)
transposedfinal = read_csv("data/transposedfinal.csv")
Transcriptselected<- select(transposedfinal, accession_id,AT3G03850)
head(Transcriptselected)
df4<-merge(Phenotypesselected, Transcriptselected, by = "accession_id", sort = TRUE)
head(df4)

###Plot1= altitude vs leaf area time potential, short days well-watered
formula <- y ~ x
p4 <- ggplot(df4, aes(x=AT3G03850 , y= Dittberner_moleco18_Delta_13C)) +
  geom_point(color='black',alpha = 0.4,  size=12)+  
  geom_smooth(method = "lm", se = T, fill="black", colour="black", size=0.5, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                         label.x.npc = "right", label.y.npc = "top",
                                                                                                         formula = formula, parse = TRUE, size = 6, colour="black") +  theme_bw()# +

p4<- p4 + theme(axis.line = element_line(size=1, colour = "black"),
                panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm")) 
p4<- p4 + 
  theme(axis.text = element_text(family = "Arial", color="black", face="bold",size=20)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=30)) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(5, 5, 5,5)
  )+ 
  #ylim(-1,1)+
  ylab(expression(bold(paste("WUE"[i]," (\u03B4"^"13", "C (\u2030))"))))+   
  xlab(expression(bold(paste(bolditalic("SAUR26")~"(FPKM)"))))
p4
png("panels/p4.png", width = 11, height = 7, units = 'in', res = 350)
p4
dev.off() 



#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#####################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

df5 = read_csv("data/Iberian_d13C.csv")
head(phenotypes)


Phenotypesselected<- select(df5, accession_id, Dittberner_moleco18_Delta_13C)
head(Phenotypesselected)
transposedfinal = read_csv("data/transposedfinal.csv")
Transcriptselected<- select(transposedfinal, accession_id,AT5G38420)
head(Transcriptselected)
df5<-merge(Phenotypesselected, Transcriptselected, by = "accession_id", sort = TRUE)
head(df5)

###Plot1= altitude vs leaf area time potential, short days well-watered
formula <- y ~ x
p5 <- ggplot(df5, aes(x=AT5G38420 , y= Dittberner_moleco18_Delta_13C)) +
  geom_point(color='black',alpha = 0.4,  size=12)+  
  geom_smooth(method = "lm", se = T, fill="black", colour="black", size=0.5, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                         label.x.npc = "right", label.y.npc = "top",
                                                                                                         formula = formula, parse = TRUE, size = 6, colour="black") +  theme_bw()# +

p5<- p5 + theme(axis.line = element_line(size=1, colour = "black"),
                panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm")) 
p5<- p5 + 
  theme(axis.text = element_text(family = "Arial", color="black", face="bold",size=20)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=30)) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(5, 5, 5,5)
  )+ 
  #ylim(-1,1)+
  ylab(expression(bold(paste("WUE"[i]," (\u03B4"^"13", "C)"))))+  
  xlab(expression(bold(paste(bolditalic("RBCS2B")~"(FPKM)"))))
  
p5
png("panels/p5.png", width = 11, height = 7, units = 'in', res = 350)
p5
dev.off() 

#######################################################################################################################################
#######################################################################################################################################


df5a = read_csv("data/Phenotypesb.csv")
head(df5a)


Phenotypesselected<- select(df5a, accession_id, iWUE)
head(Phenotypesselected)
transposedfinal = read_csv("data/transposedfinal.csv")
Transcriptselected<- select(transposedfinal, accession_id,AT3G03850)
head(Transcriptselected)
df5a<-merge(Phenotypesselected, Transcriptselected, by = "accession_id", sort = TRUE)
head(df5a)

###Plot1= altitude vs leaf area time potential, short days well-watered
formula <- y ~ x
p5a <- ggplot(df5a, aes(x=AT3G03850 , y= iWUE)) +
  geom_point(color='black',alpha = 0.4,  size=12)+  
  geom_smooth(method = "lm", se = T, fill="black", colour="black", size=0.5, alpha = 0.1) + stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label..,  sep = "~~~")), 
                                                                                                         label.x.npc = "right", label.y.npc = "top",
                                                                                                         formula = formula, parse = TRUE, size = 6, colour="black") +  theme_bw()# +

p5a<- p5a + theme(axis.line = element_line(size=1, colour = "black"),
                panel.grid.major = element_line(colour = "#d3d3d3"), panel.grid.minor = element_blank(),
                panel.border = element_blank(), panel.background = element_blank(), axis.ticks.length=unit(0.3,"cm")) 
p5a<- p5a + 
  theme(axis.text = element_text(family = "Arial", color="black", face="bold",size=20)) +
  theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=30)) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(5, 5, 5,5)
  )+ 
  #ylim(-1,1)+
  ylab(expression(bold(paste("WUE"[i]," (\u03B4"^"13", "C)"))))+  
  xlab(expression(bold(paste(bolditalic("SAUR26")~"(FPKM)"))))

p5a
png("panels/p5a.png", width = 11, height = 7, units = 'in', res = 350)
p5a
dev.off() 
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#####################################################################################################################################
####################


M1 <- image_read("panels/p1.png")
M2 <- image_read("panels/p2.png")
M3 <- image_read("panels/GO_terms_Tree.png")
M4 <- image_read("panels/p4.png")
M5 <- image_read("panels/p5.png")

plot_grid(rasterGrob(M5),rasterGrob(M4),labels=c('d',"e"), greedy = TRUE, ncol =1, nrow = 2, align = 'v', vjust=0.2, label_size=65)
png("panels/bottomright.png", width = 8, height = 10, units = 'in', res = 350)
plot_grid(rasterGrob(M5),rasterGrob(M4),labels=c('d',"e"), greedy = TRUE, ncol =1, nrow = 2, align = 'v', vjust=0.2, label_size=65)
dev.off()

bottomright <- image_read("panels/bottomright.png")

plot_grid(rasterGrob(M3),rasterGrob(bottomright),labels=c('c',""), greedy = TRUE, ncol =2, nrow = 1, align = 'v', hjust=0, label_size=55)
png("panels/bottom.png", width = 11, height = 6, units = 'in', res = 350)
plot_grid(rasterGrob(M3),rasterGrob(bottomright),labels=c('c',""), greedy = TRUE, ncol =2, nrow = 1, align = 'v', hjust=0, label_size=55)
dev.off()
#top <- image_read("top.png")
bottom <- image_read("panels/bottom.png")

plot_grid(rasterGrob(M1),rasterGrob(M2),rasterGrob(bottom),labels=c('a', 'b',""), greedy = TRUE, ncol =1, nrow = 3, align = 'v', hjust=0, label_size=35)
png("Figure 5.png", width = 10, height = 14, units = 'in', res = 350)
plot_grid(rasterGrob(M1),rasterGrob(M2),rasterGrob(bottom),labels=c('a', 'b',""), greedy = TRUE, ncol =1, nrow = 3, align = 'v', hjust=0, label_size=35)
dev.off()



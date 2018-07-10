### a ajuster pour ta machine ##
setwd("/Users/vincentcastric/Dropbox/Dossier partagé VC-NB/Thèse/1-Chapitre 1/data_script_R/figures")
################################


library(weights)
library(ggplot2)
library(RColorBrewer)
library(csvread)
library(vioplot)
library(stats)
library(lattice)
library(nlme)
library(scales)
library(lme4)
library(lmerTest)

################################

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
italic.labs<-element_text(face= "italic")



############################################
### Figure 1 (expression ~ stage) ##########
############################################
stage=read.table("stages31_08_2017.txt", h=T)
SRK=read.table("stages_SRK_31_08_2017.txt", h=T)

p1<-ggplot()+
  geom_line(aes(x=stage$x, y=stage$y, color=as.factor(stage$Allele)), size=1.15)+
  guides(colour=guide_legend(title="Allele"))+
  scale_colour_manual(values=c("S01"="coral3","S02"="chocolate1", "S03"="darkgoldenrod2", "S04"="darkolivegreen4", "S10"="LightBlue",
                               "S12"="dodgerblue3", "S13"="darkorchid3", "S20"="forestgreen", "S29"="Black"))+
  labs(title=("SCR"))+
  xlab("developmental stage")+
  ylab("Relative expression across stages")+
  theme(panel.background = element_rect(fill = "white"),axis.ticks.x=element_blank(), legend.key=element_rect(fill="white"), legend.text=element_text(size=12), 
        axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),legend.position="none",
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        legend.title=element_text(size=14),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_reverse(breaks=c(1,2,3,4), labels=c("A","B","C","D"))
p1

## normalisé par la médiane des valeurs au stade avec niveau d'expression le plus élévé
pt2<-ggplot(stage, aes(x=x, y=y3, fill=as.factor(x)))+

  geom_boxplot()+
  scale_fill_manual(values=c("Grey","Grey","Grey","Grey"))+
  labs(title=("SCR"))+
  xlab("flower development stage (mm)")+
  ylab("Relative expression across stages")+
  theme(panel.background = element_rect(fill = "white"),axis.ticks.x=element_blank(), legend.key=element_rect(fill="white"), legend.text=element_text(size=12), 
        axis.title.x=element_text(size=15, color="Black"), axis.title.y=element_text(size=15, color="Black"),legend.position="none",
        axis.text.x=element_text(size=15, color="Black"), axis.text.y=element_text(size=15, color="Black"),
        legend.title=element_text(size=15, color="Black"),
        axis.line.x=element_line(colour = "Black", size = 1),
        plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_reverse(breaks=c(1,2,3,4), labels=c("open flower",">1","0.5-1","<0.5"))

pt2


## Stage SRK ##
p2<-ggplot()+
  geom_line(aes(x=SRK$x, y=SRK$y, color=as.factor(SRK$Allele)), size=1.15)+
  labs(title=("SRK"))+
  xlab("developmental stage")+
  ylab("Relative expression across stages")+
  scale_colour_manual(values=c("S01"="coral3", "S03"="darkgoldenrod2", "S04"="darkolivegreen4", "S10"="LightBlue",
                               "S12"="dodgerblue3", "S29"="Black"))+
  
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), 
        legend.position="none",
        axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_reverse(breaks=c(1,2,3,4), labels=c("A","B","C","D"))
p2

##normalisé par la médianne des valeurs au stage  avec niveau d'expression le plus élevé
p2t<-ggplot(SRK, aes(x=x, y=y4, fill=as.factor(x)))+
  geom_boxplot()+
  labs(title=("SRK"))+
  xlab("flower development stage (mm)")+
  ylab("Relative expression")+
  scale_fill_manual(values=c("Grey","Grey","Grey","Grey"))+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), 
        legend.position="none",
        plot.title=italic.labs,
        axis.title.x=element_text(size=15, color="Black"), axis.title.y=element_text(size=15, color="Black"),
        axis.text.x=element_text(size=15, color="Black"), axis.text.y=element_text(size=15, color="Black"),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_reverse(breaks=c(1,2,3,4), labels=c("open flower",">1","0.5-1","<0.5"))
p2t 



multiplot(pt2, p2t, cols=1)

############################################
### Figure 2 (expression SCR) ##############
############################################
exp<-read.table("expression_C&D_modif2.txt", h=T)
italic.labs<-element_text(face= "italic")
## subset par allèle ##
s1<-subset(exp, allele=="SCR01")
s2<-subset(exp, allele=="SCR02")
s3<-subset(exp, allele=="SCR03")
s4<-subset(exp, allele=="SCR04")
s10<-subset(exp, allele=="SCR10")
s12<-subset(exp, allele=="SCR12")
s13<-subset(exp, allele=="SCR13")
s20<-subset(exp, allele=="SCR20")
s29<-subset(exp, allele=="SCR29")


## dans le code, il est possible de changer y=ct par y=ct2 ou ct2, tout dépend
## représentation souhaité (normalisé par la valeur d'expression plus haute par allele (ct2)
##ou par la moyenne la plus haute par allele (ct3), médianne la plus élévée (ct4)
## ou pas normalisé). Dans l'un ou l'autre, penser à ajuster dans geom_label la valeur
## de y, pour que les annotation s'intègre bien sur le graphique.

##SCR01##

p<- ggplot(s1, aes(x=other_allele, y=ct4, fill=as.factor(etat)))+
  geom_boxplot()+
  geom_label(label=s1$info, y=7.5,label.size = NA, fill="white")+ ##0.03
  scale_fill_manual(values=c("Grey","Grey","Grey"))+
  guides(fill=FALSE)+
  labs(title=("SCR01"))+
  xlab("")+
  ylab("Relative Expression")+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        plot.title = italic.labs,
        axis.text.x=element_text(size=9, angle =45, hjust=1),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_discrete(limits=c("Hom","SCR03","SCR29","SCR02","SCR10",
                            "SCR04","SCR12","SCR13","SCR20"), 
                   labels= c("SCR01(Hom)","SCR03","SCR29","SCR02","SCR10",
                             "SCR04","SCR12","SCR13","SCR20"))
#p
##SCR02##

p2<- ggplot(s2, aes(x=other_allele, y=ct4, fill=as.factor(etat)))+
  geom_boxplot()+
  geom_label(label=s2$info, y=1.6,label.size = NA, fill="white")+ ##0.4
  scale_fill_manual(values=c("Grey","Grey","Grey"))+
  labs(title=("SCR02"))+
  guides(fill=FALSE)+
  xlab("")+
  ylab("Relative Expression")+
  geom_text(x=6, y=0.03, label="Na")+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        plot.title = italic.labs,
        axis.text.x=element_text(size=9, angle =45, hjust=1),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_discrete(limits=c("SCR01","SCR03","SCR29","Hom","SCR10",
                            "SCR04","SCR12","SCR13","SCR20"), 
                   labels= c("SCR01","SCR03","SCR29","SCR02 (Hom)","SCR10",
                             "SCR04","SCR12","SCR13","SCR20"))
#p2
##SCR03##

p3<- ggplot(s3, aes(x=other_allele, y=ct4, fill=as.factor(etat)))+
  geom_boxplot()+
  geom_label(label=s3$info, y=2.3,label.size = NA, fill="white")+ ##0.06
  scale_fill_manual(values=c("Grey","Grey","Grey"))+
  guides(fill=FALSE)+
  xlab("")+
  ylab("Relative Expression")+
  labs(title=("SCR03"))+
  geom_text(x=5, y=0.005, label="Na")+
  geom_text(x=8, y=0.005, label="Na")+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        plot.title = italic.labs,
        axis.text.x=element_text(size=9, angle =45, hjust=1),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_discrete(limits=c("SCR01","Hom","SCR29","SCR02","SCR10",
                            "SCR04","SCR12","SCR13","SCR20"), 
                   labels= c("SCR01","SCR03(Hom)","SCR29","SCR02","SCR10",
                             "SCR04","SCR12","SCR13","SCR20"))
#p3
##SCR04##

p4<- ggplot(s4, aes(x=other_allele, y=ct4, fill=as.factor(etat)))+
  geom_boxplot()+
  geom_label(label=s4$info, y=2.1 ,label.size = NA, fill="white")+ ##1
  scale_fill_manual(values=c("Grey","Grey","Grey"))+
  guides(fill=FALSE)+
  xlab("")+
  ylab("Relative Expression")+
  labs(title=("SCR04"))+
  geom_text(x=4, y=0.07, label="Na")+
  geom_text(x=6, y=0.07, label="Na")+
  theme(panel.background = element_rect(fill = "white"),
        plot.title = italic.labs,
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        axis.text.x=element_text(size=9, angle =45, hjust=1),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_discrete(limits=c("SCR01","SCR03","SCR29","SCR02","SCR10",
                            "Hom","SCR12","SCR13","SCR20"), 
                   labels= c("SCR01","SCR03","SCR29","SCR02","SCR10",
                             "SCR04(Hom)","SCR12","SCR13","SCR20"))

#p4
##SCR10##

p5<- ggplot(s10, aes(x=other_allele, y=ct4, fill=as.factor(etat)))+
  geom_boxplot()+
  geom_label(label=s10$info, y=2.9,label.size = NA, fill="white")+ ##0.015
  scale_fill_manual(values=c("Grey","Grey","Grey"))+
  guides(fill=FALSE)+
  xlab("")+
  ylab("Relative Expression")+
  labs(title=("SCR10"))+
  geom_text(x=2, y=0.1, label="Na")+
  geom_text(x=5, y=0.1, label="Na")+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        plot.title = italic.labs,
        axis.text.x=element_text(size=9, angle =45, hjust=1),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_discrete(limits=c("SCR01","SCR03","SCR29","SCR02","Hom",
                            "SCR04","SCR12","SCR13","SCR20"), 
                   labels= c("SCR01","SCR03","SCR29","SCR02","SCR10(Hom)",
                             "SCR04","SCR12","SCR13","SCR20"))

#p5
##SCR12##

p6<- ggplot(s12, aes(x=other_allele, y=ct4, fill=as.factor(etat)))+
  geom_boxplot()+
  geom_label(label=s12$info, y=2.15,label.size = NA, fill="white")+##0.35
  scale_fill_manual(values=c("Grey","Grey","Grey"))+
  guides(fill=FALSE)+
  xlab("")+
  ylab("Relative Expression")+
  labs(title=("SCR12"))+
  geom_text(x=7, y=0.04, label="Na")+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        plot.title = italic.labs,
        axis.text.x=element_text(size=9, angle =45, hjust=1),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_discrete(limits=c("SCR01","SCR03","SCR02","SCR04","SCR29",
                            "SCR10","Hom","SCR13","SCR20"), 
                   labels= c("SCR01","SCR03","SCR02","SCR04","SCR29",
                             "SCR10","SCR12(Hom)","SCR13","SCR20"))

#p6
##SCR13##

p7<- ggplot(s13, aes(x=other_allele, y=ct4, fill=as.factor(etat)))+
  geom_boxplot()+
  geom_label(label=s13$info, y=1.25,label.size = NA, fill="white")+ ##0.085
  scale_fill_manual(values=c("Grey","Grey","Grey"))+
  guides(fill=FALSE)+
  xlab("")+
  ylab("Relative Expression")+
  labs(title=("SCR13"))+
  geom_text(x=2, y=0.03, label="Na")+
  geom_text(x=8, y=0.03, label="Na")+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        plot.title = italic.labs,
        axis.text.x=element_text(size=9, angle =45, hjust=1),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_discrete(limits=c("SCR01","SCR03","SCR02","SCR04","SCR29",
                            "SCR10","SCR12","Hom","SCR20"), 
                   labels= c("SCR01","SCR03","SCR02","SCR04","SCR29",
                             "SCR10","SCR12","SCR13(Hom)","SCR20"))

#p7
##SCR20##

p8<- ggplot(s20, aes(x=other_allele, y=ct4, fill=as.factor(etat)))+
  geom_boxplot()+
  geom_label(label=s20$info, y=5.7,label.size = NA, fill="white")+##0.6
  scale_fill_manual(values=c("Grey","Grey","Grey"))+
  guides(fill=FALSE)+
  xlab("")+
  ylab("Relative Expression")+
  labs(title=("SCR20"))+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        plot.title = italic.labs,
        axis.text.x=element_text(size=9, angle =45, hjust=1),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_discrete(limits=c("SCR01","SCR03","SCR02","SCR04","SCR29",
                            "SCR10","SCR12","SCR13","Hom"), 
                   labels= c("SCR01","SCR03","SCR02","SCR04","SCR29",
                             "SCR10","SCR12","SCR13","SCR20(Hom)"))

#p8
##SCR29##

p9<- ggplot(s29, aes(x=other_allele, y=ct4, fill=as.factor(etat)))+
  geom_boxplot()+
  geom_label(label=s29$info, y=2,label.size = NA, fill="white")+
  scale_fill_manual(values=c("Grey","Grey","Grey"))+
  #scale_fill_manual(values=c("Dom"="cornflowerblue","Rec"="Red","Unk"="White"))+
  guides(fill=FALSE)+
  xlab("")+
  ylab("Relative Expression")+
  labs(title=("SCR29"))+
  geom_text(x=3, y=0.02, label="Na")+
  theme(panel.background = element_rect(fill = "white"),
        plot.title = italic.labs,
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        axis.text.x=element_text(size=9, angle =45, hjust=1),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_discrete(limits=c("SCR01","SCR03","Hom","SCR02","SCR10",
                            "SCR04","SCR12","SCR13","SCR20"), 
                   labels= c("SCR01","SCR03","SCR29(Hom)","SCR02","SCR10",
                             "SCR04","SCR12","SCR13","SCR20"))



#p9

multiplot(p, p2, p6, p3, p5, p7, p9, p4, p8, cols=3)

############################################
### Figure 3 (expression pour SRK) #########
############################################
SRK<-read.table("expression_SRK2.txt", h=T)
srk1<-subset(SRK, allele=="SRK01")
srk3<-subset(SRK, allele=="SRK03")
srk4<-subset(SRK, allele=="SRK04")
srk10<-subset(SRK, allele=="SRK10")
srk12<-subset(SRK, allele=="SRK12")
srk29<-subset(SRK, allele=="SRK29")

q1<- ggplot(srk1, aes(x=other_allele, y=ct2, fill=as.factor(etat)))+
  geom_boxplot()+
  geom_label(label=srk1$info, y=5,label.size = NA, fill="white", size =5)+
  scale_fill_manual(values=c("Grey","Grey","Grey"))+
  #scale_fill_manual(values=c("Dom"="cornflowerblue","Rec"="Red","Unk"="White"))+
  guides(fill=FALSE)+
  labs(title=("SRK01"))+
  xlab("")+
  ylab("Relative Expression")+
  geom_text(x=6, y=0.04, label="Na")+
  geom_text(x=8, y=0.04, label="Na")+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        plot.title = italic.labs,
        axis.title =element_text(size=12),
        axis.text.x=element_text(size=11, angle =45, hjust=1),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_discrete(limits=c("Hom","SRK03","SRK29","SRK02","SRK10",
                            "SRK04","SRK12","SRK13","SRK20"), 
                   labels= c("SRK01(Hom)","SRK03","SRK29","SRK02","SRK10",
                             "SRK04","SRK12","SRK13","SRK20"))
q1

q3<- ggplot(srk3, aes(x=other_allele, y=ct2, fill=as.factor(etat)))+
  geom_boxplot()+
  geom_label(label=srk3$info, y=1.7,label.size = NA, fill="white", size =5)+
  scale_fill_manual(values=c("Grey","Grey","Grey"))+
  guides(fill=FALSE)+
  labs(title=("SRK03"))+
  xlab("")+
  ylab("Relative Expression")+
  geom_text(x=5, y=0.2, label="Na")+
  geom_text(x=8, y=0.2, label="Na")+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        plot.title = italic.labs,
        axis.title =element_text(size=12),
        axis.text.x=element_text(size=11, angle =45, hjust=1),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_discrete(limits=c("SRK01","Hom","SRK29","SRK02","SRK10",
                            "SRK04","SRK12","SRK13","SRK20"), 
                   labels= c("SRK01","SRK03(Hom)","SRK29","SRK02","SRK10",
                             "SRK04","SRK12","SRK13","SRK20"))
q3

q29<- ggplot(srk29, aes(x=other_allele, y=ct2, fill=as.factor(etat)))+
  geom_boxplot()+
  geom_label(label=srk29$info, y=2.5,label.size = NA, fill="white", size =5)+
  scale_fill_manual(values=c("Grey","Grey","Grey"))+
  #scale_fill_manual(values=c("Dom"="cornflowerblue","Rec"="Red","Unk"="White"))+
  guides(fill=FALSE)+
  labs(title=("SRK29"))+
  xlab("")+
  ylab("Relative Expression")+
  geom_text(x=1, y=0.2, label="Na")+
  geom_text(x=3, y=0.2, label="Na")+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        plot.title = italic.labs,
        axis.title =element_text(size=12),
        axis.text.x=element_text(size=11, angle =45, hjust=1),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_discrete(limits=c("SRK01","SRK03","Hom","SRK02","SRK10",
                            "SRK04","SRK12","SRK13","SRK20"), 
                   labels= c("SRK01","SRK03","SRK29(Hom)","SRK02","SRK10",
                             "SRK04","SRK12","SRK13","SRK20"))
q29

q10<- ggplot(srk10, aes(x=other_allele, y=ct2, fill=as.factor(etat)))+
  geom_boxplot()+
  geom_label(label=srk10$info, y=3,label.size = NA, fill="white", size =5)+
  scale_fill_manual(values=c("Grey","Grey","Grey","Grey"))+
  guides(fill=FALSE)+
  labs(title=("SRK10"))+
  xlab("")+
  ylab("Relative Expression")+
  geom_text(x=2, y=0.1, label="Na")+
  geom_text(x=5, y=0.1, label="Na")+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        plot.title = italic.labs,
        axis.title =element_text(size=12),
        axis.text.x=element_text(size=11, angle =45, hjust=1),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_discrete(limits=c("SRK01","SRK03","SRK29","SRK02","Hom",
                            "SRK04","SRK12","SRK13","SRK20"), 
                   labels= c("SRK01","SRK03","SRK29","SRK02","SRK10(Hom)",
                             "SRK04","SRK12","SRK13","SRK20"))
q10

q4<- ggplot(srk4, aes(x=other_allele, y=ct2, fill=as.factor(etat)))+
  geom_boxplot()+
  geom_label(label=srk4$info, y=3.2,label.size = NA, fill="white", size =5)+
  scale_fill_manual(values=c("Grey","Grey","Grey"))+
  guides(fill=FALSE)+
  labs(title=("SRK4"))+
  xlab("")+
  ylab("Relative Expression")+
  geom_text(x=1, y=0.2, label="Na")+
  geom_text(x=4, y=0.2, label="Na")+
  geom_text(x=6, y=0.2, label="Na")+
  geom_text(x=9, y=0.2, label="Na")+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        plot.title = italic.labs,
        axis.title =element_text(size=12),
        axis.text.x=element_text(size=11, angle =45, hjust=1),
        axis.text.y=element_text(size=11),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_discrete(limits=c("SRK01","SRK03","SRK29","SRK02","SRK10",
                            "Hom","SRK12","SRK13","SRK20"), 
                   labels= c("SRK01","SRK03","SRK29","SRK02","SRK10",
                             "SRK04(Hom)","SRK12","SRK13","SRK20)"))
q4

q12<- ggplot(srk12, aes(x=other_allele, y=ct2, fill=as.factor(etat)))+
  geom_boxplot()+
  geom_label(label=srk12$info, y=0.23,label.size = NA, fill="white", size =5)+
  scale_fill_manual(values=c("Grey","Grey","Grey"))+
  #scale_fill_manual(values=c("Dom"="cornflowerblue","Rec"="Red","Unk"="White"))+
  guides(fill=FALSE)+
  labs(title=("SRK12"))+
  xlab("")+
  ylab("Relative Expression")+
  geom_text(x=4, y=0.03, label="Na")+
  geom_text(x=5, y=0.03, label="Na")+
  geom_text(x=7, y=0.03, label="Na")+
  geom_text(x=8, y=0.03, label="Na")+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        plot.title = italic.labs,
        axis.text.x=element_text(size=11, angle =45, hjust=1),
        axis.text.y=element_text(size=11),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_discrete(limits=c("SRK01","SRK03","SRK29","SRK02","SRK10",
                            "SRK04","Hom","SRK13","SRK20"), 
                   labels= c("SRK01","SRK03","SRK29","SRK02","SRK10",
                             "SRK04","SRK12(Hom)","SRK13","SRK20"))
q12

multiplot(q1, q29 , q4, q3, q10, cols=2)


############################################
### Figure 4 (expression ~ score et AIC) ###
############################################
s=read.table("score.txt", h=T)

ggplot()+
  geom_point(data=s, aes(x=score, y=y, colour=as.factor(etat), shape=as.factor(etat), size = as.factor(etat), fill=as.factor(etat)))+
  geom_smooth(data=s, aes(x=score, y=y), span=0.5, colour="Black", size=1.1)+
  scale_colour_manual(values=c("Dominant"="black", "Recessive"="black", "Unknown"="Grey"))+
  scale_fill_manual(values=c("black", "White", "Black"))+
  scale_shape_manual(values=c(16, 21, 16))+
  scale_size_manual(values=c(3, 3, 3))+
  guides(fill=FALSE)+
  labs(title=(""))+
  xlab("target score")+
  ylab("relative expression")+
  theme(panel.background = element_rect(fill = "white"),
        legend.background= element_rect(fill="white"),
        legend.text=element_text(size=15, colour="Black"),
        legend.key=element_rect(fill="white"),
        legend.title=element_blank(),
        axis.title.x = element_text(size=20, colour="black"),
        axis.title.y = element_text(size=20,colour="black"),
        axis.text.x=element_text(size=20, colour="black"),
        axis.text.y=element_text(size=20, colour="black"),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_continuous(breaks = c(13, 14, 15, 16, 17, 18, 19, 20, 21, 22))

s3<-read.table("AIC_score.txt", h=T)

ggplot(s3, aes(x=score, y=AIC, fill=as.factor(score)))+
  geom_col()+
  scale_fill_manual(values=c("Black","Black","Black","Black","Black","Black","Black","Black","Black"))+
  scale_x_continuous(breaks=c(14,15,16,17,18,19,20,21,22))+
  coord_cartesian(ylim=c(4000,4070))+
  theme(panel.background = element_rect(fill = "white"),
        legend.background= element_rect(fill="white"),
        legend.text=element_text(size=15, colour="Black"),
        legend.key=element_rect(fill="white"),
        legend.title=element_blank(),
        legend.position="none",
        axis.title.x = element_text(size=20, colour="black"),
        axis.title.y = element_text(size=20,colour="black"),
        axis.text.x=element_text(size=15, colour="black"),
        axis.text.y=element_text(size=15, colour="black"),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))

############################################
### Figure S1 (controle négatif) ###########
############################################
exp2<-read.table("plq2.txt", h=T)


ggplot(exp2, aes(x=Etat, y=Cp, fill=as.factor(Etat)))+
  geom_boxplot()+
  scale_fill_brewer(type = "seq", palette = "Greys")+
  guides(fill=FALSE)+
  labs(title=(""))+
  xlab("")+
  ylab("Ct SCR")+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        axis.title.x = element_text(size=20, colour="black"),
        axis.title.y = element_text(size=20,colour="black"),
        axis.text.y=element_text(size=15, colour="black"),
        axis.text.x=element_text(size=15, colour="black"),
        axis.line.x=element_line(colour = "Black", size = 1),
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_x_discrete(limits=c("positif", "sample","negatif"),label=c("Positive Control","Cross Amplification", "Water"))

############################################
### Figure S2 (Dilution) ###################
############################################

#fonction pour transformer l'axe des x en log et inverse
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

#Fonction pour ajouter l'équation sur le graph, en faisant appelle à annotate (x,y,label=lm_eqn(df), parse=TRUE)
#où df=d1,d2 etc
lm_eqn <- function(df){
  m <- lm(ct~log10(dilution), df);
  eq <- substitute(~~italic(r)^2~"="~r2,
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

lm_eqn2 <- function(df){
  m <- lm(ct2~log10(dilution), df);
  eq <- substitute(~italic(r)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}


##exp(("b"-d1d$ct/"a"))/500
## ** = a*ln(x)+b
## b aussi donnée par la formule coef...

d1=read.table("SCR01_control.txt", h=T)
coef(lm(-d1$ct~log10(d1$dilution)))
coef(lm(-d1$ct2~log10(d1$dilution)))

s1<-ggplot()+
  geom_point(data=d1,aes(x=d1$dilution,y=d1$ct), size=1.15)+
  labs(title=("SCR01"))+
  geom_abline(intercept = -16.524348, slope = -3.678699)+
  geom_abline(intercept = -15.692534, slope = -4.296618, linetype=2)+
  annotate("text",x=1e-05, y= 20, label=lm_eqn(d1), parse=TRUE)+
  annotate("text",x=1e-05, y= 22, label=lm_eqn2(d1), parse=TRUE)+
  annotate("segment", x=1e-04, xend=0.5e-04, y=20.2, yend=20.2)+
  annotate("segment", x=1e-04, xend=0.5e-04, y=22.2, yend=22.2, linetype=2)+
  xlab("dilution")+
  ylab("ct")+
  theme(panel.background = element_rect(fill = "white"),axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.line.x=element_line(colour = "Black", size = 1),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_y_reverse()+
  scale_x_continuous(trans=reverselog_trans(10), breaks=c(1,1e-02,1e-04,1e-06))

s1 
###SCR02###

d2=read.table("SCR02_control.txt", h=T)
coef(lm(-d2$ct~log10(d2$dilution)))
coef(lm(-d2$ct2~log10(d2$dilution)))

s3<-ggplot()+
  geom_point(data=d2,aes(x=d2$dilution,y=d2$ct), size=1.15)+
  labs(title=("SCR02"))+
  geom_abline(intercept = -18.04188, slope = -4.69425)+
  geom_abline(intercept = -17.540278, slope = -5.101806, linetype=2)+
  annotate("text",x=1e-05, y= 20, label=lm_eqn(d2), parse=TRUE)+
  annotate("text",x=1e-05, y= 22, label=lm_eqn2(d2), parse=TRUE)+
  annotate("segment", x=1e-04, xend=0.5e-04, y=20.2, yend=20.2)+
  annotate("segment", x=1e-04, xend=0.5e-04, y=22.2, yend=22.2, linetype=2)+
  xlab("dilution")+
  ylab("ct")+
  theme(panel.background = element_rect(fill = "white"),axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.line.x=element_line(colour = "Black", size = 1),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_y_reverse()+
  scale_x_continuous(trans=reverselog_trans(10), breaks=c(1,1e-02,1e-04,1e-06))

s3
###SCR03###

d3=read.table("SCR03_control.txt", h=T)
coef(lm(-d3$ct~log10(d3$dilution)))
coef(lm(-d3$ct2~log10(d3$dilution)))



s2<-ggplot()+
  geom_point(data=d3,aes(x=d3$dilution,y=d3$ct), size=1.15)+
  labs(title=("SCR03"))+
  geom_abline(intercept = -18.410500, slope = -2.787438)+
  geom_abline(intercept = -17.961250, slope = -3.124375, linetype=2)+
  annotate("text",x=1e-05, y= 20, label=lm_eqn(d3), parse=TRUE)+
  annotate("text",x=1e-05, y= 22, label=lm_eqn2(d3), parse=TRUE)+
  annotate("segment", x=1e-04, xend=0.5e-04, y=20.2, yend=20.2)+
  annotate("segment", x=1e-04, xend=0.5e-04, y=22.2, yend=22.2, linetype=2)+
  xlab("dilution")+
  ylab("ct")+
  theme(panel.background = element_rect(fill = "white"),axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.line.x=element_line(colour = "Black", size = 1),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_y_reverse()+
  scale_x_continuous(trans=reverselog_trans(10), breaks=c(1,1e-02,1e-04,1e-06))


s2
###SCR04### 

d4=read.table("SCR04_control.txt", h=T)
coef(lm(-d4$ct~log10(d4$dilution)))
coef(lm(-d4$ct2~log10(d4$dilution)))


s6<-ggplot()+
  geom_point(data=d4,aes(x=d4$dilution,y=d4$ct), size=1.15)+
  labs(title=("SCR04"))+
  geom_abline(intercept = -14.795788, slope = -3.458854)+
  geom_abline(intercept = -13.867222, slope = -4.155278, linetype=2)+
  annotate("text",x=1e-05, y= 20, label=lm_eqn(d4), parse=TRUE)+
  annotate("text",x=1e-05, y= 22, label=lm_eqn2(d4), parse=TRUE)+
  annotate("segment", x=1e-04, xend=0.5e-04, y=20.2, yend=20.2)+
  annotate("segment", x=1e-04, xend=0.5e-04, y=22.2, yend=22.2, linetype=2)+
  xlab("dilution")+
  ylab("ct")+
  theme(panel.background = element_rect(fill = "white"),axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.line.x=element_line(colour = "Black", size = 1),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_y_reverse()+
  scale_x_continuous(trans=reverselog_trans(10), breaks=c(1,1e-02,1e-04,1e-06))

s6
###SCR10###

d5=read.table("SCR10_control.txt", h=T)
coef(lm(-d5$ct~log10(d5$dilution)))
coef(lm(-d5$ct2~log10(d5$dilution)))


s5<-ggplot()+
  geom_point(data=d5,aes(x=d5$dilution,y=d5$ct), size=1.15)+
  labs(title=("SCR10"))+
  geom_abline(intercept = -17.097087, slope = -3.444977)+
  geom_abline(intercept = -16.792642, slope = -3.687583, linetype=2)+
  annotate("text",x=1e-05, y= 20, label=lm_eqn(d5), parse=TRUE)+
  annotate("text",x=1e-05, y= 22, label=lm_eqn2(d5), parse=TRUE)+
  annotate("segment", x=1e-04, xend=0.5e-04, y=20.2, yend=20.2)+
  annotate("segment", x=1e-04, xend=0.5e-04, y=22.2, yend=22.2, linetype=2)+
  xlab("dilution")+
  ylab("ct")+
  theme(panel.background = element_rect(fill = "white"),axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.line.x=element_line(colour = "Black", size = 1),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_y_reverse()+
  scale_x_continuous(trans=reverselog_trans(10), breaks=c(1,1e-02,1e-04,1e-06))

s5
###SCR12###

d6=read.table("SCR12_control.txt", h=T)
coef(lm(-d6$ct~log10(d6$dilution)))
coef(lm(-d6$ct2~log10(d6$dilution)))


s7<-ggplot()+
  geom_point(data=d6,aes(x=d6$dilution,y=d6$ct), size=1.15)+
  labs(title=("SCR12"))+
  geom_abline(intercept = -14.63338, slope = -2.88450)+
  geom_abline(intercept = -14.037292, slope = -3.331562, linetype=2)+
  annotate("text",x=1e-05, y= 20, label=lm_eqn(d6), parse=TRUE)+
  annotate("text",x=1e-05, y= 22, label=lm_eqn2(d6), parse=TRUE)+
  annotate("segment", x=1e-04, xend=0.5e-04, y=20.2, yend=20.2)+
  annotate("segment", x=1e-04, xend=0.5e-04, y=22.2, yend=22.2, linetype=2)+
  xlab("dilution")+
  ylab("ct")+
  theme(panel.background = element_rect(fill = "white"),axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.line.x=element_line(colour = "Black", size = 1),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_y_reverse()+
  scale_x_continuous(trans=reverselog_trans(10), breaks=c(1,1e-02,1e-04,1e-06))

###SCR13###

d7=read.table("SCR13_control.txt", h=T)
coef(lm(-d7$ct~log10(d7$dilution)))
coef(lm(-d7$ct2~log10(d7$dilution)))


s8<-ggplot()+
  geom_point(data=d7,aes(x=d7$dilution,y=d7$ct), size=1.15)+
  labs(title=("SCR13"))+
  geom_abline(intercept = -17.535750, slope = -2.583812)+
  annotate("text",x=1e-05, y= 20, label=lm_eqn(d7), parse=TRUE)+
  annotate("text",x=1e-05, y= 22, label=lm_eqn2(d7), parse=TRUE)+
  annotate("segment", x=1e-04, xend=0.5e-04, y=20.2, yend=20.2)+
  annotate("segment", x=1e-04, xend=0.5e-04, y=22.2, yend=22.2, linetype=2)+
  xlab("dilution")+
  ylab("ct")+
  theme(panel.background = element_rect(fill = "white"),axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.line.x=element_line(colour = "Black", size = 1),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_y_reverse()+
  scale_x_continuous(trans=reverselog_trans(10), breaks=c(1,1e-02,1e-04,1e-06))

s8
###SCR20###

d8=read.table("SCR20_control.txt", h=T)
coef(lm(-d8$ct~log10(d8$dilution)))
coef(lm(-d8$ct2~log10(d8$dilution)))


s9<-ggplot()+
  geom_point(data=d8,aes(x=d8$dilution,y=d8$ct), size=1.15)+
  labs(title=("SCR20"))+
  geom_abline(intercept = -17.101563, slope = -3.307621)+
  geom_abline(intercept = -16.498978, slope = -3.786949, linetype=2)+
  annotate("text",x=1e-05, y= 20, label=lm_eqn(d8), parse=TRUE)+
  annotate("text",x=1e-05, y= 23, label=lm_eqn2(d8), parse=TRUE)+
  annotate("segment", x=1e-04, xend=0.5e-04, y=20.2, yend=20.2)+
  annotate("segment", x=1e-04, xend=0.5e-04, y=23.2, yend=23.2, linetype=2)+
  xlab("dilution")+
  ylab("ct")+
  theme(panel.background = element_rect(fill = "white"),axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.line.x=element_line(colour = "Black", size = 1),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_y_reverse()+
  scale_x_continuous(trans=reverselog_trans(10), breaks=c(1,1e-02,1e-04,1e-06))

s9
###SCR29### 

d9=read.table("SCR29_control.txt", h=T)
coef(lm(-d9$ct~log10(d9$dilution)))
coef(lm(-d9$ct2~log10(d9$dilution)))


s4<-ggplot()+
  geom_point(data=d9,aes(x=d9$dilution,y=d9$ct), size=1.15)+
  labs(title=("SCR29"))+
  geom_abline(intercept = -16.65238, slope = -2.75150)+
  geom_abline(intercept = -15.549792, slope = -3.578437, linetype=2)+
  annotate("text",x=1e-05, y= 20, label=lm_eqn(d9), parse=TRUE)+
  annotate("text",x=1e-05, y= 22, label=lm_eqn2(d9), parse=TRUE)+
  annotate("segment", x=1e-04, xend=0.5e-04, y=20.2, yend=20.2)+
  annotate("segment", x=1e-04, xend=0.5e-04, y=22.2, yend=22.2, linetype=2)+
  xlab("dilution")+
  ylab("ct")+
  theme(panel.background = element_rect(fill = "white"),axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.line.x=element_line(colour = "Black", size = 1),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_y_reverse()+
  scale_x_continuous(trans=reverselog_trans(10), breaks=c(1,1e-02,1e-04,1e-06))

s4
multiplot(s1,s4,s7,s2,s5,s8,s3,s6,s9, cols=3)


### SRK01 ###

p1=read.table("SRK01_control.txt", h=T)
coef(lm(-p1$ct~log10(p1$dilution)))

s01<-ggplot()+
  geom_point(data=p1,aes(x=p1$dilution,y=p1$ct), size=1.15)+
  labs(title=("SRK01"))+
  geom_abline(intercept = -17.126473, slope = -2.140134)+
  annotate("text",x=1e-07, y= 20, label=lm_eqn(p1), parse=TRUE)+
  xlab("dilution")+
  ylab("ct")+
  theme(panel.background = element_rect(fill = "white"),axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.line.x=element_line(colour = "Black", size = 1),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_y_reverse()+
  scale_x_continuous(trans=reverselog_trans(10), breaks=c(1,1e-02,1e-04,1e-06, 1e-08))
s01


### SRK03 ###

p3=read.table("SRK03_control.txt", h=T)
coef(lm(-p3$ct~log10(p3$dilution)))


s03<-ggplot()+
  geom_point(data=p3,aes(x=p3$dilution,y=p3$ct), size=1.15)+
  labs(title=("SRK03"))+
  geom_abline(intercept = -11.045341, slope = -3.122928)+
  annotate("text",x=1e-05, y= 20, label=lm_eqn(p3), parse=TRUE)+
  xlab("dilution")+
  ylab("ct")+
  theme(panel.background = element_rect(fill = "white"),axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.line.x=element_line(colour = "Black", size = 1),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_y_reverse()+
  scale_x_continuous(trans=reverselog_trans(10), breaks=c(1e-02,1e-03,1e-04,1e-05,1e-06))
s03
### SRK04 ###

p4=read.table("SRK04_control.txt", h=T)
coef(lm(-p4$ct~log10(p4$dilution)))


s04<-ggplot()+
  geom_point(data=p4,aes(x=p4$dilution,y=p4$ct), size=1.15)+
  labs(title=("SRK04"))+
  geom_abline(intercept = -10.911993, slope = -3.078229)+
  annotate("text",x=2e-06, y= 20, label=lm_eqn(p4), parse=TRUE)+
  xlab("dilution")+
  ylab("ct")+
  theme(panel.background = element_rect(fill = "white"),axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.line.x=element_line(colour = "Black", size = 1),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_y_reverse()+
  scale_x_continuous(trans=reverselog_trans(10), breaks=c(1e-04,1e-05,1e-06))
s04



###SRK10###
p10=read.table("SRK10_control.txt", h=T)
coef(lm(-p10$ct~log10(p10$dilution)))



s10<-ggplot()+
  geom_point(data=p10,aes(x=p10$dilution,y=p10$ct), size=1.15)+
  labs(title=("SRK10"))+
  geom_abline(intercept = -9.648462, slope = -4.385960)+
  annotate("text",x=2e-5, y= 20, label=lm_eqn(p10), parse=TRUE)+
  xlab("dilution")+
  ylab("ct")+
  theme(panel.background = element_rect(fill = "white"),axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.line.x=element_line(colour = "Black", size = 1),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_y_reverse()+
  scale_x_continuous(trans=reverselog_trans(10), breaks=c(1e-03,1e-04,1e-05,1e-06))
s10

###SRK12###
p12=read.table("SRK12_control.txt", h=T)
coef(lm(-p12$ct~log10(p12$dilution)))


s12<-ggplot()+
  geom_point(data=p12,aes(x=p12$dilution,y=p12$ct), size=1.15)+
  labs(title=("SRK12"))+
  geom_abline(intercept = -9.926250, slope = -3.675625)+
  annotate("text",x=1e-05, y= 10, label=lm_eqn(p12), parse=TRUE)+
  xlab("dilution")+
  ylab("ct")+
  theme(panel.background = element_rect(fill = "white"),axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.line.x=element_line(colour = "Black", size = 1),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_y_reverse()+
  scale_x_continuous(trans=reverselog_trans(10), breaks=c(1,1e-02,1e-04,1e-06))
s12
###SRK29###

p29=read.table("SRK29_control.txt", h=T)
coef(lm(-p29$ct~log10(p29$dilution)))


s29<-ggplot()+
  geom_point(data=p29,aes(x=p29$dilution,y=p29$ct), size=1.15)+
  labs(title=("SRK29"))+
  geom_abline(intercept = -5.859934, slope = -5.139583)+
  annotate("text",x=1e-04, y= 15, label=lm_eqn(p29), parse=TRUE)+
  xlab("dilution")+
  ylab("ct")+
  theme(panel.background = element_rect(fill = "white"),axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12), axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12), axis.text.y=element_text(size=12),
        axis.line.x=element_line(colour = "Black", size = 1),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))+
  scale_y_reverse()+
  scale_x_continuous(trans=reverselog_trans(10), breaks=c(1e-01,1e-02,1e-03,1e-04,1e-05))
s29

multiplot(s01, s10, s03, s04, s29, s12, cols=2)


############################################
### Figure S3 (variance) ###################
############################################
exp<-read.table("tech_biol_total2.txt", h=T)
s1<-subset(exp, allele_measured=="SCR01")
s2<-subset(exp, allele_measured=="SCR02")
s3<-subset(exp, allele_measured=="SCR03")
s4<-subset(exp, allele_measured=="SCR04")
s10<-subset(exp, allele_measured=="SCR10")
s12<-subset(exp, allele_measured=="SCR12")
s13<-subset(exp, allele_measured=="SCR13")
s20<-subset(exp, allele_measured=="SCR20")
s29<-subset(exp, allele_measured=="SCR29")


p1r<- ggplot(s1, aes(x=Other_allele_inGenot, y=Ct2, fill=as.factor(recc_tech)))+
  geom_boxplot(position=position_dodge(1))+
  scale_fill_brewer(type="div",palette="Greys")+
  guides(fill=FALSE)+
  xlab("")+
  ylab("ct SCR/Actin")+
  labs(title=("SCR01"))+
  facet_grid( . ~ Other_allele_inGenot, scales="free_x",switch="x")+
  geom_text(aes(label=s1$info, y=0.05),position= position_dodge(1),size=2)+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        strip.text = element_text(size=8),
        axis.text.x = element_blank(),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))
p1r

p2r<- ggplot(s2, aes(x=Other_allele_inGenot, y=Ct2, fill=as.factor(recc_tech)))+
  geom_boxplot(position=position_dodge(1))+
  scale_fill_brewer(type="div",palette="Greys")+
  guides(fill=FALSE)+
  xlab("")+
  ylab("ct SCR/Actin")+
  labs(title=("SCR02"))+
  facet_grid( . ~ Other_allele_inGenot, scales="free_x",switch="x")+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        strip.text = element_text(size=8),
        axis.text.x = element_blank(),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))
p2r

p3r<- ggplot(s3, aes(x=Other_allele_inGenot, y=Ct2, fill=as.factor(recc_tech)))+
  geom_boxplot(position=position_dodge(1))+
  scale_fill_brewer(type="div",palette="Greys")+
  guides(fill=FALSE)+
  xlab("")+
  ylab("ct SCR/Actin")+
  labs(title=("SCR03"))+
  facet_grid( . ~ Other_allele_inGenot, scales="free_x",switch="x")+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        strip.text = element_text(size=8),
        axis.text.x = element_blank(),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))
p3r

p4r<- ggplot(s4, aes(x=Other_allele_inGenot, y=Ct2, fill=as.factor(recc_tech)))+
  geom_boxplot(position=position_dodge(1))+
  scale_fill_brewer(type="div",palette="Greys")+
  guides(fill=FALSE)+
  xlab("")+
  ylab("ct SCR/Actin")+
  labs(title=("SCR04"))+
  facet_grid( . ~ Other_allele_inGenot, scales="free_x",switch="x")+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        strip.text = element_text(size=8),
        axis.text.x = element_blank(),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))
p4r

p5r<- ggplot(s10, aes(x=Other_allele_inGenot, y=Ct2, fill=as.factor(recc_tech)))+
  geom_boxplot(position=position_dodge(1))+
  scale_fill_brewer(type="div",palette="Greys")+
  guides(fill=FALSE)+
  xlab("")+
  ylab("ct SCR/Actin")+
  labs(title=("SCR10"))+
  facet_grid( . ~ Other_allele_inGenot, scales="free_x",switch="x")+
  facet_grid( . ~ Other_allele_inGenot, scales="free_x",switch="x")+
  geom_text(aes(label=s10$info, y=0.15),position= position_dodge(1),size=2)+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        strip.text = element_text(size=8),
        axis.text.x = element_blank(),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))
p5r

p6r<- ggplot(s12, aes(x=Other_allele_inGenot, y=Ct2, fill=as.factor(recc_tech)))+
  geom_boxplot(position=position_dodge(1))+
  scale_fill_brewer(type="div",palette="Greys")+
  guides(fill=FALSE)+
  xlab("")+
  ylab("ct SCR/Actin")+
  labs(title=("SCR12"))+
  facet_grid( . ~ Other_allele_inGenot, scales="free_x",switch="x")+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        strip.text = element_text(size=8),
        axis.text.x = element_blank(),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))
p6r

p7r<- ggplot(s13, aes(x=Other_allele_inGenot, y=Ct2, fill=as.factor(recc_tech)))+
  geom_boxplot(position=position_dodge(1))+
  scale_fill_brewer(type="div",palette="Greys")+
  guides(fill=FALSE)+
  xlab("")+
  ylab("ct SCR/Actin")+
  labs(title=("SCR13"))+
  facet_grid( . ~ Other_allele_inGenot, scales="free_x",switch="x")+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        strip.text = element_text(size=8),
        axis.text.x = element_blank(),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))
p7r

p8r<- ggplot(s20, aes(x=Other_allele_inGenot, y=Ct2, fill=as.factor(recc_tech)))+
  geom_boxplot(position=position_dodge(1))+
  scale_fill_brewer(type="div",palette="Greys")+
  guides(fill=FALSE)+
  xlab("")+
  ylab("ct SCR/Actin")+
  labs(title=("SCR20"))+
  facet_grid( . ~ Other_allele_inGenot, scales="free_x",switch="x")+
  geom_text(aes(label=s20$info, y=1.6),position= position_dodge(1),size=2)+
  theme(panel.background = element_rect(fill = "white"),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        strip.text = element_text(size=8),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))
p8r

p9r<- ggplot(s29, aes(x=Other_allele_inGenot, y=Ct2, fill=as.factor(recc_tech)))+
  geom_boxplot(position=position_dodge(1))+
  scale_fill_brewer(type="div",palette="Greys")+
  guides(fill=FALSE)+
  xlab("")+
  ylab("ct SCR/Actin")+
  labs(title=("SCR29"))+
  facet_grid( . ~ Other_allele_inGenot, scales="free_x",switch="x")+
  theme(panel.background = element_rect(fill = "white"),
        axis.ticks.x=element_blank(), legend.background= element_rect(fill="white"),
        legend.key=element_rect(fill="white"),
        strip.text = element_text(size=8),
        axis.text.x = element_blank(),plot.title = italic.labs,
        axis.line.y=element_line(colour = "Black", size = 1))
p9r


multiplot(p1r, p2r, p6r, p3r, p5r, p7r, p9r, p4r, p8r, cols=3)


############################################
### Figure S4 (normalité) ##################
############################################
s=read.table("modele_9_10_17_scr.txt", h=T)
mfamtest3<-lmer(log(s1$Ct_SCR.actine) ~ (1|allele_measured:stade)+stade*dom_phenotype+
                  (1|replicatBiol_genotype/replicat_Techclone) , 
                data =s, na.action=na.omit)
E7<- resid(mfamtest3)
hist(E7, xlab ="residuals", main="", breaks=150)
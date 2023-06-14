library(vegan)
library(ggplot2)
library(readr)
library(emmeans)
NGSdata<-read_csv("OTU Frequencies and Taxa.csv")
NGSdata<-NGSdata[1:80, 1:62]
Metadata<-read_csv("Sequencing Metadata.csv")
Taxonomy<-read_csv("OTU Taxonomy.csv")
library(tidyverse)
library(ggplot2)
mine.long <- pivot_longer(data = NGSdata, 
                          cols = -Abbreviation,
                          names_to = "OTU",
                          values_to = "Abundance")

rel_abund<-inner_join(Metadata, mine.long, by="Abbreviation") %>% inner_join(., Taxonomy, by="OTU") %>% group_by(Abbreviation) %>% 
  mutate(relative_abundance = Abundance) %>% ungroup()
library(qualpalr)
pal = qualpal(61, colorspace=list(h=c(0,360), s=c(0.3,1), l=c(0.2,0.8)))

rel_abund2 <-rel_abund %>% group_by(State, Invasion, Abbreviation, Taxonomy) %>% summarize(relative_abundance=sum(relative_abundance), .groups="drop") %>%
  group_by(State, Invasion, Taxonomy) %>% 
    summarize(mean_rel_abund = 100*mean(relative_abundance), .groups="drop")

ggplot(aes(x=Invasion, y=mean_rel_abund, fill=fct_reorder(Taxonomy, mean_rel_abund)),data=rel_abund2) + geom_bar(stat="identity", position = position_fill(reverse = TRUE)) + facet_wrap(~State) + labs(x=NULL, y="Mean Relative Abundance (%)") +theme_classic() +theme(legend.position="right") + scale_fill_manual(values=pal$hex) + theme(legend.title = element_blank()) + guides(fill=guide_legend(ncol=1, reverse=TRUE)) + theme(legend.key.width=unit(2, 'mm')) + theme(legend.key.height=unit(1, 'mm'))+ theme(legend.text = element_text(size=8))


rel_abund3<-rel_abund %>% group_by(Abbreviation, Family) %>% summarize(relative_abundance=sum(relative_abundance), .groups="drop") %>%
  group_by(Family) %>% 
  summarize(mean_rel_abund = 100*mean(relative_abundance), .groups="drop")

# Relative abundances of families by State
rel_abund4<-rel_abund %>% select(Abbreviation, State, Family, relative_abundance) %>% group_by(State, Family)

CG<-rel_abund4 %>% filter(Family=="Clarideoglomeraceae")
kruskal.test(relative_abundance~State, data=CG)

GL<-rel_abund4 %>% filter(Family=="Glomeraceae")
kruskal.test(relative_abundance~State, data=GL)

AR<-rel_abund4 %>% filter(Family=="Archaeosporaceae")
kruskal.test(relative_abundance~State, data=AR)

AM<-rel_abund4 %>% filter(Family=="Ambisporaceae")
kruskal.test(relative_abundance~State, data=AM)

DI<-rel_abund4 %>% filter(Family=="Diversisporaceae")
kruskal.test(relative_abundance~State, data=DI)

PA<-rel_abund4 %>% filter(Family=="Paraglomeraceae")
kruskal.test(relative_abundance~State, data=PA)

# Relative Abundance of Families by Invasion
rel_abund5<-rel_abund %>% select(Abbreviation, Invasion, Family, relative_abundance) %>% group_by(Invasion, Family)
CG1<-rel_abund5 %>% filter(Family=="Clarideoglomeraceae")
kruskal.test(relative_abundance~Invasion, data=CG1)

GL1<-rel_abund5 %>% filter(Family=="Glomeraceae")
kruskal.test(relative_abundance~Invasion, data=GL1)
# Sig.

AR1<-rel_abund5 %>% filter(Family=="Archaeosporaceae")
kruskal.test(relative_abundance~Invasion, data=AR1)

AM1<-rel_abund5 %>% filter(Family=="Ambisporaceae")
kruskal.test(relative_abundance~Invasion, data=AM1)

DI1<-rel_abund5 %>% filter(Family=="Diversisporaceae")
kruskal.test(relative_abundance~Invasion, data=DI1)

PA1<-rel_abund5 %>% filter(Family=="Paraglomeraceae")
kruskal.test(relative_abundance~Invasion, data=PA1)

# OTU richness by State and Invasion

rel_abund6 <- rel_abund %>% group_by(State, Taxonomy) %>% summarize(relative_abundance=sum(relative_abundance), .groups="drop") %>% 
  group_by(State)
rel_abund6$presence<- ifelse(rel_abund6$relative_abundance>0, 1,0)
rel_abund6 %>% group_by(State) %>% summarize(richness=sum(presence), .groups="drop") %>%
  group_by(State)
  
rel_abund7 <-rel_abund %>% group_by(State, Invasion, Abbreviation, Taxonomy) %>% summarize(relative_abundance=sum(relative_abundance), .groups="drop") %>%
  group_by(State, Invasion, Taxonomy) %>% 
  summarize(mean_rel_abund = 100*mean(relative_abundance), .groups="drop")
rel_abund7$presence<- ifelse(rel_abund7$mean_rel_abund>0, 1,0)

rel_abund8 <-rel_abund7 %>% group_by(State, Invasion) %>% summarize(richness=sum(presence), .groups="drop") %>%
  group_by(State, Invasion)

glm1<-glm(rel_abund7$presence~ rel_abund7$State*rel_abund7$Invasion, family = "binomial")
library(car)
Anova(glm1, type=3)
glm1pairs<-emmeans(glm1, pairwise~State*Invasion)
glm1pairs
confint(glm1pairs)

# Family richness by State and Invasion
rel_abund9 <-rel_abund %>% group_by(State, Invasion, Abbreviation, Family) %>% summarize(relative_abundance=sum(relative_abundance), .groups="drop") %>%
  group_by(State, Invasion, Family) %>% 
  summarize(mean_rel_abund = 100*mean(relative_abundance), .groups="drop")
rel_abund8$presence<- ifelse(rel_abund8$mean_rel_abund>0, 1,0)

rel_abund10 <- rel_abund %>%
  group_by(State, Family) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(relative_abundance),
    sd=sd(relative_abundance)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

rel_abund11 <- rel_abund %>%
  group_by(Invasion, Family) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(relative_abundance),
    sd=sd(relative_abundance)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

library (RColorBrewer)
library(ggpubr)
myPalette<-c(IL="grey", CO="White")
myPalette2<-c(Invaded="grey", Uninvaded="White")
ggplot(rel_abund10, aes(x=Family, y=mean, fill=State)) + geom_bar(stat="identity", color="black",position=position_dodge2(preserve="single", padding=0)) + geom_errorbar(aes(ymin=mean, ymax=mean+se), width=.2,position=position_dodge(.9)) + theme_bw(base_size = 15) + theme(panel.grid = element_blank(), legend.position = c(0.85,0.9)) + labs(x= "Family", y ="Mean Relative Abundance", pattern = "Site") + scale_fill_manual(values=myPalette) + theme(axis.title.x = element_text(size=10, colour="black"), axis.title.y = element_text(size= 10, colour = "black"), legend.title = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(), legend.text = element_text(size =10, colour = "black"), legend.key.size=unit(4,"mm"), legend.background = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_signif(y_position=c(0.03, 0.08), xmin=c(1.75, 3.75), xmax=c(2.25, 4.25), annotation=c("***", "*"), fill=NULL, tip_length=0.01)
ggplot(rel_abund11, aes(x=Family, y=mean, fill=Invasion)) + geom_bar(stat="identity", color="black",position=position_dodge2(preserve="single", padding=0)) + geom_errorbar(aes(ymin=mean, ymax=mean+se), width=.2,position=position_dodge(.9)) + theme_bw(base_size = 15) + theme(panel.grid = element_blank(), legend.position = c(0.85,0.9)) + labs(x= "Family", y ="Mean Relative Abundance", pattern = "Site") + scale_fill_manual(values=myPalette2) + theme(axis.title.x = element_text(size=10, colour="black"), axis.title.y = element_text(size= 10, colour = "black"), legend.title = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank(), legend.text = element_text(size =10, colour = "black"), legend.key.size=unit(4,"mm"), legend.background = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_signif(y_position=0.03, xmin=4.75, xmax=5.25, annotation="*", fill=NULL, tip_length=0.01)

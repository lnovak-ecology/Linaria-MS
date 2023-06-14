biomass <- read.csv("~/Desktop/MS Thesis/Biomass Data.csv")
names(biomass)
library(tidyverse)
library(ggplot2)
library(bbmle)
names(biomass)

biomass$Total.Root.Mass.1<- as.numeric(as.character(biomass$Total.Root.Mass.1))
biomass$Total.Root.Mass.2<- as.numeric(as.character(biomass$Total.Root.Mass.2))
biomass$`Shoot.Mass.1`[is.na(biomass$`Shoot.Mass.1`)]<-0
biomass$`Total.Root.Mass.1`[is.na(biomass$`Total.Root.Mass.1`)]<-0

biomass$Grown[biomass$Grown=="Native and Linaria"]<-NA
names(biomass)
biomass2<-biomass[complete.cases(biomass[, 2]),]

biomass2 <- biomass2 %>% 
  ungroup() %>% 
  dplyr::mutate(ID = row_number()) %>% 
  group_by(ID) %>% 
  dplyr::mutate(AllBiomass = (sum(`Shoot.Mass.1`, `Total.Root.Mass.1`, na.rm = T)))


biomassCO<-biomass2 %>% filter(State == "CO")
biomassIL<-biomass2 %>% filter(State == "IL")
biomassCO$presence <- ifelse(biomassCO$AllBiomass>0, 1,0)
biomassIL$presence <- ifelse(biomassIL$AllBiomass>0, 1,0)

# CO Presence
library(tidyverse)
library(car)
library(lme4)

# no t value because glm, report just P for post-hoc
glm7<-glmer(biomassCO$presence~biomassCO$Invasion*biomassCO$Grown + (1|biomassCO$Site), family = "binomial")
plot(glm7)
summary(glm7)
Anova(glm7, type=2)
# Interaction very sig.
confint(glm7)
library(emmeans)
glm7pairs<-emmeans(glm7, specs=pairwise~Invasion*Grown)
glm7pairs
glm7pairs$contrasts %>% rbind(adjust="sidak")
confint(glm7pairs)
aggregate(presence~Invasion*Grown, data=biomassCO, FUN=sum)

# IL Presence - not sig.
glm25<-glmer(biomassIL$presence~biomassIL$Invasion*biomassIL$Grown + (1|biomassIL$Site), family = "binomial")
library(emmeans)
Anova(glm25, type=2)
glm25pairs<-emmeans(glm25, pairwise~Invasion:Grown)
glm25pairs
glm25pairs$contrasts %>% rbind(adjust="sidak")
confint(glm25pairs)

library(ggplot2)

# Calculate mean and SE
library(doBy)
ILdata <- summaryBy(presence ~ Grown + Invasion + Treatment, data=biomassIL, FUN=c(length,mean,sd))
ILdata
ILdata$presence.mean<-ILdata$presence.mean*100
names(ILdata)[names(ILdata)=="presence.length"] <- "N"
ILdata$presence.se <- ILdata$presence.sd / sqrt(ILdata$N)
ILdata$presence.se<-ILdata$presence.se*100

COdata <- summaryBy(presence ~ Grown + Invasion + Treatment, data=biomassCO, FUN=c(length,mean,sd))
COdata$presence.mean<-COdata$presence.mean*100
names(COdata)[names(COdata)=="presence.length"] <- "N"
COdata$presence.se <- COdata$presence.sd / sqrt(COdata$N)
COdata$presence.se<-COdata$presence.se*100

pd <- position_dodge(1)
library (RColorBrewer)
cbPalette <- c(brewer.pal(4, "YlGnBu"))
ggplot(ILdata, aes(x = Invasion, y = presence.mean, fill = Treatment)) + geom_bar(stat = "identity", position = "dodge", color="black") + geom_errorbar(aes(ymin = presence.mean, ymax = presence.mean + presence.se), width = 0.5, position = pd) + theme(panel.grid = element_blank(),legend.position = "bottom") + xlab("Soil Type") + ylab("Percent of Pots with Plant Growth") + facet_wrap(~ILdata$Grown, ncol=1, scales = "free_y") +
scale_fill_manual(values= cbPalette)+theme(axis.text.x = element_text(size=10, colour= "black", angle = 0),axis.text.y = element_text(size=10, colour="black"),strip.text = element_text(size = 10))+theme(axis.title.x = element_text(size=12, colour="black"), axis.title.y = element_text(size= 12, colour = "black"))+ theme(legend.title = element_blank())+theme(axis.ticks.x = element_blank()) + theme(legend.title = element_blank())+theme(legend.text = element_text(size = 10, colour = "black"))+theme(legend.background = element_blank())+ylim(0,100)

       
ggplot(COdata, aes(x = Invasion, y = presence.mean, fill = Treatment)) + geom_bar(stat = "identity", position = "dodge", color="black") + geom_errorbar(aes(ymin = presence.mean, ymax = presence.mean + presence.se), width = 0.5, position = pd) + theme(panel.grid = element_blank(),legend.position = "bottom") + xlab("Soil Type") + ylab("Percent of Pots with Plant Growth") + facet_wrap(~COdata$Grown, ncol=1, scales = "free_y") +
  scale_fill_manual(values= cbPalette)+theme(axis.text.x = element_text(size=10, colour= "black", angle = 0),axis.text.y = element_text(size=10, colour="black"),strip.text = element_text(size = 10))+theme(axis.title.x = element_text(size=12, colour="black"), axis.title.y = element_text(size= 12, colour = "black"))+ theme(legend.title = element_blank())+theme(axis.ticks.x = element_blank())+ theme(legend.title = element_blank())+theme(legend.text = element_text(size = 10, colour = "black"))+theme(legend.background = element_blank()) + geom_signif(data = data.frame(Grown = "Linaria"), aes(y_position=101, xmin=1.0, xmax=2.0, annotations="**", fill=NULL), tip_length=0.01, manual=T) + facet_grid(~ Grown, scales = "free")



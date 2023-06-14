library(readr)
NGSdata<-read_csv("OTU Frequencies and Taxa.csv")
NGSdata<-NGSdata[1:80, 1:62]
Metadata<-read_csv("Sequencing Metadata.csv")
Taxonomy<-read_csv("OTU Taxonomy.csv")
Taxonomy1<-as.data.frame(Taxonomy)
rownames(Taxonomy1)<-Taxonomy1[,1]
Taxonomy1[,1] <- NULL
SoilNGS<-read_csv("Nutrient Analyses For Sequences.csv")
names(SoilNGS)
SoilNGS$pH <-as.numeric(as.character(SoilNGS$pH))
SoilNGS$`Ammonium ug/g` <-as.numeric(as.character(SoilNGS$`Ammonium ug/g`))
SoilNGS$`Nitrate u/g` <-as.numeric(as.character(SoilNGS$`Nitrate u/g`))
SoilNGS$`Phosphate u/g` <-as.numeric(as.character(SoilNGS$`Phosphate u/g`))
SoilNGS<-SoilNGS[-c(73,74),]
metadata<-SoilNGS %>% dplyr::select(Abbreviation, State, Invasion, Site, `Ammonium ug/g`, `Nitrate u/g`, `Phosphate u/g`, pH, Sand, Clay, Silt)
metadata1<-SoilNGS %>% group_by(State, Invasion) %>% summarize(`Ammonium ug/g`=mean(`Ammonium ug/g`), `Nitrate u/g`=mean(`Nitrate u/g`), `Phosphate u/g`=mean(`Phosphate u/g`),  pH=mean(pH), Sand=mean(Sand), Clay=mean(Clay),  Silt=mean(Silt), .groups = "drop")
metadata2<-metadata %>% dplyr:: select(Abbreviation, State, Invasion)
library(tidyverse)
library(ggplot2)
mine.long <- pivot_longer(data = NGSdata, 
                          cols = -Abbreviation,
                          names_to = "OTU",
                          values_to = "Abundance")
rel_abund<-inner_join(Metadata, mine.long, by="Abbreviation") %>% inner_join(., Taxonomy, by="OTU") %>% group_by(Abbreviation) %>% 
  mutate(relative_abundance = Abundance) %>% ungroup()

rel_abund$presence<- ifelse(rel_abund$relative_abundance>0, 1,0)

# OTU Richness 
OTUsum<-rel_abund %>% dplyr::group_by(Abbreviation, State, Invasion, Site) %>% dplyr::summarize(sum=sum(presence)) %>% dplyr::ungroup()
OTUsum<-OTUsum[-c(67:68),]
OTUSoil<-inner_join(metadata, OTUsum, by="Abbreviation")

# Shows that 1 predictor is best for model
library(car)
library(lme4)
aov1<-lmer(sum~State*Invasion + (1|Site), data=OTUsum)
summary(aov1)
Anova(aov1, type=2)
# Not sig.

# Not Sig.
lm1<-lmer(OTUSoil$sum~OTUSoil$`Ammonium ug/g`+ (1|OTUSoil$Site))
summary(lm1)
Anova(lm1)
confint(lm1)
plot(lm1)

# Not Sig.
lm2<-lmer(OTUSoil$sum~OTUSoil$`Nitrate u/g` + (1|OTUSoil$Site))
summary(lm2)
Anova(lm2)
confint(lm2)
plot(lm2)

# Not sig.
lm3<-lmer(OTUSoil$sum~OTUSoil$`Phosphate u/g`+ (1|OTUSoil$Site))
Anova(lm3)

# Not sig.
lm4<-lmer(OTUSoil$sum~OTUSoil$pH + (1|OTUSoil$Site))
Anova(lm4)
confint(lm4)

# Not sig.
lm5<-lmer(OTUSoil$sum~OTUSoil$Sand+ (1|OTUSoil$Site))
Anova(lm5, type=2)
confint(lm5)

# Clay is sig.
lm6<-lmer(OTUSoil$sum~OTUSoil$Clay + (1|OTUSoil$Site))
Anova(lm6, type=2)
confint(lm6)

# Silt is not sig.
lm7<-lmer(OTUSoil$sum~OTUSoil$Silt + (1|OTUSoil$Site))
Anova(lm7, type=2)
confint(lm7)

library(bbmle)
AICctab(lm1, lm2, lm3, lm4, lm5, lm6, lm7)

# Family richness
rel_abundFamily<-rel_abund %>% dplyr::group_by(Abbreviation, State, Invasion, Site, Family) %>% dplyr::summarize(relative_abundance=sum(relative_abundance), .groups="drop")
rel_abundFamily$presence<- ifelse(rel_abundFamily$relative_abundance>0, 1,0)
Familysum<-rel_abundFamily %>% dplyr::group_by(Abbreviation, State, Invasion, Site) %>% dplyr::summarize(sum=sum(presence)) %>% dplyr::ungroup()
FamilySoil<-inner_join(Familysum, metadata, by="Abbreviation")
aov2<-lmer(sum~State*Invasion + (1|Site), data=Familysum)
Anova(aov2, type=2)
summary(aov2)
Familypairs<-emmeans(aov2, pairwise~State:Invasion) # Nothing Significant

# Not sig.
lm10<-lmer(FamilySoil$sum~FamilySoil$`Ammonium ug/g`+ (1|FamilySoil$Site))
Anova(lm10, type=2)
confint(lm10)

# Not sig.
lm11<-lmer(FamilySoil$sum~FamilySoil$`Nitrate u/g` + (1|FamilySoil$Site))
Anova(lm11, type=2)
confint(lm11)

# Not sig.
lm12<-lmer(FamilySoil$sum~FamilySoil$`Phosphate u/g`+ (1|FamilySoil$Site))
Anova(lm12)
confint(lm12)

# Not sig.
lm13<-lmer(FamilySoil$sum~FamilySoil$pH + (1|FamilySoil$Site))
Anova(lm13, type=2)
confint(lm13)

# Not sig.
lm14<-lmer(FamilySoil$sum~FamilySoil$Sand + (1|FamilySoil$Site))
Anova(lm14, type=2)
confint(lm14)

# Not sig.
lm15<-lmer(FamilySoil$sum~FamilySoil$Clay + (1|FamilySoil$Site))
Anova(lm15, type=2)
confint(lm15)

lm16<-lmer(FamilySoil$sum~FamilySoil$Silt + (1|FamilySoil$Site))
Anova(lm16, type=2)
confint(lm16)

AICctab(lm10, lm11, lm12, lm13, lm14, lm15, lm16)

# Family Relative Abundances by Site
familyabund<-rel_abund %>% dplyr::group_by(Abbreviation, State, Invasion, Site, Family) %>% dplyr::summarize(mean_abund=mean(relative_abundance), .groups="drop")
ambisporaceae<-familyabund %>% filter(Family=="Ambisporaceae")
ambisporaceae %>% group_by(State) %>% summarise(mean=mean(mean_abund))
archaeosporaceae<-familyabund %>% filter(Family=="Archaeosporaceae")
archaeosporaceae %>% group_by(State) %>% summarise(mean=mean(mean_abund))
clarideoglomeraceae<-familyabund %>% filter(Family=="Clarideoglomeraceae")
clarideoglomeraceae %>% group_by(State) %>% summarise(mean=mean(mean_abund))
diversisporaceae<-familyabund %>% filter(Family=="Diversisporaceae")
diversisporaceae %>% group_by(State) %>% summarise(mean=mean(mean_abund))
glomeraceae<-familyabund %>% filter(Family=="Glomeraceae")
glomeraceae %>% group_by(State) %>% summarise(mean=mean(mean_abund))
paraglomeraceae<-familyabund %>% filter(Family=="Paraglomeraceae")



hist(ambisporaceae$mean_abund)
qqPlot(ambisporaceae$mean_abund)
kruskal.test(ambisporaceae$mean_abund~ambisporaceae$State)
kruskal.test(ambisporaceae$mean_abund~ambisporaceae$Invasion)

hist(archaeosporaceae$mean_abund)
kruskal.test(archaeosporaceae$mean_abund~archaeosporaceae$State)
kruskal.test(archaeosporaceae$mean_abund~archaeosporaceae$Invasion)

hist(clarideoglomeraceae$mean_abund)
kruskal.test(clarideoglomeraceae$mean_abund~clarideoglomeraceae$State)
kruskal.test(clarideoglomeraceae$mean_abund~clarideoglomeraceae$Invasion)

hist(diversisporaceae$mean_abund)
kruskal.test(diversisporaceae$mean_abund~diversisporaceae$State)
kruskal.test(diversisporaceae$mean_abund~diversisporaceae$Invasion)

hist(glomeraceae$mean_abund)
kruskal.test(glomeraceae$mean_abund~glomeraceae$State)
kruskal.test(glomeraceae$mean_abund~glomeraceae$Invasion)

hist(paraglomeraceae$mean_abund)
kruskal.test(paraglomeraceae$mean_abund~paraglomeraceae$State)
kruskal.test(paraglomeraceae$mean_abund~paraglomeraceae$Invasion)

# Number of reads, not relative abundances
library(vegan)
absolutedata<-read_csv("OTU Absolute Abundances.csv")
colnames(absolutedata)[1]<-"Abbreviation"
absolutedata1<-absolutedata[,2:62]
absolutedata1<-absolutedata[rowSums(is.na(absolutedata)) != ncol(absolutedata),]
absolutedata1<-as.data.frame(absolutedata1)
rownames(absolutedata1)<-absolutedata$Abbreviation
absolutedata1<-absolutedata1[,2:62]
absolutedata1<-absolutedata1[-c(67:68),]

# Rarefaction Curve
S <- specnumber(absolutedata1) # observed number of species
(raremax <- min(rowSums(absolutedata1)))
Srare <- rarefy(absolutedata1, sample=raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(absolutedata1, step = 10, sample = raremax, col = "blue", label=FALSE)
set.seed(12345)

# Shannon and Simpson Diversity Indexes
OTUsum$Simpson<-diversity(absolutedata1, index="simpson")
OTUsum$Shannon<-diversity(absolutedata1)
hist(OTUsum$Simpson)
qqPlot(OTUsum$Simpson)
hist(OTUsum$Shannon)
qqPlot(OTUsum$Shannon)

# Chao1 Richness
richness<-estimateR(absolutedata1)
richness1<-data.frame(t(richness))
OTUsum$Chao1<-richness1$S.chao1
OTUsum$ChaoSE<-richness1$se.chao1
hist(OTUsum$Chao1)
qqPlot(OTUsum$Chao1)
lm3<-lmer(OTUsum$Chao1~OTUsum$State*OTUsum$Invasion + (1|OTUsum$Site))
summary(lm3)
Anova(lm3)

# Richness
OTUsum$Richness<-specnumber(absolutedata1)
hist(OTUsum$Richness)
qqPlot(OTUsum$Richness)
lm19<-lmer(OTUsum$Richness~OTUsum$State*OTUsum$Invasion + (1|OTUsum$Site))
summary(lm19)
Anova(lm19)

# Pielou's J
OTUsum$J<-OTUsum$Shannon/log(specnumber(absolutedata1))
hist(OTUsum$J)
qqPlot(OTUsum$J)
lm20<-lmer(OTUsum$J~OTUsum$State*OTUsum$Invasion + (1|OTUsum$Site))
summary(lm20)
aov20<-aov(lm20)
Anova(lm20)
library(emmeans)
lm20pair<-emmeans(lm20, specs=pairwise~State*Invasion)
lm20pair
lm20pair$contrasts %>% rbind(adjust="sidak")




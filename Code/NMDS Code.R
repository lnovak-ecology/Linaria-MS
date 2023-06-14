library(vegan)
library(ggplot2)
library(dplyr)
library(readr)
library(plyr)
set.seed(123456)
NGSdata<-read_csv("OTU Relative Frequencies.csv")
Metadata<-read_csv("Sequencing Metadata.csv")
names(Metadata)
Metadata<-dplyr::select(Metadata, -c(Site))
Soildata<-read_csv("Sequence Soil.csv")
Soildata$`Ammonium ug/g`<-as.numeric(Soildata$`Ammonium ug/g`)
Soildata$`Phosphate u/g`<-as.numeric(Soildata$`Phosphate u/g`)
Soildata$`Nitrate u/g`<-as.numeric(Soildata$`Nitrate u/g`)
Soildata<- na.omit(Soildata)
Soildata1<-dplyr::select(Soildata, -c(Abbreviation, Invasion, State, Site, Sample, `Mass (g)`))
Metadata<-Metadata[-c(67,68),]
Metadata$Abbreviation<-as.factor(Metadata$Abbreviation)
Metadata$State<-as.factor(Metadata$State)
Metadata$Invasion<-as.factor(Metadata$Invasion)

# Both states
NGSdata1<-NGSdata[-c(67,68),]
NGSdata2<-NGSdata1[,2:62]
NGSdata2<-as.matrix(NGSdata2)
grp<-Metadata$State
grp
Clay<-Soildata1$Clay
Silt<-Soildata1$Silt
Sand<-Soildata1$Sand
sol <- metaMDS(NGSdata2)
en = envfit(sol, Soildata1, permutations = 999, na.rm = TRUE)
en_coord_cat = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
plot(sol)
plot(en)
NMDS = data.frame(MDS1 = sol$points[,1], MDS2 = sol$points[,2])
data.scores <- as.data.frame(scores(NMDS))
data.scores$grp <- grp
data.scores$Clay <- Clay
data.scores$Silt <- Silt
data.scores$Sand <- Sand
head(data.scores)

species.scores <- as.data.frame(scores(sol, "species"))
species.scores$species <- rownames(species.scores)
head(species.scores)

ds <- data.scores
find_hull2 <- function(ds) ds[chull(ds$MDS1, ds$MDS2), ]

hulls2 <- ddply(ds, "grp", find_hull2)
hulls2

library (RColorBrewer)
library(ggrepel)
myPalette <- c(CO="gray", IL="black")
ggplot() + geom_point(data=ds,aes(x=MDS1,y=MDS2,colour=grp, shape=grp),size=2) +
  geom_polygon(data=hulls2,aes(x=MDS1,y=MDS2,fill=grp,group=grp), color = "black", alpha=0.30, size=0.1) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.9,0.8),
        legend.box.background = element_rect()) + scale_color_manual(values=myPalette) +
        geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cat, size =1, alpha = 0.5, colour = "grey30") +
        geom_text_repel(data = en_coord_cat, aes(x = NMDS1, y = NMDS2), 
            label = row.names(en_coord_cat), colour = "black")

# Colorado Only
NGSdataCO<-NGSdata2[1:40,]
MetadataCO<-Metadata[Metadata$State=="CO",]
SoildataCO<-Soildata[Soildata$State=="CO",]
SoildataCO2<-SoildataCO %>% select(Abbreviation, Invasion, pH, `Ammonium ug/g`, `Nitrate u/g`, `Phosphate u/g`, Clay, Silt, Sand)
grp1<-MetadataCO$Invasion
grp1

sol1<-metaMDS(NGSdataCO)
plot(sol1)
NMDS1 = data.frame(NMDS1 = sol1$points[,1], NMDS2 = sol1$points[,2])
data.scores1 <- as.data.frame(scores(sol1, display="sites"))
data.scores1$site <- rownames(data.scores1)
data.scores1$grp <- grp1
head(data.scores1) 
species.scores1 <- as.data.frame(scores(sol1, "species"))
species.scores1$species <- rownames(species.scores1)
head(species.scores1)
en1 = envfit(sol1, SoildataCO2, permutations = 999, na.rm = TRUE)
en_coord_cat1 = as.data.frame(scores(en1, "vectors")) * ordiArrowMul(en1)

df1 <- data.scores1
find_hull1 <- function(df1) df1[chull(df1$NMDS1, df1$NMDS2), ]
hulls1 <- ddply(df1, "grp", find_hull1)
hulls1
geom_polygon(data=hulls1,aes(x=NMDS1,y=NMDS2,fill=grp,group=grp),alpha=0.30)

cbPalette <- c(Uninvaded = "green", Invaded = "yellow")
ggplot() + geom_point(data=data.scores1,aes(x=NMDS1,y=NMDS2, shape=grp),color="black",size=2) +
  geom_polygon(data=hulls1,aes(x=NMDS1,y=NMDS2,fill=grp,group=grp), color = "black", alpha=0.30, size=0.1) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_blank(),  # remove x-axis text
                   axis.text.y = element_blank(), # remove y-axis text
                   axis.ticks = element_blank(),  # remove axis ticks
                   axis.title.x = element_text(size=18), # remove x-axis labels
                   axis.title.y = element_text(size=18), # remove y-axis labels
                   panel.background = element_blank(), 
                   panel.grid.major = element_blank(),  #remove major-grid labels
                   panel.grid.minor = element_blank(),  #remove minor-grid labels
                   plot.background = element_blank(),
                  legend.title = element_blank(),
                  legend.position = c(0.2,0.8),
                  legend.box.background = element_rect()) + 
                  scale_color_manual(values=cbPalette) + scale_fill_manual(values = cbPalette) +
                  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
                  data = en_coord_cat1, size =1, alpha = 0.5, colour = "black") +
                  geom_text_repel(data = en_coord_cat1, aes(x = NMDS1, y = NMDS2), 
                  label = row.names(en_coord_cat1), colour = "black")


# Illinois Only
NGSdataIL<-NGSdata2[41:78,]
MetadataIL<-Metadata[Metadata$State=="IL",]
SoildataIL<-Soildata[Soildata$State=="IL",]
names(SoildataIL)
SoildataIL2<-SoildataIL %>% select(Abbreviation, State, Invasion, pH, `Ammonium ug/g`, `Nitrate u/g`, `Phosphate u/g`, Clay, Silt, Sand)

grp2<-MetadataIL$Invasion
grp2
sol2<-metaMDS(NGSdataIL)
NMDS2 = data.frame(MDS1 = sol2$points[,1], MDS2 = sol2$points[,2])
data.scores2 <- as.data.frame(scores(sol2, display="sites"))
data.scores2$site <- rownames(data.scores2)
data.scores2$grp <- grp2
head(data.scores2)
species.scores2 <- as.data.frame(scores(sol2, "species"))
species.scores2$species <- rownames(species.scores2)
head(species.scores2)
en2 = envfit(sol2, SoildataIL2, permutations = 999, na.rm = TRUE)
en_coord_cat2 = as.data.frame(scores(en2, "vectors")) * ordiArrowMul(en2)

df2 <- data.scores2
find_hull2 <- function(df2) df2[chull(df2$NMDS1, df2$NMDS2), ]
hulls2 <- ddply(df2, "grp", find_hull2)
hulls2
geom_polygon(data=hulls2,aes(x=NMDS1,y=NMDS2,fill=grp2,group=grp2),alpha=0.30)

ggplot() + geom_point(data=data.scores2,aes(x=NMDS1,y=NMDS2, shape=grp), colour="black", size=2) +
  geom_polygon(data=hulls2,aes(x=NMDS1,y=NMDS2,fill=grp,group=grp),color="black", alpha=0.30, size=0.1) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=18), # remove x-axis labels
        axis.title.y = element_text(size=18), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.2,0.2),
        legend.box.background = element_rect()) + scale_color_manual(values=cbPalette) + scale_fill_manual(values = cbPalette) + 
        geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
        data = en_coord_cat2, size =1, alpha = 0.5, colour = "grey30") +
        geom_text_repel(data = en_coord_cat2, aes(x = NMDS1, y = NMDS2), 
        label = row.names(en_coord_cat2), colour = "black")  

# Kruskal Wallis - State
names(NGSdata2)
NGSmerge<-inner_join(Metadata, NGSdata1, by="Abbreviation")
means<-summaryBy(OTU25~State, data=NGSmerge, FUN=c(mean, length, sd))
names(means)[names(means)=="OTU25.length"] <- "N"
means$se <- means$OTU25.sd / sqrt(means$N)
means

kruskal.test(NGSmerge$OTU20~NGSmerge$State) # Sig
kruskal.test(NGSdata2$OTU6~Metadata$State)
kruskal.test(NGSdata2$OTU48~Metadata$State) # Sig
kruskal.test(NGSdata2$OTU2~Metadata$State) # Sig
kruskal.test(NGSdata2$OTU16~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU45~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU22~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU34~Metadata$State)
kruskal.test(NGSdata2$OTU36~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU15~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU39~Metadata$State)
kruskal.test(NGSdata2$OTU37~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU59~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU38~Metadata$State)
kruskal.test(NGSdata2$OTU52~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU58~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU21~Metadata$State) #Sig

boxplot(NGSdata2$OTU20~Metadata$State, ylab="% Relative Abundance", main="OTU20")
kruskal.test(NGSdata2$OTU57~Metadata$State)
kruskal.test(NGSdata2$OTU49~Metadata$State)
kruskal.test(NGSdata2$OTU29~Metadata$State) # Sig
kruskal.test(NGSdata2$OTU14~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU17~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU18~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU12~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU11~Metadata$State)
kruskal.test(NGSdata2$OTU1~Metadata$State)
kruskal.test(NGSdata2$OTU42~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU61~Metadata$State)
kruskal.test(NGSdata2$OTU30~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU33~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU31~Metadata$State)
kruskal.test(NGSdata2$OTU35~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU50~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU13~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU27~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU46~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU7~Metadata$State)
kruskal.test(NGSdata2$OTU26~Metadata$State)
kruskal.test(NGSdata2$OTU23~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU4~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU5~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU28~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU2~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU3~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU55~Metadata$State)#Sig
kruskal.test(NGSdata2$OTU47~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU56~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU32~Metadata$State)
kruskal.test(NGSdata2$OTU29~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU44~Metadata$State)
kruskal.test(NGSdata2$OTU9~Metadata$State)
kruskal.test(NGSdata2$OTU40~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU24~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU10~Metadata$State)
kruskal.test(NGSdata2$OTU51~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU25~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU41~Metadata$State)
kruskal.test(NGSdata2$OTU43~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU8~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU53~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU19~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU54~Metadata$State) #Sig
kruskal.test(NGSdata2$OTU60~Metadata$State) #Sig
library(doBy)


# Kruskal Wallis - Invasion
kruskal.test(NGSdata2$OTU20~Metadata$Invasion)
kruskal.test(NGSdata2$OTU6~Metadata$Invasion)
kruskal.test(NGSdata2$OTU48~Metadata$Invasion)
kruskal.test(NGSdata2$OTU2~Metadata$Invasion)
kruskal.test(NGSdata2$OTU16~Metadata$Invasion)
kruskal.test(NGSdata2$OTU45~Metadata$Invasion)
kruskal.test(NGSdata2$OTU22~Metadata$Invasion)
kruskal.test(NGSdata2$OTU34~Metadata$Invasion) #Sig
kruskal.test(NGSdata2$OTU36~Metadata$Invasion) #Sig
kruskal.test(NGSdata2$OTU15~Metadata$Invasion)
kruskal.test(NGSdata2$OTU39~Metadata$Invasion)
kruskal.test(NGSdata2$OTU37~Metadata$Invasion)
kruskal.test(NGSdata2$OTU59~Metadata$Invasion)
kruskal.test(NGSdata2$OTU38~Metadata$Invasion) #Sig
kruskal.test(NGSdata2$OTU52~Metadata$Invasion)
kruskal.test(NGSdata2$OTU58~Metadata$Invasion)
kruskal.test(NGSdata2$OTU21~Metadata$Invasion)
kruskal.test(NGSdata2$OTU57~Metadata$Invasion) #Sig
kruskal.test(NGSdata2$OTU49~Metadata$Invasion)
kruskal.test(NGSdata2$OTU29~Metadata$Invasion)
kruskal.test(NGSdata2$OTU14~Metadata$Invasion) #Sig
kruskal.test(NGSdata2$OTU17~Metadata$Invasion)
kruskal.test(NGSdata2$OTU18~Metadata$Invasion)
kruskal.test(NGSdata2$OTU12~Metadata$Invasion)
kruskal.test(NGSdata2$OTU11~Metadata$Invasion)
kruskal.test(NGSdata2$OTU1~Metadata$Invasion) #Sig
kruskal.test(NGSdata2$OTU42~Metadata$Invasion)
kruskal.test(NGSdata2$OTU61~Metadata$Invasion)
kruskal.test(NGSdata2$OTU30~Metadata$Invasion)
kruskal.test(NGSdata2$OTU33~Metadata$Invasion)
kruskal.test(NGSdata2$OTU31~Metadata$Invasion)
kruskal.test(NGSdata2$OTU35~Metadata$Invasion)
kruskal.test(NGSdata2$OTU50~Metadata$Invasion)
kruskal.test(NGSdata2$OTU13~Metadata$Invasion)
kruskal.test(NGSdata2$OTU27~Metadata$Invasion)
kruskal.test(NGSdata2$OTU46~Metadata$Invasion)
kruskal.test(NGSdata2$OTU7~Metadata$Invasion)
kruskal.test(NGSdata2$OTU26~Metadata$Invasion)
kruskal.test(NGSdata2$OTU23~Metadata$Invasion)
kruskal.test(NGSdata2$OTU4~Metadata$Invasion)
kruskal.test(NGSdata2$OTU5~Metadata$Invasion)
kruskal.test(NGSdata2$OTU28~Metadata$Invasion)
kruskal.test(NGSdata2$OTU2~Metadata$Invasion)
kruskal.test(NGSdata2$OTU3~Metadata$Invasion)
kruskal.test(NGSdata2$OTU55~Metadata$Invasion)
kruskal.test(NGSdata2$OTU47~Metadata$Invasion)
kruskal.test(NGSdata2$OTU56~Metadata$Invasion)
kruskal.test(NGSdata2$OTU32~Metadata$Invasion)
kruskal.test(NGSdata2$OTU29~Metadata$Invasion)
kruskal.test(NGSdata2$OTU44~Metadata$Invasion)
kruskal.test(NGSdata2$OTU9~Metadata$Invasion)
kruskal.test(NGSdata2$OTU40~Metadata$Invasion)
kruskal.test(NGSdata2$OTU24~Metadata$Invasion)
kruskal.test(NGSdata2$OTU10~Metadata$Invasion)
kruskal.test(NGSdata2$OTU51~Metadata$Invasion)
kruskal.test(NGSdata2$OTU25~Metadata$Invasion)
kruskal.test(NGSdata2$OTU41~Metadata$Invasion)
kruskal.test(NGSdata2$OTU43~Metadata$Invasion)
kruskal.test(NGSdata2$OTU8~Metadata$Invasion)
kruskal.test(NGSdata2$OTU53~Metadata$Invasion)
kruskal.test(NGSdata2$OTU19~Metadata$Invasion)
kruskal.test(NGSdata2$OTU54~Metadata$Invasion) #Sig
kruskal.test(NGSdata2$OTU60~Metadata$Invasion)


# Indicator Species
library(indicspecies)
OTUs<-NGSdataCO
COinvade<-SoildataCO2$Invasion
inv=multipatt(OTUs, COinvade, func = "r.g", control=how(nperm=999))
summary(inv)

mergeCO<-NGSmerge %>% filter(State=="CO")
means<-summaryBy(OTU36~Invasion, data=mergeCO, FUN=c(mean, length, sd))
names(means)[names(means)=="OTU36.length"] <- "N"
means$se <- means$OTU36.sd / sqrt(means$N)
means

OTUs2<-NGSdataIL
ILinvade<-SoildataIL2$Invasion
inv2=multipatt(OTUs2, ILinvade, func = "r.g", control=how(nperm=999))
summary(inv2)

mergeIL<-NGSmerge %>% filter(State=="IL")
means<-summaryBy(OTU61~Invasion, data=mergeIL, FUN=c(mean, length, sd))
names(means)[names(means)=="OTU61.length"] <- "N"
means$se <- means$OTU61.sd / sqrt(means$N)
means

BothOTUs<-NGSdata2
states<-Soildata$State
inv3=multipatt(BothOTUs, states, func = "r.g", control=how(nperm=999))
summary(inv3)

invasion<-Soildata$Invasion
inv4=multipatt(BothOTUs, invasion, func= "r.g", control=how(nperm=999))
summary(inv4)




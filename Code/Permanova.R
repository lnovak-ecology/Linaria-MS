library(vegan)
library(readr)
NGSdata<-read_csv("OTU Relative Frequencies.csv")
names(NGSdata)
NGSdata<-NGSdata[-c(67,68),]
NGSdata[,2:62] <- lapply(NGSdata[,2:62], as.numeric)
SoilNGS<-read_csv("Nutrient Analyses For Sequences.csv")
names(SoilNGS)
SoilNGS$pH <-as.numeric(as.character(SoilNGS$pH))
SoilNGS$`Ammonium ug/g` <-as.numeric(as.character(SoilNGS$`Ammonium ug/g`))
SoilNGS$`Nitrate u/g` <-as.numeric(as.character(SoilNGS$`Nitrate u/g`))
SoilNGS$`Phosphate u/g` <-as.numeric(as.character(SoilNGS$`Phosphate u/g`))
SoilNGS<-SoilNGS[-c(73,74),]
library(dplyr)

metadata<-SoilNGS %>% dplyr::select(State, Invasion, Site, `Ammonium ug/g`, `Nitrate u/g`, `Phosphate u/g`, pH, Sand, Clay, Silt)
metadata$State<-as.factor(metadata$State)
metadata$Invasion<-as.factor(metadata$Invasion)
metadata$Site<-as.factor(metadata$Site)

set.seed(123)
NGSmerge<-inner_join(NGSdata, SoilNGS, by="Abbreviation")
NGSmerge<-as.data.frame(NGSmerge)
names(NGSmerge)
NGS1<-NGSmerge[,2:62]
names(NGS1)
NGS1[is.na(NGS1)]<-0


matrix<-as.matrix(NGS1)
sppr<-specnumber(matrix)
sppr_aov<-aov(sppr~State*Invasion, data=metadata)
summary(sppr_aov)
sum2<-as.data.frame(rowSums(matrix))
NGSdist<-vegdist(matrix, method = "bray")
NGS.div<-adonis2(NGSdist~State*Invasion, data=metadata, permutations = 999, method="bray")
NGS.div

library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis//pairwiseAdonis", force=T)
library(pairwiseAdonis)
factors<-interaction(metadata$State, metadata$Invasion)
names(NGS1)
pairwise.adonis(x=NGS1, factors=factors,sim.function='vegdist', sim.method='bray',p.adjust.m='bonferroni')

NGS.div0<-adonis2(NGSdist~State*Invasion*pH, data=metadata, permutations = 999, method="bray")
NGS.div0

NGS.div1<-adonis2(NGSdist~State*Invasion*`Ammonium ug/g`, data=metadata, permutations = 999, method="bray")
NGS.div1

NGS.div2<-adonis2(NGSdist~State*Invasion*`Nitrate u/g`, data=metadata, permutations = 999, method="bray")
NGS.div2

NGS.div3<-adonis2(NGSdist~State*Invasion*`Phosphate u/g`, data=metadata, permutations = 999, method="bray")
NGS.div3

NGS.div4<-adonis2(NGSdist~State*Invasion*Sand, data=metadata, permutations = 999, method="bray")
NGS.div4

NGS.div5<-adonis2(NGSdist~State*Invasion*Clay, data=metadata, permutations = 999, method="bray")
NGS.div5

NGS.div6<-adonis2(NGSdist~State*Invasion*Silt, data=metadata, permutations = 999, method="bray")
NGS.div6

# CO Locations
CONGS<-NGSdata[-c(41:78),]
CONGS[,2:62] <- lapply(CONGS[,2:62], as.numeric)
CONGS<-CONGS[,2:62]
COmetadata<-metadata[-c(41:78),]
matrix1<-as.matrix(CONGS)
CONGSdist<-vegdist(matrix1, method = "bray")

CONGS.div<-adonis2(CONGSdist~Invasion, data=COmetadata, permutations = 999, method="bray")
CONGS.div

CONGS.div0<-adonis2(CONGSdist~pH, data=COmetadata, permutations = 999, method="bray")
CONGS.div0

CONGS.div1<-adonis2(CONGSdist~`Ammonium ug/g`, data=COmetadata, permutations = 999, method="bray")
CONGS.div1

CONGS.div2<-adonis2(CONGSdist~`Nitrate u/g`, data=COmetadata, permutations = 999, method="bray")
CONGS.div2

CONGS.div3<-adonis2(CONGSdist~`Phosphate u/g`, data=COmetadata, permutations = 999, method="bray")
CONGS.div3

CONGS.div4<-adonis2(CONGSdist~Sand, data=COmetadata, permutations = 999, method="bray")
CONGS.div4

CONGS.div5<-adonis2(CONGSdist~Clay, data=COmetadata, permutations = 999, method="bray")
CONGS.div5

CONGS.div6<-adonis2(CONGSdist~Silt, data=COmetadata, permutations = 999, method="bray")
CONGS.div6


# IL Locations
ILNGS<-NGSdata[-c(1:40),]
ILNGS[,2:62] <- lapply(ILNGS[,2:62], as.numeric)
ILNGS<-ILNGS[,2:62]
ILmetadata<-metadata[-c(1:40),]
matrix2<-as.matrix(ILNGS)
ILNGSdist<-vegdist(matrix2, method = "bray")

ILNGS.div<-adonis2(ILNGSdist~Invasion, data=ILmetadata, permutations = 999, method="bray")
ILNGS.div

ILNGS.div0<-adonis2(ILNGSdist~pH, data=ILmetadata, permutations = 999, method="bray")
ILNGS.div0

ILNGS.div1<-adonis2(ILNGSdist~`Ammonium ug/g`, data=ILmetadata, permutations = 999, method="bray")
ILNGS.div1

ILNGS.div2<-adonis2(ILNGSdist~`Nitrate u/g`, data=ILmetadata, permutations = 999, method="bray")
ILNGS.div2

ILNGS.div3<-adonis2(ILNGSdist~`Phosphate u/g`, data=ILmetadata, permutations = 999, method="bray")
ILNGS.div3

ILNGS.div4<-adonis2(ILNGSdist~Sand, data=ILmetadata, permutations = 999, method="bray")
ILNGS.div4

ILNGS.div5<-adonis2(ILNGSdist~Clay, data=ILmetadata, permutations = 999, method="bray")
ILNGS.div5

ILNGS.div6<-adonis2(ILNGSdist~Silt, data=ILmetadata, permutations = 999, method="bray")
ILNGS.div6

# CO Invaded
CONGSINV<-CONGS[-c(3, 4, 7, 8, 11, 13, 16, 17, 20, 21, 23, 25, 26, 29, 30, 33, 34, 35, 38, 39),]
matrix3<-as.matrix(CONGSINV)
CONGSINV[,1:61] <- lapply(CONGSINV[,1:61], as.numeric)
CONGSINVdist<-vegdist(matrix3, method = "bray")
COInv<-COmetadata[COmetadata$Invasion=="Invaded",]

COInv.div<-adonis2(CONGSINVdist~pH, data=COInv, permutations = 999, method="bray")
COInv.div

COInv.div1<-adonis2(CONGSINVdist~`Ammonium ug/g`, data=COInv, permutations = 999, method="bray")
COInv.div1

COInv.div2<-adonis2(CONGSINVdist~`Nitrate u/g`, data=COInv, permutations = 999, method="bray")
COInv.div2

COInv.div3<-adonis2(CONGSINVdist~`Phosphate u/g`, data=COInv, permutations = 999, method="bray")
COInv.div3

COInv.div4<-adonis2(CONGSINVdist~Sand, data=COInv, permutations = 999, method="bray")
COInv.div4

COInv.div5<-adonis2(CONGSINVdist~Clay, data=COInv, permutations = 999, method="bray")
COInv.div5

COInv.div6<-adonis2(CONGSINVdist~Silt, data=COInv, permutations = 999, method="bray")
COInv.div6

# CO Uninvaded
CONGSUninvaded<-CONGS[-c(1,2,5,6,9,10,12,14,15,18,19,22,24,27,28,31,32,36,37,40),]
CONGSUninvaded[,1:61] <- lapply(CONGSUninvaded[,1:61], as.numeric)
matrix4<-as.matrix(CONGSUninvaded)
CONGSUNdist<-vegdist(matrix4, method = "bray")
COUninvaded<-COmetadata[COmetadata$Invasion=="Uninvaded",]

COUnv.div<-adonis2(CONGSUNdist~pH, data=COUninvaded, permutations = 999, method="bray")
COUnv.div

COUnv.div1<-adonis2(CONGSUNdist~`Ammonium ug/g`, data=COUninvaded, permutations = 999, method="bray")
COUnv.div1

COUnv.div2<-adonis2(CONGSUNdist~`Nitrate u/g`, data=COUninvaded, permutations = 999, method="bray")
COInv.div2

COUnv.div3<-adonis2(CONGSUNdist~`Phosphate u/g`, data=COUninvaded, permutations = 999, method="bray")
COUnv.div3

COUnv.div4<-adonis2(CONGSUNdist~Sand, data=COUninvaded, permutations = 999, method="bray")
COUnv.div4

COUnv.div5<-adonis2(CONGSUNdist~Clay, data=COUninvaded, permutations = 999, method="bray")
COUnv.div5

COUnv.div6<-adonis2(CONGSUNdist~Silt, data=COUninvaded, permutations = 999, method="bray")
COUnv.div6

# IL Invaded
ILNGSINV<-ILNGS[-c(3, 4, 7, 8, 11, 13, 16, 17, 20, 21, 23, 25, 26, 27,28, 31,32, 33, 36, 37),]
matrix5<-as.matrix(ILNGSINV)
ILNGSINV[,1:61] <- lapply(ILNGSINV[,1:61], as.numeric)
ILNGSINVdist<-vegdist(matrix5, method = "bray")
ILInv<-ILmetadata[ILmetadata$Invasion=="Invaded",]

ILInv.div<-adonis2(ILNGSINVdist~pH, data=ILInv, permutations = 999, method="bray")
ILInv.div

ILInv.div1<-adonis2(ILNGSINVdist~`Ammonium ug/g`, data=ILInv, permutations = 999, method="bray")
ILInv.div1

ILInv.div2<-adonis2(ILNGSINVdist~`Nitrate u/g`, data=ILInv, permutations = 999, method="bray")
ILInv.div2

ILInv.div3<-adonis2(ILNGSINVdist~`Phosphate u/g`, data=ILInv, permutations = 999, method="bray")
ILInv.div3

ILInv.div4<-adonis2(ILNGSINVdist~Sand, data=ILInv, permutations = 999, method="bray")
ILInv.div4

ILInv.div5<-adonis2(ILNGSINVdist~Clay, data=ILInv, permutations = 999, method="bray")
ILInv.div5

ILInv.div6<-adonis2(ILNGSINVdist~Silt, data=ILInv, permutations = 999, method="bray")
ILInv.div6

# IL Uninvaded
ILNGSUNV<-ILNGS[-c(1,2,5,6,9,10,12,14,15,18,19,22,24,29,30,34,35,38),]
matrix6<-as.matrix(ILNGSUNV)
ILNGSUNV[,1:61] <- lapply(ILNGSUNV[,1:61], as.numeric)
ILNGSUNVdist<-vegdist(matrix6, method = "bray")
ILUnv<-ILmetadata[ILmetadata$Invasion=="Uninvaded",]

ILUnv.div<-adonis2(ILNGSUNVdist~pH, data=ILUnv, permutations = 999, method="bray")
ILUnv.div

ILUnv.div1<-adonis2(ILNGSUNVdist~`Ammonium ug/g`, data=ILUnv, permutations = 999, method="bray")
ILUnv.div1

ILUnv.div2<-adonis2(ILNGSUNVdist~`Nitrate u/g`, data=ILUnv, permutations = 999, method="bray")
ILUnv.div2

ILUnv.div3<-adonis2(ILNGSUNVdist~`Phosphate u/g`, data=ILUnv, permutations = 999, method="bray")
ILUnv.div3

ILUnv.div4<-adonis2(ILNGSUNVdist~Sand, data=ILUnv, permutations = 999, method="bray")
ILUnv.div4

ILUnv.div5<-adonis2(ILNGSUNVdist~Clay, data=ILUnv, permutations = 999, method="bray")
ILUnv.div5

ILUnv.div6<-adonis2(ILNGSUNVdist~Silt, data=ILUnv, permutations = 999, method="bray")
ILUnv.div6


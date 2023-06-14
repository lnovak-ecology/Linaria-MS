library(readr)
isodata<-read_csv("Uptake for Analyses.csv")
names(isodata)
ILisodata<-isodata[isodata$State=="IL",]
COisodata<-isodata[isodata$State=="CO",]
library(car)
ILisodata$Invasion<-as.factor(ILisodata$Invasion)
COisodata$Invasion<-as.factor(COisodata$Invasion)
ILisodata$Fertilizer<-as.factor(ILisodata$Fertilizer)
COisodata$Fertilizer<-as.factor(COisodata$Fertilizer)
ILisodata$NUR<-as.integer(ILisodata$NUR)
COisodata$NUR<-as.integer(COisodata$NUR)
library(dplyr)
ILlinaria<-ILisodata[ILisodata$Grown=="Linaria",]
hist(ILlinaria$NUR)

# IL Linaria - not enough data
qqPlot(ILlinaria$NUR)

# IL Native - not sig.
ILnative<-ILisodata[ILisodata$Grown=="Native",]
hist(ILnative$NUR)
qqPlot(ILnative$NUR)
mean(ILnative$NUR)
var(ILnative$NUR)
shapiro.test(ILnative$NUR)
shapiro.test(log(ILnative$NUR+1))

lm8<-lm(log(ILnative$NUR+1)~ILnative$Invasion*ILnative$Fertilizer)
plot(lm8)
Anova(lm8)
confint(lm8)
chisq.test(ILnative$NUR, ILnative$Invasion)
chisq.test(ILnative$NUR, ILnative$Fertilizer)

# CO Linaria - not sig.
COlinaria<-COisodata[COisodata$Grown=="Linaria",]
hist(COlinaria$NUR)
qqPlot(COlinaria$NUR)
shapiro.test(COlinaria$NUR)
shapiro.test(log(COlinaria$NUR+1))
# overdispersed

lm4<-lm(log(COlinaria$NUR+1)~COlinaria$Invasion*COlinaria$Fertilizer)
plot(lm4)
Anova(lm4)

# CO Native - not sig.
COnative<-COisodata[COisodata$Grown=="Native",]
hist(COnative$NUR)
qqPlot(COnative$NUR)
shapiro.test(COnative$NUR)

lm9<-lm(COnative$NUR~COnative$Invasion*COnative$Fertilizer)
plot(lm9)
Anova(lm9)
# Best fit, nothing significant

library (RColorBrewer)
cbPalette <- c(brewer.pal(4, "YlGnBu"))
library(doBy)
COdata <- summaryBy(NUR ~ Grown + Invasion + Fertilizer, data=COisodata, FUN=c(length,mean,sd))
COdata
names(COdata)[names(COdata)=="NUR.length"] <- "N"
COdata$NUR.se <- COdata$NUR.sd / sqrt(COdata$N)

ILdata <- summaryBy(NUR ~ Grown + Invasion + Fertilizer, data=ILisodata, FUN=c(length,mean,sd))
ILdata
names(ILdata)[names(ILdata)=="NUR.length"] <- "N"
ILdata$NUR.se <- ILdata$NUR.sd / sqrt(ILdata$N)

library(ggplot2)
ggplot(COdata, aes(x=Invasion, y=NUR.mean, fill=Fertilizer))+geom_bar(stat="identity",
color="black",position=position_dodge()) + geom_errorbar(aes(ymin=NUR.mean, ymax=NUR.mean+NUR.se),
width=.2,position=position_dodge(.9)) +ylim(0, 1000)+labs(y="NUR (μg N/g shoot)",
size=14)+labs(x= "Soil Type", size =13)+scale_fill_manual(values=cbPalette)+theme(axis.text.x = element_text(size=11, colour= "black", angle = 0),axis.text.y =
element_text(size=12, colour="black"),strip.text = element_text(size = 10))+theme(axis.title.x =
element_text(size=12, colour="black"), axis.title.y = element_text(size= 12, colour = "black")) + theme(legend.title = element_blank())+theme(panel.grid = element_blank(),legend.position = "bottom")+facet_wrap(~ COdata$Grown, ncol=1,scales = "free_y")+theme(axis.ticks.x = element_blank())+theme(legend.background = element_blank()) + theme(panel.grid = element_blank())

ggplot(ILdata, aes(x=Invasion, y=NUR.mean, fill=Fertilizer))+geom_bar(stat="identity",color="black",position=position_dodge2(preserve="single", padding=0)) + geom_errorbar(aes(ymin=NUR.mean, ymax=NUR.mean+NUR.se),
width=.2,position=position_dodge(.9)) +ylim(0, 2000)+labs(y="NUR (μg N/g shoot)", size=14)+labs(x= "Soil Type", size =13)+scale_fill_manual(values=cbPalette)+theme(axis.text.x = element_text(size=12, colour= "black", angle = 0),axis.text.y =
element_text(size=12, colour="black"),strip.text = element_text(size = 10))+theme(axis.title.x = element_text(size=12, colour="black"), axis.title.y = element_text(size= 12, colour = "black")) + theme(legend.title = element_blank())+facet_wrap(~ ILdata$Grown, ncol=1,scales = "free_y")+theme(axis.ticks.x = element_blank())+theme(panel.grid = element_blank(),legend.position = "bottom")+theme(legend.background = element_blank()) + theme(panel.grid = element_blank())

#Wt% N
# IL Linaria - not enough data
# IL Native
hist(ILlinaria$`wt% N`)
hist(ILnative$`wt% N`)
qqPlot(ILlinaria$`wt% N`)
qqPlot(ILnative$`wt% N`)
mean(ILnative$`wt% N`)
var(ILnative$`wt% N`)

lm13<-lm(ILnative$`wt% N`~ILnative$Invasion+ILnative$Fertilizer)
summary(aov(lm13))
confint(lm13)
# Best model, but nothing significant

lm14<-lm(ILnative$`wt% N`~ILnative$Invasion*ILnative$Fertilizer)

AICctab(lm13, lm14)

# CO Linaria
qqPlot(COlinaria$`wt% N`)
hist(COlinaria$`wt% N`)
mean(COlinaria$`wt% N`)
var(COlinaria$`wt% N`)

lm18<-lm(COlinaria$`wt% N`~COlinaria$Invasion+COlinaria$Fertilizer)
summary(aov(lm18))
confint(lm18)
# Best model but nothing significant

lm19<-lm(COlinaria$`wt% N`~COlinaria$Invasion*COlinaria$Fertilizer)

AICctab(lm18, lm19)

# CO Native
qqPlot(COnative$`wt% N`)
hist(COnative$`wt% N`)
mean(COnative$`wt% N`)
var(COnative$`wt% N`)

lm20<-lm(COnative$`wt% N`~COnative$Invasion+COnative$Fertilizer)
summary(aov(lm20))
confint(lm20)
# Best model but nothing significant

lm21<-lm(COnative$`wt% N`~COnative$Invasion*COnative$Fertilizer)
AICctab(lm20, lm21)

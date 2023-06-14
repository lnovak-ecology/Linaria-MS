library(readr)
biomass<-read_csv("Biomass Data.csv")
names(biomass)
library(ggplot2)
library(tidyverse)
biomass$`Total Root Mass 1` <- as.numeric(as.character(biomass$`Total Root Mass 1`))
biomass$`Total Root Mass 2` <- as.numeric(as.character(biomass$`Total Root Mass 2`))
biomass$TotalRootMass1<-biomass$`Total Root Mass 1`
biomass$ShootMass1<-biomass$`Shoot Mass 1`

biomass$ShootMass1[biomass$ShootMass1==0]<-NA
biomass$TotalRootMass1[biomass$TotalRootMass1==0]<-NA
names(biomass)
biomass<-biomass[complete.cases(biomass[, 14:15]),]
biomass$Grown[biomass$Grown=="Native and Linaria"]<-NA
biomass3<-biomass[complete.cases(biomass[, 2]),]
biomassIL<-biomass3 %>% filter(State == "IL")
biomassIL$Invasion <- as.factor((biomassIL$Invasion))
biomassIL$Treatment <- as.factor((biomassIL$Treatment))
biomassIL$Grown <- as.factor((biomassIL$Grown))

biomassCO<-biomass3 %>% filter(State == "CO")
biomassCO$Invasion <- as.factor((biomassCO$Invasion))
biomassCO$Treatment <- as.factor((biomassCO$Treatment))
biomassCO$Grown <- as.factor((biomassCO$Grown))
names(biomassCO)


COLinariaShoots<-biomassCO[biomassCO$Grown=="Linaria",]
hist(COLinariaShoots$ShootMass1)
qqPlot(COLinariaShoots$ShootMass1)

CONativeShoots<-biomassCO[biomassCO$Grown=="Native",]
hist(CONativeShoots$ShootMass1)
qqPlot(CONativeShoots$ShootMass1)

hist(ILLinariaShoots$ShootMass1)
qqPlot(ILLinariaShoots$ShootMass1)

hist(ILNativeShoots$ShootMass1)
qqPlot(ILNativeShoots$ShootMass1)

# Roots
hist(COLinariaShoots$TotalRootMass1)
qqPlot(COLinariaShoots$TotalRootMass1)

hist(CONativeShoots$TotalRootMass1)
qqPlot(CONativeShoots$TotalRootMass1)

hist(ILLinariaShoots$TotalRootMass1)
qqPlot(ILLinariaShoots$TotalRootMass1)

hist(ILNativeShoots$TotalRootMass1)
qqPlot(ILNativeShoots$TotalRootMass1)

# CO Linaria Shoot Mass - Invasion significant!
library(lmer)
lm3<-lm(COLinariaShoots$ShootMass1~COLinariaShoots$Invasion*COLinariaShoots$Treatment)
plot(lm3)
shapiro.test(COLinariaShoots$`Shoot Mass 1`)
# Not normal
resid.ssq2 <- sum(residuals(lm3,type="pearson")^2)  
resid.df2 <- nrow(subset(COLinariaShoots,!is.na(Invasion) & !is.na(Treatment)))-length(coef(lm3))
resid.ssq2/resid.df2
# Not overdispersed, but heteroscedastic
lm4<-lmer((COLinariaShoots$ShootMass1)^2~COLinariaShoots$Invasion*COLinariaShoots$Treatment+(1|COLinariaShoots$Site))
plot(lm4)
shapiro.test((COLinariaShoots$`Shoot Mass 1`)^2)
summary(lm4)
Anova(lm4)

# CO Native Shoots - not sig.
# No random effect because not enough data (isSingular error)
lm9<-lm(CONativeShoots$ShootMass1~CONativeShoots$Invasion*CONativeShoots$Treatment)
plot(lm9)
summary(lm9)
Anova(lm9)


# IL Linaria Shoots - not enough data
ILLinariaShoots<-biomassIL[biomassIL$Grown=="Linaria",]
hist(ILLinariaShoots$ShootMass1)
qqPlot(ILLinariaShoots$ShootMass1)

lm14<-lm(ILLinariaShoots$ShootMass1~ILLinariaShoots$Invasion*ILLinariaShoots$Treatment)
AIC(lm14)
resid14<-resid(lm14)
plot(resid14)

AICctab(lm13, lm14)

# IL Native Shoots - not sig.
ILNativeShoots<-biomassIL[biomassIL$Grown=="Native",]
hist(ILNativeShoots$ShootMass1)
qqPlot(ILNativeShoots$ShootMass1)
shapiro.test(ILNativeShoots$`Shoot Mass 1`)

lm18<-lmer(ILNativeShoots$ShootMass1~ILNativeShoots$Invasion*ILNativeShoots$Treatment+ (1|ILNativeShoots$Site))
resid18<-resid(lm18)
plot(resid18)
Anova(lm18)
resid.ssq2 <- sum(residuals(lm18,type="pearson")^2)  
resid.df2 <- nrow(subset(ILNativeShoots,!is.na(Invasion) & !is.na(Treatment)))-length(coef(lm18))
resid.ssq2/resid.df2


# CO Native Roots - Invasion sig.
hist(CONativeShoots$TotalRootMass1)
qqPlot(CONativeShoots$TotalRootMass1)
shapiro.test(CONativeShoots$TotalRootMass1)
shapiro.test(log(CONativeShoots$`Total Root Mass 1`))


lm23<-lmer(log(CONativeShoots$TotalRootMass1+1)~ CONativeShoots$Invasion*CONativeShoots$Treatment + (1|CONativeShoots$Site))
plot(lm23)
Anova(lm23, type=2)
library(emmeans)
aggregate(TotalRootMass1~Invasion, data=CONativeShoots, FUN=mean)

# CO Linaria Roots - Invasion sig!
hist(COLinariaShoots$TotalRootMass1)
qqPlot(COLinariaShoots$TotalRootMass1)
shapiro.test(log(COLinariaShoots$`Total Root Mass 1`))

lm28<-lmer(log(COLinariaShoots$TotalRootMass1+1)~COLinariaShoots$Invasion*COLinariaShoots$Treatment + (1|COLinariaShoots$Site))
summary(lm28)
Anova(lm28)
lm28pair<-emmeans(lm28, specs=pairwise~Invasion*Treatment)
lm28pair
lm28pair$contrasts %>% rbind(adjust="sidak")
aggregate(TotalRootMass1~Invasion, data=COLinariaShoots, FUN=mean)

# IL Linaria Roots - not enough data
hist(ILLinariaShoots$TotalRootMass1)
qqPlot(ILLinariaShoots$TotalRootMass1)

lm34<-lm(ILLinariaShoots$TotalRootMass1~ILLinariaShoots$Invasion*ILLinariaShoots$Treatment)


# IL Native Roots - not sig.
hist(ILNativeShoots$TotalRootMass1)
qqPlot(ILNativeShoots$TotalRootMass1)
shapiro.test(ILNativeShoots$`Total Root Mass 1`)
shapiro.test(log(ILNativeShoots$`Total Root Mass 1`))

lm38<-lmer(log(ILNativeShoots$TotalRootMass1+1)~ILNativeShoots$Invasion*ILNativeShoots$Treatment + (1|ILNativeShoots$Site))
summary(lm38)
Anova(lm38)
plot(lm38)

ILroots <- biomassIL %>%
  group_by(Grown, Invasion, Treatment) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(TotalRootMass1),
    sd=sd(TotalRootMass1)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))
ILshoots <- biomassIL %>%
  group_by(Grown, Invasion, Treatment) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(ShootMass1),
    sd=sd(ShootMass1)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))
COroots <- biomassCO %>%
  group_by(Grown, Invasion, Treatment) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(TotalRootMass1),
    sd=sd(TotalRootMass1)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))
COshoots <- biomassCO %>%
  group_by(Grown, Invasion, Treatment) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(ShootMass1),
    sd=sd(ShootMass1)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

library(ggplot2)
library(devtools)
library(remotes)
library(RColorBrewer)
library(ggpubr)
cbPalette <- c(brewer.pal(4, "YlGnBu"))

ggplot(COroots, aes(x = Invasion, y = mean, fill = Treatment)) + geom_bar(stat="identity", color="black",position=position_dodge2(preserve="single", padding=0)) + geom_errorbar(aes(ymin=mean, ymax=mean+se), width=.2,position=position_dodge(.9)) + theme(panel.grid = element_blank(),legend.position = "bottom") + ylab("Root Mass (g)") + xlab("Soil Type") + scale_fill_manual(values= cbPalette)+theme(axis.text.x = element_text(size=11, colour= "black", angle = 0),axis.text.y = element_text(size=12, colour="black"),strip.text = element_text(size = 10))+theme(axis.title.x = element_text(size=12, colour="black"), axis.title.y = element_text(size= 12, colour = "black"))+ theme(legend.title = element_blank())+theme(axis.ticks.x = element_blank())+theme(panel.grid = element_blank(),legend.position = "bottom")+theme(legend.title = element_blank())+theme(legend.text = element_text(size = 10, colour = "black"))+theme(legend.background = element_blank())+ylim(0, 2.0) + geom_signif(data = data.frame(Grown = c("Linaria","Native")), aes(y_position=c(0.9, 1.95), xmin=c(1.0, 1.0), xmax=c(2.0, 2.0), annotations=c("*", "*"), fill=NULL), tip_length=0.01, manual=T) + facet_grid(~ Grown, scales = "free")

ggplot(ILroots, aes(x = Invasion, y = mean, fill = Treatment)) + geom_bar(stat="identity", color="black",position=position_dodge2(preserve="single", padding=0)) + geom_errorbar(aes(ymin=mean, ymax=mean+se), width=.2,position=position_dodge(.9)) + theme(panel.grid = element_blank(),legend.position = "bottom") + ylab("Root Mass (g)") + xlab("Soil Type") + scale_fill_manual(values= cbPalette)+theme(axis.text.x = element_text(size=11, colour= "black", angle = 0),axis.text.y = element_text(size=12, colour="black"),strip.text = element_text(size = 10))+theme(axis.title.x = element_text(size=12, colour="black"), axis.title.y = element_text(size= 12, colour = "black"))+ theme(legend.title = element_blank())+theme(axis.ticks.x = element_blank())+ theme(panel.grid = element_blank(), legend.position = "bottom") + theme(legend.title = element_blank())+theme(legend.text = element_text(size = 10, colour = "black"))+theme(legend.background = element_blank())+ylim(0,1.0) + facet_grid(~Grown, scales = "free_y")

ggplot(COshoots, aes(x = Invasion, y = mean, fill = Treatment)) + geom_bar(stat="identity", color="black",position=position_dodge2(preserve="single", padding=0)) + geom_errorbar(aes(ymin=mean, ymax=mean+se), width=.2,position=position_dodge(.9)) + theme_bw(base_size = 15) + theme(panel.grid = element_blank(),legend.position = "bottom") + ylab("Shoot Mass (g)") + xlab("Soil Type") + scale_fill_manual(values= cbPalette)+theme(axis.text.x = element_text(size=11, colour= "black", angle = 0),axis.text.y = element_text(size=12, colour="black"),strip.text = element_text(size = 14))+theme(axis.title.x = element_text(size=16, colour="black"), axis.title.y = element_text(size= 20, colour = "black"))+ theme(legend.title = element_blank())+theme(axis.ticks.x = element_blank())+ theme(legend.title = element_blank())+theme(legend.text = element_text(size = 12, colour = "black"))+theme(legend.background = element_blank()) + geom_signif(data = data.frame(Grown = "Linaria"), aes(y_position=1.1, xmin=1.0, xmax=2.0, annotations="*", fill=NULL), tip_length=0.01, manual=T) + facet_grid(~ Grown, scales = "free")

ggplot(ILshoots, aes(x = Invasion, y = mean, fill = Treatment)) + geom_bar(stat="identity", color="black",position=position_dodge2(preserve="single", padding=0)) + geom_errorbar(aes(ymin=mean, ymax=mean+se), width=.2,position=position_dodge(.9)) + facet_wrap(~ ILshoots$Grown, ncol=1, scales = "free_y")+ theme_bw(base_size = 15) + theme(panel.grid = element_blank(),legend.position = "bottom") + ylab("Shoot Mass (g)") + xlab("Soil Type") + scale_fill_manual(values= cbPalette)+theme(axis.text.x = element_text(size=11, colour= "black", angle = 0),axis.text.y = element_text(size=12, colour="black"),strip.text = element_text(size = 14))+theme(axis.title.x = element_text(size=16, colour="black"), axis.title.y = element_text(size= 20, colour = "black"))+ theme(legend.title = element_blank())+theme(axis.ticks.x = element_blank())+theme(legend.position = c(0.15, 0.90))+ theme(legend.title = element_blank())+theme(legend.text = element_text(size = 12, colour = "black"))+theme(legend.background = element_blank())



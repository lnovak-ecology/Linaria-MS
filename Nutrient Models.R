library(readr)
nutrients<-read_csv("Nutrient Analyses For Sequences.csv")
names(nutrients)
nutrients<-dplyr::select(nutrients, -c(Clay, Silt, Sand))
nutrients<-nutrients[complete.cases(nutrients[,6:10]),]
hist(nutrients$`Ammonium ug/g`)
qqPlot(nutrients$`Ammonium ug/g`)
mean(nutrients$`Ammonium ug/g`)
var(nutrients$`Ammonium ug/g`)
hist(nutrients$pH)
qqPlot(nutrients$pH)
library(MASS)
library(bbmle)
hist(nutrients$`Nitrate u/g`)
qqPlot(nutrients$`Nitrate u/g`)
mean(nutrients$`Nitrate u/g`)
var(nutrients$`Nitrate u/g`)
hist(nutrients$`Phosphate u/g`)
qqPlot(nutrients$`Phosphate u/g`)
mean(nutrients$`Phosphate u/g`)
var(nutrients$`Phosphate u/g`)
# overdispersed


#Both States - Nitrate
library(lme4)
lm1<-lm(log(nutrients$`Nitrate u/g`+10)~nutrients$Invasion*nutrients$State)
hist(lm1$residuals)
plot(lm1)
lm1<-lmer(log(nutrients$`Nitrate u/g`+10)~nutrients$Invasion*nutrients$State + (1|nutrients$Site))
Anova(lm1)
chisq.test(nutrients$`Nitrate u/g`, nutrients$Invasion, correct=F)
chisq.test(nutrients$`Nitrate u/g`, nutrients$State, correct=F)

# Both States- Ammonium
lm2<-lm(nutrients$`Ammonium ug/g`~nutrients$Invasion*nutrients$State)
plot(lm2)
hist(lm2$residuals)
lm2<-lmer(log(nutrients$`Ammonium ug/g`+10)~nutrients$Invasion*nutrients$State + (1|nutrients$Site))
summary(lm2)
Anova(lm2)

# Both States- phosphate
lm3<-lm(`Phosphate u/g`~Invasion*State, data=nutrients)
plot(lm3)
hist(lm3$residuals)
lm3<-lmer(log(`Phosphate u/g`+10)~Invasion*State + (1|Site), data=nutrients)
plot(lm3)
Anova(lm3)
chisq.test(nutrients$`Phosphate u/g`, nutrients$Invasion, correct=F)
chisq.test(nutrients$`Phosphate u/g`, nutrients$State, correct=F)
confint(lm3)
em<-emmeans(lm3, specs=pairwise~Invasion*State, adjust="tukey")
em$contrasts %>%
  rbind(adjust="sidak")

# Both States pH
lm4<-lm(nutrients$pH~nutrients$Invasion*nutrients$State)
plot(lm4)
hist(lm4$residuals)
lm4<-lmer(log(nutrients$pH+10)~nutrients$Invasion*nutrients$State + (1|nutrients$Site))
Anova(lm4)
confint(lm4)
chisq.test(nutrients$pH, nutrients$Invasion, correct=F)
chisq.test(nutrients$pH, nutrients$State, correct=F)


#CO
COnutrients<-nutrients[nutrients$State=="CO",]
aov5<-lm(COnutrients$`Nitrate u/g`~COnutrients$Invasion)
summary(aov5)
COnutrients$Invasion<-as.factor(COnutrients$Invasion)
CONO3mean <- COnutrients %>%
  group_by(Invasion) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(`Nitrate u/g`), sd=sd(`Nitrate u/g`)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))


aov6<-lm(COnutrients$`Ammonium ug/g`~COnutrients$Invasion)
summary(aov6)
CONO4mean <- COnutrients %>%
  group_by(Invasion) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(`Ammonium ug/g`), sd=sd(`Ammonium ug/g`)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

aov7<-aov(COnutrients$`Phosphate u/g`~COnutrients$Invasion)
summary(aov7)
COnutrients<-COnutrients[!is.na(COnutrients$`Phosphate u/g`),]
COPmean <- COnutrients %>%
  group_by(Invasion) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(`Phosphate u/g`), sd=sd(`Phosphate u/g`)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

aov8<-lm(COnutrients$pH~COnutrients$Invasion)
summary(aov8)
COphmean <- COnutrients %>%
  group_by(Invasion) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(pH), sd=sd(pH)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))




#IL
ILnutrients<-nutrients[nutrients$State=="IL",]
aov9<-lm(ILnutrients$`Nitrate u/g`~ILnutrients$Invasion)
summary(aov9)
ILNO3mean <- ILnutrients %>%
  group_by(Invasion) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(`Nitrate u/g`), sd=sd(`Nitrate u/g`)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

aov10<-lm(ILnutrients$`Ammonium ug/g`~ILnutrients$Invasion)
summary(aov10)
ILNO4mean <- ILnutrients %>%
  group_by(Invasion) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(`Ammonium ug/g`),
    sd=sd(`Ammonium ug/g`)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

aov11<-lm(ILnutrients$`Phosphate u/g`~ILnutrients$Invasion)
summary(aov11)
ILnutrients<-ILnutrients[!is.na(ILnutrients$`Phosphate u/g`),]
ILPmean <- ILnutrients %>%
  group_by(Invasion) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(`Phosphate u/g`), sd=sd(`Phosphate u/g`)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

aov12<-lm(ILnutrients$pH~ILnutrients$Invasion)
summary(aov12)
ILphmean <- ILnutrients %>%
  group_by(Invasion) %>%
  dplyr::summarise(
    n=n(),
    mean=mean(pH), sd=sd(pH)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

# Both Means
Bothphmean <- nutrients %>%
  group_by(State) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(pH), sd=sd(pH)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

BothNO3mean <- nutrients %>%
  group_by(State) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(`Nitrate u/g`), sd=sd(`Nitrate u/g`)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))
BothNO4mean <- nutrients %>%
  group_by(State) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(`Ammonium ug/g`), sd=sd(`Ammonium ug/g`)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))
BothPO4mean <- nutrients %>%
  group_by(State) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(`Phosphate u/g`), sd=sd(`Phosphate u/g`)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))
ggplot(Bothphmean, aes(x = State, y = mean, fill=State)) + geom_bar(stat="identity", color="black") + geom_errorbar(aes(ymin=mean, ymax=mean+se), width=.2, position=position_dodge(.9)) + theme_bw(base_size = 15) + theme(panel.grid = element_blank(), legend.position = "bottom") + ylab("Mean pH") + xlab("State")+theme(axis.text.x = element_text(size=11, colour= "black", angle = 0), axis.text.y = element_text(size=12, colour="black"), strip.text = element_text(size = 14))+theme(axis.title.x = element_text(size=16, colour="black"), axis.title.y = element_text(size= 20, colour = "black")) +theme(axis.ticks.x = element_blank())
ggplot(BothNO3mean, aes(x = State, y = mean, fill=State)) + geom_bar(stat="identity", color="black") + geom_errorbar(aes(ymin=mean, ymax=mean+se), width=.2, position=position_dodge(.9)) + theme_bw(base_size = 15) + theme(panel.grid = element_blank(), legend.position = "bottom") + ylab("Mean NO3") + xlab("State")+theme(axis.text.x = element_text(size=11, colour= "black", angle = 0), axis.text.y = element_text(size=12, colour="black"), strip.text = element_text(size = 14))+theme(axis.title.x = element_text(size=16, colour="black"), axis.title.y = element_text(size= 20, colour = "black")) +theme(axis.ticks.x = element_blank())
ggplot(BothNO4mean, aes(x = State, y = mean, fill=State)) + geom_bar(stat="identity", color="black") + geom_errorbar(aes(ymin=mean, ymax=mean+se), width=.2, position=position_dodge(.9)) + theme_bw(base_size = 15) + theme(panel.grid = element_blank(), legend.position = "bottom") + ylab("Mean NO4") + xlab("State")+theme(axis.text.x = element_text(size=11, colour= "black", angle = 0), axis.text.y = element_text(size=12, colour="black"), strip.text = element_text(size = 14))+theme(axis.title.x = element_text(size=16, colour="black"), axis.title.y = element_text(size= 20, colour = "black")) +theme(axis.ticks.x = element_blank())
ggplot(BothPO4mean, aes(x = State, y = mean, fill=State)) + geom_bar(stat="identity", color="black") + geom_errorbar(aes(ymin=mean, ymax=mean+se), width=.2, position=position_dodge(.9)) + theme_bw(base_size = 15) + theme(panel.grid = element_blank(), legend.position = "bottom") + ylab("Mean PO4") + xlab("State")+theme(axis.text.x = element_text(size=11, colour= "black", angle = 0), axis.text.y = element_text(size=12, colour="black"), strip.text = element_text(size = 14))+theme(axis.title.x = element_text(size=16, colour="black"), axis.title.y = element_text(size= 20, colour = "black")) +theme(axis.ticks.x = element_blank())  

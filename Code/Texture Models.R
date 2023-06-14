texture<-read_csv("Nutrient Analyses For Sequences.csv")
names(texture)
library(dplyr)
library(bbmle)
textureIL<-texture %>% filter(State == "IL")
textureCO<-texture %>% filter(State == "CO")
qqPlot(texture$Clay)
hist(texture$Clay)
qqPlot(texture$Silt)
hist(texture$Silt)
qqPlot(texture$Sand)
hist(texture$Sand)

# Try 2 variables interacting
library(lme4)
# Clay
lm1<-lm(texture$Clay ~ texture$State + texture$Sand + texture$Silt)
plot(lm1)
summary(lm1)
hist(lm1$residuals)
lm1<-lmer(textureCO$Clay ~ textureCO$Invasion + textureCO$Silt + textureCO$Sand + (1|textureCO$Site))
Anova(lm1)

# Silt
lm2<-lm(texture$Silt~ texture$State*texture$Invasion)
plot(lm2)
hist(lm2$residuals)
lm2<-lmer(Silt~State*Invasion+(1|Site), data=texture)
Anova(lm2)
chisq.test(texture$Silt, texture$State)
chisq.test(texture$Silt, texture$Invasion)

#Sand
lm3<-lm(texture$Sand~texture$State*texture$Invasion)
plot(lm3)
hist(lm3$residuals)
lm3<-lmer(Sand~State*Invasion+(1|Site), data=texture)
Anova(lm3)
chisq.test(texture$Sand, texture$State)
chisq.test(texture$Sand, texture$Invasion)

lm7<-lm(texture$Sand~texture$State+texture$Invasion)
summary(aov(lm7))

AICctab(lm6, lm7)

# Means

ILclay <- textureIL %>%
  group_by(Invasion) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(Clay), sd=sd(Clay)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

ILsilt <- textureIL %>%
  group_by(Invasion) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(Silt), sd=sd(Silt)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

ILsand <- textureIL %>%
  group_by(Invasion) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(Sand), sd=sd(Sand)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

COclay <- textureCO %>%
  group_by(Invasion) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(Clay), sd=sd(Clay)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

COsilt <- textureCO %>%
  group_by(Invasion) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(Silt), sd=sd(Silt)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

COsand <- textureCO %>%
  group_by(Invasion) %>%
  dplyr::summarise( 
    n=n(),
    mean=mean(Sand), sd=sd(Sand)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))


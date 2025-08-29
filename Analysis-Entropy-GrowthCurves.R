#SF1 analysis
remove(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)

#import data
setwd("XXXXX")

SF1<-read.csv("dataset-IMPE-filtered_SF1_m4lag1scale50_all.csv", check.names = FALSE)
metadata<-read.csv("dataset-sensors-synchronisation-dec2022.csv")
metadata<-select(metadata, c("ID", "Preterm"))

df<-merge(metadata, SF1[-1], by="ID")

df_long<-pivot_longer(df, cols=c(5:54), names_to="Scale", values_to="Entropy")%>%
  mutate(ID=as.factor(ID),
         Scale=as.numeric(Scale),
         Entropy=Entropy/log2(4*3*2*1),
         Preterm=factor(Preterm, levels=c("0", "1"), labels=c("Term", "Preterm")),
         Sensor=factor(Sensor, levels=c("Ankle-Left", "Ankle-Right", "Wrist-Left", "Wrist-Right", "Torso")),
         exclude=ifelse((ID==XXXX & Sensor=="Wrist-Left") | (ID==XXXX & Sensor=="Wrist-Right"), 1, 0)) %>%
  filter(exclude==0)


#remove XXXX Wrist L and Wrist R data, 
#data reduction to timescales of interest
timescales=c(1, 4, 8, 13, 25, 50)
df_reduced<-filter(df_long, Scale %in% timescales)%>% 
  select(-exclude)%>%
  mutate(Scale_f=as.factor(Scale))

df_subset<-filter(df_long, ID==XXXX | ID==XXXX | ID==XXXX | ID== XXXX  | ID==XXXX | ID==XXXX | ID==XXXX | ID == XXXX | ID==XXXX)

#plot data, show that the curve retains its global shape even after data reduction
ggplot(df_long, aes(Scale, Entropy, colour=Preterm)) + geom_path(aes(group=ID), alpha=0.3) + facet_wrap(~Sensor) +
  stat_summary(fun=mean, geom="line") +
  theme_bw(base_size = 16)

ggplot(df_reduced, aes(Scale, Entropy, colour=Preterm)) + 
  geom_path(aes(group=ID), alpha=0.2) + facet_wrap(~Sensor) +
  stat_summary(fun=mean, geom="line") +
  theme_bw(base_size = 16)+
  facet_wrap(~Sensor, 
             nrow = 3,
             ncol = 2) 

ggplot(df_long, aes(Scale, Entropy, colour=ID)) + 
  geom_point(alpha=0.3)+
  geom_path(aes(group=ID), alpha=0.3) + facet_wrap(~Sensor)+
  theme_bw(base_size=18) +
  facet_wrap(~Sensor, 
             nrow = 3,
             ncol = 2) +
  ylab("Permutation Entropy (bit)")+
  xlab("Scale factor")

#build model
library(lme4)
library(sjPlot)
library(lmerTest)
set.seed(060921)
M1<-lmer(Entropy~Scale_f*Preterm*Sensor + (1|ID), df_reduced, REML=FALSE)
M2<-lmer(Entropy~Scale_f*Preterm + Scale_f*Sensor + Sensor*Preterm + (1|ID), df_reduced, REML=FALSE)
M3<-lmer(Entropy~Scale_f*Preterm + Scale_f*Sensor + (1|ID), df_reduced, REML=FALSE)
M4<-lmer(Entropy~Scale_f*Preterm + Sensor + (1|ID), df_reduced, REML=FALSE)

M1_M2<-anova(M1, M2)%>%broom::tidy()
M2_M3<-anova(M2, M3)%>%broom::tidy()
M3_M4<-anova(M3, M4)%>%broom::tidy()

M_modcompare<-rbind(M1_M2, M2_M3, M3_M4)

#obtain anova table
M<-lmer(Entropy~Scale_f*Preterm + Scale_f*Sensor + Preterm*Scale_f*Sensor + (1|ID), df_reduced, REML=TRUE)
M_anova<-anova(M)


#obtain estimated marginal effects and conduct contrasts
#useful websites:https://cran.r-project.org/web/packages/emmeans/vignettes/interactions.html
#https://cran.r-project.org/web/packages/emmeans/vignettes/sophisticated.html
#useful for lmer + how to obtain tidy output: https://benwhalley.github.io/just-enough-r/contrasts-lmer.html
library(emmeans)
emm.M<-emmeans(M, "Preterm", by=c("Scale_f"), lmer.df="satterthwaite", weights="equal")
emm.contrasts<-pairs(emm.M, adjust="none")
emm.contrasts.df <-
  emm.contrasts %>%
  broom::tidy()


emm.all<-emmeans(M, "Preterm", by=c("Scale_f", "Sensor"), lmer.df="satterthwaite", weights="equal")

emm.df <-
  emm.all %>%
  broom::tidy()%>%
  mutate(Scale=as.numeric(Scale_f),
         ci.low=estimate-1.96*std.error, 
         ci.high=estimate+1.96*std.error,
         Sensor=factor(Sensor, levels=c("Ankle-Left", "Ankle-Right", "Wrist-Left", "Wrist-Right", "Torso")))

emm.contrasts.df<-pairs(emm.all, adjust="none")%>%
  broom::tidy()%>%
  mutate(ci.low=estimate-1.96*std.error, ci.high=estimate+1.96*std.error)

#plot estimated marginal effects with raw data 1200 x 750, portrait: 850 x 1000
#https://stats.idre.ucla.edu/r/faq/how-can-i-calculate-standard-errors-for-variance-components-from-mixed-models/
ggplot() + geom_point(df_long, mapping=aes(Scale, Entropy, colour=Preterm, group=ID), size=0.7, alpha=0.1) + 
  geom_line(emm.df, mapping=aes(Scale, estimate, colour=Preterm)) + 
  geom_point(emm.df, mapping=aes(Scale, estimate, colour=Preterm)) + 
  geom_errorbar(emm.df, mapping=aes(Scale, ymin=ci.low, ymax=ci.high, colour=Preterm), alpha=0.7)+
  facet_wrap(~Sensor, 
             nrow = 3,
             ncol = 2) +
  theme_bw(base_size=18) +
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  ylab("Permutation Entropy (bit)")+
  xlab("Scale factor")+
  labs(colour="Group")

#save output
setwd("XXXX")
tab_model(M, show.ci=0.95, df.method="satterthwaite", p.style ="numeric_stars", collapse.ci = TRUE, show.aic=TRUE)
#write.csv(emm.contrasts.df, "contrasts-MSE-PT.csv")

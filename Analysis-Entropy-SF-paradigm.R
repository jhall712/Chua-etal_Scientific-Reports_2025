#SF phases
remove(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)

#import data
setwd("XXXXX")

P<-read.csv("dataset-IMPE-filtered_P_m4lag1scale50_all.csv", check.names = FALSE)
SF1<-read.csv("dataset-IMPE-filtered_SF1_m4lag1scale50_all.csv", check.names = FALSE)
R1<-read.csv("dataset-IMPE-filtered_R1_m4lag1scale50_all.csv", check.names = FALSE)
SF2<-read.csv("dataset-IMPE-filtered_SF2_m4lag1scale50_all.csv", check.names = FALSE)
R2<-read.csv("dataset-IMPE-filtered_R2_m4lag1scale50_all.csv", check.names = FALSE)

metadata<-read.csv("dataset-sensors-synchronisation-dec2022.csv")
metadata<-select(metadata, c("ID", "Preterm"))

exclude.P<-c(XXXX, XXXX)
exclude.R1<-c(XXXX, XXXX)
exclude.WLWR<-c(XXXX)

#Merge data with preterm information
SF<-rbind(P, SF1, R1, SF2, R2)
SF<-select(SF, -1)
df<-merge(metadata, SF, by=c("ID"))

df<-df%>%mutate(
  exclude=ifelse((ID==XXXX & Sensor=="Wrist-Left") | (ID==XXXX & Sensor=="Wrist-Right"), 1,
                 ifelse(ID %in% exclude.P & Phase=="P", 1, 
                        ifelse(ID %in% exclude.R1 & Phase=="R1", 1, 0)))) %>%
  subset(exclude==0)%>%
  select(-exclude)

library(tidyr) 
#Compute complexity index (all) and at each frequency band
df_CI_bands<-pivot_longer(df, cols=c(5:54), names_to="Scale", values_to="Entropy")%>%
  mutate(
    Scale=as.numeric(Scale), 
    freq_band=ifelse(Scale>=1 & Scale<4, "Gamma", 
                     ifelse(Scale>=4 & Scale<8, "Beta",
                            ifelse(Scale>=8 & Scale<13, "Alpha", 
                                   ifelse(Scale>=13 & Scale<25, "Theta", "Delta")))))%>%
  group_by(ID, Preterm, Phase, Sensor, freq_band)%>%
  summarise(CI=sum(Entropy))

df_CI_bands<-df_CI_bands%>%
  mutate(freq_band=factor(freq_band, levels=c("Gamma", "Beta", "Alpha", "Theta", "Delta"),
                                     labels=c("Gamma (30.5-45Hz)", 
                                              "Beta (14-30Hz)", 
                                              "Alpha (8-13.5Hz)", 
                                              "Theta (4.5-7.5Hz)", 
                                              "Delta (0.5-4Hz)")),
         Phase=factor(Phase, levels=c("P", "SF1", "R1", "SF2", "R2"), labels=c("Play", "SF1", "R1", "SF2", "R2"), ordered=FALSE),
         Preterm=factor(Preterm, levels=c(0, 1), labels=c("Term", "Preterm")),
         Sensor=factor(Sensor, levels=c("Torso", "Wrist-Left", "Wrist-Right", "Ankle-Left", "Ankle-Right")),
         ID=factor(ID))
  
df_CI<-df_CI_bands%>%
  group_by(ID, Preterm, Phase, Sensor)%>%
  summarise(CI=sum(CI))%>%
  mutate(Phase=factor(Phase, levels=c("Play", "SF1", "R1", "SF2", "R2"),labels=c("Play", "SF1", "R1", "SF2", "R2"), ordered=FALSE),
         Preterm=factor(Preterm, levels=c("Term", "Preterm"), labels=c("Term", "Preterm")),
         Sensor=factor(Sensor, levels=c("Torso", "Wrist-Left", "Wrist-Right", "Ankle-Left", "Ankle-Right")),
         ID=factor(ID))

#plot data, 
ggplot(df_CI_bands, aes(Phase, CI, colour=freq_band)) + 
  geom_point(aes(group=interaction(ID, freq_band)), alpha=0.3)+
  geom_line(aes(group=interaction(ID, freq_band)), alpha=0.3) + facet_wrap(~Sensor) +
  stat_summary(aes(group=freq_band), fun=mean, geom="line") +
  theme_bw(base_size = 16)

df_gamma<-df_CI_bands%>%filter(freq_band=="Gamma (30.5-45Hz)")
ggplot(df_gamma, aes(Phase, CI, colour=Preterm)) + 
  geom_point(aes(group=ID), alpha=0.3)+
  geom_line(aes(group=ID), alpha=0.3) + facet_wrap(~Sensor) +
  stat_summary(aes(group=Preterm, colour=Preterm), fun=mean, geom="line") +
  theme_bw(base_size = 16)

df_beta<-df_CI_bands%>%filter(freq_band=="Beta (14-30Hz)")
df_beta%>%
  ggplot(aes(Phase, CI, colour=Preterm)) + 
  geom_point(aes(group=ID), alpha=0.3)+
  geom_line(aes(group=ID), alpha=0.3) + facet_wrap(~Sensor) +
  stat_summary(aes(group=Preterm, colour=Preterm), fun=mean, geom="line") +
  theme_bw(base_size = 16)

df_alpha<-df_CI_bands%>%filter(freq_band=="Alpha (8-13.5Hz)")
df_alpha%>%
  ggplot(aes(Phase, CI, colour=Preterm)) + 
  geom_point(aes(group=ID), alpha=0.3)+
  geom_line(aes(group=ID), alpha=0.3) + facet_wrap(~Sensor) +
  stat_summary(aes(group=Preterm, colour=Preterm), fun=mean, geom="line") +
  theme_bw(base_size = 16)


df_theta<-df_CI_bands%>%filter(freq_band=="Theta (4.5-7.5Hz)")
df_theta%>%
  ggplot(aes(Phase, CI, colour=Preterm)) + 
  geom_point(aes(group=ID), alpha=0.3)+
  geom_line(aes(group=ID), alpha=0.3) + facet_wrap(~Sensor) +
  stat_summary(aes(group=Preterm, colour=Preterm), fun=mean, geom="line") +
  theme_bw(base_size = 16)

df_delta<-df_CI_bands%>%filter(freq_band=="Delta (0.5-4Hz)")
df_delta%>%
  ggplot(aes(Phase, CI, colour=Preterm)) + 
  geom_point(aes(group=ID), alpha=0.3)+
  geom_line(aes(group=ID), alpha=0.3) + facet_wrap(~Sensor) +
  stat_summary(aes(group=Preterm, colour=Preterm), fun=mean, geom="line") +
  theme_bw(base_size = 16)

ggplot(df_CI, aes(Phase, CI, colour=ID)) + 
  geom_point(alpha=0.3)+
  geom_path(aes(group=ID), alpha=0.3) + facet_wrap(~Sensor) +
  stat_summary(aes(group=Preterm, colour=Preterm), fun=mean, geom="line") +
  theme_bw(base_size = 16)

#Fit mixed models
library(lme4)
library(lmerTest)

#CI overall
set.seed(060921)
M1<-lmer(CI~Phase*Preterm*Sensor + (1|ID), df_CI)
M2<-lmer(CI~Phase*Preterm + Phase*Sensor + Preterm*Sensor + (1|ID), df_CI)
M3<-lmer(CI~Phase*Preterm + Preterm*Sensor + (1|ID), df_CI)
M4<-lmer(CI~Sensor + Phase*Preterm + (1|ID), df_CI)

M1_M2<-anova(M1, M2)%>%broom::tidy()
M2_M3<-anova(M2, M3)%>%broom::tidy()
M3_M4<-anova(M3, M4)%>%broom::tidy()

M_modcompare<-rbind(M1_M2, M2_M3, M3_M4)
#setwd("XXXX")
#write.csv(M_modcompare, file="model-building-CI.csv")

#View final model M3
anova(M3)
summary(M3)

#CI by frequency band, model building
#gamma
set.seed(060921)
C1<-lmer(CI~Phase*Preterm*Sensor + (1|ID), df_gamma)
C2<-lmer(CI~Phase*Preterm + Phase*Sensor + Preterm*Sensor + (1|ID), df_gamma)
C3<-lmer(CI~Phase*Preterm + Preterm*Sensor + (1|ID), df_gamma)
C4<-lmer(CI~Sensor + Phase*Preterm + (1|ID), df_gamma)

C1_C2<-anova(C1, C2)%>%broom::tidy()
C2_C3<-anova(C2, C3)%>%broom::tidy()
C3_C4<-anova(C3, C4)%>%broom::tidy()

C_modcompare<-rbind(C1_C2, C2_C3, C3_C4)
#write.csv(C_modcompare, file="model-building-C.csv")
anova(C3)

#beta
set.seed(060921)
B1<-lmer(CI~Phase*Preterm*Sensor + (1|ID), df_beta)
B2<-lmer(CI~Phase*Preterm + Phase*Sensor + Preterm*Sensor + (1|ID), df_beta)
B3<-lmer(CI~Phase*Preterm + Preterm*Sensor + (1|ID), df_beta)
B4<-lmer(CI~Sensor + Phase*Preterm + (1|ID), df_beta)

B1_B2<-anova(B1, B2)%>%broom::tidy()
B2_B3<-anova(B2, B3)%>%broom::tidy()
B3_B4<-anova(B3, B4)%>%broom::tidy()

B_modcompare<-rbind(B1_B2, B2_B3, B3_B4)
anova(B3)
#write.csv(B_modcompare, file="model-building-B.csv")

#alpha
set.seed(060921)
A1<-lmer(CI~Phase*Preterm*Sensor + (1|ID), df_alpha)
A2<-lmer(CI~Phase*Preterm + Phase*Sensor + Preterm*Sensor + (1|ID), df_alpha)
A3<-lmer(CI~Phase*Preterm + Preterm*Sensor + (1|ID), df_alpha)
A4<-lmer(CI~Sensor + Phase*Preterm + (1|ID), df_alpha)

A1_A2<-anova(A1, A2)%>%broom::tidy()
A2_A3<-anova(A2, A3)%>%broom::tidy()
A3_A4<-anova(A3, A4)%>%broom::tidy()

A_modcompare<-rbind(A1_A2, A2_A3, A3_A4)
#write.csv(A_modcompare, file="model-building-A.csv")
anova(A3)

#delta
set.seed(060921)
D1<-lmer(CI~Phase*Preterm*Sensor + (1|ID), df_delta)
D2<-lmer(CI~Phase*Preterm + Phase*Sensor + Preterm*Sensor + (1|ID), df_delta)
D3<-lmer(CI~Phase*Preterm + Preterm*Sensor + (1|ID), df_delta)
D4<-lmer(CI~Sensor + Phase*Preterm + (1|ID), df_delta)

D1_D2<-anova(D1, D2)%>%broom::tidy()
D2_D3<-anova(D2, D3)%>%broom::tidy()
D3_D4<-anova(D3, D4)%>%broom::tidy()

D_modcompare<-rbind(D1_D2, D2_D3, D3_D4)
#write.csv(D_modcompare, file="model-building-D.csv")
anova(D3)

#theta
set.seed(060921)
E1<-lmer(CI~Phase*Preterm*Sensor + (1|ID), df_theta)
E2<-lmer(CI~Phase*Preterm + Phase*Sensor + Preterm*Sensor + (1|ID), df_theta)
E3<-lmer(CI~Phase*Preterm + Preterm*Sensor + (1|ID), df_theta)
E4<-lmer(CI~Sensor + Phase*Preterm + (1|ID), df_theta)

E1_E2<-anova(E1, E2)%>%broom::tidy()
E2_E3<-anova(E2, E3)%>%broom::tidy()
E3_E4<-anova(E3, E4)%>%broom::tidy()

E_modcompare<-rbind(E1_E2, E2_E3, E3_E4)
#write.csv(E_modcompare, file="model-building-E.csv")
anova(E3)

#Compute estimated marginal means and run contrasts
library(emmeans)
#Effect of phase at high frequency bands
skip_comp.emmc <- function(levels, skip = 1, reverse = FALSE) {
  if((k <- length(levels)) < skip + 1)
    stop("Need at least ", skip + 1, " levels")
  coef <- data.frame()
  coef <- as.data.frame(lapply(seq_len(k - skip - 1), function(i) {
    sgn <- ifelse(reverse, -1, 1)
    sgn * c(rep(0, i - 1), 1, rep(0, skip), -1, rep(0, k - i - skip - 1))
  }))
  names(coef) <- sapply(coef, function(x)
    paste(which(x == 1), "-", which(x == -1)))
  attr(coef, "adjust") = "none"   # fdr is default adjustment method
  coef
}

skip_comp.emmc(1:5, skip=0)

#Contrast to identify effect of Phase at gamma, alpha, beta frequency bands only
emm.gamma.P<-emmeans(C3, "Phase", lmer.df="satterthwaite", weights="equal")
contrasts.gamma.P<-contrast(emm.gamma.P, "skip_comp", skip=0, reverse=TRUE, adjust=NULL)%>%
  broom::tidy()%>%
  mutate(p.adj=p.value*5)

emm.beta.P<-emmeans(B3, "Phase", lmer.df="satterthwaite", weights="equal")
contrasts.beta.P<-contrast(emm.beta.P, "skip_comp", skip=0, reverse=TRUE, adjust=NULL) %>%
  broom::tidy()%>%
  mutate(p.adj=p.value*5)

emm.alpha.P<-emmeans(A3, "Phase", lmer.df="satterthwaite", weights="equal")
contrasts.alpha.P<-contrast(emm.alpha.P, "skip_comp", skip=0, reverse=TRUE, adjust=NULL) %>%
  broom::tidy()%>%
  mutate(p.adj=p.value*5)

#emm.delta.P<-emmeans(D3, "Phase", by=c("Sensor"), lmer.df="satterthwaite", weights="equal")
#contrasts.delta.P<-contrast(emm.delta.P, "skip_comp", skip=0, reverse=TRUE, adjust=NULL) %>%
#  broom::tidy()%>%
#  mutate(p.adj=p.value*5)

#emm.theta.P<-emmeans(E3, "Phase", by=c("Sensor"), lmer.df="satterthwaite", weights="equal")
#contrasts.theta.P<-contrast(emm.theta.P, "skip_comp", skip=0, reverse=TRUE, adjust=NULL) %>%
#  broom::tidy()%>%
#  mutate(p.adj=p.value*5)

#Contrast to identify effect of Preterm or Preterm x sensor interaction after adjustment, at all frequency bands except alpha
emm.gamma.PT<-emmeans(C3, "Preterm", by=c("Sensor"), lmer.df="satterthwaite", weights="equal")
contrasts.gamma.PT<-pairs(emm.gamma.PT, adjust="none")%>%
  broom::tidy()%>%
  mutate(p.adj=p.value*5)

emm.beta.PT<-emmeans(B3, "Preterm", by=c("Sensor"), lmer.df="satterthwaite", weights="equal")
contrasts.beta.PT<-pairs(emm.beta.PT, adjust="none")%>%
  broom::tidy()%>%
  mutate(p.adj=p.value*5)

#emm.alpha.PT<-emmeans(A3, "Preterm", by=c("Sensor"), lmer.df="satterthwaite", weights="equal")
#contrasts.alpha.PT<-pairs(emm.alpha.PT, adjust="none")%>%
#  broom::tidy()%>%
#  mutate(p.adj=p.value*5)

emm.theta.PT<-emmeans(E3, "Preterm", by=c("Sensor"), lmer.df="satterthwaite", weights="equal")
contrasts.theta.PT<-pairs(emm.theta.PT, adjust="none")%>%
  broom::tidy()%>%
  mutate(p.adj=p.value*5)

emm.delta.PT<-emmeans(D3, "Preterm", by=c("Sensor"), lmer.df="satterthwaite", weights="equal")
contrasts.delta.PT<-pairs(emm.delta.PT, adjust="none")%>%
  broom::tidy()%>%
  mutate(p.adj=p.value*5)

emm.CI.PT<-emmeans(M3, "Preterm", by=c("Sensor"), lmer.df="satterthwaite", weights="equal")
contrasts.CI.PT<-pairs(emm.CI.PT, adjust="none")%>%
  broom::tidy()


#Plot models
emm.CI.plot<-emmeans(M3, "Phase", by=c("Preterm", "Sensor"), lmer.df="satterthwaite", weights="equal")%>%
  broom::tidy()%>%
  mutate(ci.low=estimate-std.error*1.96,
         ci.high=estimate+std.error*1.96,
         Sensor=factor(Sensor, levels=c("Ankle-Left", "Ankle-Right", "Wrist-Left", "Wrist-Right", "Torso")))

emm.gamma.plot<-emmeans(C3, "Phase", by=c("Preterm", "Sensor"), lmer.df="satterthwaite", weights="equal")%>%
  broom::tidy()%>%
  mutate(ci.low=estimate-std.error*1.96,
         ci.high=estimate+std.error*1.96,
         Sensor=factor(Sensor, levels=c("Ankle-Left", "Ankle-Right", "Wrist-Left", "Wrist-Right", "Torso")))
emm.beta.plot<-emmeans(B3, "Phase", by=c("Preterm", "Sensor"), lmer.df="satterthwaite", weights="equal")%>%
  broom::tidy()%>%
  mutate(ci.low=estimate-std.error*1.96,
         ci.high=estimate+std.error*1.96,
         Sensor=factor(Sensor, levels=c("Ankle-Left", "Ankle-Right", "Wrist-Left", "Wrist-Right", "Torso")))
emm.alpha.plot<-emmeans(A3, "Phase", by=c("Preterm", "Sensor"), lmer.df="satterthwaite", weights="equal")%>%
  broom::tidy()%>%
  mutate(ci.low=estimate-std.error*1.96,
         ci.high=estimate+std.error*1.96,
         Sensor=factor(Sensor, levels=c("Ankle-Left", "Ankle-Right", "Wrist-Left", "Wrist-Right", "Torso")))
emm.delta.plot<-emmeans(D3, "Phase", by=c("Preterm", "Sensor"), lmer.df="satterthwaite", weights="equal")%>%
  broom::tidy()%>%
  mutate(ci.low=estimate-std.error*1.96,
         ci.high=estimate+std.error*1.96,
         Sensor=factor(Sensor, levels=c("Ankle-Left", "Ankle-Right", "Wrist-Left", "Wrist-Right", "Torso")))
emm.theta.plot<-emmeans(E3, "Phase", by=c("Preterm", "Sensor"), lmer.df="satterthwaite", weights="equal")%>%
  broom::tidy()%>%
  mutate(ci.low=estimate-std.error*1.96,
         ci.high=estimate+std.error*1.96,
         Sensor=factor(Sensor, levels=c("Ankle-Left", "Ankle-Right", "Wrist-Left", "Wrist-Right", "Torso")))



#plot estimated marginal effects with raw data
df_CI<-df_CI%>%mutate(Sensor=factor(Sensor, levels=c("Ankle-Left", "Ankle-Right", "Wrist-Left", "Wrist-Right", "Torso")))
ggplot() + geom_point(df_CI, mapping=aes(Phase, CI, colour=Preterm), size=1.5, alpha=0.2) + 
  geom_line(emm.CI.plot, mapping=aes(Phase, estimate, colour=Preterm, group=Preterm)) + 
  geom_point(emm.CI.plot, mapping=aes(Phase, estimate, colour=Preterm)) + 
  geom_errorbar(emm.CI.plot, mapping=aes(Phase, ymin=ci.low, ymax=ci.high, colour=Preterm), width=0.2, alpha=0.7)+
  facet_wrap(~Sensor, 
             nrow = 3,
             ncol = 2)+
  theme_bw(base_size=18) +
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")

df_gamma<-df_gamma%>%mutate(Sensor=factor(Sensor, levels=c("Ankle-Left", "Ankle-Right", "Wrist-Left", "Wrist-Right", "Torso")))
ggplot() + geom_point(df_gamma, mapping=aes(Phase, CI, colour=Preterm), size=1.5, alpha=0.2) + 
  geom_line(emm.gamma.plot, mapping=aes(Phase, estimate, colour=Preterm, group=Preterm)) + 
  geom_point(emm.gamma.plot, mapping=aes(Phase, estimate, colour=Preterm)) + 
  geom_errorbar(emm.gamma.plot, mapping=aes(Phase, ymin=ci.low, ymax=ci.high, colour=Preterm), width=0.2, alpha=0.7)+
  facet_wrap(~Sensor, 
             nrow = 3,
             ncol = 2)+
  theme_bw(base_size=18) +
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")



df_beta<-df_beta%>%mutate(Sensor=factor(Sensor, levels=c("Ankle-Left", "Ankle-Right", "Wrist-Left", "Wrist-Right", "Torso")))
ggplot() + geom_point(df_beta, mapping=aes(Phase, CI, colour=Preterm), size=1.5, alpha=0.2) + 
  geom_line(emm.beta.plot, mapping=aes(Phase, estimate, colour=Preterm, group=Preterm)) + 
  geom_point(emm.beta.plot, mapping=aes(Phase, estimate, colour=Preterm)) + 
  geom_errorbar(emm.beta.plot, mapping=aes(Phase, ymin=ci.low, ymax=ci.high, colour=Preterm), width=0.2, alpha=0.7)+
  facet_wrap(~Sensor, 
             nrow = 3,
             ncol = 2)+
  theme_bw(base_size=18) +
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")

df_alpha<-df_alpha%>%mutate(Sensor=factor(Sensor, levels=c("Ankle-Left", "Ankle-Right", "Wrist-Left", "Wrist-Right", "Torso")))
ggplot() + geom_point(df_alpha, mapping=aes(Phase, CI, colour=Preterm), size=1.5, alpha=0.2) + 
  geom_line(emm.alpha.plot, mapping=aes(Phase, estimate, colour=Preterm, group=Preterm)) + 
  geom_point(emm.alpha.plot, mapping=aes(Phase, estimate, colour=Preterm)) + 
  geom_errorbar(emm.alpha.plot, mapping=aes(Phase, ymin=ci.low, ymax=ci.high, colour=Preterm), width=0.2, alpha=0.7)+
  facet_wrap(~Sensor, 
             nrow = 3,
             ncol = 2)+
  theme_bw(base_size=18) +
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")

df_delta<-df_delta%>%mutate(Sensor=factor(Sensor, levels=c("Ankle-Left", "Ankle-Right", "Wrist-Left", "Wrist-Right", "Torso")))
ggplot() + geom_point(df_delta, mapping=aes(Phase, CI, colour=Preterm), size=1.5, alpha=0.2) + 
  geom_line(emm.delta.plot, mapping=aes(Phase, estimate, colour=Preterm, group=Preterm)) + 
  geom_point(emm.delta.plot, mapping=aes(Phase, estimate, colour=Preterm)) + 
  geom_errorbar(emm.delta.plot, mapping=aes(Phase, ymin=ci.low, ymax=ci.high, colour=Preterm), width=0.2, alpha=0.7)+
  facet_wrap(~Sensor, 
             nrow = 3,
             ncol = 2)+
  theme_bw(base_size=18) +
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")

df_theta<-df_theta%>%mutate(
  Sensor=factor(Sensor, levels=c("Ankle-Left", "Ankle-Right", "Wrist-Left", "Wrist-Right", "Torso")))
ggplot() + geom_point(df_theta, mapping=aes(Phase, CI, colour=Preterm), size=1.5, alpha=0.2) + 
  geom_line(emm.theta.plot, mapping=aes(Phase, estimate, colour=Preterm, group=Preterm)) + 
  geom_point(emm.theta.plot, mapping=aes(Phase, estimate, colour=Preterm)) + 
  geom_errorbar(emm.theta.plot, mapping=aes(Phase, ymin=ci.low, ymax=ci.high, colour=Preterm), width=0.2, alpha=0.7)+
  facet_wrap(~Sensor, 
             nrow = 3,
             ncol = 2)+
  theme_bw(base_size=18) +
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")

tab_model(C3, A3, B3, D3, E3, show.ci=0.95, df.method="satterthwaite", p.style ="numeric_stars", collapse.ci = TRUE, show.aic=TRUE)

#write.csv(emm.alpha.plot, "EMM-alpha.csv")
#write.csv(emm.beta.plot, "EMM-beta.csv")
#write.csv(emm.gamma.plot, "EMM-gamma.csv")
#write.csv(emm.delta.plot, "EMM-delta.csv")
#write.csv(emm.theta.plot, "EMM-theta.csv")

#write.csv(contrasts.alpha.P, "contrasts-Phase-alpha.csv")
#write.csv(contrasts.beta.P, "contrasts-Phase-beta.csv")
#write.csv(contrasts.gamma.P, "contrasts-Phase-gamma.csv")

#write.csv(contrasts.beta.PT, "contrasts-Preterm-beta.csv")
#write.csv(contrasts.gamma.PT, "contrasts-Preterm-gamma.csv")
#write.csv(contrasts.delta.PT, "contrasts-Preterm-delta.csv")
#write.csv(contrasts.theta.PT, "contrasts-Preterm-theta.csv")
#write.csv(contrasts.CI.PT, "contrasts-Preterm-CI.csv")

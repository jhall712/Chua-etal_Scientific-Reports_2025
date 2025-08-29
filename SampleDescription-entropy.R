remove(list=ls())

library(dplyr)
library(tidyr)
###############################################################################   
##### PREPARE DATASETS FOR ANALYSIS
#Load datasets
setwd("XXXXX")
df_anth<-read.csv("dataset-anthropometry-9m.csv")
df_demog<-read.csv("dataset-demographics.csv")
df_age<-read.csv("dataset-age-at-visit.csv")
df_motor<-read.csv("dataset-VABS-motor.csv")

setwd("XXXXX")
df_entropy<-read.csv("dataset-IMPE-filtered_SF1_m4lag1scale50_all.csv")

#obtain list of IDs with entropy data
id<-df_entropy%>%filter(Sensor=="Torso")%>%select(ID)

#obtain selection of each dataset to be merged
A<-df_anth%>%select(c(record_id, height_cm, weight_kg))
B<-df_age%>%select(c(record_id, age_9m_m, age_corrected_9m_m))
C<-df_motor%>%select(c(ID, sits_unsupported, MOT, MOT_gross, MOT_fine))

#merge datasets for analysis
df_sample<-merge(id, df_demog, by.x=c("ID"), by.y="record_id")
df_sample<-merge(df_sample, A, by.x=c("ID"), by.y="record_id")
df_sample<-merge(df_sample, B, by.x=c("ID"), by.y="record_id")
df_sample<-merge(df_sample, C, by=c("ID"))


df_sample<-
  df_sample%>% 
  mutate(ID=as.factor(ID),
         Sex=factor(sex, levels=c(1, 2), labels=c("Male", "Female")),
         Preterm=factor(preterm, levels=c(0,1), labels=c("Term", "Preterm")),
         Ethnicity=ifelse(child_s_ethnicity==1, "Any White background", 
                          ifelse(child_s_ethnicity<=5, "Any Mixed background", 
                                 ifelse(child_s_ethnicity<=9, "Any Asian background",
                                        ifelse(is.na(child_s_ethnicity)==TRUE, NA, 
                                               "Any other ethnic group")))),
         age_or_corrage=ifelse(Preterm=="Term", age_9m_m, age_corrected_9m_m)) %>%
  select(-X)

#save datasets
setwd("XXXXX")
#write.csv(df_sample, file="sensor-subsample-description-final-with-MOT.csv")

################### Obtain descriptive statistics
remove(list=ls())

library(dplyr)
library(tidyr)

#Load datasets
setwd("XXXXX")

df_sample=read.csv("sensor-subsample-description-final-with-MOT.csv")

#data wrangling
df_sample<-df_sample%>% mutate(ID=as.factor(ID),
                               Sex=factor(sex, levels=c(1, 2), labels=c("Male", "Female")),
                               Preterm=factor(preterm, levels=c(0,1), labels=c("Term", "Preterm")),
                               Ethnicity=ifelse(child_s_ethnicity==1, "Any White background", 
                                                ifelse(child_s_ethnicity<=5, "Any Mixed background", 
                                                       ifelse(child_s_ethnicity<=9, "Any Asian background",
                                                              ifelse(is.na(child_s_ethnicity)==TRUE, NA, 
                                                                     "Any other ethnic group")))),
                               age_or_corrage=ifelse(Preterm=="Term", age_9m_m, age_corrected_9m_m),
                               Singleton=ifelse(no_of_infants==1, "Singleton", "Twin"),
                               Singleton=factor(Singleton, levels=c("Singleton", "Twin"), labels=c("Singleton", "Twin")),
                               Sits_unsupported=factor(sits_unsupported, levels=c(0,1), labels=c("sits supported", "sits unsupported")))



df_sample_sf<-merge(df_sample_sf, df_sample[c("ID", "Preterm")], by.x="ID", by.y="ID")%>%
  select(-X)%>%
  mutate(SF_phase=factor(SF_phase, levels=c("Play", "SF1", "Reunion1", "SF2", "Reunion2"), labels=c("Play", "SF1", "R1", "SF2", "R2")))

#exclude the infant who took off wrist sensor
df_sample_acc=read.csv("dataset-sensors-acceleration-SF-all.csv")
df_sample_acc<-merge(df_sample_acc, df_sample[c("ID", "Preterm")], by.x="ID", by.y="ID")%>%
  mutate(SF_phase=factor(SF_phase, levels=c("P", "SF1", "R1", "SF2", "R2"), labels=c("Play", "SF1", "R1", "SF2", "R2")))%>%
  mutate(Sensor=factor(Sensor, levels=c("Torso", "Wrist-Left", "Wrist-Right", "Ankle-Left", "Ankle-Right")))

df_sample_acc<-df_sample_acc%>%
  mutate(exclude=ifelse((ID==XXXX & Sensor=="Wrist-Left") | (ID==XXXX & Sensor=="Wrist-Right") | (ID==XXXX & SF_phase=="R1"), 1, 0))%>%
  subset(exclude==0)%>%
  select(-exclude)

#convert to wide dataframe with mean acceleration for all 5 sensor locations
df_acc<-df_sample_acc%>%select(ID, Acc_mean, Preterm, Sensor, SF_phase) %>% pivot_wider(names_from=Sensor, values_from=Acc_mean)

#merge acceleration dataset with behavioural dataset
df_motor<-merge(df_acc, df_sample_sf, by.x=c("ID", "SF_phase", "Preterm"), by.y=c("ID", "SF_phase", "Preterm"), all=TRUE)

#merge complexity dataset with behavioural dataset
df_CI<-read.csv("dataset-sensors-CI.csv")%>%
  select(-X)%>%
  pivot_wider(id_cols=c("ID", "Phase"), values_from=CI, names_from=Sensor, names_prefix = "CI_")%>%
  mutate(Phase=ifelse(Phase=="P", "Play", Phase))
df_motor<-merge(df_motor, df_CI, by.x=c("ID", "SF_phase"), by.y=c("ID", "Phase"), all=TRUE) #both datasets include same infants

df_motor_SF1<-filter(df_motor, SF_phase=="SF1")
df_motor_SF1<-merge(df_sample, df_motor_SF1, by.x=c("ID", "Preterm"),  by.y=c("ID", "Preterm"))

##############################
#descriptive tables - https://cran.r-project.org/web/packages/gtsummary/vignettes/tbl_summary.html
#https://education.rstudio.com/blog/2020/07/gtsummary/#inline
#install.packages("gtsummary")
library(gtsummary)

df_sample_table<-df_sample%>%select(Preterm, age_or_corrage, birthweight_g, birthweight_z_score, gestation, Sex, Ethnicity, simd_quintile, Singleton, height_cm, weight_kg, Sits_unsupported, MOT, MOT_gross, MOT_fine)

#XXXX weight wrong - change from 81kg to 8.1kg (checked on REDCAP)
#infant characteristics: age at visit, demographics (gender, ethnicity, simd, twin birth), 9m anthropometry
table1<- tbl_summary(df_sample_table, by="Preterm",
            statistic=list(all_continuous()~c("{mean} ({sd})")),
            digits=list(everything()~1))%>%
  add_overall()
library(flextable)

#behavioural still-face data - use figures to describe data
#XXXX and XXXX video coding data nonSF phase incomplete.
term_sample_sf<-df_sample_sf%>%filter(Preterm=="Term")%>%select(-ID, -Preterm)
preterm_sample_sf<-df_sample_sf%>%filter(Preterm=="Preterm")%>%select(-ID, -Preterm)

#table2<- tbl_summary(df_sample_sf, by=c("SF_phase"),
#                      statistic=list(all_continuous()~c("{mean} ({sd})"))) %>%
#                      add_p(all_continuous() ~ "aov")
table2<-df_sample_sf%>%
  select(-ID)%>%
  tbl_strata (strata=SF_phase,
              .tbl_fun = ~.x%>%
                tbl_summary(by=Preterm, missing="no",
                            type = list(c("composite_obj", "ICEP_neg", "caregiver_moves_infant") ~ 'continuous'))
  )

neg<-ggplot(df_sample_sf, aes(SF_phase, ICEP_neg, fill=Preterm)) + geom_boxplot(width=0.5, alpha=0.7) +
  theme_bw(base_size=18) +
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  ylab("Negative Affect")+
  xlab("Phase")+
  labs(fill="Group")

obj<-ggplot(df_sample_sf, aes(SF_phase, composite_obj, fill=Preterm)) + geom_boxplot(width=0.5, alpha=0.7) +
  theme_bw(base_size=18) +
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  ylab("OBJ")+
  xlab("Phase")+
  labs(fill="Group")
                          
rme<-ggplot(df_sample_sf, aes(SF_phase, composite_rme, fill=Preterm)) + geom_boxplot(width=0.5, alpha=0.7) +
  theme_bw(base_size=18) +
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  ylab("RME")+
  xlab("Phase")+
  labs(fill="Group")

sc<-ggplot(df_sample_sf, aes(SF_phase, composite_sc, fill=Preterm)) + geom_boxplot(width=0.5, alpha=0.7) +
  theme_bw(base_size=18) +
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  ylab("SC")+
  xlab("Phase")+
  labs(fill="Group")

soc<-ggplot(df_sample_sf, aes(SF_phase, composite_soc, fill=Preterm)) + geom_boxplot(width=0.5, alpha=0.7) +
  theme_bw(base_size=18) +
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  ylab("SOC")+
  xlab("Phase")+
  labs(fill="Group")

caregiver<-ggplot(df_sample_sf, aes(SF_phase, caregiver_moves_infant, fill=Preterm)) + geom_boxplot(width=0.5, alpha=0.7) +
  theme_bw(base_size=18) +
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  ylab("Caregiver movement")+
  xlab("Phase")+
  labs(fill="Group")

library(ggpubr)
fig_sf_behaviours<-ggarrange(sc, soc, obj, rme, neg, caregiver,
                        ncol=2,
                        nrow=3,
                        align = "v", 
                        common.legend=TRUE,
                        legend="right")


#motor activity
table3<-df_acc%>%
  select("Torso", "Wrist-Left", "Wrist-Right", "Ankle-Left", "Ankle-Right", SF_phase, Preterm)%>%
  tbl_strata (strata=SF_phase,
              .tbl_fun = ~.x%>%
                tbl_summary(by=Preterm, missing="no")
  )

acc_torso<-ggplot(df_acc, aes(SF_phase, Torso, fill=Preterm)) + geom_violin(width=0.5, alpha=0.7) +
  theme_bw(base_size=18) +
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  ylab(expression(paste("Acceleration (m/", s^2, ")", sep="")))+
  xlab("Phase")+
  labs(fill="Group", title="Torso")

acc_ankleL<-ggplot(df_acc, aes(SF_phase, `Ankle-Left`, fill=Preterm)) + geom_violin(width=0.5, alpha=0.7) +
  theme_bw(base_size=18) +
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  ylab(expression(paste("Acceleration (m/", s^2, ")", sep="")))+
  xlab("Phase")+
  labs(fill="Group", title="Ankle-Left")

acc_ankleR<-ggplot(df_acc, aes(SF_phase, `Ankle-Right`, fill=Preterm)) + geom_violin(width=0.5, alpha=0.7) +
  theme_bw(base_size=18) +
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  ylab(expression(paste("Acceleration (m/", s^2, ")", sep="")))+
  xlab("Phase")+
  labs(fill="Group", title="Ankle-Right")


acc_WristL<-ggplot(df_acc, aes(SF_phase, `Wrist-Left`, fill=Preterm)) + geom_violin(width=0.5, alpha=0.7) +
  theme_bw(base_size=18) +
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  ylab(expression(paste("Acceleration (m/", s^2, ")", sep="")))+
  xlab("Phase")+
  labs(fill="Group", title="Wrist-Left")

acc_WristR<-ggplot(df_acc, aes(SF_phase, `Wrist-Right`, fill=Preterm)) + geom_violin(width=0.5, alpha=0.7) +
  theme_bw(base_size=18) +
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  ylab(expression(paste("Acceleration (m/", s^2, ")", sep="")))+
  xlab("Phase")+
  labs(fill="Group", title="Wrist-Right")

library(ggpubr)
fig_acc<-ggarrange(acc_ankleL, acc_ankleR, acc_WristL, acc_WristR, acc_torso,
                   ncol=2,
                   nrow=3,
                   common.legend=TRUE,
                  legend="right")




df_acc_long<-pivot_longer(df_acc, cols=4:8, names_to = "Sensor", values_to="Acceleration")%>%
  mutate(Sensor=factor(Sensor, levels=c("Ankle-Left", "Ankle-Right", "Wrist-Left", "Wrist-Right","Torso")))%>%
  rename(Phase=SF_phase)
library(lme4)
library(lmerTest)

set.seed(060921)
M1<-lmer(Acceleration~Phase*Preterm*Sensor + (1|ID), df_acc_long, REML=FALSE)
M2<-lmer(Acceleration~Phase*Preterm + Phase*Sensor + Preterm*Sensor + (1|ID), df_acc_long, REML=FALSE)
M3<-lmer(Acceleration~Phase*Preterm + Preterm*Sensor + (1|ID), df_acc_long, REML=FALSE)
M4<-lmer(Acceleration~Sensor + Phase*Preterm + (1|ID), df_acc_long, REML=FALSE)

M1_M2<-anova(M1, M2)%>%broom::tidy()
M2_M3<-anova(M2, M3)%>%broom::tidy()
M3_M4<-anova(M3, M4)%>%broom::tidy()

M_modcompare<-rbind(M1_M2, M2_M3, M3_M4)

library(sjPlot)
M<-lmer(Acceleration~Phase*Preterm + Preterm*Sensor + (1|ID), df_acc_long, REML=TRUE)
tab_model(M, show.ci=0.95, df.method="satterthwaite", p.style ="numeric_stars", collapse.ci = TRUE)

library(htmlTable)
htmlTable(anova(M))

setwd("XXXXX")
#model building output
#write.csv(M_modcompare, "model-building-Acc.csv")

library(emmeans)
#emms
emm.acc.plot<-emmeans(M, "Phase", by=c("Preterm", "Sensor"), lmer.df="satterthwaite", weights="equal")%>%
  broom::tidy()%>%
  mutate(ci.low=estimate-1.96*std.error, 
         ci.high=estimate+1.96*std.error/2,
         Sensor=factor(Sensor, levels=c("Ankle-Left", "Ankle-Right", "Wrist-Left", "Wrist-Right","Torso")))

library(ggplot2)
predeff_acc<-ggplot() + geom_point(df_acc_long, mapping=aes(Phase, Acceleration, colour=Preterm), size=1.5, alpha=0.2) + 
  geom_line(emm.acc.plot, mapping=aes(Phase, estimate, colour=Preterm, group=Preterm)) + 
  geom_point(emm.acc.plot, mapping=aes(Phase, estimate, colour=Preterm)) + 
  geom_errorbar(emm.acc.plot, mapping=aes(Phase, ymin=ci.low, ymax=ci.high, colour=Preterm), width=0.2, alpha=0.7)+
  facet_wrap(~Sensor, nrow=3, ncol=2) +
  theme_bw(base_size=18) +
  scale_colour_brewer(palette="Dark2")+
  scale_fill_brewer(palette="Dark2")+
  ylab(expression(paste("Acceleration (m/", s^2, ")", sep="")))  

#Effect of preterm at each sensor location
emm.acc.PT<-emmeans(M, "Preterm", by=c("Sensor"), lmer.df="satterthwaite", weights="equal")
contrasts.acc.PT<-pairs(emm.acc.PT, adjust="none")%>%
  broom::tidy()%>%
  mutate(p.adj=p.value*5)


#Effect of phase
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
  attr(coef, "adjust") = "fdr"   # default adjustment method
  coef
}

skip_comp.emmc(1:5, skip=0)

emm.acc.P<-emmeans(M, "Phase", lmer.df="satterthwaite", weights="equal")
contrasts.acc.P<-contrast(emm.acc.P, "skip_comp", skip=0, reverse=TRUE, adjust=NULL)%>%
  broom::tidy()%>%
  mutate(p.adj=p.value*5)

write.csv(emm.acc.plot, "EMM-Acc.csv")
write.csv(contrasts.acc.P, "contrasts-Phase-Acc.csv")
write.csv(contrasts.acc.PT, "contrasts-PT-Acc.csv")


#correlations
#by each phase

shapiro.test(df_motor_SF1$`Torso`)
shapiro.test(df_motor_SF1$`Ankle-Right`)


T.neg<-cor.test(df_motor_SF1$`Torso`, df_motor_SF1$ICEP_neg, method="spearman")
WR.neg<-cor.test(df_motor_SF1$`Wrist-Right`, df_motor_SF1$ICEP_neg, method="spearman")
WL.neg<-cor.test(df_motor_SF1$`Wrist-Left`, df_motor_SF1$ICEP_neg, method="spearman")
AR.neg<-cor.test(df_motor_SF1$`Ankle-Right`, df_motor_SF1$ICEP_neg, method="spearman")
AL.neg<-cor.test(df_motor_SF1$`Ankle-Left`, df_motor_SF1$ICEP_neg, method="spearman")

T.rme<-cor.test(df_motor_SF1$`Torso`, df_motor_SF1$composite_rme, method="spearman")
WR.rme<-cor.test(df_motor_SF1$`Wrist-Right`, df_motor_SF1$composite_rme, method="spearman")
WL.rme<-cor.test(df_motor_SF1$`Wrist-Left`, df_motor_SF1$composite_rme, method="spearman")
AR.rme<-cor.test(df_motor_SF1$`Ankle-Right`, df_motor_SF1$composite_rme, method="spearman")
AL.rme<-cor.test(df_motor_SF1$`Ankle-Left`, df_motor_SF1$composite_rme, method="spearman")

T.h<-cor.test(df_motor_SF1$`Torso`, df_motor_SF1$height_cm, method="spearman")
WR.h<-cor.test(df_motor_SF1$`Wrist-Right`, df_motor_SF1$height_cm, method="spearman")
WL.h<-cor.test(df_motor_SF1$`Wrist-Left`, df_motor_SF1$height_cm, method="spearman")
AR.h<-cor.test(df_motor_SF1$`Ankle-Right`, df_motor_SF1$height_cm, method="spearman")
AL.h<-cor.test(df_motor_SF1$`Ankle-Left`, df_motor_SF1$height_cm, method="spearman")

T.w<-cor.test(df_motor_SF1$`Torso`, df_motor_SF1$weight_kg, method="spearman")
WR.w<-cor.test(df_motor_SF1$`Wrist-Right`, df_motor_SF1$weight_kg, method="spearman")
WL.w<-cor.test(df_motor_SF1$`Wrist-Left`, df_motor_SF1$weight_kg, method="spearman")
AR.w<-cor.test(df_motor_SF1$`Ankle-Right`, df_motor_SF1$weight_kg, method="spearman")
AL.w<-cor.test(df_motor_SF1$`Ankle-Left`, df_motor_SF1$weight_kg, method="spearman")


T.ci.neg<-cor.test(df_motor_SF1$`CI_Torso`, df_motor_SF1$ICEP_neg, method="spearman")
WR.ci.neg<-cor.test(df_motor_SF1$`CI_Wrist-Right`, df_motor_SF1$ICEP_neg, method="spearman")
WL.ci.neg<-cor.test(df_motor_SF1$`CI_Wrist-Left`, df_motor_SF1$ICEP_neg, method="spearman")
AR.ci.neg<-cor.test(df_motor_SF1$`CI_Ankle-Right`, df_motor_SF1$ICEP_neg, method="spearman")
AL.ci.neg<-cor.test(df_motor_SF1$`CI_Ankle-Left`, df_motor_SF1$ICEP_neg, method="spearman")

T.ci.rme<-cor.test(df_motor_SF1$`CI_Torso`, df_motor_SF1$composite_rme, method="spearman")
WR.ci.rme<-cor.test(df_motor_SF1$`CI_Wrist-Right`, df_motor_SF1$composite_rme, method="spearman")
WL.ci.rme<-cor.test(df_motor_SF1$`CI_Wrist-Left`, df_motor_SF1$composite_rme, method="spearman")
AR.ci.rme<-cor.test(df_motor_SF1$`CI_Ankle-Right`, df_motor_SF1$composite_rme, method="spearman")
AL.ci.rme<-cor.test(df_motor_SF1$`CI_Ankle-Left`, df_motor_SF1$composite_rme, method="spearman")

T.ci.h<-cor.test(df_motor_SF1$`CI_Torso`, df_motor_SF1$height_cm, method="spearman")
WR.ci.h<-cor.test(df_motor_SF1$`CI_Wrist-Right`, df_motor_SF1$height_cm, method="spearman")
WL.ci.h<-cor.test(df_motor_SF1$`CI_Wrist-Left`, df_motor_SF1$height_cm, method="spearman")
AR.ci.h<-cor.test(df_motor_SF1$`CI_Ankle-Right`, df_motor_SF1$height_cm, method="spearman")
AL.ci.h<-cor.test(df_motor_SF1$`CI_Ankle-Left`, df_motor_SF1$height_cm, method="spearman")

T.ci.w<-cor.test(df_motor_SF1$`CI_Torso`, df_motor_SF1$weight_kg, method="spearman")
WR.ci.w<-cor.test(df_motor_SF1$`CI_Wrist-Right`, df_motor_SF1$weight_kg, method="spearman")
WL.ci.w<-cor.test(df_motor_SF1$`CI_Wrist-Left`, df_motor_SF1$weight_kg, method="spearman")
AR.ci.w<-cor.test(df_motor_SF1$`CI_Ankle-Right`, df_motor_SF1$weight_kg, method="spearman")
AL.ci.w<-cor.test(df_motor_SF1$`CI_Ankle-Left`, df_motor_SF1$weight_kg, method="spearman")

neg.cor<-c(T.neg$estimate, T.ci.neg$estimate, WR.neg$estimate, WR.ci.neg$estimate, WL.neg$estimate, WL.ci.neg$estimate, AR.neg$estimate, AR.ci.neg$estimate, AL.neg$estimate, AL.ci.neg$estimate)
rme.cor<-c(T.rme$estimate, T.ci.rme$estimate, WR.rme$estimate, WR.ci.rme$estimate, WL.rme$estimate, WL.ci.rme$estimate, AR.rme$estimate, AR.ci.rme$estimate, AL.rme$estimate, AL.ci.rme$estimate)
h.cor<-c(T.h$estimate, T.ci.h$estimate, WR.h$estimate, WR.ci.h$estimate, WL.h$estimate, WL.ci.h$estimate, AR.h$estimate, AR.ci.h$estimate, AL.h$estimate, AL.ci.h$estimate)
w.cor<-c(T.w$estimate, T.ci.w$estimate, WR.w$estimate, WR.ci.w$estimate, WL.w$estimate, WL.ci.w$estimate, AR.w$estimate, AR.ci.w$estimate, AL.w$estimate, AL.ci.w$estimate)

neg.p<-c(T.neg$p.value, T.ci.neg$p.value, WR.neg$p.value, WR.ci.neg$p.value, WL.neg$p.value, WL.ci.neg$p.value, AR.neg$p.value, AR.ci.neg$p.value, AL.neg$p.value, AL.ci.neg$p.value)
rme.p<-c(T.rme$p.value, T.ci.rme$p.value, WR.rme$p.value, WR.ci.rme$p.value, WL.rme$p.value, WL.ci.rme$p.value, AR.rme$p.value, AR.ci.rme$p.value, AL.rme$p.value, AL.ci.rme$p.value)
h.p<-c(T.h$p.value, T.ci.h$p.value, WR.h$p.value, WR.ci.h$p.value, WL.h$p.value, WL.ci.h$p.value, AR.h$p.value, AR.ci.h$p.value, AL.h$p.value, AL.ci.h$p.value)
w.p<-c(T.w$p.value, T.ci.w$p.value, WR.w$p.value, WR.ci.w$p.value, WL.w$p.value, WL.ci.w$p.value, AR.w$p.value, AR.ci.w$p.value, AL.w$p.value, AL.ci.w$p.value)

row<-c("Torso_acc", "Torso_CI", "WR_acc", "WR_CI", "WL_acc", "WL_CI", "AR_acc", "AR_CI", "AL_acc", "AL_CI")
table.cor<-data.frame(row, neg.cor, rme.cor, h.cor, w.cor)
table.p<-data.frame(row, neg.p, rme.p, h.p, w.p)


#save tables to word doc
setwd("XXXXX")
#flextable::save_as_docx(table1, as_flex_table(table1), path="table-demog.docx")
#flextable::save_as_docx(table1, as_flex_table(table2), path="table-SFbehaviours.docx")
#flextable::save_as_docx(table1, as_flex_table(table3), path="table-acc.docx")

#write.csv(table.cor, "correlations-acc-CI-estimates.csv")
#write.csv(table.p, "correlations-acc-CI-pvalue.csv")



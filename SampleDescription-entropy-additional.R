#Additional descriptives
remove(list=ls())

#create preterm sitting ability groups
#compare sitting ability on relevant characteristics: GA, cross motor score

setwd("/Users/trb18168/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/Files/PhD/PhD studies/Data/Processed datasets/Sensor analysis_2022/")
sample<-read.csv("sensor-subsample-description-final-with-MOT.csv")%>%
  mutate(preterm_groups=ifelse(preterm==0, 0, 
                ifelse(preterm==1 & sits_unsupported==1, 1,
                       ifelse(preterm==1 & sits_unsupported==0, 2, NA))))%>%
  select(-X)

setwd("/Users/trb18168/Library/CloudStorage/OneDrive-TheUniversityofLiverpool/Files/TEBC/Research studies/Sensors analysis/Revised")
characteristics<-read.csv("YWC_ScientificReports_TEBCdata.csv")

sample<-merge(sample, characteristics, by.x="ID", by.y="ID", all.x=TRUE)

sample<-sample%>% mutate(ID=as.factor(ID),
                         Sex=factor(sex, levels=c(1, 2), labels=c("Male", "Female")),
                         Preterm=factor(preterm, levels=c(0,1), labels=c("Term", "Preterm")),
                         Ethnicity=ifelse(child_s_ethnicity==1, "Any White background", 
                                          ifelse(child_s_ethnicity<=5, "Any Mixed background", 
                                                 ifelse(child_s_ethnicity<=9, "Any Asian background",
                                                        ifelse(is.na(child_s_ethnicity)==TRUE, NA, 
                                                               "Any other ethnic group")))),
                         age_or_corrage=ifelse(Preterm=="Term", age_9m_m, age_corrected_9m_m),
                         simd=ifelse(simd_quintile %in% c(1, 2, 3), 1, 
                                     ifelse(simd_quintile %in% c(4, 5), 2, NA)),
                         simd=factor(simd, levels=c(1, 2), labels=c("Quintile 1, 2, or 3 (more deprived)", "Quintile 4 or 5 (less deprived")),
                         Singleton=ifelse(no_of_infants==1, "Singleton", "Twin"),
                         Singleton=factor(Singleton, levels=c("Singleton", "Twin"), labels=c("Singleton", "Twin")),
                         Sits_unsupported=factor(sits_unsupported, levels=c(0,1), labels=c("sits supported", "sits unsupported")),
                         Sepsis=ifelse(Early_onset_sepsis==1 | Late_onset_sepsis==1, 1, 0),
                         Retinopathy=ifelse(retinopathy_of_prematurity==1, 1, 0),
                         Sepsis=factor(Sepsis),
                         NEC=factor(Medical_or_surgical_nec),
                         Bronchopulmonary_dysplasia=factor(bronchopulmonary_dysplasia),
                         Retinopathy=factor(Retinopathy))

sample<-sample%>%mutate(degree=ifelse(Maternal_education %in% c(1,2,3,4,5), "No", 
                                       ifelse(Maternal_education %in% c(4, 5, 6, 7), "Yes", NA)),
                        degree=factor(degree, levels=c("No", "Yes")),
                        ethnicity=factor(Child_ethnicity),
                        comorbidities=rowSums(sapply(select(., c("Sepsis", "NEC",	"Bronchopulmonary_dysplasia",	"Retinopathy")), function(x) x==1), na.rm = TRUE))

library(gtsummary)
df_ptgroups_table<-sample%>%select(preterm_groups, age_or_corrage, birthweight_g, birthweight_z_score, gestation, Sex, Ethnicity, degree, simd, Singleton, height_cm, weight_kg, MOT, MOT_gross, MOT_fine, comorbidities, Sepsis, NEC, Bronchopulmonary_dysplasia, Retinopathy)
table1s<- tbl_summary(df_ptgroups_table, by="preterm_groups",
                      statistic=list(all_continuous()~c("{mean} ({sd})")),
                      digits=list(everything()~1))


df_sample_table<-sample%>%select(Preterm, age_or_corrage, birthweight_g, birthweight_z_score, gestation, Sex, Ethnicity, degree, degree, Singleton, height_cm, weight_kg, Sits_unsupported, MOT, MOT_gross, MOT_fine, comorbidities, Sepsis, NEC, Bronchopulmonary_dysplasia, Retinopathy)

#XXXX weight wrong - change from 81kg to 8.1kg (checked on REDCAP)
#infant characteristics: age at visit, demographics (gender, ethnicity, simd, twin birth), 9m anthropometry
table1<- tbl_summary(df_sample_table, by="Preterm",
                     statistic=list(all_continuous()~c("{mean} ({sd})")),
                     digits=list(everything()~1))
table1



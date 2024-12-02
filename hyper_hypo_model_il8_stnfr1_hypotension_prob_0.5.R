#Clear R environment
rm(list = ls())

#Import RESERVE-U-1-EBB dataset
reserve <- read.csv(file.choose(), header=TRUE)

#PIDs to fixed row identifiers
rownames(reserve) = reserve$pid
reserve$pid = NULL

#Restrict to patients with qSOFA>=1
reserve <- subset(reserve, qsofa1p==1)

#Create variable for hypotension - SBP<90 or MAP<65
reserve$map2 <- reserve$dbp2 + (0.333333*(reserve$sbp2-reserve$dbp2))
reserve$hypotension = ifelse((reserve$sbp3<90 | reserve$map2<65), 1, 0)

#Model: IL-8, sTNFR1, and vasopressors (hypotension)
library(dplyr)
library(ggplot2)

#Probability threshold for phenotype assignment 
threshold <- 0.5

#Apply the 1 + natural log transformation to the continuous variables
reserve <- reserve %>%
  mutate(
    il8_log = log1p(il8),
    stnfr1_log = log1p(stnfr1))

#Define the model coefficients
intercept <- -18.4764
coef_il8 <- 1.3013
coef_stnfr1 <- 1.3367
coef_hypotension <- 2.3439

#Apply the model to calculate log-odds and predicted probabilities
reserve <- reserve %>%
  mutate(
    phenotype_prob_1 = 1 / (1 + exp(-(intercept + 
                                        coef_il8 * il8_log + 
                                        coef_stnfr1 * stnfr1_log +
                                        coef_hypotension * hypotension))),
    phenotype_binary_1 = ifelse(phenotype_prob_1 > threshold, "hyperinflammatory", "hypoinflammatory")
  )

hist(reserve$phenotype_prob_1)

#Compare clinical characteristics across phenotypes
library(gmodels)
library(dplyr)
library(skimr)
library(gtsummary)

#Create variable for qSOFA score >=3
reserve$qsofa3p <- ifelse(reserve$qsofa_score==3, 1, 0)

#Table - Clinical characteristics by phenotype 
reserve_table_df <- reserve[, c("sex", 
                                "age",
                                "illnessduration_enroll",
                                "tempmax",
                                "heartrate3", 
                                "resprate3", 
                                "sbp3",
                                "o2sat3", 
                                "ams", 
                                "qsofa2p", 
                                "qsofa3p",
                                "MEWS_score",
                                "UVA_score", 
                                "hivrdtresult", 
                                "hivstage34", 
                                "malariardtresult", 
                                "microtbdx", 
                                "urinetblamresult",
                                "influenzapcrresult", 
                                "hospoutcome", 
                                "dcperf16andgr", 
                                "death30d",
                                "phenotype_binary_1")]

theme_gtsummary_compact()

#Replace KPS of 50 in patients who were transferred to higher level facility 
reserve_table_df$dcperf16andgr[reserve_table_df$hospoutcome == "Referred/transferred to other facility"] <- 50

reserve_table_df %>%
  tbl_summary(
    by = phenotype_binary_1,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n}/{N} ({p}%)"
    ),
    type = list(sex ~ "dichotomous", 
                age ~ "continuous",
                illnessduration_enroll ~ "continuous",
                tempmax ~ "continuous",
                heartrate3 ~ "continuous", 
                resprate3 ~ "continuous", 
                sbp3 ~ "continuous",
                o2sat3 ~ "continuous", 
                ams ~ "dichotomous",   
                qsofa2p ~ "dichotomous", 
                qsofa3p ~ "dichotomous",
                MEWS_score ~ "continuous",
                UVA_score ~ "continuous", 
                hivrdtresult ~ "dichotomous", 
                hivstage34 ~ "dichotomous", 
                malariardtresult ~ "dichotomous", 
                microtbdx ~ "dichotomous", 
                urinetblamresult ~ "dichotomous",
                influenzapcrresult ~ "dichotomous", 
                hospoutcome ~ "categorical", 
                dcperf16andgr ~ "continuous",
                death30d ~ "dichotomous"),
    digits = list(all_continuous() ~ 0, all_dichotomous() ~ c(0,0,1), all_categorical() ~ c(0,0,1), tempmax ~ c(1, 1)),
    missing_text = "(Missing)",
    label = list(sex ~ "Sex", 
                 age ~ "Age, years",
                 illnessduration_enroll ~ "Illness duration prior to enrollment, days",
                 tempmax ~ "Temperature, \u00b0C",
                 heartrate3 ~ "Heart rate, beats/min", 
                 resprate3 ~ "Respiratory rate, breaths/min", 
                 sbp3 ~ "Systolic blood pressure, mmHg",
                 o2sat3 ~ "Oxygen saturation, %", 
                 ams ~ "Altered mental status",   
                 qsofa2p ~ "qSOFA score \u2265 2",   
                 qsofa3p ~ "qSOFA score \u2265 3",
                 MEWS_score ~ "Modified Early Warning Score",
                 UVA_score ~ "Universal Vital Assessment score", 
                 hivrdtresult ~ "Person living with HIV (PLWH)", 
                 hivstage34 ~ "WHO HIV clinical stage 3 or 4", 
                 malariardtresult ~ "Malaria RDT positive", 
                 microtbdx ~ "Microbiological TB positive", 
                 urinetblamresult ~ "Urine TB-LAM positive if PLWH",
                 influenzapcrresult ~ "Influenza PCR result", 
                 hospoutcome ~ "Hospital outcome", 
                 dcperf16andgr ~ "KPS at alive discharge or transfer", 
                 death30d ~ "Dead at 30 days post-discharge")) %>% add_overall() %>% add_n() %>% add_p(pvalue_fun = ~style_sigfig(., digits = 3)) %>% as_gt() 

#Probability of 30d mortality by that of hyperinflammatory phenotype assignment
library(tidyverse)
reserve_death <- reserve %>% drop_na(death30d, age, sex, illnessduration_enroll, hivrdtresult)

library(ggplot2)
phenotypemodel <- glm(death30d ~ phenotype_prob_1 + age + sex + illnessduration_enroll + hivrdtresult,
                        data = reserve_death, family = "binomial")
summary(phenotypemodel)
exp(phenotypemodel$coefficients)
exp(confint.default(phenotypemodel))

reserve_death$death_prob <- predict(phenotypemodel, newdata = reserve_death, type = "response")
reserve_death

hyperinflammatory_death_prob <- ggplot(reserve_death, aes(x=phenotype_prob_1, y=death_prob)) +
  geom_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = TRUE, color = "navy") + xlim(0,1) + ylim(0,1)
hyperinflammatory_death_prob <- hyperinflammatory_death_prob + xlab("Probability of hyperinflammatory phenotype") + ylab("Probability of death at 30 days") + theme_bw()
hyperinflammatory_death_prob <- hyperinflammatory_death_prob + theme(text=element_text(size=16)) + annotate("text", x=0.75, y=0.05, label= "aOR 5.16 (95% CI 1.45-18.33)", size = 5)
hyperinflammatory_death_prob    

#Differential protein expression by phenotype
#Import Olink dataset
reserve_olink <- read.csv(file.choose(), header=TRUE)

#PIDs to fixed row identifiers
rownames(reserve_olink) = reserve_olink$pid
reserve_olink$pid = NULL

#Restrict to patients with qSOFA>=1
reserve_olink <- subset(reserve_olink, qsofa1p==1)

#Merge phenotype assignment into Olink dataset
library(tidyverse)

reserve_phenotype <- reserve %>% select(phenotype_binary_1)

reserve_olink_phenotype <- merge(reserve_olink, reserve_phenotype,
                          by = 'row.names', all = FALSE) 

#PIDs to fixed row identifiers
rownames(reserve_olink_phenotype) = reserve_olink_phenotype$Row.names
reserve_olink_phenotype$Row.names = NULL

#Drop PID 361 - extreme outlier
reserve_olink_phenotype <- reserve_olink_phenotype[row.names(reserve_olink_phenotype) != "361",]

#Select and view the log2-transformed Olink data
biomarkers <- reserve_olink_phenotype[1:242, c(225:409)]
names(biomarkers)

#Drop biomarkers with <20% of NPX values above panel’s estimated limit of detection
library(dplyr)
biomarkers <- dplyr::select(biomarkers, -c("il1alpha",
                                           "il2", 
                                           "il33", 
                                           "il4", 
                                           "il13",
                                           "prcp",
                                           "ltbp2",
                                           "sod1",
                                           "itgam",
                                           "fap",
                                           "mfap5"))  

#Reformat column names
colnames(biomarkers) <- base::toupper(colnames(biomarkers))
colnames(biomarkers)[colnames(biomarkers) == "IL8"] = "IL-8"
colnames(biomarkers)[colnames(biomarkers) == "IL2"] = "IL-2"
colnames(biomarkers)[colnames(biomarkers) == "IL7"] = "IL-7"
colnames(biomarkers)[colnames(biomarkers) == "IL6"] = "IL-6"
colnames(biomarkers)[colnames(biomarkers) == "IL18"] = "IL-18"
colnames(biomarkers)[colnames(biomarkers) == "IL15"] = "IL-15"
colnames(biomarkers)[colnames(biomarkers) == "IL5"] = "IL-5"
colnames(biomarkers)[colnames(biomarkers) == "IL10"] = "IL-10"
colnames(biomarkers)[colnames(biomarkers) == "TNF"] = "TNF"
colnames(biomarkers)[colnames(biomarkers) == "IFNGAMMA"] = "IFN-\u03B3"
colnames(biomarkers)[colnames(biomarkers) == "IL12RB1"] = "IL-12RB1"
colnames(biomarkers)[colnames(biomarkers) == "IL12"] = "IL-12"
colnames(biomarkers)[colnames(biomarkers) == "IL7R"] = "IL-7R"
colnames(biomarkers)[colnames(biomarkers) == "PDGFSUBUNITB"] = "PDGFB"
colnames(biomarkers)[colnames(biomarkers) == "LAPTGFBETA1"] = "LAPTGFB1"
colnames(biomarkers)[colnames(biomarkers) == "PDL1"] = "PD-L1"
colnames(biomarkers)[colnames(biomarkers) == "PDL2"] = "PD-L2"

biomarkers[biomarkers$PHENOTYPE_BINARY_1=="hypoinflammatory", "PHENOTYPE_BINARY"] <- "HYPO"
biomarkers[biomarkers$PHENOTYPE_BINARY_1=="hyperinflammatory", "PHENOTYPE_BINARY"] <- "HYPER"
biomarkers$PHENOTYPE_BINARY <- as.factor(biomarkers$PHENOTYPE_BINARY)

data = biomarkers

#Remove phenotype column and transpose the data so that proteins are rows and samples are columns
expression_data <- biomarkers[1:242, c(1:173)] 
names(expression_data)
expression_data <- t(expression_data)  

group_assignments = factor(biomarkers$PHENOTYPE_BINARY, levels=c("HYPO", "HYPER"))
group_assignments

# Create the design matrix for the linear model
design <- model.matrix(~ 0 + group_assignments)
colnames(design)[1:2] <- levels(group_assignments)
colnames(design)
# run limma
library(limma)
fit <- lmFit(expression_data, design)

# Apply empirical Bayes statistics
fit <- eBayes(fit)

# Create a contrast matrix for the comparison of phenotypes
contrast.matrix <- makeContrasts(HypervsHypo = HYPER - HYPO, levels = design)

# Fit the contrast model
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get the table of top differentially expressed proteins
results <- topTable(fit2, coef="HypervsHypo", adjust.method="BH", sort.by="P", number=Inf)

#Volcano plot
library(EnhancedVolcano)
# Create a data frame for the volcano plot
results$Prot <- toupper(rownames(results))  # Create a column for gene names if not already present
names(results$logFC) <- toupper(results$Prot)
names(results$adj.P.Val) <- toupper(results$Prot)
# 

results$selectedLab <- ifelse(results$adj.P.Val < 0.05 & abs(results$logFC) > 0.5, as.character(results$Prot), NA)
EnhancedVolcano(
  results,
  lab = NA,
  x = 'logFC',
  y = 'adj.P.Val',
  ylim = c(0,10),
  xlim = c(-3,3),
  pCutoff = 0.05,
  FCcutoff = 0.5,  
  pointSize = 4,
  legendPosition = 'none',
  title = NULL,
  subtitle = NULL,
  labSize = 6,
  legendLabSize = 15.0,
  axisLabSize = 30,
  max.overlaps = 150,
  drawConnectors = TRUE
)+ geom_text_repel(
  data = subset(results, !is.na(selectedLab)),  # Only use rows where labels are not NA
  aes(label = selectedLab, x = logFC, y = -log10(adj.P.Val)),  # Ensure 'label' aesthetic is correctly mapped
  size = 6,  # Ensure size matches labSize in EnhancedVolcano
  box.padding = unit(0.3, "lines"),  # Adjust padding to avoid overlapping
  point.padding = unit(0.2, "lines"),  # Increase point padding
  max.overlaps = 150,  # Allow more overlaps
  #  segment.color = NA
)

#Performance of clinical +/- microbiological variables for phenotype discrimination
#New dataset
reserve_model <- reserve
reserve_model <- reserve_model[!is.na(reserve_model$hivrdtresult),]
reserve_model <- reserve_model[!is.na(reserve_model$malariardtresult),]
reserve_model <- reserve_model[!is.na(reserve_model$microtbdx),]

#Recode variables
reserve_model$hyper <- NA
reserve_model[reserve_model$phenotype_binary_1=="hyperinflammatory", "hyper"] <- 1
reserve_model[reserve_model$phenotype_binary_1=="hypoinflammatory", "hyper"] <- 0
reserve_model$hyper <- as.factor(reserve_model$hyper)

#Clinical + micro
clinicalmicromodel <- glm(hyper ~ age + 
                            sex + 
                            tempmax + 
                           heartrate3 + 
                           resprate3 + 
                           sbp3 + 
                           o2sat3 +
                          avpu3 +
                          microtbdx +
                          hivrdtresult +
                          malariardtresult,
                         data = reserve_model, family = "binomial")
summary(clinicalmicromodel)
exp(clinicalmicromodel$coefficients)
exp(confint.default(clinicalmicromodel))

#ROC
library(pROC)

predicted <- predict(clinicalmicromodel, reserve_model, type="response")
auc(reserve_model$hyper, predicted)

rocobj <- plot.roc(reserve_model$hyper, clinicalmicromodel$fitted.values,
                   percent=FALSE,
                   ci = TRUE,                  
                   print.auc = TRUE,
                   print.auc.x = 0.6, 
                   print.auc.y = 0.05,
                   legacy.axes = TRUE,
                   print.auc.pattern = "AUROC %.2f (%.2f-%.2f)")           
ciobj <- ci.se(rocobj, boot.n=10000,                         
               specificities = seq(0, 1, 0.05)) 
plot(ciobj, type = "shape", col = "#8da0cb")     


#Clinical
#New dataset
reserve_model_2 <- reserve

#Recode variables
reserve_model_2$hyper <- NA
reserve_model_2[reserve_model_2$phenotype_binary_1=="hyperinflammatory", "hyper"] <- 1
reserve_model_2[reserve_model_2$phenotype_binary_1=="hypoinflammatory", "hyper"] <- 0
reserve_model_2$hyper <- as.factor(reserve_model_2$hyper)

clinicalmodel <- glm(hyper ~ age +
                            sex + 
                            tempmax + 
                            heartrate3 + 
                            resprate3 + 
                            sbp3 + 
                            o2sat3 +
                            avpu3,
                          data = reserve_model_2, family = "binomial")
summary(clinicalmodel)
exp(clinicalmodel$coefficients)
exp(confint.default(clinicalmodel))

#ROC
library(pROC)

predicted <- predict(clinicalmodel, reserve_model_2, type="response")
auc(reserve_model_2$hyper, predicted)



rocobj <- plot.roc(reserve_model_2$hyper, clinicalmodel$fitted.values,
                   percent=FALSE,
                   ci = TRUE,                  
                   print.auc = TRUE,
                   print.auc.x = 0.6, 
                   print.auc.y = 0.05,
                   legacy.axes = TRUE,
                   print.auc.pattern = "AUROC %.2f (%.2f-%.2f)")           
ciobj <- ci.se(rocobj, boot.n=10000,                         
               specificities = seq(0, 1, 0.05))
plot(ciobj, type = "shape", col = "#8da0cb")     



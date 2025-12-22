##KN-DISSERTATION project. Mutation data (ARID1A and CTNNB1 genes)
#survival analysis
Sys.setenv(LANG = "en")
#libraries
library(survival)
library(survminer)
library(tidyverse)
library(grid)
library(gridExtra)
library(patchwork)
library(magick)
library(survivalROC)
library(purrr)
library(broom)
library(dplyr)
library(timeROC)
library(coxphf)
library(readxl)
#set wd for plots###########################################
setwd("C:/Users/Ieva/rprojects/outputs_all/DISS/")
#upload data##############################################################
KN_data <- read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/audiniu_mut_exprs_met20250709.xlsx")
#mutation genes
expression3 <- c("ARID1A","CTNNB1", "FBXW7")
#make smaller mution focused df
Arid1a_columns <- c("Histology", "Grupė_Ieva", expression3,
                     "ARID1A_tumor_mut","ARID1A_tumor_VUS",  "ARID1A_tumor_type",
                    "CTNNB1_tumor_mut",  "TP53_tumor_mut", "KN", "CA125", "Amžius")
#make arid1a df
ARID1A_df <- KN_data[, colnames(KN_data) %in% Arid1a_columns]
#fix mutations
cols <- c("ARID1A_tumor_mut", "CTNNB1_tumor_mut",  "TP53_tumor_mut")
ARID1A_df[cols] <- lapply(ARID1A_df[cols], function(x) ifelse(is.na(x), "Be mutacijų", "Mutacija"))
##only up to kn-95 the samples had mutation data, others should be NA
rownames(ARID1A_df) <- ARID1A_df$KN
not_mut <- c("KN-96" , "KN-97",  "KN-99",  "KN-100" ,"KN-101", "KN-103", "KN-104" ,
             "KN-105", "KN-106", "KN-107", "KN-108", "KN-109", "KN-110", "KN-111", "KN-112")
ARID1A_df[row.names(ARID1A_df) %in% not_mut, cols] <- NA
#histology to lt version
ARID1A_df$Histology <- ifelse(ARID1A_df$Histology %in% c("Endometrial", "Endometriod"), "Endometrioidinis",
                              ARID1A_df$Histology )
ARID1A_df$Histology <- ifelse(ARID1A_df$Histology == "Endometriois", "Endometriozė",
                              ARID1A_df$Histology )
ARID1A_df$Histology <- ifelse(ARID1A_df$Histology == "Cystis", "Cista",
                              ARID1A_df$Histology )
ARID1A_df$Histology <- ifelse(ARID1A_df$Histology == "Mucinous", "Mucininis",
                              ARID1A_df$Histology )
ARID1A_df$Histology <- ifelse(ARID1A_df$Histology == "Serous", "Serozinis",
                              ARID1A_df$Histology )
ARID1A_df$Histology <- ifelse(ARID1A_df$Histology == "Clear cell", "Šviesių lastelių",
                              ARID1A_df$Histology )
ARID1A_df$Histology <- ifelse(ARID1A_df$Histology == "RSS", "Riziką mažinanti operacija",
                              ARID1A_df$Histology )
ARID1A_df$Histology <- ifelse(ARID1A_df$Histology == "Granuloza", "Granulosa",
                              ARID1A_df$Histology )
ARID1A_df$Histology <- factor(ARID1A_df$Histology)
#rename arid1a VUS
ARID1A_df$ARID1A_tumor_VUS[is.na(ARID1A_df$ARID1A_tumor_VUS)] <- "Be mutaciju"
ARID1A_df$ARID1A_tumor_VUS <- ifelse(ARID1A_df$ARID1A_tumor_VUS == "Pathogenic", "Patogeninė, VUS",
                                         ARID1A_df$ARID1A_tumor_VUS )
ARID1A_df$ARID1A_tumor_VUS <- ifelse(ARID1A_df$ARID1A_tumor_VUS == "VUS, benign", "Gerybinė, VUS",
                                         ARID1A_df$ARID1A_tumor_VUS )
ARID1A_df$ARID1A_tumor_VUS
#get survival data as of 2025-09-11
SURVIVAL_KN <- openxlsx::read.xlsx("../../OTHER DATA/KN-DISSERTATION FILES/KN_MIRTIES_FAILAS_20250911.xlsx")
#make only surv df
SURV <- SURVIVAL_KN[, c(2, 3,20, 21)]
#make mutation + survival df
ARID1A_df$KN %in% SURV$KN  #all the same
MAIN_DF <- left_join(ARID1A_df, SURV, by = "KN") #join
colnames(MAIN_DF)
#Make main df OC only
MAIN_DF <- MAIN_DF[ MAIN_DF$Grupė_Ieva != "Benign", ]
#SURVIVAL WITH AR1D1A MUTATION COUNT############################################
str(MAIN_DF$ARID1A_tumor_mut)
MAIN_DF$ARID1A_tumor_mut <- factor(MAIN_DF$ARID1A_tumor_mut)
with(MAIN_DF, table(ARID1A_tumor_mut, STATUS)) #no death with mutation
#univaribale cox - wrong as as there is no deaths in the mutated group
cox_model <- coxph(Surv(OS, STATUS) ~ ARID1A_tumor_mut, data = MAIN_DF)
sum_cox <- summary(cox_model)
print(sum_cox)
#do this instead
survdiff(Surv(OS, STATUS) ~ ARID1A_tumor_mut, data = MAIN_DF)
#or this # Fit a Kaplan-Meier model
km_fit_ARID1A_mut <- survfit(Surv(OS, STATUS) ~ ARID1A_tumor_mut, data = MAIN_DF)
summary(km_fit_ARID1A_mut)$table #overall summary
summary(km_fit_ARID1A_mut)$table[, "rmean"]
#KM summary at times 12, 24, 36 months
time_points <- c(12, 36, 60)
km_summary <- summary(km_fit_ARID1A_mut, times = time_points)
km_summary

#add RMS onto the KM plot##################################
km_table <- summary(km_fit_ARID1A_mut)$table
rms_values <- km_table[, "rmean"]

# get p-value from survdiff (log-rank test)
logrank_test <- survdiff(Surv(OS, STATUS) ~ ARID1A_tumor_mut, data = MAIN_DF)
p_val <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)

# make subtitle
rms_subtitle <- paste0(
  "Vidutinis išgyvenamumo laikas, mėnesiais: ",
  "Be mutacijų = ", round(rms_values["ARID1A_tumor_mut=Be mutacijų"], 1), ", ",
  "Mutacija = ", round(rms_values["ARID1A_tumor_mut=Mutacija"], 1),
  "; Log-rank p = ", signif(p_val, 2)
)

# plot km
test_survplot_Arid <- ggsurvplot(
  km_fit_ARID1A_mut,
  data = MAIN_DF,
  pval = FALSE,              
  risk.table = TRUE,
  title = expression( italic("ARID1A") * " mutacijų sąsaja su bendru išgyvenamumu KV audinių imtyje"),
  risk.table.title = "Pacientų skaičius rizikos grupėje",
  subtitle = rms_subtitle,
  xlab = "Bendras išgyvenamumo laikas",
  ylab = "Išgyvenamumo tikimybė",
  palette = c("darkblue", "maroon"),
  legend.title = "Mutacija",
  legend.labs = c("Be mutacijų", "Mutacija")
)

test_survplot_Arid

#save
png("C:/Users/Ieva/rprojects/outputs_all/DISS/ARID1A_mut_KM_ALL_20250925.png",
    width = 1000, height = 600, res = 100) # width and height in pixels, resolution in dpi
test_survplot_Arid #
dev.off() # Close the PNG device

#time roc 
MAIN_DF$ARID1A_mut_num <- ifelse(MAIN_DF$ARID1A_tumor_mut == "Mutacija", 1, 0) #make mutation a factor
roc_obj <- timeROC(
  T = MAIN_DF$OS,
  delta = MAIN_DF$STATUS,
  marker = MAIN_DF$ARID1A_mut_num,
  cause = 1,
  weighting = "marginal",
  times = c(12, 36, 60),   # 1yr, 3yr, 5yr
  iid = TRUE
)
plot(roc_obj, time = 60)   # ROC curve at 60 months
roc_obj

#coords
# Extract index for 60 months
idx <- which(roc_obj$times == 60)

# Sensitivity (TPR) and Specificity (1 - FPR) at each cutoff
sens_60 <- roc_obj$TP[, idx]
spec_60 <- 1 - roc_obj$FP[, idx]

results <- data.frame(
  threshold = roc_obj$cutoffs,
  sensitivity = sens_60,
  specificity = spec_60
)
head(results)


#english ARID1A KM################

# make subtitle
rms_subtitleEN <- paste0(
  "Median survival, months: ",
  "No mutations = ", round(rms_values["ARID1A_tumor_mut=Be mutacijų"], 1), ", ",
  "Mutations = ", round(rms_values["ARID1A_tumor_mut=Mutacija"], 1),
  "; Log-rank p = ", signif(p_val, 2)
)

# plot km
test_survplot_AridEN <- ggsurvplot(
  km_fit_ARID1A_mut,
  data = MAIN_DF,
  pval = FALSE,              
  risk.table = TRUE,
  title = expression( italic("ARID1A") * " mutations vs survival in the ovarian cancer group"),
  #risk.table.title = "Pacientų skaičius rizikos grupėje",
  subtitle = rms_subtitleEN,
  xlab = "Overall survival, months",
  #ylab = "Išgyvenamumo tikimybė",
  palette = c("darkblue", "maroon"),
  legend.title = "Mutation",
  legend.labs = c("No mutations", "Mutation")
)

test_survplot_AridEN

#save
png("C:/Users/Ieva/rprojects/outputs_all/DISS/ARID1A_mut_KM_ALL_EN20251215.png",
    width = 1000, height = 600, res = 100) # width and height in pixels, resolution in dpi
test_survplot_AridEN #
dev.off() # Close the PNG device

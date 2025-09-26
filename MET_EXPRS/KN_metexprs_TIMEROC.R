#KN-DISSERTATION project. 10 gene panel - genes selected via statistical analysis of TCGA and GTEx data
#KN 10 gene analysis data (cv 2 %, one threshold, 65 cases!, restored KN-59) 2025-02-17
#KN-DISSERTATION project. MET_EXPRS - methylation, expression, mutation data
#survival analysis, time rocs and models, only NOTCH
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
#set wd for plots
setwd("C:/Users/Ieva/rprojects/outputs_all/DISS/")
#get one file for all data
ALL_EXPRESSION_DF <- readRDS("../../OTHER DATA/KN-DISSERTATION FILES/ALL_expresssion_Df20250709.RDS")
#get survival data as of 2025-09-11
SURVIVAL_KN <- openxlsx::read.xlsx("../../OTHER DATA/KN-DISSERTATION FILES/KN_MIRTIES_FAILAS_20250911.xlsx")
#make only surv df
SURV <- SURVIVAL_KN[, c(2, 3,20, 21)]
head(SURV)
#join with main data
ALL_SURV_EXPRESSION <- left_join(ALL_EXPRESSION_DF, SURV, by = "patient_id_aud")
#remove random columns or rename
ALL_SURV_EXPRESSION <-  ALL_SURV_EXPRESSION[, -c(1, 15, 32)]
ALL_SURV_EXPRESSION <- ALL_SURV_EXPRESSION %>%
  rename(tumor = tumor.y ,
         KN = KN.y )
rownames(ALL_SURV_EXPRESSION) <- ALL_SURV_EXPRESSION$patient_id_aud
#remove NA from new surv data
ALL_SURV_EXPRESSION <- ALL_SURV_EXPRESSION%>%
  filter(!is.na(OS), !is.na(STATUS)) #58 observations
#make levels of my gene expression data
methylation <- c("ARID1A_met", "CDX2", "ALX4", "HOPX")
genes <- c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4",
           "ARID1A", "CTNNB1", "FBXW7",
           "JAG2", "DLL1", "HES1",
           "EXO1", "RAD50", "PPT2", "LUC7L2", 
           "PKP3", "CDCA5", "ZFPL1", "VPS33B", 
           "GRB7", "TCEAL4")
genes_f <- paste0(genes, "_f")
genes10 <- c("EXO1", "RAD50", "PPT2", "LUC7L2", 
             "PKP3", "CDCA5", "ZFPL1", "VPS33B", 
             "GRB7", "TCEAL4")
genes_notch <- genes[!genes %in% genes10]  #notch genes
#in plots i use names with _f so add it
genes10_f <- paste0(genes10, "_f")
genes_notch_f <- paste0(genes_notch, "_f")

#make levels of gene expresssion data
ALL_SURV_EXPRESSION <- ALL_SURV_EXPRESSION %>%
  mutate(
    across(all_of(genes),
           ~ factor(if_else(. > median(., na.rm = TRUE), "High", "Low")),
           .names = "{.col}_f")
  )

#rename AGE (comes up later in multivariable cox)
ALL_SURV_EXPRESSION <- ALL_SURV_EXPRESSION %>%
  rename(Age = Amžius )
colnames(ALL_SURV_EXPRESSION)

#make factors
ALL_SURV_EXPRESSION$Stage2 <- factor(ALL_SURV_EXPRESSION$Stage2)
ALL_EXPRESSION_DF$Stage4 <- factor(ALL_EXPRESSION_DF$Stage4)
ALL_EXPRESSION_DF$Grade2 <- factor(ALL_EXPRESSION_DF$Grade2)
ALL_EXPRESSION_DF$CA125_f <- factor(ALL_EXPRESSION_DF$CA125_f)

#remove non-OC (benign) cases#############################
OC_SURV_EXPRESSION <- ALL_SURV_EXPRESSION %>%
  filter(Grupė_Ieva != "Benign") #55 left#

#HGSOC ONLY##########################
HGSOC_SURV_EXPRESSION <- ALL_SURV_EXPRESSION %>%
  filter(Grupė_Ieva == "HGSOC") #41 left#
table(HGSOC_SURV_EXPRESSION$STATUS, useNA = "a") #17 dead

#RISK SCORE NOTCH GENES, OC ONLY#####################################################
#gene_data2NOTCH - only OC, notch
gene_data2NOTCH <- OC_SURV_EXPRESSION[, colnames(OC_SURV_EXPRESSION) %in% genes_notch]
# Example with 10 biomarkers
cox_model_OC_notch <- coxph(
  Surv(OS, STATUS) ~ NOTCH1 + NOTCH2 + NOTCH3 + NOTCH4 + ARID1A + CTNNB1 + FBXW7 ,
  data = OC_SURV_EXPRESSION
)
summary(cox_model_OC_notch)
#get coeficients
coefs2_notch <- coef(cox_model_OC_notch)
# Risk score = sum( gene_expression * coefficient )
risk_scores_test2notch <- rowSums(sweep(gene_data2NOTCH, 2, coefs2_notch, "*"))
# View the risk scores
print(risk_scores_test2notch) # now I have some risk scores
#add risk scores to the clin_df_joined_test
OC_SURV_EXPRESSION$RiskScore_notch <- risk_scores_test2notch
#create df wih survival data
surv_df_test2notch <- OC_SURV_EXPRESSION[, colnames(OC_SURV_EXPRESSION) %in%
                                           c("OS", "STATUS", genes_notch, "RiskScore_notch", "patient_id_aud", "Age", "CA125")]

rownames(surv_df_test2notch) <- surv_df_test2notch$patient_id_aud

# Calculate the median risk score
median_risk2notch <- median(surv_df_test2notch$RiskScore_notch, na.rm = TRUE) #3.064142
# Create a new factor column based on the median value
surv_df_test2notch$RiskGroup_notch <- ifelse(surv_df_test2notch$RiskScore_notch <= median_risk2notch,
                                             "Low Risk", "High Risk")
#Create a survival object
surv_object2notch <- Surv(time = surv_df_test2notch$OS,
                          event = surv_df_test2notch$STATUS )
# Fit Univariable cox model for notch model ############################
cox_model_notch_oc <- coxph(Surv(OS, STATUS) ~RiskGroup_notch , data = surv_df_test2notch)
summary(cox_model_notch_oc)
# Fit Multivariable cox model for notch model ############################
cox_model_notch_oc2 <- coxph(Surv(OS, STATUS) ~RiskGroup_notch + Age + CA125, data = surv_df_test2notch)
summary(cox_model_notch_oc2)
# Fit a Kaplan-Meier model for notch model#################################
km_fit2notch <- survfit(surv_object2notch ~ RiskGroup_notch, data = surv_df_test2notch)
# Plot the Kaplan-Meier curve using ggsurvplot
test_survplot2_notch <- ggsurvplot(km_fit2notch, data = surv_df_test2notch, 
                                   #pval = TRUE,  # Show p-value of the log-rank test
                                   risk.table = TRUE, 
                                   risk.table.title = "Pacientų skaičius rizikos grupėje",
                                   title = bquote("A    " ~
                                                    "Didelės vs. mažos rizikos atvejai KV imtyje, Notch, Wnt ir" ~ italic("ARID1A") ~ "genų raiškos kombinacija"),
                                   xlab = "Bendras išgyvenamumo laikas",
                                   ylab = "Išgyvenamumo tikimybė",
                                   palette = c("turquoise", "deeppink"),  # Color palette for groups
                                   legend.title = "Rizikos grupė", 
                                   legend.labs = c("Mažas rizikos balas", "Didelis rizikos balas"))
# Add subtitle form cox result
test_survplot2_notch$plot <- test_survplot2_notch$plot +
  labs(subtitle = "Uni PR = 0,41  (95 % PI: 0,15–1,07); Multi PR = 0,48 (95 % PI: 0,16–1,47), Long-rank p = 0,059, n= 47")

print(test_survplot2_notch)
#save
png("KM_plot_notch_arid1a_OC_w_HR20250925.png", width = 1800, height = 1100, res = 150)
print(test_survplot2_notch)  # print the full ggsurvplot object
dev.off()

#NOTCH and METHYLATION GENES, OC ONLY#####################################################
#gene_data2NOTCH - only OC, notch
gene_data14 <- OC_SURV_EXPRESSION[, colnames(OC_SURV_EXPRESSION) %in% c(genes_notch, methylation)]
# Example with 10 biomarkers
cox_model_OC_14 <- coxph(
  Surv(OS, STATUS) ~ NOTCH1 + NOTCH2 + NOTCH3 + NOTCH4 + ARID1A + CTNNB1 + FBXW7 + JAG2 + DLL1 + HES1 +ARID1A_met + ALX4 + HOPX + CDX2 ,
  data = OC_SURV_EXPRESSION
)
summary(cox_model_OC_14)
#get coeficients
coefs2_14 <- coef(cox_model_OC_14)
# Risk score = sum( gene_expression * coefficient )
risk_scores_14 <- rowSums(sweep(gene_data14, 2, coefs2_14, "*"))
# View the risk scores
print(risk_scores_14) # now I have some risk scores
#add risk scores to the clin_df_joined_test
OC_SURV_EXPRESSION$RiskScore_14 <- risk_scores_14
#create df wih survival data
surv_df_test14<- OC_SURV_EXPRESSION[, colnames(OC_SURV_EXPRESSION) %in%
                                      c("OS", "STATUS", genes_notch, "RiskScore_14", "patient_id_aud", "Age", "CA125")]

rownames(surv_df_test14) <- surv_df_test14$patient_id_aud

# Calculate the median risk score
median_risk14<- median(surv_df_test14$RiskScore_14, na.rm = TRUE) #3.064142
# Create a new factor column based on the median value
surv_df_test14$RiskGroup_14 <- ifelse(surv_df_test14$RiskScore_14 <= median_risk14,
                                      "Low Risk", "High Risk")
#Create a survival object
surv_object14<- Surv(time = surv_df_test14$OS,
                     event = surv_df_test14$STATUS )
# Fit Univariable cox model for notch methylation model ############################
cox_model_14 <- coxph(Surv(OS, STATUS) ~RiskGroup_14 , data = surv_df_test14)
summary(cox_model_14)
# Fit Multivariable cox model for notch methylation model ############################
cox_model_14m <- coxph(Surv(OS, STATUS) ~RiskGroup_14 + Age + CA125, data = surv_df_test14)
summary(cox_model_14m)
# Fit a Kaplan-Meier model for notch methylation model#################################
km_fit14<- survfit(surv_object14 ~ RiskGroup_14, data = surv_df_test14)
# Plot the Kaplan-Meier curve using ggsurvplot
test_survplot14<- ggsurvplot(km_fit14, data = surv_df_test14, 
                             #pval = TRUE,  # Show p-value of the log-rank test
                             risk.table = TRUE,  # Add risk table below the plot
                             risk.table.title = "Pacientų skaičius rizikos grupėje",
                             title = bquote("B  " ~  
                                              "Didelės vs. mažos rizikos atvejai KV imtyje, Promotorių metilinimo, Notch, Wnt ir" ~ italic("ARID1A") ~ "genų raiškos kombinacija"),
                             xlab = "Bendras išgyvenamumo laikas",
                             ylab = "Išgyvenamumo tikimybė",
                             palette = c("turquoise", "deeppink"),  # Color palette for groups
                             legend.title = "Rizikos grupė", 
                             legend.labs = c("Mažas rizikos balas", "Didelis rizikos balas"))
# Add subtitle form cox result
test_survplot14$plot <- test_survplot14$plot +
  labs(subtitle = "Uni PR =  0.22  (95 % PI: 0,07–0,66); Multi PR = 0.22 (95 % PI: 0,06–0,77), Long-rank p = 0,0035, n= 47")

print(test_survplot14)
#save
png("KM_plot_14_OC_w_HR20250925.png", width = 1800, height = 1100, res = 150)
print(test_survplot14)  # print the full ggsurvplot object
dev.off()

#TIME ROC, NOTCH, OC #########################################
surv_df_notch_oc <- OC_SURV_EXPRESSION[, colnames(OC_SURV_EXPRESSION) %in%
                                         c("OS", "STATUS", genes_notch, methylation, "RiskScore_notch", "RiskScore_14" , "patient_id_aud")]

##time rocs for separate NOTCH biomarkers######################################
coxdf2_notch_oc <- surv_df_notch_oc[, (colnames(surv_df_notch_oc) %in% c(genes_notch, methylation) )]
dim(coxdf2_notch_oc)
#time roc
t_eval <- c(12, 36, 60)  # time points
#make surf df but only of my genes!
rez_list2_notch_OC <- apply(coxdf2_notch_oc, 2, timeROC,
                            T = surv_df_notch_oc$OS,       # Survival time from df
                            delta =  surv_df_notch_oc$STATUS,# Event indicator from df
                            #marker  # Predictor already in the df
                            cause = 1,         # Event of interest
                            times = t_eval,    # Time points for ROC
                            iid = TRUE )        # Compute confidence intervals)

auc_table_notch_oc <- map_dfr(names(rez_list2_notch_OC), function(gene) {
  roc <- rez_list2_notch_OC[[gene]]
  
  tibble(
    gene = gene,
    time = roc$times,
    cases = roc$cases,
    survivors = roc$survivors,
    censored = roc$censored,
    auc = roc$AUC,
    se = roc$inference$vect_sd_1
  )
})

View(auc_table_notch_oc)

##COORDS, seprate notch genes #####################
# Build table of AUC + Sensitivity + Specificity at 60 months
sens_spec_auc_60 <- map_dfr(names(rez_list2_notch_OC), function(gene) {
  roc <- rez_list2_notch_OC[[gene]]
  
  # Get index for t=60
  idx_60 <- which(roc$times == 60)
  
  # Sensitivity & specificity vectors across thresholds
  sens_60 <- roc$TP[, idx_60]
  spec_60 <- 1 - roc$FP[, idx_60]
  
  # Use Youden index to select optimal cutoff
  youden <- sens_60 + spec_60 - 1
  best_idx <- which.max(youden)
  
  tibble(
    gene = gene,
    time = roc$times[idx_60],
    auc = roc$AUC[idx_60]*100,
    sens = sens_60[best_idx],
    spec = spec_60[best_idx],
    cutoff = roc$cutoffs[best_idx]
  )
})

View(sens_spec_auc_60)


##time roc, all notch genes ###################
roc_result_model_notch <- timeROC(
  T = surv_df_notch_oc$OS,       # Survival time from df
  delta = surv_df_notch_oc$STATUS, # Event indicator from df
  marker = surv_df_notch_oc$RiskScore_notch, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = t_eval,    # Time points for ROC
  iid = TRUE         # Compute confidence intervals
)
roc_result_model_notch

##COORDS time roc, notch model#############################
# Sensitivity and specificity vectors at 60 months
sens_60 <- roc_result_model_notch$TP[, which(roc_result_model_notch$times == 60)]
spec_60 <- 1 - roc_result_model_notch$FP[, which(roc_result_model_notch$times == 60)]

# Youden index
youden_60 <- sens_60 + spec_60 - 1
best_index_60 <- which.max(youden_60)

# Extract best specificity
best_spec_60 <- spec_60[best_index_60]
best_sens_60 <- sens_60[best_index_60]
best_cutoff_60 <- roc_result_model_notch$cutoffs[best_index_60]

cat("At 60 months:\n")
cat("  Specificity:", round(best_spec_60, 3), "\n")
cat("  Sensitivity:", round(best_sens_60, 3), "\n")
cat("  Cutoff:", best_cutoff_60, "\n")

##time roc, all notch, methylation genes ###################
roc_result_model_notch_met<- timeROC(
  T = surv_df_notch_oc$OS,       # Survival time from df
  delta = surv_df_notch_oc$STATUS, # Event indicator from df
  marker = surv_df_notch_oc$RiskScore_14, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = t_eval,    # Time points for ROC
  iid = TRUE         # Compute confidence intervals
)
roc_result_model_notch_met

##COORDS time roc, notch methylation model#############################
# Sensitivity and specificity vectors at 60 months
sens_60 <- roc_result_model_notch_met$TP[, which(roc_result_model_notch_met$times == 60)]
spec_60 <- 1 - roc_result_model_notch_met$FP[, which(roc_result_model_notch_met$times == 60)]

# Youden index
youden_60 <- sens_60 + spec_60 - 1
best_index_60 <- which.max(youden_60)

# Extract best specificity
best_spec_60 <- spec_60[best_index_60]
best_sens_60 <- sens_60[best_index_60]
best_cutoff_60 <- roc_result_model_notch$cutoffs[best_index_60]

cat("At 60 months:\n")
cat("  Specificity:", round(best_spec_60, 3), "\n")
cat("  Sensitivity:", round(best_sens_60, 3), "\n")
cat("  Cutoff:", best_cutoff_60, "\n")

##plot at year 1  with plotting###############
# Choose target time
target_time <- 12        
time_index <- which(rez_list2_notch_OC[[1]]$times == target_time)

png("C:/Users/Ieva/rprojects/outputs_all/DISS/tissues_OC_timeROC_12_test20250924.png",
    width = 1500, height = 1200, res = 200) # width and height in pixels, resolution in dpi
# Base ROC curve plot (same as before)
par(pty="s")
plot(
  rez_list2_notch_OC[[1]]$FP[, time_index],
  rez_list2_notch_OC[[1]]$TP[, time_index],
  type = "l",
  col = 1,
  lwd = 2,
  xlab = "Specifiškumas",
  ylab = "Jautrumas",
  main = paste("Nuo laiko priklausomos ROC kreivės,\n",
               "1 metai po diagnozės audinių imtyje"),
  xlim = c(0, 1),
  ylim = c(0, 1),
  asp = 1
)

# Add ROC lines for all genes
for (i in 2:length(rez_list2_notch_OC)) {
  lines(
    rez_list2_notch_OC[[i]]$FP[, time_index],
    rez_list2_notch_OC[[i]]$TP[, time_index],
    col = i,
    lwd = 2
  )
}

# Add risk score ROC line (Notch)
lines(
  roc_result_model_notch$FP[, time_index],
  roc_result_model_notch$TP[, time_index],
  col = "maroon",
  lwd = 3,
  lty = 1
)

# Add second risk score ROC line (14 biomarkers)
lines(
  roc_result_model_notch_met$FP[, time_index],
  roc_result_model_notch_met$TP[, time_index],
  col = "darkgreen",
  lwd = 3,
  lty = 2
)

# Add diagonal reference line
abline(0, 1, lty = 2, col = "gray")

# --- Legend preparation ---

# Get AUCs for each gene at time_index
auc_list <- sapply(rez_list2_notch_OC, function(x) x$AUC[time_index])
auc_risk_notch <- roc_result_model_notch$AUC[time_index]
auc_risk_met   <- roc_result_model_notch_met$AUC[time_index]

# Build gene labels with italic names and AUCs
legend_labels <- mapply(function(name, auc) {
  paste0("italic('", name, "')~'(AUC = ", sprintf("%.3f", auc), ")'")
}, names(rez_list2_notch_OC), auc_list)

# Add risk score labels with AUC
legend_labels <- c(
  legend_labels,
  paste0("'Genų raiškos rizikos balas (AUC = ", sprintf("%.3f", auc_risk_notch), ")'"),
  paste0("'14 biožymenų rizikos balas (AUC = ", sprintf("%.3f", auc_risk_met), ")'")
)

# Add legend
legend(
  "bottomright",
  legend = parse(text = legend_labels),
  col = c(1:length(rez_list2_notch_OC), "maroon", "darkgreen"),
  lwd = c(rep(2, length(rez_list2_notch_OC)), 3, 3),
  lty = c(rep(1, length(rez_list2_notch_OC)), 1, 2),
  cex = 0.6,
  bty = "n"
)

# Optional: Add panel label "D" above title
#mtext("D", side = 3, line = 2.5, adj = -0.2, font = 2, cex = 1.5)

#run plot
dev.off() # Close the PNG device

##plot at year 3  with plotting###############
# Choose target time
target_time <- 36        
time_index <- which(rez_list2_notch_OC[[1]]$times == target_time)

png("C:/Users/Ieva/rprojects/outputs_all/DISS/tissues_OC_timeROC_36_test20250924.png",
    width = 1500, height = 1200, res = 200) # width and height in pixels, resolution in dpi
# Set up base plot with gene 1
par(pty="s")
# Set up base plot with gene 1
plot(
  rez_list2_notch_OC[[1]]$FP[, time_index],
  rez_list2_notch_OC[[1]]$TP[, time_index],
  type = "l",
  col = 1,
  lwd = 2,
  xlab = "Specifiškumas",
  ylab = "Jautrumas",
  main = paste("Nuo laiko priklausomos ROC kreivės,
3 metai po diagnozės audinių imtyje"),
  xlim = c(0, 1),
  ylim = c(0, 1),
  asp = 1
)

# Add ROC lines for all genes
for (i in 2:length(rez_list2_notch_OC)) {
  lines(
    rez_list2_notch_OC[[i]]$FP[, time_index],
    rez_list2_notch_OC[[i]]$TP[, time_index],
    col = i,
    lwd = 2
  )
}

# Add risk score ROC line in bold black
lines(
  roc_result_model_notch$FP[, time_index],
  roc_result_model_notch$TP[, time_index],
  col = "maroon",
  lwd = 3,
  lty = 1
)

# Add diagonal reference line
abline(0, 1, lty = 2, col = "gray")

# Build legend names: italic gene names + "risk score"
legend_labels <- c(
  parse(text = paste0("italic('", names(rez_list2_notch_OC), "')")),
  "Risk Score"
)

# Get AUCs for each gene at time_index
auc_list <- sapply(rez_list2_notch_OC, function(x) x$AUC[time_index])
auc_risk <- roc_result_model_notch$AUC[time_index]

# Build gene labels with italic names and AUCs
legend_labels <- mapply(function(name, auc) {
  paste0("italic('", name, "')~'(AUC = ", sprintf("%.3f", auc), ")'")
}, names(rez_list2_notch_OC), auc_list)

# Add risk score with AUC
legend_labels <- c(legend_labels,
                   paste0("'Risk Score (AUC = ", sprintf("%.3f", auc_risk), ")'"))

# Add legend
legend(
  "bottomright",
  legend = parse(text = legend_labels),
  col = c(1:length(rez_list2_notch_OC), "maroon"),
  lwd = c(rep(2, length(rez_list2_notch_OC)), 3),
  cex = 0.6,
  bty = "n"
)
#add A label during png
#mtext("D", side = 3, line = 2.5, adj = -0.2, font = 2, cex = 1.5)

#run plot
dev.off() # Close the PNG device

##plot at year 5  with plotting###############
# Choose target time
target_time <- 60        
time_index <- which(rez_list2_notch_OC[[1]]$times == target_time)

png("C:/Users/Ieva/rprojects/outputs_all/DISS/tissues_OCtimeROC_test20250924x.png",
    width = 1800, height = 1800, res = 220) # width and height in pixels, resolution in dpi
# Set up base plot with gene 1
par(pty="s")
# Set up base plot with gene 1
plot(
  rez_list2_notch_OC[[1]]$FP[, time_index],
  rez_list2_notch_OC[[1]]$TP[, time_index],
  type = "l",
  col = 1,
  lwd = 2,
  xlab = "Specifiškumas",
  ylab = "Jautrumas",
  main = paste("Nuo laiko priklausomos ROC kreivės,
5 metai po diagnozės audinių imtyje"),
  xlim = c(0, 1),
  ylim = c(0, 1),
  asp = 1
)

# Add ROC lines for all genes
for (i in 2:length(rez_list2_notch_OC)) {
  lines(
    rez_list2_notch_OC[[i]]$FP[, time_index],
    rez_list2_notch_OC[[i]]$TP[, time_index],
    col = i,
    lwd = 2
  )
}

# Add risk score ROC line in bold black
lines(
  roc_result_model_notch$FP[, time_index],
  roc_result_model_notch$TP[, time_index],
  col = "maroon",
  lwd = 3,
  lty = 1
)

# Add second risk score ROC line (14 biomarkers)
lines(
  roc_result_model_notch_met$FP[, time_index],
  roc_result_model_notch_met$TP[, time_index],
  col = "darkgreen",
  lwd = 3,
  lty = 2
)

# Add diagonal reference line
abline(0, 1, lty = 2, col = "gray")

# Get AUCs for each gene at time_index
auc_list <- sapply(rez_list2_notch_OC, function(x) x$AUC[time_index])
auc_risk_notch <- roc_result_model_notch$AUC[time_index]
auc_risk_met   <- roc_result_model_notch_met$AUC[time_index]

# Build gene labels with italic names and AUCs
legend_labels <- mapply(function(name, auc) {
  paste0("italic('", name, "')~'(AUC = ", sprintf("%.3f", auc), ")'")
}, names(rez_list2_notch_OC), auc_list)

# Add risk score labels with AUC
legend_labels <- c(
  legend_labels,
  paste0("'Genų raiškos rizikos balas (AUC = ", sprintf("%.3f", auc_risk_notch), ")'"),
  paste0("'14 biožymenų rizikos balas (AUC = ", sprintf("%.3f", auc_risk_met), ")'")
)

# Add legend
legend(
  "bottomright",
  legend = parse(text = legend_labels),
  col = c(1:length(rez_list2_notch_OC), "maroon", "darkgreen"),
  lwd = c(rep(2, length(rez_list2_notch_OC)), 3, 3),
  lty = c(rep(1, length(rez_list2_notch_OC)), 1, 2),
  cex = 0.6,
  bty = "n"
)
#run plot
dev.off() # Close the PNG device

#SAME 5 year plot, no AUCS#######################
# Choose target time
target_time <- 60        
time_index <- which(rez_list2_notch_OC[[1]]$times == target_time)

png("C:/Users/Ieva/rprojects/outputs_all/DISS/tissuesNO_AUCS_OCtimeROC_test20250924.png",
    width = 1500, height = 1500, res = 200) # width and height in pixels, resolution in dpi
# Set up base plot with gene 1
# Base ROC curve plot
par(pty="s")
plot(
  rez_list2_notch_OC[[1]]$FP[, time_index],
  rez_list2_notch_OC[[1]]$TP[, time_index],
  type = "l",
  col = 1,
  lwd = 2,
  xlab = "Specifiškumas",
  ylab = "Jautrumas",
  main = paste("Nuo laiko priklausomos ROC kreivės,\n",
               "5 metai po diagnozės audinių imtyje"),
  xlim = c(0, 1),
  ylim = c(0, 1),
  asp = 1
)

# Add ROC lines for all genes
for (i in 2:length(rez_list2_notch_OC)) {
  lines(
    rez_list2_notch_OC[[i]]$FP[, time_index],
    rez_list2_notch_OC[[i]]$TP[, time_index],
    col = i,
    lwd = 2
  )
}

# Add risk score ROC line (Notch)
lines(
  roc_result_model_notch$FP[, time_index],
  roc_result_model_notch$TP[, time_index],
  col = "maroon",
  lwd = 3,
  lty = 1
)

# Add second risk score ROC line (14 biomarkers)
lines(
  roc_result_model_notch_met$FP[, time_index],
  roc_result_model_notch_met$TP[, time_index],
  col = "darkgreen",
  lwd = 3,
  lty = 2
)

# Add diagonal reference line
abline(0, 1, lty = 2, col = "gray")

# --- Legend preparation without AUCs ---

# Gene names
legend_labels <- names(rez_list2_notch_OC)

# Add risk score labels
legend_labels <- c(
  legend_labels,
  "Genų raiškos rizikos balas",
  "14 biožymenų rizikos balas"
)

# Convert labels to italic expressions
legend_labels_italic <- lapply(legend_labels, function(x) bquote(italic(.(x))))

# Add legend
legend(
  "bottomright",
  legend = legend_labels_italic,
  col = c(1:length(rez_list2_notch_OC), "maroon", "darkgreen"),
  lwd = c(rep(2, length(rez_list2_notch_OC)), 3, 3),
  lty = c(rep(1, length(rez_list2_notch_OC)), 1, 2),
  cex = 0.6,
  bty = "n"
)

# Optional: Add panel label "D" above title
#mtext("D", side = 3, line = 2.5, adj = -0.2, font = 2, cex = 1.5)
dev.off() # Close the PNG device

#TABLE of aucs, sensitivity, specificity##########################
# Extract best sensitivity, specificity, cutoff for Notch risk score at 60 months
sens_60_risk <- roc_result_model_notch$TP[, which(roc_result_model_notch$times == 60)]
spec_60_risk <- 1 - roc_result_model_notch$FP[, which(roc_result_model_notch$times == 60)]

youden_risk <- sens_60_risk + spec_60_risk - 1
best_idx_risk <- which.max(youden_risk)

risk_score_row <- tibble(
  gene = "Genų raiškos rizikos balas",
  time = roc_result_model_notch$times[which(roc_result_model_notch$times == 60)],
  auc = roc_result_model_notch$AUC[which(roc_result_model_notch$times == 60)]*100,
  sens = sens_60_risk[best_idx_risk],
  spec = spec_60_risk[best_idx_risk],
  cutoff = roc_result_model_notch$cutoffs[best_idx_risk]
)

# Combine with gene-level table
sens_spec_auc_60_all <- bind_rows(sens_spec_auc_60, risk_score_row)

View(sens_spec_auc_60_all)

# Extract best sensitivity, specificity, cutoff for 14-biomarker risk score at 60 months
sens_60_met <- roc_result_model_notch_met$TP[, which(roc_result_model_notch_met$times == 60)]
spec_60_met <- 1 - roc_result_model_notch_met$FP[, which(roc_result_model_notch_met$times == 60)]

youden_60_met <- sens_60_met + spec_60_met - 1
best_idx_met <- which.max(youden_60_met)

risk_score_met_row <- tibble(
  gene = "14 biožymenų rizikos balas",
  time = roc_result_model_notch_met$times[which(roc_result_model_notch_met$times == 60)],
  auc = roc_result_model_notch_met$AUC[which(roc_result_model_notch_met$times == 60)]*100,
  sens = sens_60_met[best_idx_met],
  spec = spec_60_met[best_idx_met],
  cutoff = roc_result_model_notch_met$cutoffs[best_idx_met]
)

# Combine with previous table (genes + Notch risk score)
sens_spec_auc_60_all <- bind_rows(sens_spec_auc_60_all, risk_score_met_row)

View(sens_spec_auc_60_all)

#make prety table
# Prepare table: rename columns for display
roc_table_display <- sens_spec_auc_60_all %>%
  rename(
    Biožymuo = gene,
    Laikas_mėn = time,
    AUC = auc,
    Jautrumas = sens,
    Specifiškumas = spec
  )

# Create gt table
gt_table_roc_60 <- roc_table_display %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai",
    subtitle = "Prognostinis tikslumas 5 metų išgyvenimui"
  ) %>%
  fmt_number(
    columns = vars(AUC, Jautrumas, Specifiškumas),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymuo))
  )

# Show table
gt_table_roc_60

#there is no other convenient way to save gt outputs
gtsave(gt_table_roc_60,
       filename = "timeroc_table_20250924.png")

#Combine the images
roc_image <- image_read("tissuesNO_AUCS_OCtimeROC_test20250924.png")
table_image <- image_read("timeroc_table_20250924.png")

combined_image <- image_append(c(roc_image, table_image), stack = F)

# Save the combined image
image_write(combined_image, 
            "roc_notch_all_oc_with_table_20250924.png")


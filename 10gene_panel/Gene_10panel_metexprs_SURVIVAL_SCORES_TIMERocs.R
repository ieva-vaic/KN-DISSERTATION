#KN-DISSERTATION project. 10 gene panel - genes selected via statistical analysis of TCGA and GTEx data
#KN 10 gene analysis data (cv 2 %, one threshold, 65 cases!, restored KN-59) 2025-02-17
#KN-DISSERTATION project. MET_EXPRS - methylation, expression, mutation data
#survival analysis, time rocs and models
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
library(gt)
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

#ALL 10 GENES#####################################################
#gene data - all samples
gene_data <- ALL_SURV_EXPRESSION[, colnames(ALL_SURV_EXPRESSION) %in% genes10]
# Example with 10 biomarkers
cox_model_ALL <- coxph(
  Surv(OS, STATUS) ~ EXO1 + RAD50 + PPT2 + LUC7L2 + PKP3 + CDCA5 + ZFPL1 + VPS33B + GRB7 + TCEAL4,
  data = ALL_SURV_EXPRESSION
)
summary(cox_model_ALL)
#get coeficients
coefs <- coef(cox_model_ALL)
# Risk score = sum( gene_expression * coefficient )
risk_scores_test <- rowSums(sweep(gene_data, 2, coefs, "*"))
# View the risk scores
print(risk_scores_test) # now I have some risk scores
#add risk scores to the clin_df_joined_test
ALL_SURV_EXPRESSION$RiskScore <- risk_scores_test
#create df wih survival data
surv_df_test <- ALL_SURV_EXPRESSION[, colnames(ALL_SURV_EXPRESSION) %in%
                                      c("OS", "STATUS", genes10, "RiskScore", "patient_id_aud", "Age", "CA125")]

rownames(surv_df_test) <- surv_df_test$patient_id_aud

# Calculate the median risk score
median_risk <- median(surv_df_test$RiskScore, na.rm = TRUE) #-7.700461
# Create a new factor column based on the median value
surv_df_test$RiskGroup <- ifelse(surv_df_test$RiskScore <= median_risk,
                                 "Low Risk", "High Risk")
#Create a survival object
surv_object <- Surv(time = surv_df_test$OS,
                    event = surv_df_test$STATUS )

# Fit a Kaplan-Meier model
km_fit <- survfit(surv_object ~ RiskGroup, data = surv_df_test)

# Fit Univariable cox model for 10 gene model 
cox_model <- coxph(Surv(OS, STATUS) ~RiskGroup , data = surv_df_test)
summary(cox_model)
# Fit Multivariable cox model for  10 gene model
multicox_model <- coxph(Surv(OS, STATUS) ~ RiskGroup + Age + CA125, data = surv_df_test)
summary(multicox_model)

# Plot the Kaplan-Meier curve using ggsurvplot
test_survplot<- ggsurvplot(km_fit, data = surv_df_test, 
                           #pval = TRUE,  # Show p-value of the log-rank test
                           risk.table = TRUE,  # Add risk table below the plot
                           risk.table.title = "Pacienčių skaičius rizikos grupėje",
                           title = paste0("Didelės vs. mažos rizikos atvejai KV imtyje,\n",
                           "10 genų raiškos kombinacija"),
                           xlab = "Bendras išgyvenamumo laikas",
                           ylab = "Išgyvenamumo tikimybė",
                           palette = c("turquoise", "deeppink"),  # Color palette for groups
                           legend.title = "Rizikos grupė", 
                           legend.labs = c("Mažos rizikos balas", "Didelės rizikos balas"))#2026.01.20 changed from "mažas rizikos balas" and "didelis rizikos balas"
# Add subtitle form cox result
test_survplot$plot <- test_survplot$plot +
  labs(
    subtitle = paste0(
      "Uni PR = 0,26 (95% PI: 0,08–0.85)\n",
      "Multi PR = 0,24 (95% PI: 0,06–1.03)\n",
      "Log-rank p = 0,017, n = 47"
    )
  )
print(test_survplot)
#save
png("C:/Users/Ieva/rprojects/outputs_all/DISS/KM_10_gene_sig_20260120LT.png",
    width = 1100, height = 600, res = 100) # width and height in pixels, resolution in dpi
test_survplot #
dev.off() # Close the PNG device

#ALL 10 GENES, HGSOC#####################################################
#gene dataH - only HGSOC
gene_dataH <- HGSOC_SURV_EXPRESSION[, colnames(HGSOC_SURV_EXPRESSION) %in% genes10]
# Example with 10 biomarkers
cox_model_HGSOC <- coxph(
  Surv(OS, STATUS) ~ EXO1 + RAD50 + PPT2 + LUC7L2 + PKP3 + CDCA5 + ZFPL1 + VPS33B + GRB7 + TCEAL4,
  data = HGSOC_SURV_EXPRESSION
)
summary(cox_model_HGSOC)
#get coeficients
coefs <- coef(cox_model_HGSOC)
# Risk score = sum( gene_expression * coefficient )
risk_scores_test_h <- rowSums(sweep(gene_dataH, 2, coefs, "*"))
# View the risk scores
print(risk_scores_test_h) # now I have some risk scores
#add risk scores to the clin_df_joined_test
HGSOC_SURV_EXPRESSION$RiskScore <- risk_scores_test_h
#create df wih survival data
surv_df_test_h <- HGSOC_SURV_EXPRESSION[, colnames(HGSOC_SURV_EXPRESSION) %in%
                                          c("OS", "STATUS", genes10, "RiskScore", "patient_id_aud")]

rownames(surv_df_test_h) <- surv_df_test_h$patient_id_aud

# Calculate the median risk score
median_risk_h <- median(surv_df_test_h$RiskScore, na.rm = TRUE) #-7.700461
# Create a new factor column based on the median value
surv_df_test_h$RiskGroup <- ifelse(surv_df_test_h$RiskScore <= median_risk_h,
                                   "Low Risk", "High Risk")
#Create a survival object
surv_object <- Surv(time = surv_df_test_h$OS,
                    event = surv_df_test_h$STATUS )

# Fit a Kaplan-Meier model
km_fitH <- survfit(surv_object ~ RiskGroup, data = surv_df_test_h)
# Plot the Kaplan-Meier curve using ggsurvplot
test_survplotH <- ggsurvplot(km_fitH, data = surv_df_test_h, 
                             pval = TRUE,  # Show p-value of the log-rank test
                             risk.table = TRUE,  # Add risk table below the plot
                             title = "Kaplan-Meier kreivė: Didelės vs. mažos rizikos atvejai HGSOC audinių imtyje",
                             xlab = "Bendras išgyvenamumo laikas",
                             ylab = "Išgyvenamumo tikimybė",
                             palette = c("turquoise", "deeppink"),  # Color palette for groups
                             legend.title = "Rizikos grupė", 
                             legend.labs = c("Mažas rizikos balas", "Didelis rizikos balas"))
test_survplotH
#save
png("C:/Users/Ieva/rprojects/outputs_all/DISS/KM_10_gene_hgsocsig_20260121.png",
    width = 800, height = 600, res = 100) # width and height in pixels, resolution in dpi
test_survplotH #
dev.off() # Close the PNG device


#ALL 10 GENES, OC ONLY#####################################################
#gene_data2 - only OC
gene_data2 <- OC_SURV_EXPRESSION[, colnames(OC_SURV_EXPRESSION) %in% genes10]
# Example with 10 biomarkers
cox_model_OC <- coxph(
  Surv(OS, STATUS) ~ EXO1 + RAD50 + PPT2 + LUC7L2 + PKP3 + CDCA5 + ZFPL1 + VPS33B + GRB7 + TCEAL4,
  data = OC_SURV_EXPRESSION
)
summary(cox_model_OC)
#get coeficients
coefs2 <- coef(cox_model_OC)
# Risk score = sum( gene_expression * coefficient )
risk_scores_test2 <- rowSums(sweep(gene_data2, 2, coefs2, "*"))
# View the risk scores
print(risk_scores_test2) # now I have some risk scores
#add risk scores to the clin_df_joined_test
OC_SURV_EXPRESSION$RiskScore <- risk_scores_test2
#create df wih survival data
surv_df_test2 <- OC_SURV_EXPRESSION[, colnames(OC_SURV_EXPRESSION) %in%
                                      c("OS", "STATUS", genes10, "RiskScore", "patient_id_aud", "Age", "CA125")]

rownames(surv_df_test2) <- surv_df_test2$patient_id_aud

# Calculate the median risk score
median_risk2 <- median(surv_df_test2$RiskScore, na.rm = TRUE) #-7.396536
# Create a new factor column based on the median value
surv_df_test2$RiskGroup <- ifelse(surv_df_test2$RiskScore <= median_risk2,
                                  "Low Risk", "High Risk")
#Create a survival object
surv_object2 <- Surv(time = surv_df_test2$OS,
                     event = surv_df_test2$STATUS )

# Fit a Kaplan-Meier model
km_fit2 <- survfit(surv_object2 ~ RiskGroup, data = surv_df_test2)


# Fit Univariable cox model for 10 gene model 
cox_model2 <- coxph(Surv(OS, STATUS) ~RiskGroup , data = surv_df_test2)
summary(cox_model2)
# Fit Multivariable cox model for  10 gene model
multicox_model2 <- coxph(Surv(OS, STATUS) ~ RiskGroup + Age + CA125, data = surv_df_test2)
summary(multicox_model2)

# Plot the Kaplan-Meier curve using ggsurvplot
test_survplot2 <- ggsurvplot(km_fit2, data = surv_df_test2, 
                             #pval = TRUE,  # Show p-value of the log-rank test
                             risk.table = TRUE,  # Add risk table below the plot
                             risk.table.title = "Pacienčių skaičius rizikos grupėje",
                             title = paste0("Didelės vs. mažos rizikos atvejai KV imtyje,\n",
                                            "10 genų raiškos kombinacija"),
                             xlab = "Bendras išgyvenamumo laikas, mėnesiais",
                             ylab = "Išgyvenamumo tikimybė",
                             palette = c(  "deeppink", "turquoise"),  # Color palette for groups
                             legend.title = "Rizikos grupė", 
                             legend.labs = c( "Didelės rizikos balas", "Mažos rizikos balas")) #2026.01.20 changed from "didelis / mažas rizikos balas"
# Add subtitle form cox result
test_survplot2$plot <- test_survplot2$plot +
  labs(subtitle = paste0("Uni PR =  0,24  (95 % PI: 0,07–0,81);\n",
                         "Multi PR = 0,21 (95 % PI: 0,05–0,91);\n",
                         "Log-rank p = 0,013, n = 47"))

test_survplot2
#save
png("C:/Users/Ieva/rprojects/outputs_all/DISS/KM_10_gene_OCa_20260120short.png",
    width = 800, height = 620, res = 120) # width and height in pixels, resolution in dpi
test_survplot2 #
dev.off() # Close the PNG device

#TIME ROC, Risk score, OC #######################
surv_df_test_x <- surv_df_test2[!is.na(surv_df_test2$RiskScore),] #remove NA
#create features for timeroc
nobs <- NROW(surv_df_test_x)
#make surf df but only of my genes!
time <- surv_df_test_x$OS
event <- surv_df_test_x$STATUS

#time roc, risk score
t_eval <- c(12, 36, 60)  # time points
roc_result <- timeROC(
  T = time,       # Survival time from df
  delta = event, # Event indicator from df
  marker = surv_df_test_x$RiskScore, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = t_eval,    # Time points for ROC
  iid = TRUE         # Compute confidence intervals
)
roc_result

##time rocs for separate biomarkers, ALL ######################################
coxdf <- surv_df_test_x[, (colnames(surv_df_test_x) %in% genes10)]
dim(coxdf)

rez_list <- apply(coxdf, 2, timeROC,
                  T = time,       # Survival time from df
                  delta = event, # Event indicator from df
                  #marker  # Predictor already in the df
                  cause = 1,         # Event of interest
                  times = t_eval,    # Time points for ROC
                  iid = TRUE )        # Compute confidence intervals)

auc_table <- map_dfr(names(rez_list), function(gene) {
  roc <- rez_list[[gene]]
  
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

auc_table

##plot at year 5  with plotting, no AUCS in the plot###############
# Choose target time
target_time <- 60        
time_index <- which(rez_list[[1]]$times == target_time)

png("C:/Users/Ieva/rprojects/outputs_all/DISS/tissues_OCtimeROC_test20260121LT.png",
    width = 1500, height = 1200, res = 200) # width and height in pixels, resolution in dpi
# Set up base plot with gene 1
par(pty="s")
# Set up base plot with gene 1
plot(
  rez_list[[1]]$FP[, time_index],
  rez_list[[1]]$TP[, time_index],
  type = "l",
  col = 1,
  lwd = 2,
  xlab = "1 - Specifiškumas",
  ylab = "Jautrumas",
  main = paste("Nuo laiko priklausomos ROC kreivės,
5 metai po diagnozės audinių imtyje"),
  xlim = c(0, 1),
  ylim = c(0, 1),
  asp = 1
)

# Add ROC lines for all genes
for (i in 2:length(rez_list)) {
  lines(
    rez_list[[i]]$FP[, time_index],
    rez_list[[i]]$TP[, time_index],
    col = i,
    lwd = 2
  )
}

# Add risk score ROC line in bold black
lines(
  roc_result$FP[, time_index],
  roc_result$TP[, time_index],
  col = "maroon",
  lwd = 3,
  lty = 1
)

# Add diagonal reference line
abline(0, 1, lty = 2, col = "gray")

# Build legend names: italic gene names + "risk score"
legend_labels <- c(
  parse(text = paste0("italic('", names(rez_list), "')")),
  "Risk Score"
)

# Get AUCs for each gene at time_index
auc_list <- sapply(rez_list, function(x) x$AUC[time_index])
auc_risk <- roc_result$AUC[time_index]

# Build gene labels with italic names and AUCs
legend_labels <- mapply(function(name, auc) {
  paste0("italic('", name, "')~'(AUC = ", sprintf("%.3f", auc), ")'")
}, names(rez_list), auc_list)

# Add risk score with AUC
legend_labels <- c(legend_labels,
                   paste0("'Rizikos balas (AUC = ", sprintf("%.3f", auc_risk), ")'"))

# Add legend
legend(
  "bottomright",
  legend = parse(text = legend_labels),
  col = c(1:length(rez_list), "maroon"),
  lwd = c(rep(2, length(rez_list)), 3),
  cex = 0.6,
  bty = "n"
)

#run plot
dev.off() # Close the PNG device

#5 year plot with table #####################
target_time <- 60        
time_index <- which(rez_list[[1]]$times == target_time)

png("C:/Users/Ieva/rprojects/outputs_all/DISS/tissues_OCtimeROC_test20260122LT.png",
    width = 1100, height = 1100, res = 200) # width and height in pixels, resolution in dpi
# Set up base plot with gene 1
# Base ROC curve plot
par(pty="s")
plot(
  rez_list[[1]]$FP[, time_index],
  rez_list[[1]]$TP[, time_index],
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
for (i in 2:length(rez_list)) {
  lines(
    rez_list[[i]]$FP[, time_index],
    rez_list[[i]]$TP[, time_index],
    col = i,
    lwd = 2
  )
}

# Add risk score ROC line 
lines(
  roc_result$FP[, time_index],
  roc_result$TP[, time_index],
  col = "maroon",
  lwd = 3,
  lty = 1
)


# Add diagonal reference line
abline(0, 1, lty = 2, col = "gray")

# --- Legend preparation without AUCs ---

# Gene names
legend_labels <- names(rez_list)

# Add risk score labels
legend_labels <- c(
  legend_labels,
  "Rizikos balas"
)

# Convert labels to italic expressions
legend_labels_italic <- lapply(legend_labels, function(x) bquote(italic(.(x))))

# Add legend
legend(
  "bottomright",
  legend = legend_labels_italic,
  col = c(1:length(rez_list), "maroon"),
  lwd = c(rep(2, length(rez_list)), 3, 3),
  lty = c(rep(1, length(rez_list)), 1, 2),
  cex = 0.6,
  bty = "n"
)
dev.off() # Close the PNG device


#TABLE of aucs, sensitivity, specificity##########################
# Extract best sensitivity, specificity, cutoff for Notch risk score at 60 months
sens_60_risk <- roc_result$TP[, which(roc_result$times == 60)]
spec_60_risk <- 1 - roc_result$FP[, which(roc_result$times == 60)]

youden_risk <- sens_60_risk + spec_60_risk - 1
best_idx_risk <- which.max(youden_risk)

risk_score_row <- tibble(
  gene = "10 Genų raiškos rizikos balas",
  time = roc_result$times[which(roc_result$times == 60)],
  auc = roc_result$AUC[which(roc_result$times == 60)],
  sens = sens_60_risk[best_idx_risk],
  spec = spec_60_risk[best_idx_risk],
  cutoff = roc_result$cutoffs[best_idx_risk]
)

# Build table of AUC + Sensitivity + Specificity at 60 months
sens_spec_auc_60 <- map_dfr(names(rez_list), function(gene) {
  roc <- rez_list[[gene]]
  
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
    auc = roc$AUC[idx_60],
    sens = sens_60[best_idx],
    spec = spec_60[best_idx],
    cutoff = roc$cutoffs[best_idx]
  )
})

# Combine with gene-level table
sens_spec_auc_60_all <- bind_rows(sens_spec_auc_60, risk_score_row)

#make pretty table
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
       filename = "10genetimeroc_table_20260121.png")

#Combine the images
roc_image <- image_read("tissues_OCtimeROC_test20260122LT.png")
table_image <- image_read("10genetimeroc_table_20260121.png")

combined_image <- image_append(c(roc_image, table_image), stack = F)

# Save the combined image
image_write(combined_image, 
            "tissues_OCtimeROCw_table_test20260121.png")



#USE COEFS FROM TCGA FOR RISK SCORE#####################
coef_tcga <- readRDS("C:/Users/Ieva/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/coefs.RDS")

##OC only, TCGA COEFS##################################
# Risk score = sum( gene_expression * coefficient )
risk_scores_tcga_type2 <- rowSums(sweep(gene_data2, 2, coef_tcga, "*"))
# View the risk scores
print(risk_scores_tcga_type2) # now I have some risk scores
#add risk scores to the clin_df_joined_test
OC_SURV_EXPRESSION$RiskScore_tcga <- risk_scores_tcga_type2
#create df wih survival data
surv_df_tcga_type2 <- OC_SURV_EXPRESSION[, colnames(OC_SURV_EXPRESSION) %in%
                                           c("OS", "STATUS", genes10, "RiskScore_tcga", "patient_id_aud")]

rownames(surv_df_tcga_type2) <- surv_df_tcga_type2$patient_id_aud

# Calculate the median risk score
median_risktcga2 <- median(surv_df_tcga_type2$RiskScore_tcga, na.rm = TRUE) #-7.396536
# Create a new factor column based on the median value
surv_df_tcga_type2$RiskGrouptcga <- ifelse(surv_df_tcga_type2$RiskScore_tcga <= median_risktcga2,
                                           "Low Risk", "High Risk")
#Create a survival object
surv_objecttcga2 <- Surv(time = surv_df_tcga_type2$OS,
                         event = surv_df_tcga_type2$STATUS )

# Fit a Kaplan-Meier model
km_fitctga2 <- survfit(surv_objecttcga2 ~ RiskGrouptcga, data = surv_df_tcga_type2)
# Plot the Kaplan-Meier curve using ggsurvplot
test_survplotctga2 <- ggsurvplot(km_fitctga2, data = surv_df_tcga_type2, 
                                 pval = TRUE,  # Show p-value of the log-rank test
                                 risk.table = TRUE,  # Add risk table below the plot
                                 title = "Kaplan-Meier kreivė: Didelės vs. mažos rizikos atvejai KV audinių imtyje",
                                 xlab = "Bendras išgyvenamumo laikas",
                                 ylab = "Išgyvenamumo tikimybė",
                                 palette = c("darkblue", "maroon"),  # Color palette for groups
                                 legend.title = "Rizikos grupė", 
                                 legend.labs = c("Mažas rizikos balas", "Didelis rizikos balas"))
test_survplotctga2

#save
png("C:/Users/Ieva/rprojects/outputs_all/DISS/KM_10_gene_TCGA_COEF_OCa_20150915.png",
    width = 800, height = 600, res = 100) # width and height in pixels, resolution in dpi
test_survplotctga2 #
dev.off() # Close the PNG device

##ALL CASES, TCGA COEFS#############################
# Risk score = sum( gene_expression * coefficient )
risk_scores_tcga_type1 <- rowSums(sweep(gene_data, 2, coef_tcga, "*"))
# View the risk scores
print(risk_scores_tcga_type1) # now I have some risk scores
#add risk scores to the clin_df_joined_test
ALL_SURV_EXPRESSION$RiskScore_tcga <- risk_scores_tcga_type1
#create df wih survival data
surv_df_tcga_type1 <- ALL_SURV_EXPRESSION[, colnames(OC_SURV_EXPRESSION) %in%
                                            c("OS", "STATUS", genes10, "RiskScore_tcga", "patient_id_aud")]

rownames(surv_df_tcga_type1) <- surv_df_tcga_type1$patient_id_aud

# Calculate the median risk score
median_risktcga1 <- median(surv_df_tcga_type1$RiskScore_tcga, na.rm = TRUE) #-7.396536
# Create a new factor column based on the median value
surv_df_tcga_type1$RiskGrouptcga <- ifelse(surv_df_tcga_type1$RiskScore_tcga <= median_risktcga1,
                                           "Low Risk", "High Risk")
#Create a survival object
surv_objecttcga1 <- Surv(time = surv_df_tcga_type1$OS,
                         event = surv_df_tcga_type1$STATUS )

# Fit a Kaplan-Meier model
km_fitctga1 <- survfit(surv_objecttcga1 ~ RiskGrouptcga, data = surv_df_tcga_type1)
# Plot the Kaplan-Meier curve using ggsurvplot
test_survplotctga1 <- ggsurvplot(km_fitctga1, data = surv_df_tcga_type1, 
                                 pval = TRUE,  # Show p-value of the log-rank test
                                 risk.table = TRUE,  # Add risk table below the plot
                                 title = "Kaplan-Meier kreivė: Didelės vs. mažos rizikos atvejai KV audinių imtyje",
                                 xlab = "Bendras išgyvenamumo laikas",
                                 ylab = "Išgyvenamumo tikimybė",
                                 palette = c("darkblue", "maroon"),  # Color palette for groups
                                 legend.title = "Rizikos grupė", 
                                 legend.labs = c("Mažas rizikos balas", "Didelis rizikos balas"))
test_survplotctga1

#save
png("C:/Users/Ieva/rprojects/outputs_all/DISS/KM_10_gene_TCGA_COEF_ALL_20150915.png",
    width = 800, height = 600, res = 100) # width and height in pixels, resolution in dpi
test_survplotctga1 #
dev.off() # Close the PNG device

##HGSOC CASES, TCGA COEFS#############################
# Risk score = sum( gene_expression * coefficient )
risk_scores_tcga_type3 <- rowSums(sweep(gene_dataH, 2, coef_tcga, "*"))
# View the risk scores
print(risk_scores_tcga_type3) # now I have some risk scores
#add risk scores to the clin_df_joined_test
HGSOC_SURV_EXPRESSION$RiskScore_tcga <- risk_scores_tcga_type3
#create df wih survival data
surv_df_tcga_type3 <- HGSOC_SURV_EXPRESSION[, colnames(OC_SURV_EXPRESSION) %in%
                                              c("OS", "STATUS", genes10, "RiskScore_tcga", "patient_id_aud")]

rownames(surv_df_tcga_type3) <- surv_df_tcga_type3$patient_id_aud

# Calculate the median risk score
median_risktcga3 <- median(surv_df_tcga_type3$RiskScore_tcga, na.rm = TRUE) #-7.396536
# Create a new factor column based on the median value
surv_df_tcga_type3$RiskGrouptcga <- ifelse(surv_df_tcga_type3$RiskScore_tcga <= median_risktcga3,
                                           "Low Risk", "High Risk")
#Create a survival object
surv_objecttcga3 <- Surv(time = surv_df_tcga_type3$OS,
                         event = surv_df_tcga_type3$STATUS )

# Fit a Kaplan-Meier model
km_fitctga3 <- survfit(surv_objecttcga3 ~ RiskGrouptcga, data = surv_df_tcga_type3)
# Plot the Kaplan-Meier curve using ggsurvplot
test_survplotctga3 <- ggsurvplot(km_fitctga3, data = surv_df_tcga_type3, 
                                 pval = TRUE,  # Show p-value of the log-rank test
                                 risk.table = TRUE,  # Add risk table below the plot
                                 title = "Kaplan-Meier kreivė: Didelės vs. mažos rizikos atvejai KV audinių imtyje",
                                 xlab = "Bendras išgyvenamumo laikas",
                                 ylab = "Išgyvenamumo tikimybė",
                                 palette = c("darkblue", "maroon"),  # Color palette for groups
                                 legend.title = "Rizikos grupė", 
                                 legend.labs = c("Mažas rizikos balas", "Didelis rizikos balas"))
test_survplotctga3

#save
png("C:/Users/Ieva/rprojects/outputs_all/DISS/KM_10_gene_TCGA_COEF_hgsoc_20150915.png",
    width = 800, height = 600, res = 100) # width and height in pixels, resolution in dpi
test_survplotctga3 #
dev.off() # Close the PNG device

#ENGLISH KM PLOT ###############################
test_survplot2EN <- ggsurvplot(km_fit2, data = surv_df_test2, 
                               #pval = TRUE,  # Show p-value of the log-rank test
                               risk.table = TRUE,  # Add risk table below the plot
                               title = "Risk score association with overall survival, OC group",
                               palette = c("turquoise", "deeppink"),  # Color palette for groups
                               xlab = "Overall survival, months",
                               legend.title = "Risk score", 
                               legend.labs = c( "High risk score", "Low risk score"))
# Add subtitle form cox result
test_survplot2EN$plot <- test_survplot2EN$plot +
  labs(subtitle = paste0("Uni HR =  0.24  (95% CI: 0.07–0.81);\n",
                         "Multi HR = 0.21 (95% CI: 0.05–0.91),\n",
                         "Long-rank p = 0.013, n = 47") )

test_survplot2EN
#save
png("C:/Users/Ieva/rprojects/outputs_all/DISS/KM_10_gene_OCa_20260120ENshort.png",
    width = 800, height = 600, res = 120) # width and height in pixels, resolution in dpi
test_survplot2EN #
dev.off() # Close the PNG device

#ENGLISH plot at year 5  with plotting, no AUCS in the plot###############
# Choose target time (5 years)
target_time <- 60
time_index <- which(rez_list[[1]]$times == target_time)

png("C:/Users/Ieva/rprojects/outputs_all/DISS/tissues_OCtimeROC_test20260121EN.png",
    width = 1200, height = 1200, res = 200)

par(pty = "s")

# Base plot (gene 1)
plot(
  rez_list[[1]]$FP[, time_index],
  rez_list[[1]]$TP[, time_index],
  type = "l",
  col = 1,
  lwd = 2,
  xlab = "1 - Specificity",
  ylab = "Sensitivity",
  main = "ROC curves, 5 years from diagnosis, OC group",
  xlim = c(0, 1),
  ylim = c(0, 1),
  asp = 1
)

# Add ROC curves for remaining genes
for (i in 2:length(rez_list)) {
  lines(
    rez_list[[i]]$FP[, time_index],
    rez_list[[i]]$TP[, time_index],
    col = i,
    lwd = 2
  )
}

# Add risk score ROC (bold maroon)
lines(
  roc_result$FP[, time_index],
  roc_result$TP[, time_index],
  col = "maroon",
  lwd = 3
)

# Diagonal reference line
abline(0, 1, lty = 2, col = "gray")

# Legend labels: italic gene names + Risk score (NO AUC)
legend_labels <- c(
  parse(text = paste0("italic('", names(rez_list), "')")),
  "Risk score"
)

legend(
  "bottomright",
  legend = legend_labels,
  col = c(1:length(rez_list), "maroon"),
  lwd = c(rep(2, length(rez_list)), 3),
  cex = 0.6,
  bty = "n"
)

dev.off()

#TABLE of aucs, sensitivity, specificity##########################
# Extract best sensitivity, specificity, cutoff for Notch risk score at 60 months
sens_60_risk <- roc_result$TP[, which(roc_result$times == 60)]
spec_60_risk <- 1 - roc_result$FP[, which(roc_result$times == 60)]

youden_risk <- sens_60_risk + spec_60_risk - 1
best_idx_risk <- which.max(youden_risk)

risk_score_row <- tibble(
  gene = "Risk score",
  time = roc_result$times[which(roc_result$times == 60)],
  auc = roc_result$AUC[which(roc_result$times == 60)]*100,
  sens = sens_60_risk[best_idx_risk],
  spec = spec_60_risk[best_idx_risk],
  cutoff = roc_result$cutoffs[best_idx_risk]
)

# Build table of AUC + Sensitivity + Specificity at 60 months
sens_spec_auc_60 <- map_dfr(names(rez_list), function(gene) {
  roc <- rez_list[[gene]]
  
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


# Combine with gene-level table
sens_spec_auc_60_all <- bind_rows(sens_spec_auc_60, risk_score_row)

print(sens_spec_auc_60_all)

#make prety table
# Prepare table: rename columns for display
roc_table_display <- sens_spec_auc_60_all %>%
  rename(
    Biomarker = gene,
    Time_months = time,
    AUC = auc,
    Sensitivity = sens,
    Specificity = spec
  )

# Create gt table
gt_table_roc_60 <- roc_table_display %>%
  gt() %>%
  tab_header(
    title = "ROC metrics",
    subtitle = "Prognostic criteria, 5 years from diagnosis"
  ) %>%
  fmt_number(
    columns = vars(AUC, Sensitivity, Specificity),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biomarker))
  )

# Show table
gt_table_roc_60

#there is no other convenient way to save gt outputs
gtsave(gt_table_roc_60,
       filename = "10genetimeroc_table_20260121EN.png")

#Combine the images
roc_image <- image_read("tissues_OCtimeROC_test20260121EN.png")
table_image <- image_read("10genetimeroc_table_20260121EN.png")

combined_image <- image_append(c(roc_image, table_image), stack = F)

# Save the combined image
image_write(combined_image, 
            "tissues_OCtimeROCw_table_test20260121EN.png")


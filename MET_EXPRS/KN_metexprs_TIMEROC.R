#KN-DISSERTATION project. 10 gene panel - genes selected via statistical analysis of TCGA and GTEx data
#KN 10 gene analysis data (cv 2 %, one threshold, 65 cases!, restored KN-59) 2026-01-20
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
  Surv(OS, STATUS) ~ NOTCH1 + NOTCH2 + NOTCH3 + NOTCH4 + ARID1A + CTNNB1 + FBXW7 + HES1 + DLL1 + JAG2,
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
                                   risk.table.title = "Pacienčių skaičius rizikos grupėje",
                                   title = expression(
                                     bold("A ") * "Notch, Wnt ir " * italic("ARID1A") * " genų raiškos kombinacija"
                                   ),
                                   xlab = "Bendras išgyvenamumo laikas, mėnesiais",
                                   ylab = "Išgyvenamumo tikimybė",
                                   palette = c("deeppink", "turquoise"),  # Color palette for groups
                                   legend.title = "Rizikos grupė", 
                                   legend.labs = c( "Didelės rizikos balas", "Mažos rizikos balas")) #changed form "mažas/ didelis rizikos balas"
# Add subtitle form cox result
test_survplot2_notch$plot <- test_survplot2_notch$plot +
  labs(subtitle = 
paste0("Uni PR = 0,32 (95 % PI: 0,12–0,92);\n",
"Multi PR = 0,32 (95 % PI: 0,10–1,06);\n",
"Log-rank p = 0,03, n = 47") )

print(test_survplot2_notch)
#save
png("KM_plot_notch_arid1a_OC_w_HR20260130short.png",
    width = 15, height = 15, res = 200, units = "cm")
print(test_survplot2_notch)  # print the full ggsurvplot object
dev.off()

# EN A Kaplan-Meier curve############################
test_survplot2_notchEN <- ggsurvplot(km_fit2notch, data = surv_df_test2notch, 
                                     #pval = TRUE,  # Show p-value of the log-rank test
                                     risk.table = TRUE, 
                                     # risk.table.title = "Pacienčių skaičius rizikos grupėje",
                                     title = expression(
                                       bold("A ") * "Notch, Wnt and " * italic("ARID1A") * " gene expression combination"
                                     ),
                                     xlab = "Overall survival, months",
                                     #ylab = "Išgyvenamumo tikimybė",
                                     palette = c("deeppink", "turquoise"),  # Color palette for groups
                                     #legend.title = "Rizikos grupė", 
                                     font.main = 14, 
                                     legend.labs = c("High risk score", "Low risk score"))
# Add subtitle form cox result
test_survplot2_notchEN$plot <- test_survplot2_notchEN$plot +
  labs(subtitle = 
         "Uni HR = 0,32 (95 % CI: 0.12–0.92);
Multi HR = 0.32 (95 % CI: 0.10–1.06);
Log-rank p = 0.03, n = 47")

print(test_survplot2_notchEN)

#save
png("KM_plot_notch_arid1a_OC_w_HR_EN20260130short.png",
    width = 15, height = 15, res = 200, units = "cm")
print(test_survplot2_notchEN)  # print the full ggsurvplot object
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
                             risk.table.title = "Pacienčių skaičius rizikos grupėje",
                             title = expression(
                               bold("B ") * "Promotorių metilinimo, Notch, Wnt ir" * italic(" ARID1A") * " genų raiškos kombinacija"
                             ),
                             xlab = "Bendras išgyvenamumo laikas, mėnesiais",
                             ylab = "Išgyvenamumo tikimybė",
                             palette = c( "deeppink", "turquoise"),  # Color palette for groups
                             legend.title = "Rizikos grupė", 
                             font.main = 11, 
                             legend.labs = c( "Didelės rizikos balas", "Mažos rizikos balas")) #changed form "mažs/didelės rizikos balas"
# Add subtitle form cox result
test_survplot14$plot <- test_survplot14$plot +
  labs(subtitle = paste0("Uni PR =  0.22  (95 % PI: 0,07–0,66);\n",
"Multi PR = 0.22 (95 % PI: 0,06–0,77);\n",
"Log-rank p = 0,003, n = 47") )

print(test_survplot14)
#save
png("KM_plot_14_OC_w_HR20260130short.png",
    width = 15, height = 15, res = 200, units = "cm")
print(test_survplot14)  # print the full ggsurvplot object
dev.off()

# EN B Kaplan-Meier curve########################
test_survplot14EN<- ggsurvplot(km_fit14, data = surv_df_test14, 
                               #pval = TRUE,  # Show p-value of the log-rank test
                               risk.table = TRUE,  # Add risk table below the plot
                               #risk.table.title = "Pacienčių skaičius rizikos grupėje",
                               title = expression(
                                 bold("B ") * "Promoter methylation, Notch, Wnt and" * italic(" ARID1A") * " gene expression combination"
                               ),
                               xlab = "Overall survival, months",
                               #ylab = "Išgyvenamumo tikimybė",
                               palette = c( "deeppink", "turquoise"),  # Color palette for groups
                               #legend.title = "Rizikos grupė", 
                               font.main = 10, 
                               legend.labs = c("High risk score", "Low risk score")) #changed places
# Add subtitle form cox result
test_survplot14EN$plot <- test_survplot14EN$plot +
  labs(subtitle = paste0("Uni HR =  0.22  (95 % CI: 0.07–0.66);\n",
"Multi HR = 0.22 (95 % CI: 0.06–0.77);\n",
"Log-rank p = 0.003, n = 47") )

print(test_survplot14EN)
#save
png("KM_plot_14_OC_w_HR_en20260130short.png",
    width = 15, height = 15, res = 200, units = "cm")
print(test_survplot14EN)  # print the full ggsurvplot object
dev.off()

#NOTCH2 and HES1 risk score####################################
#gene_dataNOTCH_HES1 - only OC, notch and hes
gene_dataNOTCH_HES1 <- OC_SURV_EXPRESSION[, colnames(OC_SURV_EXPRESSION) %in% c("HES1", "NOTCH3")]
# Example with 10 biomarkers
cox_model_OC_notch_hes <- coxph(
  Surv(OS, STATUS) ~ NOTCH3 + HES1 ,
  data = OC_SURV_EXPRESSION
)
summary(cox_model_OC_notch_hes)

#get coeficients
coefs_notch_hes <- coef(cox_model_OC_notch_hes)
# Risk score = sum( gene_expression * coefficient )
risk_scores_notch_hes <- rowSums(sweep(gene_dataNOTCH_HES1, 2, coefs_notch_hes, "*"))
# View the risk scores
print(risk_scores_notch_hes) # now I have some risk scores
#add risk scores to the clin_df_joined_test
OC_SURV_EXPRESSION$RiskScore_notch_hes <- risk_scores_notch_hes
#create df wih survival data
surv_df_testnotch_hes <- OC_SURV_EXPRESSION[, colnames(OC_SURV_EXPRESSION) %in%
                                              c("OS", "STATUS", genes_notch, "RiskScore_notch_hes", "patient_id_aud", "Age", "CA125")]

rownames(surv_df_testnotch_hes) <- surv_df_testnotch_hes$patient_id_aud


# Calculate the median risk score
median_risk_notch_hes <- median(surv_df_testnotch_hes$RiskScore_notch_hes, na.rm = TRUE) #-0.08001295
# Create a new factor column based on the median value
surv_df_testnotch_hes$RiskScore_notch_hes <- ifelse(surv_df_testnotch_hes$RiskScore_notch_hes <= median_risk_notch_hes,
                                                    "Low Risk", "High Risk")
#Create a survival object
surv_object_notch_hes<- Surv(time = surv_df_testnotch_hes$OS,
                             event = surv_df_testnotch_hes$STATUS )
# Fit Univariable cox model for notch3 hes model ############################
cox_model_notch_hes <- coxph(Surv(OS, STATUS) ~ RiskScore_notch_hes , data = surv_df_testnotch_hes)
summary(cox_model_notch_hes)
# Fit Multivariable cox model for notch3 hes model ############################
cox_model_notch_hes2 <- coxph(Surv(OS, STATUS) ~RiskScore_notch_hes + Age + CA125, data = surv_df_testnotch_hes)
summary(cox_model_notch_hes2)
# Fit a Kaplan-Meier model for notch3 hes model#################################
km_fitnotch_hes <- survfit(surv_object_notch_hes ~ RiskScore_notch_hes, data = surv_df_testnotch_hes)
# Plot the Kaplan-Meier curve using ggsurvplot
test_survplot_notch_hes <- ggsurvplot(km_fitnotch_hes, data = surv_df_testnotch_hes, 
                                      #pval = TRUE,  # Show p-value of the log-rank test
                                      risk.table = TRUE, 
                                      risk.table.title = "Pacienčių skaičius rizikos grupėje",
                                      title =  expression(
                                        bold("C ") * italic("NOTCH3") * " ir " * italic(" HES1") * " genų raiškos kombinacija"
                                      ),
                                      xlab = "Bendras išgyvenamumo laikas, mėnesiais",
                                      ylab = "Išgyvenamumo tikimybė",
                                      palette = c( "deeppink", "turquoise"),  # Color palette for groups
                                      legend.title = "Rizikos grupė", 
                                      legend.labs = c( "Didelės rizikos balas", "Mažos rizikos balas")) #changed places form "mažas / didelis rizikos balas"
# Add subtitle form cox result
test_survplot_notch_hes$plot <- test_survplot_notch_hes$plot +
  labs(subtitle = paste0("Uni PR = 0,37  (95 % PI: 0,14–0,98);\n",
"Multi PR = 0,50 (95 % PI: 0,17–1,47);\n",
"Log-rank p = 0,04, n = 48") )

print(test_survplot_notch_hes)

#save
png("KM_plot_notch_HES_OC_w_HR20260130short.png",
    width = 15, height = 15, res = 200, units = "cm")
print(test_survplot_notch_hes)  # print the full ggsurvplot object
dev.off()
#EN  C Kaplan-Meier #############################
test_survplot_notch_hesEN <- ggsurvplot(km_fitnotch_hes, data = surv_df_testnotch_hes, 
                                        #pval = TRUE,  # Show p-value of the log-rank test
                                        risk.table = TRUE, 
                                        #risk.table.title = "Pacienčių skaičius rizikos grupėje",
                                        title =  expression(
                                          bold("C ") * italic("NOTCH3") * " and " * italic(" HES1") * " gene expression combination"
                                        ),
                                        xlab = "Overall survival, months",
                                        #ylab = "Išgyvenamumo tikimybė",
                                        palette = c( "deeppink", "turquoise"),  # Color palette for groups
                                        #legend.title = "Rizikos grupė", 
                                        legend.labs = c("High risk score", "Low risk score")) #changed places
# Add subtitle form cox result
test_survplot_notch_hesEN$plot <- test_survplot_notch_hesEN$plot +
  labs(subtitle = paste0("Uni HR = 0.37  (95 % CI: 0.14–0.98);\n",
"Multi HR = 0.50 (95 % CI: 0.17–1.47);\n",
"Log-rank p = 0.04, n = 48") )

print(test_survplot_notch_hesEN)

#save
png("KM_plot_notch_HES_OC_w_HREN20260130short.png",
    width = 15, height = 15, res = 200, units = "cm")
print(test_survplot_notch_hesEN)  # print the full ggsurvplot object
dev.off()
#COMBINE IMAGES##############################################
img1 <- image_read("KM_plot_notch_arid1a_OC_w_HR20260130short.png")
img2 <- image_read("KM_plot_14_OC_w_HR20260130short.png")
img3 <- image_read("KM_plot_notch_HES_OC_w_HR20260130short.png")

# Match heights to avoid misalignment
target_height <- min(
  image_info(img1)$height,
  image_info(img2)$height,
  image_info(img3)$height
)

img1 <- image_resize(img1, paste0("x", target_height))
img2 <- image_resize(img2, paste0("x", target_height))
img3 <- image_resize(img3, paste0("x", target_height))

combined_km <- image_append(c(img1, img2, img3), stack = FALSE)

image_write(
  combined_km,
  "KM_plots_combined_horizontal_20260130.png"
)

#COMBINE STAKED IMAGE###############################
#save
png("KM_plot_notch_arid1a_OC_w_HR20260131short.png",
    width = 25, height = 15, res = 200, units = "cm")
print(test_survplot2_notch)  # print the full ggsurvplot object
dev.off()
test_survplot14<- ggsurvplot(km_fit14, data = surv_df_test14, 
                             #pval = TRUE,  # Show p-value of the log-rank test
                             risk.table = TRUE,  # Add risk table below the plot
                             risk.table.title = "Pacienčių skaičius rizikos grupėje",
                             title = expression(
                               bold("B ") * "Promotorių metilinimo, Notch, Wnt ir" * italic(" ARID1A") * " genų raiškos kombinacija"
                             ),
                             xlab = "Bendras išgyvenamumo laikas, mėnesiais",
                             ylab = "Išgyvenamumo tikimybė",
                             palette = c( "deeppink", "turquoise"),  # Color palette for groups
                             legend.title = "Rizikos grupė", 
                             font.main = 15, 
                             legend.labs = c( "Didelės rizikos balas", "Mažos rizikos balas")) #changed form "mažs/didelės rizikos balas"
# Add subtitle form cox result
test_survplot14$plot <- test_survplot14$plot +
  labs(subtitle = paste0("Uni PR =  0.22  (95 % PI: 0,07–0,66);\n",
                         "Multi PR = 0.22 (95 % PI: 0,06–0,77);\n",
                         "Log-rank p = 0,003, n = 47") )

print(test_survplot14)
#save
png("KM_plot_14_OC_w_HR20260131short.png",
    width = 25, height = 15, res = 200, units = "cm")
print(test_survplot14)  # print the full ggsurvplot object
dev.off()
#save
png("KM_plot_notch_HES_OC_w_HR20260131short.png",
    width = 25, height = 15, res = 200, units = "cm")
print(test_survplot_notch_hes)  # print the full ggsurvplot object
dev.off()

# Read images
img1 <- image_read("KM_plot_notch_arid1a_OC_w_HR20260131short.png")
img2 <- image_read("KM_plot_14_OC_w_HR20260131short.png")
img3 <- image_read("KM_plot_notch_HES_OC_w_HR20260131short.png")

# Match widths to avoid misalignment
target_width <- min(
  image_info(img1)$width,
  image_info(img2)$width,
  image_info(img3)$width
)

img1 <- image_resize(img1, paste0(target_width, "x"))
img2 <- image_resize(img2, paste0(target_width, "x"))
img3 <- image_resize(img3, paste0(target_width, "x"))

# Stack vertically
combined_km_vertical <- image_append(c(img1, img2, img3), stack = TRUE)

# Save the stacked image
image_write(
  combined_km_vertical,
  "KM_plots_combined_vertical_20260130.png"
)

##combine EN image##################
img1en <- image_read("KM_plot_notch_arid1a_OC_w_HR_EN20260130short.png")
img2en <- image_read("KM_plot_14_OC_w_HR_en20260130short.png")
img3en <- image_read("KM_plot_notch_HES_OC_w_HREN20260130short.png")

# Match heights to avoid misalignment
target_height <- min(
  image_info(img1en)$height,
  image_info(img2en)$height,
  image_info(img3en)$height
)

img1en <- image_resize(img1en, paste0("x", target_height))
img2en <- image_resize(img2en, paste0("x", target_height))
img3en <- image_resize(img3en, paste0("x", target_height))

combined_kmen <- image_append(c(img1en, img2en, img3en), stack = FALSE)

image_write(
  combined_kmen,
  "KM_plots_combined_horizontal_EN20260130.png"
)

#staked en image
# Match heights for all panels
target_height <- min(
  image_info(img1en)$height,
  image_info(img2en)$height,
  image_info(img3en)$height
)

img1en <- image_resize(img1en, paste0("x", target_height))
img2en <- image_resize(img2en, paste0("x", target_height))
img3en <- image_resize(img3en, paste0("x", target_height))

# Create blank image (same size as one panel)
panel_width <- image_info(img1en)$width
blank_panel <- image_blank(
  width  = panel_width,
  height = target_height,
  color  = "white"
)

# Row 1: two images
row1 <- image_append(c(img1en, img2en), stack = FALSE)

# Row 2: image + blank space
row2 <- image_append(c(img3en, blank_panel), stack = FALSE)

# Stack rows
combined_kmen <- image_append(c(row1, row2), stack = TRUE)

# Save
image_write(
  combined_kmen,
  "KM_plots_combined_2rows_with_space_EN20260130.png"
)
#TIME ROC, NOTCH, OC #########################################
surv_df_notch_oc <- OC_SURV_EXPRESSION[, colnames(OC_SURV_EXPRESSION) %in%
                                         c("OS", "STATUS", genes_notch, methylation, "RiskScore_notch", "RiskScore_14", "RiskScore_notch_hes", "patient_id_aud")]

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

print(auc_table_notch_oc)

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
    auc = roc$AUC[idx_60],
    sens = sens_60[best_idx],
    spec = spec_60[best_idx],
    cutoff = roc$cutoffs[best_idx]
  )
})

print(sens_spec_auc_60)


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

##time roc, notch3 and hes1 genes ###################
roc_result_model_notch_hes <- timeROC(
  T = surv_df_notch_oc$OS,       # Survival time from df
  delta = surv_df_notch_oc$STATUS, # Event indicator from df
  marker = surv_df_notch_oc$RiskScore_notch_hes, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = t_eval,    # Time points for ROC
  iid = TRUE         # Compute confidence intervals
)
roc_result_model_notch_hes

##COORDS time roc, notch3 and hes1#############################
# Sensitivity and specificity vectors at 60 months
sens_60 <- roc_result_model_notch_hes$TP[, which(roc_result_model_notch_hes$times == 60)]
spec_60 <- 1 - roc_result_model_notch_hes$FP[, which(roc_result_model_notch_hes$times == 60)]

# Youden index
youden_60 <- sens_60 + spec_60 - 1
best_index_60 <- which.max(youden_60)

# Extract best specificity
best_spec_60 <- spec_60[best_index_60]
best_sens_60 <- sens_60[best_index_60]
best_cutoff_60 <- roc_result_model_notch_hes$cutoffs[best_index_60]

cat("At 60 months:\n")
cat("  Specificity:", round(best_spec_60, 3), "\n")
cat("  Sensitivity:", round(best_sens_60, 3), "\n")
cat("  Cutoff:", best_cutoff_60, "\n")

##plot at year 1  with plotting###############
# Choose target time
target_time <- 12        
time_index <- which(rez_list2_notch_OC[[1]]$times == target_time)

png("C:/Users/Ieva/rprojects/outputs_all/DISS/tissues_OC_timeROC_12_test20260121.png",
    width = 15, height = 15, res = 510, units = "cm") 
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

png("C:/Users/Ieva/rprojects/outputs_all/DISS/tissues_OC_timeROC_36_test20260121.png",
    width = 15, height = 15, res = 510, units = "cm") # width and height in pixels, resolution in dpi
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

png("C:/Users/Ieva/rprojects/outputs_all/DISS/tissues_OCtimeROC_test20260121_5yr.png",
    width = 15, height = 15, res = 510, units = "cm") # width and height in pixels, resolution in dpi
# Set up base plot with gene 1
par(pty="s")
# Set up base plot with gene 1
plot(
  rez_list2_notch_OC[[1]]$FP[, time_index],
  rez_list2_notch_OC[[1]]$TP[, time_index],
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

# FINAL 5 year plot, no AUCs #######################
# Choose target time
target_time <- 60        
time_index <- which(rez_list2_notch_OC[[1]]$times == target_time)

png("C:/Users/Ieva/rprojects/outputs_all/DISS/tissuesNO_AUCS_OCtimeROC_test20202601215yr.png",
    width = 1500, height = 1500, res = 200) # width and height in pixels, resolution in dpi

par(pty="s")

# --- Base ROC curve plot with first gene ---
plot(
  rez_list2_notch_OC[[1]]$FP[, time_index],
  rez_list2_notch_OC[[1]]$TP[, time_index],
  type = "l",
  col = 1,
  lwd = 2,
  xlab = "1 - Specifiškumas",
  ylab = "Jautrumas",
  main = paste("Nuo laiko priklausomos ROC kreivės,\n",
               "5 metai po diagnozės audinių imtyje"),
  xlim = c(0, 1),
  ylim = c(0, 1),
  asp = 1
)

# --- Add ROC lines for all genes ---
for (i in 2:length(rez_list2_notch_OC)) {
  lines(
    rez_list2_notch_OC[[i]]$FP[, time_index],
    rez_list2_notch_OC[[i]]$TP[, time_index],
    col = i,
    lwd = 2
  )
}


# --- Add risk score ROC curves ---
lines(
  roc_result_model_notch$FP[, time_index],
  roc_result_model_notch$TP[, time_index],
  col = "maroon",
  lwd = 3,
  lty = 2
)

lines(
  roc_result_model_notch_met$FP[, time_index],
  roc_result_model_notch_met$TP[, time_index],
  col = "darkgreen",
  lwd = 3,
  lty = 2
)


lines(
  roc_result_model_notch_hes$FP[, time_index],
  roc_result_model_notch_hes$TP[, time_index],
  col = "darkblue",
  lwd = 3,
  lty = 2
)

# --- Diagonal reference line ---
abline(0, 1, lty = 2, col = "gray")

# --- Legend preparation (no AUCs) ---
legend_labels <- c(
  names(rez_list2_notch_OC),
  "HES1 ir NOTCH3 raiškos rizikos balas",
  "Genų raiškos rizikos balas",
  "14 biožymenų rizikos balas"
)

# italicize gene names only, leave risk scores normal
legend_labels_expr <- c(
  lapply(names(rez_list2_notch_OC), function(x) bquote(italic(.(x)))),
  legend_labels[(length(rez_list2_notch_OC)+1):length(legend_labels)]
)

legend(
  "bottomright",
  legend = legend_labels_expr,
  col   = c(1:length(rez_list2_notch_OC),  "maroon", "darkgreen", "darkblue"),
  lwd   = c(rep(2, length(rez_list2_notch_OC)), 3, 3, 3),
  lty   = c(rep(1, length(rez_list2_notch_OC)), 2, 2, 2),
  cex   = 0.6,
  bty   = "n"
)

dev.off()


#TABLE of aucs, sensitivity, specificity##########################
# Extract best sensitivity, specificity, cutoff for Notch risk score at 60 months
sens_60_risk <- roc_result_model_notch$TP[, which(roc_result_model_notch$times == 60)]
spec_60_risk <- 1 - roc_result_model_notch$FP[, which(roc_result_model_notch$times == 60)]

youden_risk <- sens_60_risk + spec_60_risk - 1
best_idx_risk <- which.max(youden_risk)

risk_score_row <- tibble(
  gene = "Genų raiškos rizikos balas",
  time = roc_result_model_notch$times[which(roc_result_model_notch$times == 60)],
  auc = roc_result_model_notch$AUC[which(roc_result_model_notch$times == 60)],
  sens = sens_60_risk[best_idx_risk],
  spec = spec_60_risk[best_idx_risk],
  cutoff = roc_result_model_notch$cutoffs[best_idx_risk]
)

# Combine with gene-level table
sens_spec_auc_60_all <- bind_rows(sens_spec_auc_60, risk_score_row)

# Extract best sensitivity, specificity, cutoff for 14-biomarker risk score at 60 months
sens_60_met <- roc_result_model_notch_met$TP[, which(roc_result_model_notch_met$times == 60)]
spec_60_met <- 1 - roc_result_model_notch_met$FP[, which(roc_result_model_notch_met$times == 60)]

youden_60_met <- sens_60_met + spec_60_met - 1
best_idx_met <- which.max(youden_60_met)

risk_score_met_row <- tibble(
  gene = "14 biožymenų rizikos balas",
  time = roc_result_model_notch_met$times[which(roc_result_model_notch_met$times == 60)],
  auc = roc_result_model_notch_met$AUC[which(roc_result_model_notch_met$times == 60)],
  sens = sens_60_met[best_idx_met],
  spec = spec_60_met[best_idx_met],
  cutoff = roc_result_model_notch_met$cutoffs[best_idx_met]
)

# Combine with previous table (genes + Notch risk score)
sens_spec_auc_60_all <- bind_rows(sens_spec_auc_60_all, risk_score_met_row)

# Extract best sensitivity, specificity, cutoff for notch3 and hes risk score at 60 months
sens_60_met <- roc_result_model_notch_hes$TP[, which(roc_result_model_notch_hes$times == 60)]
spec_60_met <- 1 - roc_result_model_notch_hes$FP[, which(roc_result_model_notch_hes$times == 60)]

youden_60_met <- sens_60_met + spec_60_met - 1
best_idx_met <- which.max(youden_60_met)

risk_score_notch_hes_row <- tibble(
  gene = "NOTCH3 ir HES1 rizikos balas",
  time = roc_result_model_notch_hes$times[which(roc_result_model_notch_hes$times == 60)],
  auc = roc_result_model_notch_hes$AUC[which(roc_result_model_notch_hes$times == 60)],
  sens = sens_60_met[best_idx_met],
  spec = spec_60_met[best_idx_met],
  cutoff = roc_result_model_notch_hes$cutoffs[best_idx_met]
)

# Combine with previous table (genes + Notch risk score)
sens_spec_auc_60_all <- bind_rows(sens_spec_auc_60_all, risk_score_notch_hes_row)

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
    subtitle = "Prognostinis tikslumas 5 metų išgyvenamumui"
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
       filename = "timeroc_table_20260121.png")

#Combine the images
#save images together
roc_image2   <- image_read("tissuesNO_AUCS_OCtimeROC_test20202601215yr.png")
table_image2 <- image_read("timeroc_table_20260121.png")

# Get height of ROC image
roc_height <- image_info(roc_image2)$height

# Resize table to same height (keeps aspect ratio)
table_image2_resized <- image_resize(
  table_image2,
  geometry = paste0("x", roc_height)
)

combined_image2 <- image_append(
  c(roc_image2, table_image2_resized),
  stack = FALSE
)

image_write(
  combined_image2,
  "roc_notch_all_oc_with_table_20260130.png"
)

# ENGLSIH FINAL 5 year plot, no AUCs #######################
# Choose target time
target_time <- 60        
time_index <- which(rez_list2_notch_OC[[1]]$times == target_time)
#rename ARID1A_met
names(rez_list2_notch_OC)[names(rez_list2_notch_OC) == "ARID1A_met"] <- 
  "ARID1A methylation"
png("C:/Users/Ieva/rprojects/outputs_all/DISS/tissuesNO_AUCS_OCtimeROC_testEN20260121.png",
    width = 20, height = 20, res = 510, units = "cm") # width and height in pixels, resolution in dpi

par(pty="s")

# --- Base ROC curve plot with first gene ---
plot(
  rez_list2_notch_OC[[1]]$FP[, time_index],
  rez_list2_notch_OC[[1]]$TP[, time_index],
  type = "l",
  col = 1,
  lwd = 2,
  xlab = "1 - Specificity",
  ylab = "Sensitivity",
  cex.main = 1.5,  
  main = paste("ROC curves, 5 years form diagnosis, OC group"),
  xlim = c(0, 1),
  ylim = c(0, 1),
  asp = 1
)

# --- Add ROC lines for all genes ---
for (i in 2:length(rez_list2_notch_OC)) {
  lines(
    rez_list2_notch_OC[[i]]$FP[, time_index],
    rez_list2_notch_OC[[i]]$TP[, time_index],
    col = i,
    lwd = 2
  )
}


# --- Add risk score ROC curves ---
lines(
  roc_result_model_notch$FP[, time_index],
  roc_result_model_notch$TP[, time_index],
  col = "maroon",
  lwd = 3,
  lty = 2
)

lines(
  roc_result_model_notch_met$FP[, time_index],
  roc_result_model_notch_met$TP[, time_index],
  col = "darkgreen",
  lwd = 3,
  lty = 2
)


lines(
  roc_result_model_notch_hes$FP[, time_index],
  roc_result_model_notch_hes$TP[, time_index],
  col = "darkblue",
  lwd = 3,
  lty = 2
)

# --- Diagonal reference line ---
abline(0, 1, lty = 2, col = "gray")

# --- Legend preparation (no AUCs) ---
legend_labels <- c(
  names(rez_list2_notch_OC),
  "HES1 and NOTCH3 gene expression combination",
  "Notch, Wnt and ARID1A gene expression combination",
  "14 biomarker combination"
)

# # italicize gene names only, leave risk scores normal
# legend_labels_expr <- c(
#   lapply(names(rez_list2_notch_OC), function(x) bquote(italic(.(x)))),
#   legend_labels[(length(rez_list2_notch_OC)+1):length(legend_labels)]
# )
legend_labels_expr <- c(
  # gene-specific ROC curves (italic)
  lapply(names(rez_list2_notch_OC), function(x) bquote(italic(.(x)))),
  
  # additional ROC curves (italic)
  bquote(italic("HES1 and NOTCH3 gene expression combination")),
  bquote(italic("Notch, Wnt and ARID1A gene expression combination")),
  bquote(italic("14 biomarker combination"))
)

legend(
  "bottomright",
  legend = legend_labels_expr,
  col   = c(1:length(rez_list2_notch_OC),  "maroon", "darkgreen", "darkblue"),
  lwd   = c(rep(2, length(rez_list2_notch_OC)), 3, 3, 3),
  lty   = c(rep(1, length(rez_list2_notch_OC)), 2, 2, 2),
  cex   = 0.6,
  bty   = "n"
)

dev.off()
#make prety table
# Prepare table: rename columns for display
roc_table_display2 <- sens_spec_auc_60_all %>%
  rename(
    Predictor = gene,
    Time_months = time,
    AUC = auc,
    Sensitivity = sens,
    Specificity = spec
  )
roc_table_display2$Predictor <- c(  "NOTCH1",
                                    "NOTCH2",
                                    "NOTCH3",
                                    "NOTCH4",
                                    "ARID1A",
                                    "CTNNB1" ,
                                    "FBXW7",
                                    "JAG2",
                                    "DLL1",
                                    "HES1",
                                    "HOPX",
                                    "ALX4" ,
                                    "CDX2",
                                    "ARID1A methylation",
                                    "Notch, Wnt and ARID1A gene expression combination"  ,
                                    "14 biomarker combination" ,
                                    "HES1 and NOTCH3 gene expression combination")
# Create gt table
gt_table_roc_60EN <- roc_table_display2 %>%
  gt() %>%
  tab_header(
    title = "ROC metrics",
    subtitle = "Prognostic citeria, 5 years from diagnosis "
  ) %>%
  fmt_number(
    columns = vars(AUC, Sensitivity, Specificity),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Predictor))
  )

# Show table
gt_table_roc_60EN


#there is no other convenient way to save gt outputs
gtsave(gt_table_roc_60EN,vwidth = 1000, 
       filename = "timeroc_table_EN20260121.png")

#save images together
roc_image2   <- image_read("tissuesNO_AUCS_OCtimeROC_testEN20260121.png")
table_image2 <- image_read("timeroc_table_EN20260121.png")

# Get height of ROC image
roc_height <- image_info(roc_image2)$height

# Resize table to same height (keeps aspect ratio)
table_image2_resized <- image_resize(
  table_image2,
  geometry = paste0("x", roc_height)
)

combined_image2 <- image_append(
  c(roc_image2, table_image2_resized),
  stack = FALSE
)

image_write(
  combined_image2,
  "roc_notch_all_oc_with_table_en20260130.png"
)


#compare time rocs ##################################################
#separate
rez_list2_notch_OC_60 <- apply(surv_df_notch_oc[, c(2:15)], 2, timeROC,
                               T = surv_df_notch_oc$OS,       # Survival time from df
                               delta =  surv_df_notch_oc$STATUS,# Event indicator from df
                               #marker  # Predictor already in the df
                               cause = 1,         # Event of interest
                               times = 60,    # Time points for ROC
                               iid = TRUE )        # Compute confidence intervals)
#gene expresion model
roc_result_model_notch_60 <- timeROC(
  T = surv_df_notch_oc$OS,       # Survival time from df
  delta = surv_df_notch_oc$STATUS, # Event indicator from df
  marker = surv_df_notch_oc$RiskScore_notch, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = 60,    # Time points for ROC
  iid = TRUE         # Compute confidence intervals
)
#gene expression + methylation model
roc_result_model_notch_met_60 <- timeROC(
  T = surv_df_notch_oc$OS,       # Survival time from df
  delta = surv_df_notch_oc$STATUS, # Event indicator from df
  marker = surv_df_notch_oc$RiskScore_14, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = 60,    # Time points for ROC
  iid = TRUE      )   # Compute confidence intervals
#compare time roc models
compare(roc_result_model_notch_met_60, roc_result_model_notch_60, adjusted = FALSE) #0.5482305 


#compare all time roc models for separate genes to the 14 biomarker model
#compare full model to arid1a_met
#remove case KN-100 because it is not compatable with comaprisons (empty gene expression)
surv_df_notch_oc100 <- surv_df_notch_oc[!c(surv_df_notch_oc$patient_id_aud %in% "KN-100"),]

#notch3 and hes 1 model
roc_result_model_notch_hes_60 <- timeROC(
  T = surv_df_notch_oc100$OS,       # Survival time from df
  delta = surv_df_notch_oc100$STATUS, # Event indicator from df
  marker = surv_df_notch_oc100$RiskScore_notch_hes, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = 60,    # Time points for ROC
  iid = TRUE      )   # Compute confidence intervals

#compare time roc models
compare(roc_result_model_notch_hes_60, roc_result_model_notch_60, adjusted = FALSE) #0.7189816  
compare(roc_result_model_notch_met_60, roc_result_model_notch_hes_60, adjusted = FALSE) #0.8340145 

#ARID1A MET
roc_result_model_ARID1A_met_60 <- timeROC(
  T = surv_df_notch_oc100$OS,       # Survival time from df
  delta = surv_df_notch_oc100$STATUS, # Event indicator from df
  marker = surv_df_notch_oc100$ARID1A_met, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = 60,    # Time points for ROC
  iid = TRUE      )   # Compute confidence intervals

compare(roc_result_model_notch_met_60, roc_result_model_ARID1A_met_60, adjusted = F) #0.223373 

#HOPX MET
roc_result_model_HOPX_met_60 <- timeROC(
  T = surv_df_notch_oc100$OS,       # Survival time from df
  delta = surv_df_notch_oc100$STATUS, # Event indicator from df
  marker = surv_df_notch_oc100$HOPX, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = 60,    # Time points for ROC
  iid = TRUE      )   # Compute confidence intervals
compare(roc_result_model_notch_met_60, roc_result_model_HOPX_met_60, adjusted = F) #0.0006880568  

#ALX4 MET
roc_result_model_ALX4_met_60 <- timeROC(
  T = surv_df_notch_oc100$OS,       # Survival time from df
  delta = surv_df_notch_oc100$STATUS, # Event indicator from df
  marker = surv_df_notch_oc100$ALX4, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = 60,    # Time points for ROC
  iid = TRUE      )   # Compute confidence intervals

compare(roc_result_model_notch_met_60, roc_result_model_ALX4_met_60, adjusted = F) #0.02423082  

#cdx2 MET
roc_result_model_CDX2_met_60 <- timeROC(
  T = surv_df_notch_oc100$OS,       # Survival time from df
  delta = surv_df_notch_oc100$STATUS, # Event indicator from df
  marker = surv_df_notch_oc100$CDX2, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = 60,    # Time points for ROC
  iid = TRUE      )   # Compute confidence intervals

compare(roc_result_model_notch_met_60, roc_result_model_CDX2_met_60, adjusted = F) #0.006333155    
#hes1
roc_result_model_HES1_60 <- timeROC(
  T = surv_df_notch_oc100$OS,       # Survival time from df
  delta = surv_df_notch_oc100$STATUS, # Event indicator from df
  marker = surv_df_notch_oc100$HES1, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = 60,    # Time points for ROC
  iid = TRUE      )   # Compute confidence intervals

compare(roc_result_model_notch_met_60, roc_result_model_HES1_60, adjusted = F) #0.001634476 

#jag2
roc_result_model_JAG2_60 <- timeROC(
  T = surv_df_notch_oc100$OS,       # Survival time from df
  delta = surv_df_notch_oc100$STATUS, # Event indicator from df
  marker = surv_df_notch_oc100$JAG2, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = 60,    # Time points for ROC
  iid = TRUE      )   # Compute confidence intervals

compare(roc_result_model_notch_met_60, roc_result_model_JAG2_60, adjusted = F) #0.04426828    

#dll1
roc_result_model_DLL1_60 <- timeROC(
  T = surv_df_notch_oc100$OS,       # Survival time from df
  delta = surv_df_notch_oc100$STATUS, # Event indicator from df
  marker = surv_df_notch_oc100$DLL1, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = 60,    # Time points for ROC
  iid = TRUE      )   # Compute confidence intervals

compare(roc_result_model_notch_met_60, roc_result_model_DLL1_60, adjusted = F) #0.0158335    
#notch1
roc_result_model_NOTCH1_60 <- timeROC(
  T = surv_df_notch_oc100$OS,       # Survival time from df
  delta = surv_df_notch_oc100$STATUS, # Event indicator from df
  marker = surv_df_notch_oc100$NOTCH1, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = 60,    # Time points for ROC
  iid = TRUE      )   # Compute confidence intervals

compare(roc_result_model_notch_met_60, roc_result_model_NOTCH1_60, adjusted = F) #0.252605  

#notch2
roc_result_model_NOTCH2_60 <- timeROC(
  T = surv_df_notch_oc100$OS,       # Survival time from df
  delta = surv_df_notch_oc100$STATUS, # Event indicator from df
  marker = surv_df_notch_oc100$NOTCH2, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = 60,    # Time points for ROC
  iid = TRUE      )   # Compute confidence intervals

compare(roc_result_model_notch_met_60, roc_result_model_NOTCH2_60, adjusted = F) #0.0256091 

#notch3
roc_result_model_NOTCH3_60 <- timeROC(
  T = surv_df_notch_oc100$OS,       # Survival time from df
  delta = surv_df_notch_oc100$STATUS, # Event indicator from df
  marker = surv_df_notch_oc100$NOTCH3, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = 60,    # Time points for ROC
  iid = TRUE      )   # Compute confidence intervals

compare(roc_result_model_notch_met_60, roc_result_model_NOTCH3_60, adjusted = F) #0.2970238  

#notch4
roc_result_model_NOTCH4_60 <- timeROC(
  T = surv_df_notch_oc100$OS,       # Survival time from df
  delta = surv_df_notch_oc100$STATUS, # Event indicator from df
  marker = surv_df_notch_oc100$NOTCH4, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = 60,    # Time points for ROC
  iid = TRUE      )   # Compute confidence intervals

compare(roc_result_model_notch_met_60, roc_result_model_NOTCH4_60, adjusted = F) #0.01927686 

#ctnnb1
roc_result_model_CTNNB1_60 <- timeROC(
  T = surv_df_notch_oc100$OS,       # Survival time from df
  delta = surv_df_notch_oc100$STATUS, # Event indicator from df
  marker = surv_df_notch_oc100$CTNNB1, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = 60,    # Time points for ROC
  iid = TRUE      )   # Compute confidence intervals

compare(roc_result_model_notch_met_60, roc_result_model_CTNNB1_60, adjusted = F) #0.003989765     
#fbxw7
roc_result_model_FBXW7_60 <- timeROC(
  T = surv_df_notch_oc100$OS,       # Survival time from df
  delta = surv_df_notch_oc100$STATUS, # Event indicator from df
  marker = surv_df_notch_oc100$FBXW7, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = 60,    # Time points for ROC
  iid = TRUE      )   # Compute confidence intervals

compare(roc_result_model_notch_met_60, roc_result_model_FBXW7_60, adjusted = F) #0.04926349   

#ARID1A EXPRESSION
roc_result_model_ARID1A_60 <- timeROC(
  T = surv_df_notch_oc100$OS,       # Survival time from df
  delta = surv_df_notch_oc100$STATUS, # Event indicator from df
  marker = surv_df_notch_oc100$ARID1A, # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = 60,    # Time points for ROC
  iid = TRUE      )   # Compute confidence intervals

compare(roc_result_model_notch_met_60, roc_result_model_ARID1A_60, adjusted = F) #0.1225441     


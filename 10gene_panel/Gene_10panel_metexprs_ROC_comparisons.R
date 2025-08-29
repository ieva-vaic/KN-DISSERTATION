#KN-DISSERTATION project. 10 gene panel - genes selected via statistical analysis of TCGA and GTEx data
#ROC tests compare 10 gene panel with NOTCH, WNT and ARID1A gene panel
Sys.setenv(LANG = "en")
#libraries
library(tidyverse)
library(pROC)
library(glmnet)
library(gtsummary)
library(gt)
library(grid)
library(brglm2)
library(reshape2)
library(rstatix) 
library(ggprism)
library(gridExtra)
library(png)
library(htmlwidgets)
library(webshot)
library(magick)
library(openxlsx)
#set directory of outputs
setwd("C:/Users/Ieva/rprojects/outputs_all/DISS/")
#data METXPRS###############################################
KN_data <- readRDS("../../OTHER DATA/KN-DISSERTATION FILES/KN_data1114_essential.rds")
colnames(KN_data)
#biomarker groups
#expression
raiska <- colnames(KN_data[18:27])
#methylation
metilinimas <- colnames(KN_data[28:31])
biomarkers <- c(raiska, metilinimas)
#change tumor
KN_data$tumor <- recode(KN_data$tumor, OvCa = "OC", Benign ="Benign")
#methylation data must be 0 1
KN_data$HOPX <- as.numeric(KN_data$HOPX) -1
KN_data$CDX2 <- as.numeric(KN_data$CDX2) -1
KN_data$ALX4 <- as.numeric(KN_data$ALX4) -1
KN_data$ARID1A_met <- as.numeric(KN_data$ARID1A_met) -1

#data 10 genes###################################################################
OC_full <- readRDS("../../OTHER DATA/KN-DISSERTATION FILES/OC_10_genes_clean_2025_02_14.RDS")
#expression genes list
expression <- c("EXO1", "RAD50","PPT2", "LUC7L2","PKP3", "CDCA5","ZFPL1","VPS33B", "GRB7","TCEAL4")
#make HGSOC vs BENIGN df
OC_HGSOC_BENIGN<- OC_full[c(OC_full$type != "Other"),] #51 cases left
table(OC_HGSOC_BENIGN$type) #51 left 42 vs 9 OK
OC_HGSOC_BENIGN$tumor <- relevel(factor(OC_HGSOC_BENIGN$type), ref = "Benign")
#make HGSOC vs OTHERS df
OC_HGSOC_OTHERS<- OC_full[c(OC_full$type != "Benign"),] #56 cases left
OC_HGSOC_OTHERS$tumor <- relevel(factor(OC_HGSOC_OTHERS$type), ref = "Other")
table(OC_HGSOC_OTHERS$tumor) #51 left 42 vs 13 OK
table(OC_HGSOC_OTHERS$tumor, OC_HGSOC_OTHERS$CA125_f)

#combine expression datasets#####################################
#fix ids
KN_data$KN <- sub("KN-0*([1-9][0-9]*)$", "KN-\\1", KN_data$patient_id_aud)
KN_data$KN %in% OC_full$KN
#remove clinical from 10 gene
OC_data_expression <- OC_full[, c(1, 3:13)]
#join
ALL_EXPRESSION_DF <- left_join(KN_data, OC_data_expression, by = "KN")
#save rds for later
saveRDS(ALL_EXPRESSION_DF, "../../OTHER DATA/KN-DISSERTATION FILES/ALL_expresssion_Df20250709.RDS")
#save exel for later
saveRDS(ALL_EXPRESSION_DF, "../../OTHER DATA/KN-DISSERTATION FILES/ALL_expresssion_Df20250709.xlsx")
#ROC HGSOC vs BENIGN######################################################################
# HGSOC vs benign df
KN_BENIGN_HGSOC <- ALL_EXPRESSION_DF[ALL_EXPRESSION_DF$Grupė_Ieva != "Other", ] 
KN_BENIGN_HGSOC$Grupė_Ieva <- droplevels(KN_BENIGN_HGSOC$Grupė_Ieva)
table(KN_BENIGN_HGSOC$Grupė_Ieva) #51 left 42 vs 9
#tumor must be a factor with reference: in this command first  level is a reference
KN_BENIGN_HGSOC$Grupė_Ieva <- factor(KN_BENIGN_HGSOC$Grupė_Ieva, levels = c("Benign", "HGSOC"))
#all biomarers in one place:
biomarkers2 <- c(biomarkers, expression)
#ROC: HGSOC vs BENIGN
roc_results_tumor_bh<- lapply(biomarkers2, function(col) {
  roc(response = KN_BENIGN_HGSOC$Grupė_Ieva, predictor = KN_BENIGN_HGSOC[[col]])
})
names(roc_results_tumor_bh) <- biomarkers2
roc_results_tumor_bh
#extract the aucs
auc_values_tumor_bh <- sapply(roc_results_tumor_bh, function(roc_obj) {auc(roc_obj)})
auc_values_tumor_bh #extracted aucs
##test best rocs for diferences###################################
#best biomarkers together
roc.test(roc_results_tumor_bh[["CTNNB1"]], roc_results_tumor_bh[["GRB7"]], method=c("delong"))# 0.5498 
roc.test(roc_results_tumor_bh[["CTNNB1"]], roc_results_tumor_bh[["TCEAL4"]])#  0.4598 
#all other METEXPRS expression biomarkers with GRB7
roc.test(roc_results_tumor_bh[["NOTCH4"]], roc_results_tumor_bh[["GRB7"]])#  0.1944 
roc.test(roc_results_tumor_bh[["NOTCH2"]], roc_results_tumor_bh[["GRB7"]])#  0.09918 
roc.test(roc_results_tumor_bh[["NOTCH1"]], roc_results_tumor_bh[["GRB7"]])#  0.01049* 1
roc.test(roc_results_tumor_bh[["NOTCH3"]], roc_results_tumor_bh[["GRB7"]])# 0.02127* 2
roc.test(roc_results_tumor_bh[["ARID1A"]], roc_results_tumor_bh[["GRB7"]])#0.0004514* 3
roc.test(roc_results_tumor_bh[["FBXW7"]], roc_results_tumor_bh[["GRB7"]])#0.3658
roc.test(roc_results_tumor_bh[["JAG2"]], roc_results_tumor_bh[["GRB7"]])#1.871e-08*** 4
roc.test(roc_results_tumor_bh[["DLL1"]], roc_results_tumor_bh[["GRB7"]])#0.04927* 5
roc.test(roc_results_tumor_bh[["HES1"]], roc_results_tumor_bh[["GRB7"]])#0.2699
#all METEXPRS methylation biomarkers with GRB7
roc.test(roc_results_tumor_bh[["HOPX"]], roc_results_tumor_bh[["GRB7"]])#4.537e-12**** 6
roc.test(roc_results_tumor_bh[["ALX4"]], roc_results_tumor_bh[["GRB7"]])#0.0007431*** 7 
roc.test(roc_results_tumor_bh[["ARID1A_met"]], roc_results_tumor_bh[["GRB7"]])#2.459e-06**** 8 
roc.test(roc_results_tumor_bh[["CDX2"]], roc_results_tumor_bh[["GRB7"]])#6.512e-05**** 9
#all 10 gene panel bioamarkers with GRB7
roc.test(roc_results_tumor_bh[["RAD50"]], roc_results_tumor_bh[["GRB7"]])#  0.02611* 1
roc.test(roc_results_tumor_bh[["EXO1"]], roc_results_tumor_bh[["GRB7"]])#  0.05898
roc.test(roc_results_tumor_bh[["PPT2"]], roc_results_tumor_bh[["GRB7"]])#  0.09991
roc.test(roc_results_tumor_bh[["LUC7L2"]], roc_results_tumor_bh[["GRB7"]])# 0.03828* 2
roc.test(roc_results_tumor_bh[["PKP3"]], roc_results_tumor_bh[["GRB7"]])# 0.494
roc.test(roc_results_tumor_bh[["CDCA5"]], roc_results_tumor_bh[["GRB7"]])#0.1037
roc.test(roc_results_tumor_bh[["ZFPL1"]], roc_results_tumor_bh[["GRB7"]])#0.0219* 3
roc.test(roc_results_tumor_bh[["VPS33B"]], roc_results_tumor_bh[["GRB7"]])#0.05743

# HGSOC vs OTHERS ###################################################################
#make df HGSOC vs OTHRS
KN_OTHERS_HGSOC <- ALL_EXPRESSION_DF[ALL_EXPRESSION_DF$Grupė_Ieva != "Benign", ] 
KN_OTHERS_HGSOC$Grupė_Ieva <- droplevels(KN_OTHERS_HGSOC$Grupė_Ieva)
table(KN_OTHERS_HGSOC$Grupė_Ieva) 
#tumor must be a factor with reference: in this command first  level is a reference
KN_OTHERS_HGSOC$Grupė_Ieva <- factor(KN_OTHERS_HGSOC$Grupė_Ieva, levels = c("Other", "HGSOC"))
#ROC: HGSOC vs OTHERS
roc_results_tumor_oh<- lapply(biomarkers2, function(col) {
  roc(response = KN_OTHERS_HGSOC$Grupė_Ieva, predictor = KN_OTHERS_HGSOC[[col]])
})
names(roc_results_tumor_oh) <- biomarkers2
roc_results_tumor_oh
#extract the aucs
auc_values_tumor_oh <- sapply(roc_results_tumor_oh, function(roc_obj) {auc(roc_obj)})
auc_values_tumor_oh #extracted aucs
##test best rocs for diferences###################################
#Best METEXPRS biomarker HES1 vs others in this panel
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["NOTCH1"]])# 0.4564553
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["NOTCH2"]])#  0.5102932
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["NOTCH3"]])#  0.05582
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["NOTCH4"]])#  0.1577
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["ARID1A"]])#  0.01608* 1
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["FBXW7"]])#0.0555
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["CTNNB1"]])#0.03293* 2
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["JAG2"]])#0.0424* 3
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["DLL1"]])#0.04305* 4
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["CDX2"]])#0.0182* 
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["HOPX"]])#0.003272** 
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["ALX4"]])#0.0007158*** 
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["ARID1A_met"]])#0.02057* 
#Best METEXPRS biomarker HES1 vs 10 gene panel
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["EXO1"]])# 0.8341 
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["RAD50"]])#  0.0002023* 5
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["PPT2"]])#  0.02664* 6
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["LUC7L2"]])#  0.02664* 7
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["PKP3"]])#  0.004837** 8
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["CDCA5"]])#0.4053
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["ZFPL1"]])#0.0004771*** 9
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["VPS33B"]])#0.0336* 10
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["GRB7"]])#0.0336* 11
roc.test(roc_results_tumor_oh[["HES1"]], roc_results_tumor_oh[["TCEAL4"]])#0.7452

#KN-DISSERTATION project. MET_EXPRS - methylation, expression, mutation data
#ROC analisies, compare ROCs
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
#setwd for plots
setwd("C:/Users/Ieva/rprojects/outputs_all/DISS/")
#DATA - FROM THE PUBLICATION MET_EXPRS
KN_data <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/KN_data1114_essential.rds")
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
#ADD ROC TESTS#######################################################
#form other scrips
load("roc_list_separate_biomarkers20250826.RData")
load("roc_list20250826.RData")
#COMPARE TESTS#####################################################
#best least biomarker in HGSOC - NOTCH combination #roc_curve2.3
##best single biomarker combination with separate markers ##################
roc.test(roc_curve2.3, roc_results_tumor_bh[["CTNNB1"]], method = "delong")
roc.test(roc_curve2, roc_results_tumor_bh[["CTNNB1"]], method = "delong")
roc.test(roc_curve_CA2, roc_results_tumor_bh[["CTNNB1"]], method = "delong")#0.08346 
roc.test(roc_curve_CA2, roc_results_tumor_bh[["ALX4"]], method = "delong")#0.4556 #methyl
roc.test(roc_results_tumor_bh[["CTNNB1"]], roc_results_tumor_bh[["ALX4"]], method = "delong")#8.256e-05
##OVCa vs benign comparison########################Ę
set.seed(1)
roc.test(roc_results_tumor[["CTNNB1"]], roc_curve1)#10 genes vs ctnnb1
roc.test(roc_curve1, roc_results_tumor[["NOTCH1"]])# 0.00213
roc.test(roc_curve1, roc_results_tumor[["NOTCH2"]])#0.07008
roc.test(roc_curve1, roc_results_tumor[["NOTCH3"]])#0.00849
roc.test(roc_curve1, roc_results_tumor[["NOTCH4"]])#0.01045
roc.test(roc_curve1, roc_results_tumor[["ARID1A"]])#0.000133
roc.test(roc_curve1, roc_results_tumor[["FBXW7"]])#0.008917
roc.test(roc_curve1, roc_results_tumor[["JAG2"]])#8.028e-08
roc.test(roc_curve1, roc_results_tumor[["DLL1"]])#0.005312
roc.test(roc_curve1, roc_results_tumor[["HES1"]])#0.000753
roc.test(roc_curve1, roc_results_tumor[["HOPX"]])#2.344e-13
roc.test(roc_curve1, roc_results_tumor[["ALX4"]])#1.416e-05
roc.test(roc_curve1, roc_results_tumor[["ARID1A_met"]])#1.021e-06
roc.test(roc_curve1, roc_results_tumor[["CDX2"]])#6.021e-05

##methylation combination with separate markers
roc.test(roc_curve1.2, roc_results_tumor[["HOPX"]])#0.0009919
roc.test(roc_curve1.2, roc_results_tumor[["ALX4"]])#0.01731
roc.test(roc_curve1.2, roc_results_tumor[["ARID1A_met"]])#0.0002838
roc.test(roc_curve1.2, roc_results_tumor[["CDX2"]])#0.03982

##hgsoc vs benign comaprisons #############################
roc.test(roc_results_tumor_bh[["CTNNB1"]], roc_curve2)#10 genes vs ctnnb1 0.1445
roc.test(roc_curve2, roc_results_tumor_bh[["NOTCH1"]])#  0.003041
roc.test(roc_curve2, roc_results_tumor_bh[["NOTCH2"]])#0.04894
roc.test(roc_curve2, roc_results_tumor_bh[["NOTCH3"]])#0.003567
roc.test(roc_curve2, roc_results_tumor_bh[["NOTCH4"]])#0.03884
roc.test(roc_curve2, roc_results_tumor_bh[["ARID1A"]])#0.0003837
roc.test(roc_curve2, roc_results_tumor_bh[["FBXW7"]])#0.09913 ###
roc.test(roc_curve2, roc_results_tumor_bh[["JAG2"]])#8.776e-08
roc.test(roc_curve2, roc_results_tumor_bh[["DLL1"]])# 0.02657
roc.test(roc_curve2, roc_results_tumor_bh[["HES1"]])#0.09506 ####
roc.test(roc_curve2, roc_results_tumor_bh[["HOPX"]])#3.284e-16
roc.test(roc_curve2, roc_results_tumor_bh[["ALX4"]])#6.799e-06
roc.test(roc_curve2, roc_results_tumor_bh[["ARID1A_met"]])#5.825e-08
roc.test(roc_curve2, roc_results_tumor_bh[["CDX2"]])#6.799e-06

##methylation combination with separate markers
roc.test(roc_curve2.2, roc_results_tumor_bh[["HOPX"]])#0.000614
roc.test(roc_curve2.2, roc_results_tumor_bh[["ALX4"]])#0.001225
roc.test(roc_curve2.2, roc_results_tumor_bh[["ARID1A_met"]])#0.001099
roc.test(roc_curve2.2, roc_results_tumor_bh[["CDX2"]])#0.004731

#HGSOC vs OTHERS comparisons ######################################
roc.test(roc_results_tumor_oh[["CTNNB1"]], roc_curvex)#10 genes vs ctnnb1 0.005265
roc.test(roc_curvex, roc_results_tumor_oh[["NOTCH1"]])#  1.709e-05
roc.test(roc_curvex, roc_results_tumor_oh[["NOTCH2"]])#4.271e-06
roc.test(roc_curvex, roc_results_tumor_oh[["NOTCH3"]])#0.0004781
roc.test(roc_curvex, roc_results_tumor_oh[["NOTCH4"]])#0.003541
roc.test(roc_curvex, roc_results_tumor_oh[["ARID1A"]])# 1.033e-05
roc.test(roc_curvex, roc_results_tumor_oh[["FBXW7"]])#0.00937 
roc.test(roc_curvex, roc_results_tumor_oh[["JAG2"]])#0.0001005
roc.test(roc_curvex, roc_results_tumor_oh[["DLL1"]])# 0.01397
roc.test(roc_curvex, roc_results_tumor_oh[["HES1"]])#0.1886 ####
roc.test(roc_curvex, roc_results_tumor_oh[["HOPX"]])#7.36e-07
roc.test(roc_curvex, roc_results_tumor_oh[["ALX4"]])#2.88e-06
roc.test(roc_curvex, roc_results_tumor_oh[["ARID1A_met"]])#3.826e-06
roc.test(roc_curvex, roc_results_tumor_oh[["CDX2"]])#4.584e-05

#with all biomarkers#######################################################
roc.test(roc_results_tumor_oh[["CTNNB1"]], roc_curvex.3)#10 genes vs ctnnb1  0.002442
roc.test(roc_curvex.3, roc_results_tumor_oh[["NOTCH1"]])#  2.05e-06
roc.test(roc_curvex.3, roc_results_tumor_oh[["NOTCH2"]])#2.409e-06
roc.test(roc_curvex.3, roc_results_tumor_oh[["NOTCH3"]])#0.0001596
roc.test(roc_curvex.3, roc_results_tumor_oh[["NOTCH4"]])#0.001745
roc.test(roc_curvex.3, roc_results_tumor_oh[["ARID1A"]])# 3.194e-06
roc.test(roc_curvex.3, roc_results_tumor_oh[["FBXW7"]])#0.002974 
roc.test(roc_curvex.3, roc_results_tumor_oh[["JAG2"]])#3.384e-05
roc.test(roc_curvex.3, roc_results_tumor_oh[["DLL1"]])# 0.004916
roc.test(roc_curvex.3, roc_results_tumor_oh[["HES1"]])#0.1006 ####
roc.test(roc_curvex.3, roc_results_tumor_oh[["HOPX"]])#1.486e-07
roc.test(roc_curvex.3, roc_results_tumor_oh[["ALX4"]])#6.277e-07
roc.test(roc_curvex.3, roc_results_tumor_oh[["ARID1A_met"]])#4.452e-07
roc.test(roc_curvex.3, roc_results_tumor_oh[["CDX2"]])#2.6e-06


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
##MAKE PRESENTATION PLOTS#############################################
#ROC COMBINATIONS HGSOC###################

##PLOT ROC COMBINATIONS HGSOC ##########################
#HGSOC vs benign models
roc_plot_4 <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_curve2, print.auc = F, col = "#911eb4", lty = 2,
           cex.main=0.8, main ="Gerybinių pakitimų atskyrimas nuo HGSOC atvejų",
           xlab = "Specifiškumas",   # Custom x-axis label (e.g., in Lithuanian)
           ylab = "Jautrumas") #7
  lines(roc_curve2.2, col = "#dcbeff", lwd =2 ) #6
  lines(roc_curve2.3, col ="#fabed4", lwd =2, lty = 4) #8
  lines(roc_curve2.4, col ="darkred", lwd =2, lty = 3) 
  lines(roc_curve2.5, col ="deeppink", lwd =2) 
  lines(roc_results_tumor_bh$CTNNB1, col ="darkgreen", lwd =2)  #ctnnb1
  lines(roc_results_tumor_bh$FBXW7, col ="darkblue", lwd =2) #FBXW7
  lines(roc_results_tumor_bh$HES1, col ="lightblue", lwd =2) #NOTCH2
  lines(roc_curve_CA2, col = "grey", lwd = 2)
  
  # Add legend
  legend("bottomright", legend = c( expression(italic("Genų raiškos žymenų kombinacija")), 
                                    expression(italic("Promotorių metilinimo statuso kombinacija")),
                                    expression(italic("Genų raiškos žymenų + promotorių metilinimo statuso kombinacija ")), 
                                    expression(italic("NOTCH genų raiškos žymenų kombinacija")),
                                    expression(italic("HOX promotorių metilinimo statuso kombinacija")),
                                    
                                    expression(italic("CTNNB1 raiška")),
                                    expression(italic("FBXW7 raiška")),
                                    expression(italic("HES1 raiška")),
                                    
                                    expression(italic("Serumo CA125 biožymens statusas")))
         ,
         col = c("#911eb4","#dcbeff", "#fabed4", "darkred",
                 "deeppink","darkgreen", "darkblue","lightblue", "grey" ), lty = 1, 
         cex = 0.8, lwd =3)
}
# Save the plot as a PNG file
png("metexprs_roc_HGSOC_output20251020.png", width = 1000, height = 1000, res = 150)
roc_plot_4()
dev.off()

##TABLE ROC COMBINATIONS HGSOC #####################
#make coords first for separate biomarkers
##Table for ROC metrics, HGSOC vs benign#################################
#extract the aucs
auc_values_tumor_bh <- sapply(roc_results_tumor_bh, function(roc_obj) {auc(roc_obj)})
auc_values_tumor_bh #extracted aucs
coords_results_tumor_bh <- lapply(roc_results_tumor_bh, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv", "tpr", "fpr"), transpose = FALSE)
})
coords_results_tumor_bh
coords_results_tumor_bh$ARID1A_met #arid1a methylation have infinate thresholds 
# Create a dataframe combining AUC values and coordinates results
results_tumor_bh<- data.frame(
  Predictor = c(raiska, "HOPX", "ALX4","CDX2" ),
  AUC = auc_values_tumor_bh[1:13],
  do.call(rbind, coords_results_tumor_bh[1:13]) 
)
rownames(results_tumor_bh) <- NULL
results_tumor_bh

#make combination coords
#ca125
coords_ca2 <- coords(roc_curve_CA2, "best", 
                     ret=c("threshold", "accuracy", "sensitivity", "specificity",
                           "precision", "npv", "tpr", "fpr"), transpose = FALSE)

#coords2
coords2 <- coords(roc_curve2, "best", 
                    ret=c("threshold", "accuracy", "sensitivity", "specificity",
                          "precision", "npv", "tpr", "fpr"), transpose = FALSE)
#coords2.2
coords2.2  <- coords(roc_curve2.2, "best", 
                    ret=c("threshold", "accuracy", "sensitivity", "specificity",
                          "precision", "npv", "tpr", "fpr"), transpose = FALSE)
#coords2.3
coords2.3 <- coords(roc_curve2.3, "best", 
                      ret=c("threshold", "accuracy", "sensitivity", "specificity",
                            "precision", "npv", "tpr", "fpr"), transpose = FALSE)
#coords2.4
coords2.4 <- coords(roc_curve2.4, "best", 
                      ret=c("threshold", "accuracy", "sensitivity", "specificity",
                            "precision", "npv", "tpr", "fpr"), transpose = FALSE)

#coords2.5
coords2.5 <- coords(roc_curve2.5, "best", 
                      ret=c("threshold", "accuracy", "sensitivity", "specificity",
                            "precision", "npv", "tpr", "fpr"), transpose = FALSE)
#MAKE ONE DF
results_roc2<- data.frame(
  Biožymenys = c("Genų raiškos žymenų kombinacija", 
                 "Promotorių metilinimo statuso kombinacija",
                 "Genų raiškos žymenų + promotorių metilinimo statuso kombinacija", 
                 "NOTCH genų raiškos žymenų kombinacija",
                 "HOX promotorių metilinimo statuso kombinacija",
                 "CTNNB1 raiška",
                 "FBXW7 raiška",
                 "HES1 raiška",
                 "Serumo CA125 biožymens statusas"),
  `plotas po kreive` = c(roc_curve2$auc, roc_curve2.2$auc, roc_curve2.3$auc, roc_curve2.4$auc, roc_curve2.5$auc,
                         roc_results_tumor_bh[["CTNNB1"]]$auc, roc_results_tumor_bh[["FBXW7"]]$auc,roc_results_tumor_bh[["HES1"]]$auc,
                         roc_curve_CA2$auc), 
  `slenkstinė vertė` = c(coords2$threshold, coords2.2$threshold , coords2.3$threshold,
                         coords2.4$threshold, coords2.5$threshold, 
                         coords_results_tumor_bh[["CTNNB1"]]$threshold,
                         coords_results_tumor_bh[["FBXW7"]]$threshold,
                         coords_results_tumor_bh[["HES1"]]$threshold,
                         coords_ca2$threshold ),
  tikslumas = c(coords2$accuracy, coords2.2$accuracy , coords2.3$accuracy,
                coords2.4$accuracy, coords2.5$accuracy,
                coords_results_tumor_bh[["CTNNB1"]]$accuracy,
                coords_results_tumor_bh[["FBXW7"]]$accuracy,
                coords_results_tumor_bh[["HES1"]]$accuracy,
                coords_ca2$accuracy ),
  jautrumas = c(coords2$sensitivity, coords2.2$sensitivity, coords2.3$sensitivity,
                coords2.4$sensitivity, coords2.5$sensitivity,
                coords_results_tumor_bh[["CTNNB1"]]$sensitivity,
                coords_results_tumor_bh[["FBXW7"]]$sensitivity,
                coords_results_tumor_bh[["HES1"]]$sensitivity,
                coords_ca2$sensitivity),
  specifiškumas = c(coords2$specificity, coords2.2$specificity, coords2.3$specificity,
                    coords2.4$specificity, coords2.5$specificity,
                    coords_results_tumor_bh[["CTNNB1"]]$specificity,
                    coords_results_tumor_bh[["FBXW7"]]$specificity,
                    coords_results_tumor_bh[["HES1"]]$specificity,
                    coords_ca2$specificity),
  ppv  = c(coords2$precision, coords2.2$precision, coords2.3$precision,
           coords2.4$precision,coords2.5$precision,
           coords_results_tumor_bh[["CTNNB1"]]$precision,
           coords_results_tumor_bh[["FBXW7"]]$precision,
           coords_results_tumor_bh[["HES1"]]$precision,
           coords_ca2$precision ),
  npv  = c(coords2$npv, coords2.2$npv, coords2.3$npv,
           coords2.4$npv, coords2.5$npv, 
           coords_results_tumor_bh[["CTNNB1"]]$npv,
           coords_results_tumor_bh[["FBXW7"]]$npv,
           coords_results_tumor_bh[["HES1"]]$npv,
           coords_ca2$npv),
  tpr  = c(coords2$tpr, coords2.2$tpr, coords2.3$tpr,
           coords2.4$tpr,coords2.5$tpr,
           coords_results_tumor_bh[["CTNNB1"]]$tpr,
           coords_results_tumor_bh[["FBXW7"]]$tpr,
           coords_results_tumor_bh[["HES1"]]$tpr,
           coords_ca2$tpr),
  fpr  = c(coords2$fpr, coords2.2$fpr, coords2.3$fpr,
           coords2.4$fpr, coords2.5$fpr,
           coords_results_tumor_bh[["CTNNB1"]]$fpr,
           coords_results_tumor_bh[["FBXW7"]]$fpr,
           coords_results_tumor_bh[["HES1"]]$fpr,
           coords_ca2$fpr),
  check.names = FALSE
)
rownames(results_roc2) <- NULL
results_roc2

gt_table2 <- results_roc2 %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai", 
    subtitle = "Gerybinių pakitimų atskyrimas nuo HGSOC atvejų",
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymenys))
  )
#show
gt_table2

#there is no other convenient way to save gt outputs
gtsave(gt_table2,vwidth = 800,
       filename = "metexprs_table_HGSOC_output20251020.png")

#Combine the images
roc_image2<- image_read("metexprs_roc_HGSOC_output20251020.png")
table_image2 <- image_read("metexprs_table_HGSOC_output20251020.png")

# Now append vertically
combined_image2 <- image_append(c(roc_image2, table_image2), stack = FALSE)

# Save the combined image
image_write(combined_image2, 
            "metexprs_Roctable_HGSOC_output20251020.png")


#PLOT ROC COMBINATIONS HGSOC vs OTHERS##########################
roc_plot_5 <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_curvex, print.auc = F, col = "#911eb4", 
           cex.main=0.8, main ="HGSOC atskyrimas nuo kitų KV atvejų",
           xlab = "Specifiškumas",   # Custom x-axis label (e.g., in Lithuanian)
           ylab = "Jautrumas") #7
  lines(roc_curvex.1, col = "#dcbeff", lwd =2 ) #6
  lines(roc_curvex.3, col ="#fabed4", lwd =2) #8
  lines(roc_curvex.4, col ="darkred", lwd =2 ) 
  lines(roc_curveX.5, col ="deeppink", lwd =2) 
  lines(roc_results_tumor_oh[["HES1"]], col = "#808000", lwd =2) #2-3
  lines(roc_curve_CA2X, col = "grey", lwd = 2)
  
  # Add legend
  legend("bottomright", legend = c( expression(italic("Genų raiškos žymenų kombinacija")), 
                                    expression(italic("Promotorių metilinimo statuso kombinacija")),
                                    expression(italic("Genų raiškos žymenų + promotorių metilinimo statuso kombinacija ")), 
                                    expression(italic("NOTCH genų raiškos žymenų kombinacija")),
                                    expression(italic("HOX promotorių metilinimo statuso kombinacija")),
                                    expression(italic("HES1 raiška")),
                                    expression(italic("Serumo CA125 biožymens statusas")))
         ,
         col = c("#911eb4","#dcbeff", "#fabed4", "darkred",
                 "deeppink", "#808000", "grey" ), lty = 1, 
         cex = 0.73, lwd =3)
}

# Save the plot as a PNG file
png("metexprs_roc_HGSOC_OTHERS_MODELS_output20251013.png", 
    width = 1000, height = 1000, res = 150)
roc_plot_5()
dev.off()

##make cooords hgsoc vs others #########################
coordsx <- coords(roc_curvex, "best", 
                  ret=c("threshold", "accuracy", "sensitivity", "specificity",
                        "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coordsx.1 <- coords(roc_curvex.1, "best", 
                   ret=c("threshold", "accuracy", "sensitivity", "specificity",
                         "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coordsx.3 <- coords(roc_curvex.3, "best", 
                   ret=c("threshold", "accuracy", "sensitivity", "specificity",
                         "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coordsx.4 <- coords(roc_curvex.4, "best", 
                   ret=c("threshold", "accuracy", "sensitivity", "specificity",
                         "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coordsX.5 <- coords(roc_curveX.5, "best", 
                    ret=c("threshold", "accuracy", "sensitivity", "specificity",
                          "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords_ca2X <- coords(roc_curve_CA2X, "best", 
                      ret=c("threshold", "accuracy", "sensitivity", "specificity",
                            "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords_HES <- coords(roc_results_tumor_oh[["HES1"]], "best", 
                      ret=c("threshold", "accuracy", "sensitivity", "specificity",
                            "precision", "npv", "tpr", "fpr"), transpose = FALSE)

##TABLE ROC COMBINATIONS HGSOC vs OTHERS#####################
results_rocx<- data.frame(
  Biožymenys = c("Genų raiškos žymenų kombinacija", 
                 "Promotorių metilinimo statuso kombinacija",
                 "Genų raiškos žymenų + promotorių metilinimo statuso kombinacija", 
                 "NOTCH genų raiškos žymenų kombinacija",
                 "HOX promotorių metilinimo statuso kombinacija",
                 "HES1 raiška",
                 "Serumo CA125 biožymens statusas"),
  `plotas po kreive` = c(roc_curvex$auc, roc_curvex.1$auc, roc_curvex.3$auc, 
                         roc_curvex.4$auc, roc_curveX.5$auc,roc_results_tumor_oh[["HES1"]]$auc,
                         roc_curve_CA2X$auc), 
  `slenkstinė vertė` = c(coordsx$threshold, coordsx.1$threshold , coordsx.3$threshold,
                         coordsx.4$threshold, coordsX.5$threshold,
                         coords_HES$threshold,
                         coords_ca2X$threshold ),
  tikslumas = c(coordsx$accuracy, coordsx.1$accuracy , coordsx.3$accuracy,
                coordsx.4$accuracy, coordsX.5$accuracy,
                coords_HES$accuracy,
                coords_ca2X$accuracy ),
  jautrumas = c(coordsx$sensitivity, coordsx.1$sensitivity, coordsx.3$sensitivity,
                coordsx.4$sensitivity, coordsX.5$sensitivity,coords_HES$sensitivity, coords_ca2X$sensitivity),
  specifiškumas = c(coordsx$specificity, coordsx.1$specificity, coordsx.3$specificity,
                    coordsx.4$specificity, coordsX.5$specificity,coords_HES$specificity, coords_ca2X$specificity),
  ppv  = c(coordsx$precision, coordsx.1$precision, coordsx.3$precision,
           coordsx.4$precision,coordsX.5$precision,coords_HES$precision, coords_ca2X$precision ),
  npv  = c(coordsx$npv, coordsx.1$npv, coordsx.3$npv,
           coordsx.4$npv, coordsX.5$npv, coords_HES$npv, coords_ca2X$npv),
  tpr  = c(coordsx$tpr, coordsx.1$tpr, coordsx.3$tpr,
           coordsx.4$tpr,coordsX.5$tpr,coords_HES$tpr, coords_ca2X$tpr),
  fpr  = c(coordsx$fpr, coordsx.1$fpr,  coordsx.3$fpr,
           coordsx.4$fpr, coordsX.5$fpr, coords_HES$fpr, coords_ca2X$fpr),
  check.names = FALSE
)
rownames(results_rocx) <- NULL
results_rocx

gt_tablex <- results_rocx %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai",
    subtitle = "HGSOC atskyrimas nuo kitų KV atvejų",
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymenys))
  )
#show
gt_tablex

#there is no other convieneat way to save gt outputs
gtsave(gt_tablex,vwidth = 800,
       filename = "metexprs_table_HGSOC_OTHERS_MODELS_output20251013.png")

#Combine the images
roc_image2<- image_read("metexprs_roc_HGSOC_OTHERS_MODELS_output20251013.png")
table_image2 <- image_read("metexprs_table_HGSOC_OTHERS_MODELS_output20251013.png")

# Now append horizontaly
combined_image2 <- image_append(c(roc_image2, table_image2), stack = F)

# Save the combined image
image_write(combined_image2, 
            "metexprs_tableroc_HGSOC_OTHERS_MODELS_output202501013.png")
#STATISTICAL MODEL GENES #################################
OC_full <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES//OC_10_genes_clean_2025_02_14.RDS")
OC_full <- OC_full[c(OC_full$KN != "KN-100"), ]
#expression genes list
expression <- c("EXO1", "RAD50","PPT2", "LUC7L2","PKP3", "CDCA5","ZFPL1","VPS33B", "GRB7","TCEAL4")
#HGSOC vs BENIGN df
OC_HGSOC_BENIGN<- OC_full[c(OC_full$type != "Other"),] #51 cases left
OC_HGSOC_BENIGN$tumor <- relevel(factor(OC_HGSOC_BENIGN$type), ref = "Benign")
#HGSOC vs OTHERS df
OC_HGSOC_OTHERS<- OC_full[c(OC_full$type != "Benign"),] #56 cases left
OC_HGSOC_OTHERS$tumor <- relevel(factor(OC_HGSOC_OTHERS$type), ref = "Other")
table(OC_HGSOC_OTHERS$tumor) #51 left 42 vs 14
table(OC_HGSOC_OTHERS$tumor, OC_HGSOC_OTHERS$CA125_f) #chek CA125
##all 10 HGSOC VS BENIGN combinations ###############################################
expr_tumor10 <- OC_HGSOC_BENIGN[colnames(OC_HGSOC_BENIGN) %in% expression]
brglm.model_10 <- glm(OC_HGSOC_BENIGN$tumor ~ ., data = expr_tumor10,
                      family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_10 <- predict.glm(brglm.model_10, type='response') 
pred_data10 <- OC_HGSOC_BENIGN[(rownames(OC_HGSOC_BENIGN) #remove incoplete rows
                                %in% names(predicted_probs_10)), ]
roc_curve10 <- roc(pred_data10$tumor, predicted_probs_10)
AUC10 <- auc(roc_curve10)
coords10 <- coords(roc_curve10,
                   "best", ret=c("threshold", "accuracy", "sensitivity",
                                 "specificity", "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords10$AUC <- AUC10
coords10
#separate markers
roc_results_tumor<- lapply(expression, function(col) {
  roc(response = OC_HGSOC_BENIGN$tumor, predictor = OC_HGSOC_BENIGN[[col]])})
names(roc_results_tumor) <- expression
roc_results_tumor
#extract the aucs
auc_values_tumor <- sapply(roc_results_tumor, function(roc_obj) {auc(roc_obj)})
auc_values_tumor #extracted aucs
#get roc features
coords_results_tumor <- lapply(roc_results_tumor, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity", "precision", "npv", "tpr", "fpr"),
         transpose = FALSE)
})
coords_results_tumor

##FIG: best 3 HGSOC vs benign PLOT ##########################
#HGSOC vs benign models plot:
roc_plot_custom <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_curve2, print.auc = F, col = "deeppink", 
           cex.main=0.8, main ="Gerybinių pakitimų atskyrimas HGSOC atvejų",
           xlab = "Specifiškumas",   # Custom x-axis label 
           ylab = "Jautrumas") #7
  lines(roc_curve10, col = "blue", lwd =2, lty = 2 ) #6
  lines(roc_results_tumor[["GRB7"]], col = "#469990", lwd =2) 
  lines(roc_results_tumor[["TCEAL4"]], col = "#808000", lwd =2)
  lines(roc_curve_CA2X, col = "grey", lwd = 2)
  # Add legend
  legend("bottomright", legend = c( expression(italic("GRB7 + TCEAL4")), 
                                    expression(italic("10 genų raiškos kombinacija")),
                                    expression(italic("GRB7 raiška")),
                                    expression(italic("TCEAL4 raiška")),
                                    expression(italic("Serumo CA125 biožymens statusas")))
         ,
         col = c("deeppink","blue","#469990", "#808000", "grey" ), lty = 1, 
         cex = 1, lwd =3)
}
#show plot
roc_plot_custom()
# Save the plot as a PNG file
png("FIG_best3_HSGOC_BENIGN20251020.png", width = 1000, height = 1000, res = 150)
roc_plot_custom()
#mtext("B", side = 3, adj = 0, line = 2.5, cex = 1.5, font = 2)
dev.off()

#FIG: best 3 HGSOC vs benign PLOT TABLE #####################
results_roc_custom <- data.frame(
  Biožymenys = c("GRB7 + TCEAL4", 
                 "10 genų raiškos kombinacija",
                 "GRB7 raiška",
                 "TCEAL4 raiška",
                 "Serumo CA125 biožymens statusas"),
  `plotas po kreive` = c(roc_curve2$auc, roc_curve10$auc, roc_results_tumor[["GRB7"]]$auc,
                         roc_results_tumor[["TCEAL4"]]$auc, roc_curve_CA2X$auc), 
  `slenkstinė vertė` = c(coords2$threshold, coords10$threshold ,
                         coords_results_tumor[["GRB7"]]$threshold,
                         coords_results_tumor[["TCEAL4"]]$threshold,
                         coords_ca2X$threshold),
  tikslumas = c(coords2$accuracy, coords10$accuracy ,
                coords_results_tumor[["GRB7"]]$accuracy,
                coords_results_tumor[["TCEAL4"]]$accuracy, coords_ca2X$accuracy),
  jautrumas = c(coords2$sensitivity, coords10$sensitivity,
                coords_results_tumor[["GRB7"]]$sensitivity,
                coords_results_tumor[["TCEAL4"]]$sensitivity, coords_ca2X$sensitivity),
  specifiškumas = c(coords2$specificity, coords10$specificity,
                    coords_results_tumor[["GRB7"]]$specificity,
                    coords_results_tumor[["TCEAL4"]]$specificity,
                    coords_ca2X$specificity),
  ppv  = c(coords2$precision, coords10$precision, 
           coords_results_tumor[["GRB7"]]$precision,
           coords_results_tumor[["TCEAL4"]]$precision, coords_ca2X$precision),
  npv  = c(coords2$npv, coords10$npv,
           coords_results_tumor[["GRB7"]]$npv,
           coords_results_tumor[["TCEAL4"]]$npv, coords_ca2X$npv),
  tpr  = c(coords2$tpr, coords10$tpr, 
           coords_results_tumor[["GRB7"]]$tpr,
           coords_results_tumor[["TCEAL4"]]$tpr, coords_ca2X$tpr),
  fpr  = c(coords2$fpr, coords10$fpr,
           coords_results_tumor[["GRB7"]]$fpr,
           coords_results_tumor[["TCEAL4"]]$fpr, coords_ca2X$fpr),
  check.names = FALSE
)
rownames(results_roc_custom) <- NULL
results_roc_custom

gt_table_cut <- results_roc_custom %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai",
    subtitle = "Gerybinių pakitimų atskyrimas nuo HGSOC atvejų",
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymenys))
  )
#show
gt_table_cut

#there is no other convieneat way to save gt outputs
gtsave(gt_table_cut,vwidth = 800,
       filename = "FIG_tabbest3_HGSOC_BENIGN20251020.png")

#Combine the images
roc_image2<- image_read("FIG_best3_HSGOC_BENIGN20251020.png")
table_image2 <- image_read("FIG_tabbest3_HGSOC_BENIGN20251020.png")


# Now append horizontaly
combined_image2 <- image_append(c(roc_image2, table_image2), stack = F)

# Save the combined image
image_write(combined_image2, 
            "FIG_COMBINED_best3_HGSOC_BENIGN20251020.png")

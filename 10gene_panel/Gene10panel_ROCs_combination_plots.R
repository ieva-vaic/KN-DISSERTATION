#KN-DISSERTATION project. 10 gene panel - genes selected via statistical analysis of TCGA and GTEx data
#KN 10 gene analysis data (cv 2 %, one threshold, 65 cases!, restored KN-59) 2025-03-05
#ROCs: combination plots and comparisons
#libraries:
library(tidyverse)
library(pROC)
library(glmnet)
library(gtsummary)
library(gt)
library(brglm2)
library(reshape2)
library(rstatix) 
library(ggprism)
library(gridExtra)
library(scales)
library(htmlwidgets)
library(webshot)
library(magick)
#set wd for plot saving
setwd("C:/Users/Ieva/rprojects/outputs_all/DISS/")
#data###################################################################
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
#HGSOC vs BENIGN MODELS###################
##combination of the best two HGSOC VS BENIGN#####
genes2 <- c( "GRB7", "TCEAL4")
expr_tumor2 <- OC_HGSOC_BENIGN[colnames(OC_HGSOC_BENIGN) %in% genes2]
brglm.model_2 <- glm(OC_HGSOC_BENIGN$tumor ~ ., data = expr_tumor2, family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_2 <- predict.glm(brglm.model_2, type='response') 
pred_data2 <- OC_HGSOC_BENIGN[(rownames(OC_HGSOC_BENIGN) #remove incomplete rows
                               %in% names(predicted_probs_2)), ]
roc_curve2 <- roc(pred_data2$tumor, predicted_probs_2)
AUC2 <- auc(roc_curve2)
coords2 <- coords(roc_curve2,"best",
        ret=c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
coords2$AUC <- AUC2
coords2

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
                  "best", ret=c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
coords10$AUC <- AUC10
coords10

##CA125 #using only HGSOC vs BENIGN data######################
KN_CA2X <- OC_HGSOC_BENIGN[!is.na(OC_HGSOC_BENIGN$CA125_f), ] #remove empty
table(KN_CA2X$tumor)
KN_CA2X$CA125_fN <- as.numeric(factor(KN_CA2X$CA125_f))- 1
roc_curve_CA2X <- roc(KN_CA2X$tumor, KN_CA2X$CA125_fN , direction = ">")
plot(roc_curve_CA2X) #auc = 0.765
auc(roc_curve_CA2X)
coords_ca2X <- coords(roc_curve_CA2X, "best", ret=c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv",
                                                    "tpr", "fpr"), transpose = FALSE)
coords_ca2X

##compare HGSOC vs BENIGN MODELS to ca125 ##################################
roc.test(roc_curve2, roc_curve_CA2X)
roc.test(roc_curve10, roc_curve_CA2X)

##compare HGSOC vs BENIGN MODELS to tecal4##########################
roc_curve_tceal4 <- roc(OC_HGSOC_BENIGN$tumor, OC_HGSOC_BENIGN$TCEAL4)
roc.test(roc_curve_tceal4, roc_curve2)#0.2758 

##FIG: best 3 HGSOC vs benign PLOT ##########################
#HGSOC vs benign models plot:
roc_plot_custom <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_curve2, print.auc = F, col = "deeppink", 
           cex.main=0.8, main ="Gerybinių pokyčių atskyrimas HGSOC atvejų",
           xlab = "Specifiškumas",   # Custom x-axis label 
           ylab = "Jautrumas") #7
  lines(roc_curve10, col = "blue", lwd =2, lty = 2 ) #6
  lines(roc_curve_CA2X, col = "grey", lwd = 2)
  # Add legend
  legend("bottomright", legend = c( expression(italic("GRB7 + TCEAL4")), 
                                    expression(italic("10 genų raiškos kombinacija")),
                                    expression(italic("Serumo CA125 biožymens statusas")))
         ,
         col = c("deeppink","blue", "grey" ), lty = 1, 
         cex = 1, lwd =3)
}
#show plot
roc_plot_custom()
# Save the plot as a PNG file
png("FIG_best3_HSGOC_BENIGN20250718.png", width = 1000, height = 1000, res = 150)
roc_plot_custom()
mtext("B", side = 3, adj = 0, line = 2.5, cex = 1.5, font = 2)
dev.off()

#FIG: best 3 HGSOC vs benign PLOT TABLE #####################
results_roc_custom <- data.frame(
  Biožymenys = c("GRB7 + TCEAL4", 
                 "10 genų raiškos kombinacija",
                 "Serumo CA125 biožymens statusas"),
  `plotas po kreive` = c(roc_curve2$auc, roc_curve10$auc, roc_curve_CA2X$auc), 
  `slenkstinė vertė` = c(coords2$threshold, coords10$threshold , coords_ca2X$threshold),
  tikslumas = c(coords2$accuracy, coords10$accuracy , coords_ca2X$accuracy),
  jautrumas = c(coords2$sensitivity, coords10$sensitivity, coords_ca2X$sensitivity),
  specifiškumas = c(coords2$specificity, coords10$specificity, coords_ca2X$specificity),
  ppv  = c(coords2$precision, coords10$precision, coords_ca2X$precision),
  npv  = c(coords2$npv, coords10$npv, coords_ca2X$npv),
  tpr  = c(coords2$tpr, coords10$tpr, coords_ca2X$tpr),
  fpr  = c(coords2$fpr, coords10$fpr, coords_ca2X$fpr),
  check.names = FALSE
)
rownames(results_roc_custom) <- NULL
results_roc_custom

gt_table_cut <- results_roc_custom %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai",
    subtitle = "Gerybinių pokyčių atskyrimas nuo HGSOC atvejų",
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
  filename = "FIG_tabbest3_HGSOC_BENIGN20250624.png")

#Combine the images
roc_image2<- image_read("FIG_best3_HSGOC_BENIGN20250718.png")
table_image2 <- image_read("FIG_tabbest3_HGSOC_BENIGN20250624.png")

# Find the max width to align both
roc_info <- image_info(roc_image2)
table_info <- image_info(table_image2)
max_width <- max(roc_info$width, table_info$width)

# Pad each image to the max width
roc_image2_padded <- image_extent(roc_image2, geometry = geometry_area(max_width, roc_info$height), gravity = "center", color = "white")
table_image2_padded <- image_extent(table_image2, geometry = geometry_area(max_width, table_info$height), gravity = "center", color = "white")

# Now append vertically
combined_image2 <- image_append(c(roc_image2_padded, table_image2_padded), stack = TRUE)

# Save the combined image
image_write(combined_image2, 
            "FIG_COMBINED_best3_HGSOC_BENIGN20250909.png")

#OVCa vs BENIGN MODELS################################################################
##manually make combination of the best two OVCa vs benign#####
genes2 <- c( "GRB7", "TCEAL4")
expr_tumor2x <- OC_full[colnames(OC_full) %in% genes2]
brglm.model_2x <- glm(OC_full$tumor ~ ., data = expr_tumor2x, 
                     family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_2x <- predict.glm(brglm.model_2x, type='response') 
pred_data2x <- OC_full[(rownames(OC_full) #remove incoplete rows
                               %in% names(predicted_probs_2x)), ]
roc_curve2x <- roc(pred_data2x$tumor, predicted_probs_2x)
AUC2x <- auc(roc_curve2x)
coords2x <- coords(roc_curve2x,
                  "best", ret=c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
coords2x$AUC <- AUC2x
coords2x
#all 10 combinations ##################
expr_tumor10x <- OC_full[colnames(OC_full) %in% expression]
brglm.model_10x <- glm(OC_full$tumor ~ ., data = expr_tumor10x,
                      family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_10x <- predict.glm(brglm.model_10x, type='response') 
pred_data10x <- OC_full[(rownames(OC_full) #remove incomplete rows
                                %in% names(predicted_probs_10x)), ]
roc_curve10x <- roc(pred_data10x$tumor, predicted_probs_10x)
AUC10x <- auc(roc_curve10x)
coords10x <- coords(roc_curve10x,
                   "best", ret=c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
coords10x$AUC <- AUC10x
coords10x

#CA125 #using only HGSIC BENIGN data
KN_CA2Xx <- OC_full[!is.na(OC_full$CA125_f), ] #remove empty
table(KN_CA2Xx$tumor)
KN_CA2Xx$CA125_fN <- as.numeric(factor(KN_CA2Xx$CA125_f))- 1
roc_curve_CA2Xx <- roc(KN_CA2Xx$tumor, KN_CA2Xx$CA125_fN , direction = ">")
plot(roc_curve_CA2Xx) #auc = 0.7751
auc(roc_curve_CA2Xx)
coords_ca2Xx <- coords(roc_curve_CA2Xx, "best", ret=c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv",
                                                    "tpr", "fpr"), transpose = FALSE)
coords_ca2Xx

##FIG best 3 OVCa vs benign PLOT ##########################
#OVCa vs benign models plot:
roc_plot_customx <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_curve2, print.auc = F, col = "lightpink", 
           cex.main=0.8, main ="Gerybinių pokyčių atskyrimas KV atvejų",
           xlab = "Specifiškumas",   # Custom x-axis label 
           ylab = "Jautrumas") #7
  lines(roc_curve10, col = "darkgreen", lwd =2, lty = 2 ) #6
  lines(roc_curve_CA2X, col = "grey", lwd = 2)
  
  # Add legend
  legend("bottomright", legend = c( expression(italic("GRB7 + TCEAL4")), 
                                    expression(italic("10 genų raiškos kombinacija")),
                                    expression(italic("Serumo CA125 biožymens statusas")))
         ,
         col = c("lightpink","darkgreen", "grey" ), lty = 1, 
         cex = 1, lwd =3)
}
# show plot
roc_plot_customx()

# Save the plot as a PNG file
png("FIG_best3ocfull_for_genes20250818.png", width = 1000, height = 1000, res = 150)
roc_plot_customx()
mtext("A", side = 3, adj = 0, line = 2.5, cex = 1.5, font = 2)
dev.off()
## FIG  best 3 OVCa vs benign TABLE #####################
results_roc_customx <- data.frame(
  Biožymenys = c("GRB7 + TCEAL4", 
                 "10 genų raiškos kombinacija",
                 "Serumo CA125 biožymens statusas"),
  `plotas po kreive` = c(roc_curve2$auc, roc_curve10$auc, roc_curve_CA2X$auc), 
  `slenkstinė vertė` = c(coords2$threshold, coords10$threshold , coords_ca2X$threshold),
  tikslumas = c(coords2$accuracy, coords10$accuracy , coords_ca2X$accuracy),
  jautrumas = c(coords2$sensitivity, coords10$sensitivity, coords_ca2X$sensitivity),
  specifiškumas = c(coords2$specificity, coords10$specificity, coords_ca2X$specificity),
  ppv  = c(coords2$precision, coords10$precision, coords_ca2X$precision),
  npv  = c(coords2$npv, coords10$npv, coords_ca2X$npv),
  tpr  = c(coords2$tpr, coords10$tpr, coords_ca2X$tpr),
  fpr  = c(coords2$fpr, coords10$fpr, coords_ca2X$fpr),
  check.names = FALSE
)
rownames(results_roc_customx) <- NULL
results_roc_customx

gt_table_cutx <- results_roc_customx %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai",
    subtitle = "Gerybinių pokyčių atskyrimas nuo KV atvejų",
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
gt_table_cutx

#there is no other convieneat way to save gt outputs
gtsave(gt_table_cutx,vwidth = 800, 
 filename = "FIG_tabbest3ocfull_for_genes20250624.png")

#Combine the images
roc_image2x<- image_read("FIG_best3ocfull_for_genes20250818.png")
table_image2x <- image_read("FIG_tabbest3ocfull_for_genes20250624.png")

# Find the max width to align both
roc_infox <- image_info(roc_image2x)
table_infox <- image_info(table_image2x)
max_widthx <- max(roc_infox$width, table_infox$width)

# Pad each image to the max width
roc_image2_paddedx <- image_extent(roc_image2x, geometry = geometry_area(max_width, roc_infox$height), gravity = "center", color = "white")
table_image2_paddedx <- image_extent(table_image2x, geometry = geometry_area(max_width, table_infox$height), gravity = "center", color = "white")

# Now append vertically
combined_image2x <- image_append(c(roc_image2_paddedx, table_image2_paddedx), stack = TRUE)

# Save the combined image
image_write(combined_image2x, 
            "FIG_COMBINED_forOC_full_genes20250909.png")

#HGSOC vs OTHERS MODELS###############################
##combination of the best 6 HGSOC vs others#####
genes6 <- c( "RAD50",	"PKP3",	"CDCA5"	,"ZFPL1",	"VPS33B",	"TCEAL4")
expr_tumor6 <- OC_HGSOC_OTHERS[colnames(OC_HGSOC_OTHERS) %in% genes6]
brglm.model_6 <- glm(OC_HGSOC_OTHERS$tumor ~ ., data = expr_tumor6, family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_6 <- predict.glm(brglm.model_6, type='response') 
pred_data6 <- OC_HGSOC_OTHERS[(rownames(OC_HGSOC_OTHERS) #remove incoplete rows
                               %in% names(predicted_probs_6)), ]
roc_curve6 <- roc(pred_data6$tumor, predicted_probs_6)
AUC6 <- auc(roc_curve6)
coords6 <- coords(roc_curve6,
                  "best", ret=c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv",
                                "tpr", "fpr"), transpose = FALSE)
coords6$AUC <- AUC6
coords6
## all 10 gene comb HGSOC vs others#####
expr_tumor10o <- OC_HGSOC_OTHERS[colnames(OC_HGSOC_OTHERS) %in% expression]
brglm.model_10o <- glm(OC_HGSOC_OTHERS$tumor ~ ., data = expr_tumor10o,
                      family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_10o<- predict.glm(brglm.model_10o, type='response') 
pred_data10o <- OC_HGSOC_OTHERS[(rownames(OC_HGSOC_OTHERS) #remove incomplete rows
                                %in% names(predicted_probs_10o)), ]
roc_curve10o <- roc(pred_data10o$tumor, predicted_probs_10o)
AUC10o <- auc(roc_curve10o)
coords10o <- coords(roc_curve10o,
                   "best", ret=c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv",
                                 "tpr", "fpr"), transpose = FALSE)
coords10o$AUC <- AUC10o
coords10o
##CA125 HGSOC vs others#####
KN_CA2Xo <- OC_HGSOC_OTHERS[!is.na(OC_HGSOC_OTHERS$CA125_f), ] #remove empty
table(KN_CA2Xo$tumor)
KN_CA2Xo$CA125_fN <- as.numeric(factor(KN_CA2Xo$CA125_f))- 1
roc_curve_CA2Xo <- roc(KN_CA2Xo$tumor, KN_CA2Xo$CA125_fN , direction = ">")
plot(roc_curve_CA2Xo) 
auc(roc_curve_CA2Xo)
coords_ca2Xo <- coords(roc_curve_CA2Xo, "best",
              ret=c("threshold", "accuracy", "sensitivity", "specificity",
                    "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords_ca2Xo
##best 2 combination HGSOC vs others#####
genes22 <- c( "EXO1", "TCEAL4")
expr_tumor22 <- OC_HGSOC_OTHERS[colnames(OC_HGSOC_OTHERS) %in% genes22]
brglm.model_22 <- glm(OC_HGSOC_OTHERS$tumor ~ ., data = expr_tumor22,
                      family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_22 <- predict.glm(brglm.model_22, type='response') 
pred_data22 <- OC_HGSOC_OTHERS[(rownames(OC_HGSOC_OTHERS) #remove incoplete rows
                               %in% names(predicted_probs_22)), ]
roc_curve22 <- roc(pred_data22$tumor, predicted_probs_22)
AUC22 <- auc(roc_curve22)
coords22 <- coords(roc_curve22,
                  "best", ret=c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv",
                                "tpr", "fpr"), transpose = FALSE)
coords22$AUC <- AUC22
coords22

##FIG: HGSOC vs OTHERS PLOT####################
roc_plot_customo <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_curve6, print.auc = F, col = "lightpink", 
           cex.main=0.8, main ="HGSOC atskyrimas kitų KV atvejų",
           xlab = "Specifiškumas",   # Custom x-axis label 
           ylab = "Jautrumas") #7
  lines(roc_curve22, col = "darkgreen", lwd =2)
  lines(roc_curve10o, col = "darkblue", lwd =2)
  lines(roc_curve_CA2Xo, col = "grey", lwd = 2)
  
  # Add legend
  legend("bottomright", legend = c( expression(italic("6 genų raiškos kombinacija")), 
                                    expression(italic("EXO1 + TCEAL4")),
                                    expression(italic("10 genų raiškos kombinacija")),
                                    expression(italic("Serumo CA125 biožymens statusas")))
         ,
         col = c("lightpink","darkgreen","darkblue" , "grey" ), lty = 1, 
         cex = 1, lwd =3)
}
# Save the plot to a variable
roc_plot_customo()

# Save the plot as a PNG file
png("FIG_bestHGSOC_VS_others20250624.png", width = 1000, height = 1000, res = 150)
roc_plot_customo()
dev.off()

##FIG: HGSOC vs OTHERS TABLE####################
results_roc_customo <- data.frame(
  Biožymenys = c("6 genų raiškos kombinacija",
                 "EXO1 + TCEAL4",
                 "10 genų raiškos kombinacija",
                 "Serumo CA125 biožymens statusas"),
  `plotas po kreive` = c(roc_curve6$auc, roc_curve22$auc, roc_curve10o$auc, roc_curve_CA2Xo$auc), 
  `slenkstinė vertė` = c(coords6$threshold, coords22$threshold ,  coords10o$threshold ,coords_ca2Xo$threshold),
  tikslumas = c(coords6$accuracy, coords22$accuracy ,coords10o$accuracy , coords_ca2Xo$accuracy),
  jautrumas = c(coords6$sensitivity, coords22$sensitivity, coords10o$sensitivity, coords_ca2Xo$sensitivity),
  specifiškumas = c(coords6$specificity, coords22$specificity,  coords10$specificity,coords_ca2X$specificity),
  ppv  = c(coords6$precision, coords22$precision,coords10o$precision, coords_ca2Xo$precision),
  npv  = c(coords6$npv, coords22$npv,coords10o$npv, coords_ca2Xo$npv),
  tpr  = c(coords6$tpr, coords22$tpr, coords10o$tpr, coords_ca2Xo$tpr),
  fpr  = c(coords6$fpr, coords22$fpr, coords10o$fpr, coords_ca2Xo$fpr),
  check.names = FALSE
)
rownames(results_roc_customo) <- NULL
results_roc_customo

gt_table_cuto <- results_roc_customo %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai",
    subtitle = "HGSOC atskyrimas kitų KV atvejų",
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
gt_table_cuto

#there is no other convieneat way to save gt outputs
gtsave(gt_table_cuto,vwidth = 800, 
    filename = "FIG_besttableHGSOC_VS_others20250624.png")

#Combine the images
roc_image2o<- image_read("FIG_bestHGSOC_VS_others20250624.png")
table_image2o <- image_read("FIG_besttableHGSOC_VS_others20250624.png")

# Find the max width to align both
roc_infoo <- image_info(roc_image2o)
table_infoo <- image_info(table_image2o)
max_widtho <- max(roc_infoo$width, table_infoo$width)

# Pad each image to the max width
roc_image2_paddedo <- image_extent(roc_image2o, geometry = geometry_area(max_width, roc_infoo$height), gravity = "center", color = "white")
table_image2_paddedo <- image_extent(table_image2o, geometry = geometry_area(max_width, table_infoo$height), gravity = "center", color = "white")

# Now append vertically
combined_image2o <- image_append(c(roc_image2_paddedo, table_image2_paddedo), stack = TRUE)

# Save the combined image
image_write(combined_image2o, 
            "FIG_COMBINED_forHGSOC_others_genes20250625.png")

#horizontal version
# Now append vertically
combined_image2o2 <- image_append(c(roc_image2o, table_image2o), stack = F)

# Save the combined image
image_write(combined_image2o2, 
            "FIG_COMBINED_forHGSOC_others_genes20250923.png")

##compare HGSOC vs OTHERs ###################################
roc.test(roc_curve10, roc_curve22, method = "delong")
roc.test(roc_curve6, roc_curve10, method = "delong")
roc.test(roc_curve6, roc_curve_CA2X)
roc_curve_tceal4 <- roc(OC_HGSOC_OTHERS$tumor, OC_HGSOC_OTHERS$TCEAL4)
roc.test(roc_curve_tceal4, roc_curve6)
roc.test(roc_curve_tceal4, roc_curve22)

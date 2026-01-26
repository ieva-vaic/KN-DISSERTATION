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
           cex.main=0.8, main ="Gerybinių pakitimų atskyrimas HGSOC atvejų",
           xlab = "1 - Specifiškumas",   # Custom x-axis label 
           ylab = "Jautrumas",
           legacy.axes = T) #7
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
png("FIG_best3_HSGOC_BENIGN20260121.png",width = 15, height = 15, res = 510, units = "cm")
roc_plot_custom()
mtext("B", side = 3, adj = 0, line = 2.5, cex = 1.5, font = 2)
dev.off()

##TABLE: best 3 HGSOC vs benign #####################
results_roc_custom <- data.frame(
  Biožymenys = c("GRB7 + TCEAL4", 
                 "10 genų raiškos kombinacija",
                 "Serumo CA125 biožymens statusas"),
  `plotas po kreive` = c(roc_curve2$auc, roc_curve10$auc, roc_curve_CA2X$auc), 
  `slenkstinė vertė` = c(coords2$threshold, coords10$threshold , coords_ca2X$threshold),
  tikslumas = c(coords2$accuracy, coords10$accuracy , coords_ca2X$accuracy),
  jautrumas = c(coords2$sensitivity, coords10$sensitivity, coords_ca2X$sensitivity),
  specifiškumas = c(coords2$specificity, coords10$specificity, coords_ca2X$specificity),
 # ppv  = c(coords2$precision, coords10$precision, coords_ca2X$precision),
 # npv  = c(coords2$npv, coords10$npv, coords_ca2X$npv),
#  tpr  = c(coords2$tpr, coords10$tpr, coords_ca2X$tpr),
 # fpr  = c(coords2$fpr, coords10$fpr, coords_ca2X$fpr),
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
gtsave(gt_table_cut,vwidth = 1000,
  filename = "FIG_tabbest3_HGSOC_BENIGN20260121.png")

#Combine the images
roc_image2<- image_read("FIG_best3_HSGOC_BENIGN20260121.png")
table_image2 <- image_read("FIG_tabbest3_HGSOC_BENIGN20260121.png")

# Match widths (use the larger width to avoid downscaling)
max_width <- max(image_info(roc_image2)$width,
                 image_info(table_image2)$width)

roc_image2   <- image_resize(roc_image2,   paste0(max_width, "x"))
table_image2 <- image_resize(table_image2, paste0(max_width, "x"))

# Append vertically
combined_image <- image_append(c(roc_image2, table_image2), stack = TRUE)

# Save with high quality
image_write(
  combined_image,
  path = "FIG_best3_HGSOC_BENIGN_combined20260121.png",
  format = "png"
)


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
##all 10 combinations ##################
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

##FIG: best 3 OVCa vs benign PLOT ##########################
#OVCa vs benign models plot:
roc_plot_customx <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_curve2, print.auc = F, col = "lightpink", 
           cex.main=0.8, main ="Gerybinių pakitimų atskyrimas KV atvejų",
           xlab = "1 - Specifiškumas",   # Custom x-axis label 
           ylab = "Jautrumas", 
           legacy.axes = T) #7
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
png("FIG_best3ocfull_for_genes20260121.png", width = 15, height = 15, res = 510, units = "cm")
roc_plot_customx()
mtext("A", side = 3, adj = 0, line = 2.5, cex = 1.5, font = 2)
dev.off()

##TABLE:  best 3 OVCa vs benign #####################
results_roc_customx <- data.frame(
  Biožymenys = c("GRB7 + TCEAL4", 
                 "10 genų raiškos kombinacija",
                 "Serumo CA125 biožymens statusas"),
  `plotas po kreive` = c(roc_curve2$auc, roc_curve10$auc, roc_curve_CA2X$auc), 
  `slenkstinė vertė` = c(coords2$threshold, coords10$threshold , coords_ca2X$threshold),
  tikslumas = c(coords2$accuracy, coords10$accuracy , coords_ca2X$accuracy),
  jautrumas = c(coords2$sensitivity, coords10$sensitivity, coords_ca2X$sensitivity),
  specifiškumas = c(coords2$specificity, coords10$specificity, coords_ca2X$specificity),
  # ppv  = c(coords2$precision, coords10$precision, coords_ca2X$precision),
  # npv  = c(coords2$npv, coords10$npv, coords_ca2X$npv),
  # tpr  = c(coords2$tpr, coords10$tpr, coords_ca2X$tpr),
  # fpr  = c(coords2$fpr, coords10$fpr, coords_ca2X$fpr),
  check.names = FALSE
)
rownames(results_roc_customx) <- NULL
results_roc_customx

gt_table_cutx <- results_roc_customx %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai",
    subtitle = "Gerybinių pakitimų atskyrimas nuo KV atvejų",
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
 filename = "FIG_tabbest3ocfull_for_genes20260121.png")


# Save the combined image
roc_image2<- image_read("FIG_best3ocfull_for_genes20260121.png")
table_image2 <- image_read("FIG_tabbest3ocfull_for_genes20260121.png")

# Match widths (use the larger width to avoid downscaling)
max_width <- max(image_info(roc_image2)$width,
                 image_info(table_image2)$width)

roc_image2   <- image_resize(roc_image2,   paste0(max_width, "x"))
table_image2 <- image_resize(table_image2, paste0(max_width, "x"))

# Append vertically
combined_image <- image_append(c(roc_image2, table_image2), stack = TRUE)

# Save with high quality
image_write(
  combined_image,
  path = "FIG_COMBINED_forOC_full_genes20260121.png",
  format = "png"
)


##COMBINE fig A and B#############################################
img2 <- image_read("FIG_best3_HGSOC_BENIGN_combined20260121.png")
img1 <- image_read("FIG_COMBINED_forOC_full_genes20260121.png")

# Match heights (use the larger height to preserve resolution)
max_height <- max(image_info(img1)$height,
                  image_info(img2)$height)

img1 <- image_resize(img1, paste0("x", max_height))
img2 <- image_resize(img2, paste0("x", max_height))

# Append horizontally
combined_horizontal <- image_append(c(img1, img2), stack = FALSE)

# Save at high quality
image_write(
  combined_horizontal,
  path = "FIG_MODELS_horizontal_combined_20260121.png",
  format = "png"
)
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
           cex.main=0.8, main ="HGSOC atskyrimas nuo kitų KV atvejų",
           xlab = "1 - Specifiškumas",   # Custom x-axis label 
           ylab = "Jautrumas", 
           legacy.axes = T) #7
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
png("FIG_bestHGSOC_VS_others20260121.png", width = 1000, height = 1000, res = 150)
roc_plot_customo()
dev.off()

##TABLE: HGSOC vs OTHERS####################
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
  # ppv  = c(coords6$precision, coords22$precision,coords10o$precision, coords_ca2Xo$precision),
  # npv  = c(coords6$npv, coords22$npv,coords10o$npv, coords_ca2Xo$npv),
  # tpr  = c(coords6$tpr, coords22$tpr, coords10o$tpr, coords_ca2Xo$tpr),
  # fpr  = c(coords6$fpr, coords22$fpr, coords10o$fpr, coords_ca2Xo$fpr),
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
    filename = "FIG_besttableHGSOC_VS_others20260121.png")

roc_image2<- image_read("FIG_bestHGSOC_VS_others20260121.png")
table_image2 <- image_read("FIG_besttableHGSOC_VS_others20260121.png")

# Match widths (use the larger width to avoid downscaling)
max_width <- max(image_info(roc_image2)$width,
                 image_info(table_image2)$width)

roc_image2   <- image_resize(roc_image2,   paste0(max_width, "x"))
table_image2 <- image_resize(table_image2, paste0(max_width, "x"))

# Append vertically
combined_image <- image_append(c(roc_image2, table_image2), stack = TRUE)

# Save with high quality
image_write(
  combined_image,
  path = "FIG_COMBINED_forHGSOC_OTHERS_20260121.png",
  format = "png"
)

##compare HGSOC vs OTHERs ###################################
roc.test(roc_curve10, roc_curve22, method = "delong")
roc.test(roc_curve6, roc_curve10, method = "delong")
roc.test(roc_curve6, roc_curve_CA2X)
roc_curve_tceal4 <- roc(OC_HGSOC_OTHERS$tumor, OC_HGSOC_OTHERS$TCEAL4)
roc.test(roc_curve_tceal4, roc_curve6)
roc.test(roc_curve_tceal4, roc_curve22)

#ENGLISH PLOTS #########################################
##EN FIG: best 3 HGSOC vs benign PLOT ##########################
#HGSOC vs benign models plot:
roc_plot_customEN <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_curve2, print.auc = F, col = "deeppink", 
           cex.main=0.8, main ="Separation of benign ovarian tumors form HGSOC",
           #xlab = "1 - Specifiškumas",   # Custom x-axis label 
           #ylab = "Jautrumas",
           legacy.axes = T) #7
  lines(roc_curve10, col = "blue", lwd =2, lty = 2 ) #6
  lines(roc_curve_CA2X, col = "grey", lwd = 2)
  # Add legend
  legend("bottomright", legend = c( expression(italic("GRB7 + TCEAL4")), 
                                    expression(italic("10 gene expression combination")),
                                    expression(italic("Serum CA125 status")))
         ,
         col = c("deeppink","blue", "grey" ), lty = 1, 
         cex = 1, lwd =3)
}
#show plot
roc_plot_customEN()
# Save the plot as a PNG file
png("FIG_best3_HSGOC_BENIGN20251218EN.png", width = 800, height = 800, res = 150)
roc_plot_customEN()
mtext("B", side = 3, adj = 0, line = 2.5, cex = 1.5, font = 2)
dev.off()

##EN TABLE: best 3 HGSOC vs benign #####################
results_roc_customEN <- data.frame(
  Predictor = c("GRB7 + TCEAL4", 
                "10 gene expression combination",
                "Serum CA125 status"),
  AUC = c(roc_curve2$auc, roc_curve10$auc, roc_curve_CA2X$auc), 
  threshold = c(coords2$threshold, coords10$threshold , coords_ca2X$threshold),
  accuracy = c(coords2$accuracy, coords10$accuracy , coords_ca2X$accuracy),
  sensitivity = c(coords2$sensitivity, coords10$sensitivity, coords_ca2X$sensitivity),
  specificity = c(coords2$specificity, coords10$specificity, coords_ca2X$specificity),
  # ppv  = c(coords2$precision, coords10$precision, coords_ca2X$precision),
  # npv  = c(coords2$npv, coords10$npv, coords_ca2X$npv),
  #  tpr  = c(coords2$tpr, coords10$tpr, coords_ca2X$tpr),
  # fpr  = c(coords2$fpr, coords10$fpr, coords_ca2X$fpr),
  check.names = FALSE
)
rownames(results_roc_customEN) <- NULL
results_roc_customEN

gt_table_cutEN <- results_roc_customEN %>%
  gt() %>%
  tab_header(
    title = "ROC metrics",
    subtitle = "Separation of benign ovarian tumors form HGSOC",
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Predictor))
  )
#show
gt_table_cutEN

#there is no other convieneat way to save gt outputs
gtsave(gt_table_cutEN,vwidth = 600,
       filename = "FIG_tabbest3_HGSOC_BENIGN20251218EN.png")

#Combine the images
roc_image2EN<- image_read("FIG_best3_HSGOC_BENIGN20251218EN.png")
table_image2EN <- image_read("FIG_tabbest3_HGSOC_BENIGN20251218EN.png")

# Find the max width to align both
roc_infoEN <- image_info(roc_image2EN)
table_infoEN <- image_info(table_image2EN)
max_widthEN <- max(roc_infoEN$width, table_infoEN$width)

# Pad each image to the max width
roc_image2_paddedEN <- image_extent(roc_image2EN, 
                                    geometry = geometry_area(max_widthEN, roc_infoEN$height),
                                    gravity = "center", color = "white")
table_image2_paddedEN <- image_extent(table_image2EN,
                                      geometry = geometry_area(max_widthEN, table_infoEN$height),
                                      gravity = "center", color = "white")

# Now append vertically
combined_image2EN <- image_append(c(roc_image2_paddedEN, table_image2_paddedEN), stack = TRUE)

# Save the combined image
image_write(combined_image2EN, 
            "FIG_COMBINED_best3_HGSOC_BENIGN20251218EN.png")

##EN FIG: best 3 OVCa vs benign PLOT ##########################
#OVCa vs benign models plot:
roc_plot_customxEN <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_curve2, print.auc = F, col = "lightpink", 
           cex.main=0.8, main ="Separation of benign ovarian tumors form ovarian cancer",
           # xlab = "1 - Specifiškumas", 
           # ylab = "Jautrumas", 
           legacy.axes = T) #7
  lines(roc_curve10, col = "darkgreen", lwd =2, lty = 2 ) #6
  lines(roc_curve_CA2X, col = "grey", lwd = 2)
  
  # Add legend
  legend("bottomright", legend = c( expression(italic("GRB7 + TCEAL4")), 
                                    expression(italic("10 gene expression combination")),
                                    expression(italic("Serum CA125 status")))
         ,
         col = c("lightpink","darkgreen", "grey" ), lty = 1, 
         cex = 1, lwd =3)
}
# show plot
roc_plot_customxEN()

# Save the plot as a PNG file
png("FIG_best3ocfull_for_genes20251218EN.png", width = 1000, height = 1000, res = 150)
roc_plot_customxEN()
mtext("A", side = 3, adj = 0, line = 2.5, cex = 1.5, font = 2)
dev.off()

##EN TABLE:  best 3 OVCa vs benign #####################
results_roc_customxEN <- data.frame(
  Predictor = c("GRB7 + TCEAL4", 
                "10 gene expression combination",
                "Serum CA125 status"),
  AUC = c(roc_curve2$auc, roc_curve10$auc, roc_curve_CA2X$auc), 
  threshold = c(coords2$threshold, coords10$threshold , coords_ca2X$threshold),
  accuracy = c(coords2$accuracy, coords10$accuracy , coords_ca2X$accuracy),
  sensitivity = c(coords2$sensitivity, coords10$sensitivity, coords_ca2X$sensitivity),
  specificity = c(coords2$specificity, coords10$specificity, coords_ca2X$specificity),
  # ppv  = c(coords2$precision, coords10$precision, coords_ca2X$precision),
  # npv  = c(coords2$npv, coords10$npv, coords_ca2X$npv),
  # tpr  = c(coords2$tpr, coords10$tpr, coords_ca2X$tpr),
  # fpr  = c(coords2$fpr, coords10$fpr, coords_ca2X$fpr),
  check.names = FALSE
)
rownames(results_roc_customxEN) <- NULL
results_roc_customxEN

gt_table_cutxEN <- results_roc_customxEN %>%
  gt() %>%
  tab_header(
    title = "ROC metrics",
    subtitle = "Separation of benign ovarian tumors form ovarian cancer",
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Predictor))
  )
#show
gt_table_cutxEN

#there is no other convieneat way to save gt outputs
gtsave(gt_table_cutxEN,vwidth = 800, 
       filename = "FIG_tabbest3ocfull_for_genes20251218EN.png")

#Combine the images
roc_image2xEN<- image_read("FIG_best3ocfull_for_genes20251218EN.png")
table_image2xEN <- image_read("FIG_tabbest3ocfull_for_genes20251218EN.png")

# Find the max width to align both
roc_infoxEN <- image_info(roc_image2xEN)
table_infoxEN <- image_info(table_image2xEN)
max_widthxEN <- max(roc_infoxEN$width, table_infoxEN$width)

# Pad each image to the max width
roc_image2_paddedxEN <- image_extent(roc_image2xEN,
                                     geometry = geometry_area(max_widthEN, roc_infoxEN$height),
                                     gravity = "center", color = "white")
table_image2_paddedxEN <- image_extent(table_image2xEN,
                                       geometry = geometry_area(max_widthEN, table_infoxEN$height),
                                       gravity = "center", color = "white")

# Now append vertically
combined_image2xEN <- image_append(c(roc_image2_paddedxEN, table_image2_paddedxEN), stack = TRUE)

# Save the combined image
image_write(combined_image2xEN, 
            "FIG_COMBINED_forOC_full_genes20251218EN.png")

##EN FIG: HGSOC vs OTHERS PLOT####################
roc_plot_customoEN <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_curve6, print.auc = F, col = "lightpink", 
           cex.main=0.8, main ="Separation of HGSOC from other ovarian cancer cases",
           # xlab = "1 - Specifiškumas",   # Custom x-axis label 
           # ylab = "Jautrumas", 
           legacy.axes = T) #7
  lines(roc_curve22, col = "darkgreen", lwd =2)
  lines(roc_curve10o, col = "darkblue", lwd =2)
  lines(roc_curve_CA2Xo, col = "grey", lwd = 2)
  
  # Add legend
  legend("bottomright", legend = c( expression(italic("6 gene expression combination")), 
                                    expression(italic("EXO1 + TCEAL4")),
                                    expression(italic("10 gene expression combination")),
                                    expression(italic("Serum CA125 status")))
         ,
         col = c("lightpink","darkgreen","darkblue" , "grey" ), lty = 1, 
         cex = 1, lwd =3)
}
# Save the plot to a variable
roc_plot_customoEN()

# Save the plot as a PNG file
png("FIG_bestHGSOC_VS_others20251218EN.png", width = 1000, height = 1000, res = 150)
roc_plot_customoEN()
dev.off()

##EN TABLE: HGSOC vs OTHERS####################
results_roc_customoEN <- data.frame(
  Predictor = c("6 gene expression combination",
                "EXO1 + TCEAL4",
                "10 gene expression combination",
                "Serum CA125 status"),
  AUC = c(roc_curve6$auc, roc_curve22$auc, roc_curve10o$auc, roc_curve_CA2Xo$auc), 
  threshold = c(coords6$threshold, coords22$threshold ,  coords10o$threshold ,coords_ca2Xo$threshold),
  accuracy = c(coords6$accuracy, coords22$accuracy ,coords10o$accuracy , coords_ca2Xo$accuracy),
  sensitivity = c(coords6$sensitivity, coords22$sensitivity, coords10o$sensitivity, coords_ca2Xo$sensitivity),
  specificity = c(coords6$specificity, coords22$specificity,  coords10$specificity,coords_ca2X$specificity),
  # ppv  = c(coords6$precision, coords22$precision,coords10o$precision, coords_ca2Xo$precision),
  # npv  = c(coords6$npv, coords22$npv,coords10o$npv, coords_ca2Xo$npv),
  # tpr  = c(coords6$tpr, coords22$tpr, coords10o$tpr, coords_ca2Xo$tpr),
  # fpr  = c(coords6$fpr, coords22$fpr, coords10o$fpr, coords_ca2Xo$fpr),
  check.names = FALSE
)
rownames(results_roc_customoEN) <- NULL
results_roc_customoEN

gt_table_cutoEN <- results_roc_customoEN %>%
  gt() %>%
  tab_header(
    title = "ROC metrics",
    subtitle = "Separation of HGSOC from other ovarian cancer cases",
  ) %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Predictor))
  )
#show
gt_table_cutoEN

#there is no other convieneat way to save gt outputs
gtsave(gt_table_cutoEN,vwidth = 1000, 
       filename = "FIG_besttableHGSOC_VS_others20251218EN.png")

#Combine the images
roc_image2oEN<- image_read("FIG_bestHGSOC_VS_others20251218EN.png")
table_image2oEN <- image_read("FIG_besttableHGSOC_VS_others20251218EN.png")

# Find the max width to align both
roc_infooEN <- image_info(roc_image2oEN)
table_infooEN <- image_info(table_image2oEN)
max_widthoEN <- max(roc_infooEN$width, table_infooEN$width)

# Pad each image to the max width
roc_image2_paddedoEN <- image_extent(roc_image2oEN,
                                     geometry = geometry_area(max_widthEN, roc_infooEN$height),
                                     gravity = "center", color = "white")
table_image2_paddedoEN <- image_extent(table_image2oEN,
                                       geometry = geometry_area(max_widthEN, table_infooEN$height),
                                       gravity = "center", color = "white")

# Now append vertically
combined_image2oEN <- image_append(c(roc_image2_paddedoEN, table_image2_paddedoEN), stack = TRUE)

# Save the combined image
image_write(combined_image2oEN, 
            "FIG_COMBINED_forHGSOC_others_genes2026122EN.png")

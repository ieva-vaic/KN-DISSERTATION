#KN-DISSERTATION project. 10 gene panel - genes selected via statistical analysis of TCGA and GTEx data
#KN 10 gene analysis data (cv 2 %, one threshold, 65 cases!, restored KN-59) 2025-02-17
#ROCs
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
#set wd for plots
setwd("C:/Users/Ieva/rprojects/outputs_all/DISS/")
#data###################################################################
OC_full <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/OC_10_genes_clean_2025_02_14.RDS")
#expression genes list
expression <- c("EXO1", "RAD50","PPT2", "LUC7L2","PKP3", "CDCA5","ZFPL1","VPS33B", "GRB7","TCEAL4")
#make HGSOC vs BENIGN df
OC_HGSOC_BENIGN<- OC_full[c(OC_full$type != "Other"),] #51 cases left
OC_HGSOC_BENIGN$tumor <- relevel(factor(OC_HGSOC_BENIGN$type), ref = "Benign")

#ROC HGSOC vs BENIGN########################################################
roc_results_tumor<- lapply(expression, function(col) {
  roc(response = OC_HGSOC_BENIGN$tumor, predictor = OC_HGSOC_BENIGN[[col]])})
names(roc_results_tumor) <- expression
roc_results_tumor
#extract the aucs
auc_values_tumor <- sapply(roc_results_tumor, function(roc_obj) {auc(roc_obj)})
auc_values_tumor #extracted aucs

##PLOT roc HGSOC vs BENIGN #############################################
roc_plot <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_tumor[["EXO1"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="Gerybinių patologijų atskyrimas nuo HGSOC",
           xlab = "Specifiškumas", 
           ylab = "Jautrumas") #title
  lines(roc_results_tumor[["RAD50"]], col = "#911eb4", lwd =2) 
  lines(roc_results_tumor[["PPT2"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_tumor[["LUC7L2"]], col = "#42d4f4", lwd =2) 
  lines(roc_results_tumor[["PKP3"]], col = "#fabed4", lwd =2) 
  lines(roc_results_tumor[["CDCA5"]], col = "#f032e6", lwd =2) 
  lines(roc_results_tumor[["ZFPL1"]], col = "#f58231", lwd =2) 
  lines(roc_results_tumor[["VPS33B"]], col = "#a9a9a9", lwd =2) 
  lines(roc_results_tumor[["GRB7"]], col = "#469990", lwd =2) 
  lines(roc_results_tumor[["TCEAL4"]], col = "#808000", lwd =2) 
  # Add legend
  #SUDĖTA NE PAGAL PAVEIKLĄ BET PAGAL AUC DYDI
  legend("bottomright", legend = c( expression(italic("EXO1")),
                                    expression(italic("RAD50")),
                                    expression(italic("PPT2")), 
                                    expression(italic("LUC7L2")), 
                                    expression(italic("PKP3")),
                                    expression(italic("CDCA5")),
                                    expression(italic("ZFPL1")), 
                                    expression(italic("VPS33B")),
                                    expression(italic("GRB7")),
                                    expression(italic("TCEAL4"))
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4", "#fabed4",
          "#f032e6", "#f58231", "#a9a9a9", "#469990", "#808000"), lty = 1, 
  cex = 0.7, lwd =3)
}

roc_plot()
## Save the plot as a PNG file###################################
png("10_genes_roc_hgsoc_benign_20250623.png",
    width = 1000, height = 1000, res = 200)
roc_plot()
dev.off()

##compare ROC HGSOC vs BENIGN ###################################
roc.test(roc_results_tumor[["EXO1"]], roc_results_tumor[["TCEAL4"]],  method=c("delong"))# 0.017 10 genes
##ROC table HGSOC vs BENIGN################################
#get roc features
coords_results_tumor <- lapply(roc_results_tumor, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity", "precision", "npv", "tpr", "fpr"),
         transpose = FALSE)
})
coords_results_tumor
#create df
results_tumor<- data.frame(
  Predictor = expression,
  AUC = auc_values_tumor,
  do.call(rbind,coords_results_tumor) 
)
#lithuanize it 
colnames(results_tumor) <- c("Biožymuo", "plotas po kreive", "slenkstinė vertė", 
                             "tikslumas", "jautrumas", "specifiškumas", 
                             "ppv", "npv", "tpr", "fpr")
#nice formating of the Table metrics for ROC OC
gt_table_tumor <- results_tumor %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai", 
    subtitle = "Gerybinių patologijų atskyrimas nuo HGSOC") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymuo))
  )
#show
gt_table_tumor

#there is no other convenient way to save gt outputs
gtsave(gt_table_tumor,
       filename = "10_genes_roc_table_hgsoc_benign_20250623.png")

#Combine the images
roc_image1<- image_read("10_genes_roc_hgsoc_benign_20250623.png")
table_image1 <- image_read("10_genes_roc_table_hgsoc_benign_20250623.png")

combined_image1 <- image_append(c(roc_image1, table_image1), stack = F)

# Find the max width to align both
roc_info <- image_info(roc_image1)
table_info <- image_info(table_image1)
max_width <- max(roc_info$width, table_info$width)

# Pad each image to the max width
roc_image1_padded <- image_extent(roc_image1, geometry = geometry_area(max_width, roc_info$height), gravity = "center", color = "white")
table_image1_padded <- image_extent(table_image1, geometry = geometry_area(max_width, table_info$height), gravity = "center", color = "white")
#combine
combined_image1 <- image_append(c(roc_image1_padded, table_image1_padded), stack = T)
# Save the combined image
image_write(combined_image1, 
            "10_genes_roc_combined_hgsoc_benign_202507.png")


##CA125 ROC HGSOC vs benign for comparison###################################
KN_CA2X <- OC_HGSOC_BENIGN[!is.na(OC_HGSOC_BENIGN$CA125_f), ] #remove empty
table(KN_CA2X$tumor)
KN_CA2X$CA125_fN <- as.numeric(factor(KN_CA2X$CA125_f))- 1
roc_curve_CA2X <- roc(KN_CA2X$tumor, KN_CA2X$CA125_fN , direction = ">")
plot(roc_curve_CA2X) #auc = 0.765
auc(roc_curve_CA2X)
coords_ca2X <- coords(roc_curve_CA2X, "best", ret=c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv",
                                                    "tpr", "fpr"), transpose = FALSE)
coords_ca2X
#roc tests
roc.test(roc_curve_CA2X, roc_results_tumor[["EXO1"]]) # 0.9377
roc.test(roc_curve_CA2X, roc_results_tumor[["RAD50"]])#0.6012
roc.test(roc_curve_CA2X, roc_results_tumor[["PPT2"]])#0.222
roc.test(roc_curve_CA2X, roc_results_tumor[["LUC7L2"]])#0.4331
roc.test(roc_curve_CA2X, roc_results_tumor[["CDCA5"]]) #0.3579
roc.test(roc_curve_CA2X, roc_results_tumor[["ZFPL1"]]) #0.6487
roc.test(roc_curve_CA2X, roc_results_tumor[["VPS33B"]]) #0.5084
roc.test(roc_curve_CA2X, roc_results_tumor[["GRB7"]]) #0.02912
roc.test(roc_curve_CA2X, roc_results_tumor[["TCEAL4"]]) # 0.05818

#ROC HGSOC vs OTHERS#########################################################
#make sure the levels are correct 
OC_HGSOC_OTHERS<- OC_full[c(OC_full$type != "Benign"),] #56 cases left
OC_HGSOC_OTHERS$tumor <- relevel(factor(OC_HGSOC_OTHERS$type), ref = "Other")
#roc HGSOC vs others
roc_results_tumor_others <- lapply(expression, function(col) {
  roc(response = OC_HGSOC_OTHERS$tumor, predictor = OC_HGSOC_OTHERS[[col]])})
names(roc_results_tumor_others) <- expression
roc_results_tumor_others
#extract the aucs
auc_values_tumor_others <- sapply(roc_results_tumor_others, function(roc_obj) {auc(roc_obj)})
auc_values_tumor_others #extracted aucs

##PLOT roc HGSOC vs others #############################################
roc_plot2 <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_tumor_others[["EXO1"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="HGSOC navikų atskyrimas nuo kitų KV atvejų",
           xlab = "Specifiškumas", 
           ylab = "Jautrumas") #title
  lines(roc_results_tumor_others[["RAD50"]], col = "#911eb4", lwd =2) 
  lines(roc_results_tumor_others[["PPT2"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_tumor_others[["LUC7L2"]], col = "#42d4f4", lwd =2) 
  lines(roc_results_tumor_others[["PKP3"]], col = "#fabed4", lwd =2) 
  lines(roc_results_tumor_others[["CDCA5"]], col = "#f032e6", lwd =2) 
  lines(roc_results_tumor_others[["ZFPL1"]], col = "#f58231", lwd =2) 
  lines(roc_results_tumor_others[["VPS33B"]], col = "#a9a9a9", lwd =2) 
  lines(roc_results_tumor_others[["GRB7"]], col = "#469990", lwd =2) 
  lines(roc_results_tumor_others[["TCEAL4"]], col = "#808000", lwd =2) 
  # Add legend
  #SUDĖTA NE PAGAL PAVEIKLĄ BET PAGAL AUC DYDI
  legend("bottomright", legend = c( expression(italic("EXO1")),
                                    expression(italic("RAD50")),
                                    expression(italic("PPT2")), 
                                    expression(italic("LUC7L2")), 
                                    expression(italic("PKP3")),
                                    expression(italic("CDCA5")),
                                    expression(italic("ZFPL1")), 
                                    expression(italic("VPS33B")),
                                    expression(italic("GRB7")),
                                    expression(italic("TCEAL4"))
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4", "#fabed4",
          "#f032e6", "#f58231", "#a9a9a9", "#469990", "#808000"), lty = 1, 
  cex = 0.7, lwd =3)
}

roc_plot2()
## Save the plot as a PNG file#######################
png("10_genes_roc_hgsoc_others_20250623.png",
    width = 1000, height = 1000, res = 200)
roc_plot2()
dev.off()

##compare to roc of 0,5###################################
set.seed(123)  # For reproducibility
random_pred <- runif(length(OC_HGSOC_OTHERS$tumor))  # Random scores between 0 and 1

roc_random <- roc(OC_HGSOC_OTHERS$tumor, random_pred)

# DeLong test
roc.test(roc_results_tumor_others[["EXO1"]], roc_random)#0.0588
roc.test(roc_results_tumor_others[["RAD50"]], roc_random)
roc.test(roc_results_tumor_others[["PPT2"]], roc_random)
roc.test(roc_results_tumor_others[["LUC7L2"]], roc_random)
roc.test(roc_results_tumor_others[["PKP3"]], roc_random)
roc.test(roc_results_tumor_others[["CDCA5"]], roc_random)
roc.test(roc_results_tumor_others[["ZFPL1"]], roc_random)
roc.test(roc_results_tumor_others[["VPS33B"]], roc_random)
roc.test(roc_results_tumor_others[["GRB7"]], roc_random)
roc.test(roc_results_tumor_others[["TCEAL4"]], roc_random)

##compare to CA125#########################################
#ca125 roc hgsoc vs other
KN_CA2X <- OC_HGSOC_OTHERS[!is.na(OC_HGSOC_OTHERS$CA125_f), ] #remove empty
table(KN_CA2X$tumor)
KN_CA2X$CA125_fN <- as.numeric(factor(KN_CA2X$CA125_f))- 1
roc_curve_CA2X <- roc(KN_CA2X$tumor, KN_CA2X$CA125_fN , direction = ">")
plot(roc_curve_CA2X) #auc = 0.765
auc(roc_curve_CA2X)
coords_ca2X <- coords(roc_curve_CA2X, "best", 
        ret=c("threshold", "accuracy", "sensitivity", "specificity",
              "precision", "npv","tpr", "fpr"), transpose = FALSE)
coords_ca2X
#roc tests
roc.test(roc_curve_CA2X, roc_results_tumor_others[["EXO1"]]) # 0.04754
roc.test(roc_curve_CA2X, roc_results_tumor_others[["RAD50"]])#0.5565
roc.test(roc_curve_CA2X, roc_results_tumor_others[["PPT2"]])#0.2999
roc.test(roc_curve_CA2X, roc_results_tumor_others[["PKP3"]])#0.624
roc.test(roc_curve_CA2X, roc_results_tumor_others[["LUC7L2"]])#0.6134
roc.test(roc_curve_CA2X, roc_results_tumor_others[["CDCA5"]]) #0.1022
roc.test(roc_curve_CA2X, roc_results_tumor_others[["ZFPL1"]]) #0.5239
roc.test(roc_curve_CA2X, roc_results_tumor_others[["VPS33B"]]) #0.4877
roc.test(roc_curve_CA2X, roc_results_tumor_others[["GRB7"]]) #0.3351
roc.test(roc_curve_CA2X, roc_results_tumor_others[["TCEAL4"]]) # 0.004414
##ROC table HGSOC vs others ################################
#get roc features
coords_results_tumor_others <- lapply(roc_results_tumor_others, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity", "precision", "npv", "tpr", "fpr"),
         transpose = FALSE)
})
coords_results_tumor_others
#create df
results_tumor_others<- data.frame(
  Predictor = expression,
  AUC = auc_values_tumor_others,
  do.call(rbind,coords_results_tumor_others) 
)
#make lt
colnames(results_tumor_others) <- c("Biožymuo", "plotas po kreive", "slenkstinė vertė", 
                             "tikslumas", "jautrumas", "specifiškumas", 
                             "ppv", "npv", "tpr", "fpr")
#nice formating of the Table metrics for ROC OC
gt_table_tumor_others <- results_tumor_others %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai",
    subtitle = "HGSOC navikų atskyrimas nuo kitų KV atvejų") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymuo))
  )
#show
gt_table_tumor_others

#there is no other convenient way to save gt outputs
gtsave(gt_table_tumor_others,
       filename = "10_genes_roc_table_hgsoc_other_20250623.png")

#Combine the images
roc_image2<- image_read("10_genes_roc_hgsoc_others_20250623.png")
table_image2 <- image_read("10_genes_roc_table_hgsoc_other_20250623.png")

combined_image2 <- image_append(c(roc_image2, table_image2), stack = F)

# Save the combined image
image_write(combined_image2, 
            "10_genes_roc_combined_hgsoc_others_20250623.png")
#ROC OVCa vs benign############################################
#ROC w OVCa vs benign
roc_results_tumor_OC<- lapply(expression, function(col) {
  roc(response = OC_full$tumor, predictor = OC_full[[col]])})
names(roc_results_tumor_OC) <- expression
roc_results_tumor_OC
#extract the aucs
auc_values_tumor_OC <- sapply(roc_results_tumor_OC, function(roc_obj) {auc(roc_obj)})
auc_values_tumor_OC #extracted aucs

##PLOT roc OVCa vs benign#############################################
roc_plot3 <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_tumor_OC[["EXO1"]], print.auc = F, col = "#dcbeff",
           cex.main=0.8, 
           main ="Gerybinių patologijų atskyrimas nuo KV atvejų", 
           xlab = "Specifiškumas", 
           ylab = "Jautrumas") #title
  lines(roc_results_tumor_OC[["RAD50"]], col = "#911eb4", lwd =2) 
  lines(roc_results_tumor_OC[["PPT2"]], col ="#ffd8b1", lwd =2) 
  lines(roc_results_tumor_OC[["LUC7L2"]], col = "#42d4f4", lwd =2) 
  lines(roc_results_tumor_OC[["PKP3"]], col = "#fabed4", lwd =2) 
  lines(roc_results_tumor_OC[["CDCA5"]], col = "#f032e6", lwd =2) 
  lines(roc_results_tumor_OC[["ZFPL1"]], col = "#f58231", lwd =2) 
  lines(roc_results_tumor_OC[["VPS33B"]], col = "#a9a9a9", lwd =2) 
  lines(roc_results_tumor_OC[["GRB7"]], col = "#469990", lwd =2) 
  lines(roc_results_tumor_OC[["TCEAL4"]], col = "#808000", lwd =2) 
  # Add legend
  #SUDĖTA NE PAGAL PAVEIKLĄ BET PAGAL AUC DYDI
  legend("bottomright", legend = c( expression(italic("EXO1")),
                                    expression(italic("RAD50")),
                                    expression(italic("PPT2")), 
                                    expression(italic("LUC7L2")), 
                                    expression(italic("PKP3")),
                                    expression(italic("CDCA5")),
                                    expression(italic("ZFPL1")), 
                                    expression(italic("VPS33B")),
                                    expression(italic("GRB7")),
                                    expression(italic("TCEAL4"))
  ),
  
  col = c("#dcbeff", "#911eb4", "#ffd8b1", "#42d4f4", "#fabed4",
          "#f032e6", "#f58231", "#a9a9a9", "#469990", "#808000"), lty = 1, 
  cex = 0.7, lwd =3)
}

roc_plot3()
## Save the plot as a PNG file#########################
png("10_genes_roc_oc_20250623.png",
    width = 1000, height = 1000, res = 200)
roc_plot3()
dev.off()

##compare to CA125#############################
#CA125 OVCa vs benign
KN_CA2X <- OC_full[!is.na(OC_full$CA125_f), ] #remove empty
table(KN_CA2X$tumor)
KN_CA2X$CA125_fN <- as.numeric(factor(KN_CA2X$CA125_f))- 1
roc_curve_CA2X <- roc(KN_CA2X$tumor, KN_CA2X$CA125_fN , direction = ">")
plot(roc_curve_CA2X) #auc = 0.765
auc(roc_curve_CA2X)
coords_ca2X <- coords(roc_curve_CA2X, "best",
    ret=c("threshold", "accuracy", "sensitivity", "specificity",
          "precision", "npv","tpr", "fpr"), transpose = FALSE)
coords_ca2X
#roc tests
roc.test(roc_curve_CA2X, roc_results_tumor_OC[["EXO1"]]) # 0.8994
roc.test(roc_curve_CA2X, roc_results_tumor_OC[["RAD50"]])#0.521
roc.test(roc_curve_CA2X, roc_results_tumor_OC[["PPT2"]])#0.3368
roc.test(roc_curve_CA2X, roc_results_tumor_OC[["LUC7L2"]])#0.5118
roc.test(roc_curve_CA2X, roc_results_tumor_OC[["CDCA5"]]) #0.547
roc.test(roc_curve_CA2X, roc_results_tumor_OC[["ZFPL1"]]) #0.4945
roc.test(roc_curve_CA2X, roc_results_tumor_OC[["VPS33B"]]) #0.6404
roc.test(roc_curve_CA2X, roc_results_tumor_OC[["GRB7"]]) #0.03597
roc.test(roc_curve_CA2X, roc_results_tumor_OC[["TCEAL4"]]) # 0.07026

#comapare together
roc.test(roc_results_tumor_OC[["EXO1"]], roc_results_tumor_OC[["TCEAL4"]],  method=c("delong"))# 0.07 10 genes
##ROC table OVCa vs benign ################################
#get roc features
coords_results_tumor_OC <- lapply(roc_results_tumor_OC, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity",
                                  "specificity", "precision", "npv", "tpr", "fpr"),
         transpose = FALSE)
})
coords_results_tumor_OC
#create df
results_tumor_OC<- data.frame(
  Predictor = expression,
  AUC = auc_values_tumor_OC,
  do.call(rbind,coords_results_tumor_OC) 
)
#make lt
colnames(results_tumor_OC) <- c("Biožymuo", "plotas po kreive", "slenkstinė vertė", 
                                    "tikslumas", "jautrumas", "specifiškumas", 
                                    "ppv", "npv", "tpr", "fpr")
#nice formating of the Table metrics for ROC OC
gt_table_tumor_OC <- results_tumor_OC %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai",
    subtitle = "Gerybinių patologijų atskyrimas nuo KV atvejų") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymuo))
  )
#show
gt_table_tumor_OC

#there is no other convenient way to save gt outputs
gtsave(gt_table_tumor_OC,
       filename = "10_genes_roc_table_oc_20250623.png")

#Combine the images
roc_image3 <- image_read("10_genes_roc_oc_20250623.png")
table_image3 <- image_read("10_genes_roc_table_oc_20250623.png")

combined_image3 <- image_append(c(roc_image3, table_image3), stack = F)

# Save the combined image
image_write(combined_image3, 
"10_genes_roc_combined_OC_20250623.png")

#combine HGSOC vs benign ####################################################
#2 best genes is GRB7 and TCEAL4
#logistic reg . with 2 genes:
genes2 <- c( "GRB7", "TCEAL4")
expr_tumor2 <- OC_HGSOC_BENIGN[colnames(OC_HGSOC_BENIGN) %in% genes2]
logistic.model_2<- glm(OC_HGSOC_BENIGN$tumor ~ ., data = expr_tumor2, family = "binomial") #will not converge
#to make it converge
brglm.model_2 <- glm(OC_HGSOC_BENIGN$tumor ~ ., data = expr_tumor2, family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_2 <- predict.glm(brglm.model_2, type='response') #43 cases
pred_data2 <- OC_HGSOC_BENIGN[(rownames(OC_HGSOC_BENIGN)
                         %in% names(predicted_probs_2)), ]
dim(pred_data2)#45 cases
#roc 
roc_curve2 <- roc(pred_data2$tumor, predicted_probs_2)
AUC2 <- auc(roc_curve2)
AUC2
roc_plot_2 <- function() {
  par(pty = "s") #sets square
  plot(roc_curve2, main = "GRB7, TCEAL4, HGSOC vs benign") #auc =  1
}
roc_plot_2() 
#table
coords2 <- coords(roc_curve2,
"best", ret=c("threshold", "accuracy", "sensitivity", "specificity",
              "precision", "npv","tpr", "fpr"), transpose = FALSE)
coords2
#df of table
coords2 <- as.data.frame(coords2)
# Format coordinates as text
coord_text <- paste0( 
  "AUC: ", paste(AUC2, collapse = ", "),"\n",
  "Sensitivity: ", paste(round(coords2$sensitivity, 2), collapse = ", "),"\n",
  "Specificity: ", paste(round(coords2$specificity, 2), collapse = ", ")
)
# Add text below the plot
roc_plot_2 <- function() {
  par(pty = "s") #sets square
  plot(roc_curve2, main = "GRB7, TCEAL4 HGSOC vs BENIGN") #auc =  1
  
  mtext(coord_text, side = 1, line = 4,adj = 0.9, cex = 0.7)
}
roc_plot_2()

# Save the plot as a PNG file
png("10_genes_roc_model2_hgsoc_benign_2025_02_18.png",
    width = 500, height = 500, res = 100)
roc_plot_2()
dev.off()

#combine HGSOC others ###############################
#3best genes is "EXO1", "CDCA5", "TCEAL4"
#get only the expresion df
genes3 <- c( "EXO1", "CDCA5", "TCEAL4")
expr_tumor3 <- OC_HGSOC_OTHERS[colnames(OC_HGSOC_OTHERS) %in% genes3]
logistic.model_3<- glm(OC_HGSOC_OTHERS$tumor ~ ., data = expr_tumor3, family = "binomial") #will not converge
#to make it converge
#brglm.model_3 <- glm(OC_HGSOC_OTHERS$tumor ~ ., data = expr_tumor3, family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_3 <- predict.glm(logistic.model_3, type='response') #43 cases
pred_data3 <- OC_HGSOC_OTHERS[(rownames(OC_HGSOC_OTHERS)
                               %in% names(predicted_probs_3)), ]
dim(pred_data3)#43 cases
#roc 
roc_curve3 <- roc(pred_data3$tumor, predicted_probs_3)
AUC3 <- auc(roc_curve3)
AUC3
roc_plot_3 <- function() {
  par(pty = "s") #sets square
  plot(roc_curve3, main = "EXO1, CDCA5, TCEAL4, HGSOC vs others") #auc =  1
}
#table
coords3 <- coords(roc_curve3,
                  "best", ret=c("threshold", "accuracy", "sensitivity", "specificity",
                                "precision", "npv","tpr", "fpr"), transpose = FALSE)
coords3
# Convert to data frame
coords3 <- as.data.frame(coords3)

# Format coordinates as text
coord_text <- paste0( 
  "AUC: ", paste(AUC3, collapse = ", "),"\n",
  "Sensitivity: ", paste(round(coords3$sensitivity, 2), collapse = ", "),"\n",
  "Specificity: ", paste(round(coords3$specificity, 2), collapse = ", ")
)

# Add text below the plot
roc_plot_3 <- function() {
  par(pty = "s") #sets square
  plot(roc_curve3, main = "EXO1, CDCA5, TCEAL4 together, HGSOC vs others") 
  
  mtext(coord_text, side = 1, line = 4,adj = 0.9, cex = 0.7)
}
roc_plot_3()
# Save the plot as a PNG file
png("10_genes_roc_model3_hgsoc_others_2025_02_18.png",
    width = 500, height = 500, res = 100)
roc_plot_3()
dev.off()

#combine all 10 HGSOC vs others################################################
#get only the expresion df, all 10 genes
expr_tumor10 <- OC_HGSOC_OTHERS[colnames(OC_HGSOC_OTHERS) %in% expression]
logistic.model_10<- glm(OC_HGSOC_OTHERS$tumor ~ ., data = expr_tumor10, family = "binomial") 
predicted_probs_10 <- predict.glm(logistic.model_10, type='response')
pred_data10 <- OC_HGSOC_OTHERS[(rownames(OC_HGSOC_OTHERS)
                               %in% names(predicted_probs_10)), ]
dim(pred_data10)#43 cases
#roc 
roc_curve10 <- roc(pred_data10$tumor, predicted_probs_10)
AUC10 <- auc(roc_curve10)
AUC10
roc_plot_10 <- function() {
  par(pty = "s") #sets square
  plot(roc_curve10, main = "All 10 genes together HGSOC vs others") 
}
#table
coords10 <- coords(roc_curve10,
"best", ret=c("threshold", "accuracy", "sensitivity", "specificity",
"precision", "npv","tpr", "fpr"), transpose = FALSE)
coords10
# Convert to data frame
coords10 <- as.data.frame(coords10)

# Format coordinates as text
coord_text <- paste0( 
  "AUC: ", paste(AUC10, collapse = ", "),"\n",
  "Sensitivity: ", paste(round(coords10$sensitivity, 2), collapse = ", "),"\n",
  "Specificity: ", paste(round(coords10$specificity, 2), collapse = ", ")
)

# Add text below the plot
roc_plot_10 <- function() {
  par(pty = "s") #sets square
  plot(roc_curve10, main = "All 10 genes together, HGSOC vs others") 
  
  mtext(coord_text, side = 1, line = 4,adj = 0.9, cex = 0.7)
}
roc_plot_10()
# Save the plot as a PNG file
png("10_genes_roc_model10_hgsoc_others_2025_02_18.png",
    width = 500, height = 500, res = 100)
roc_plot_10()
dev.off()

#combine all 10 OC vs benign###################################
#get only the expresion df, all 10 genes
expr_tumor10_oc <- OC_full[colnames(OC_full) %in% expression]
logistic.model_10_OC <- glm(OC_full$tumor ~ ., data = expr_tumor10_oc, family = "binomial") 
#maxed out need converging
#to make it converge
brglm.model_10oc <- glm(OC_full$tumor ~ ., data = expr_tumor10_oc,
                        family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_2_oc <- predict.glm(brglm.model_10oc, type='response') 
pred_data2oc <- OC_full[(rownames(OC_full)
                               %in% names(predicted_probs_2_oc)), ]
dim(pred_data2oc)#41
table(pred_data2oc$tumor)

#roc 
roc_curve10oc <- roc(pred_data2oc$tumor, predicted_probs_2_oc)
AUC10oc <- auc(roc_curve10oc)
AUC10oc
roc_plot_10oc <- function() {
  par(pty = "s") #sets square
  plot(roc_curve10oc, main = "All 10 genes together OC vs benign") 
}
#table
coords10oc <- coords(roc_curve10oc,
                   "best", ret=c("threshold", "accuracy", "sensitivity", "specificity",
                                 "precision", "npv","tpr", "fpr"), transpose = FALSE)
coords10oc
# Convert to data frame
coords10oc <- as.data.frame(coords10oc)

# Format coordinates as text
coord_text_oc <- paste0( 
  "AUC: ", paste(AUC10oc, collapse = ", "),"\n",
  "Sensitivity: ", paste(round(coords10oc$sensitivity, 2), collapse = ", "),"\n",
  "Specificity: ", paste(round(coords10oc$specificity, 2), collapse = ", ")
)

# Add text below the plot
roc_plot_10oc <- function() {
  par(pty = "s") #sets square
  plot(roc_curve10oc, main = "All 10 genes together, OC vs benign") 
  
  mtext(coord_text_oc, side = 1, line = 4,adj = 0.9, cex = 0.7)
}
roc_plot_10oc()

# Save the plot as a PNG file
png("10_genes_roc_model10_oc_2025_0305.png",
    width = 500, height = 500, res = 100)
roc_plot_10oc()
dev.off()
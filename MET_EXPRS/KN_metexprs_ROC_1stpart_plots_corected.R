#KN-DISSERTATION project. MET_EXPRS - methylation, expression, mutation data
##redo of the ROC plots for METEXPRS, no extra metrics
#seprate biomarkers
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

#ROC OVCa ###########################
#OVCa - ovarian cancer, all cases vs benign cases
#tumor must be a factor with reference: in this command first  level is a reference
KN_data$tumor <- factor(KN_data$tumor, levels = c("Benign", "OC"))
#make sure the levels are correct 
KN_data$tumor <- relevel(KN_data$tumor, ref = "Benign")
#roc
roc_results_tumor<- lapply(biomarkers, function(col) {
  roc(response = KN_data$tumor, predictor = KN_data[[col]])})
names(roc_results_tumor) <- biomarkers
roc_results_tumor
#extract the aucs
auc_values_tumor <- sapply(roc_results_tumor, function(roc_obj) {auc(roc_obj)})
auc_values_tumor #extracted aucs

#In general I want for expression lower values to be in OC (direction benign > OC)
#and for methylation lower values to be in benign (benign < OC)
#do expression biomarkers separately
roc_results_tumor_expr<- lapply(raiska, function(col) {
  roc(response = KN_data$tumor, predictor = KN_data[[col]], 
      direction = ">")
}) #direction is needed to show that controls have higher values than cases
names(roc_results_tumor_expr) <- raiska
auc_values_tumor_expr <- sapply(roc_results_tumor_expr, function(roc_obj) {auc(roc_obj)})
auc_values_tumor_expr #extracted aucs for expression
#this has slightly fixed dll1, not gonna use these values

#do methylation biomarkers separately
roc_results_tumor_met <- lapply(metilinimas, function(col) {
  roc(response = KN_data$tumor, predictor = KN_data[[col]], 
      direction = "<")
}) #direction is needed to show that controls have higher values than cases
names(roc_results_tumor_met) <- metilinimas
auc_values_tumor_met <- sapply(roc_results_tumor_met, function(roc_obj) {auc(roc_obj)})

##roc figure OVCa######################################
roc_plot <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_tumor[["NOTCH1"]], print.auc = F, col = "#dcbeff",
           cex.main=1,
           main ="Gerybinių pakitimų atskyrimas nuo visų KV atvejų", 
           xlab = "1 - Specifiškumas",   # Custom x-axis label (e.g., in Lithuanian)
           ylab = "Jautrumas", 
           legacy.axes = T) #7
  lines(roc_results_tumor[["NOTCH2"]], col = "#911eb4", lwd =2) #6
  lines(roc_results_tumor[["NOTCH3"]], col ="#ffd8b1", lwd =2) #8
  lines(roc_results_tumor[["NOTCH4"]], col = "#42d4f4", lwd =2) #4-5
  lines(roc_results_tumor[["ARID1A"]], col = "#fabed4", lwd =2) #9
  lines(roc_results_tumor[["CTNNB1"]], col = "#f032e6", lwd =3.5) #1
  lines(roc_results_tumor[["FBXW7"]], col = "#f58231", lwd =2) #2-3
  lines(roc_results_tumor[["JAG2"]], col = "#a9a9a9", lwd =2) #2-3
  lines(roc_results_tumor[["DLL1"]], col = "#469990", lwd =2) #4-5
  lines(roc_results_tumor[["HES1"]], col = "#808000", lwd =2) #2-3
  lines(roc_results_tumor[["HOPX"]], col = "black", lwd =2) #2-3
  lines(roc_results_tumor[["ALX4"]], col = "red", lwd =2) #2-3
  lines(roc_results_tumor[["CDX2"]], col = "darkred", lwd =2) #2-3
  lines(roc_results_tumor[["ARID1A_met"]], col = "deeppink", lwd =2) #2-3
  # Add legend
  legend("bottomright", legend = c( expression(italic("CTNNB1")), #f032e6
                                    expression(italic("FBXW7")), #f58231
                                    expression(italic("HES1")), #808000
                                    expression(italic("DLL1")), #469990
                                    expression(italic("NOTCH4")),#42d4f4
                                    expression(italic("NOTCH2")),#911eb4
                                    expression(italic("NOTCH1")), #dcbeff
                                    expression(italic("NOTCH3")),#ffd8b1
                                    expression(italic("ARID1A")),#fabed4
                                    expression(italic("JAG2")),
                                    #a9a9a9
                                    expression(italic("HOPX metilinimas")),
                                    expression(italic("ALX4 metilinimas")),
                                    expression(italic("CDX2 metilinimas")),
                                    expression(italic("ARID1A metilinimas"))
                                    
  ),
  col = c("#f032e6", "#f58231","#808000", "#469990", "#42d4f4",
          "#911eb4", "#dcbeff", "#ffd8b1", "#fabed4", "#a9a9a9",
          "black", "red", "darkred", "deeppink"), lty = 1, 
  cex = 0.6, lwd =3)
}
roc_plot()
#save
png("met_exprs_roc_OC_output20260130.png",
    width = 13, height = 13, res = 400, units = "cm")
roc_plot()
dev.off()

##Table for ROC metrics, OVCa##########
coords_results_tumor <- lapply(roc_results_tumor, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity"), transpose = FALSE)
})
coords_results_tumor
coords_results_tumor$ARID1A_met #arid1a methylation have infinate thresholds 
# Create a dataframe combining AUC values and coordinates results
results_tumor<- data.frame(
  Predictor = c(raiska, "HOPX", "ALX4","CDX2" ),
  AUC = auc_values_tumor[1:13],
  do.call(rbind, coords_results_tumor[1:13]) 
)
rownames(results_tumor) <- NULL
results_tumor

#fix ARID1A
arid1a_metrics_tumor <- coords(roc_results_tumor[["ARID1A_met"]],
                               best.method="closest.toplef", "best",
                               ret=c("threshold", "accuracy", "sensitivity",
                                     "specificity"), 
                               transpose = FALSE) #jei INF
arid1a_metrics_tumor
arid1a_metrics_tumor$Predictor <- "ARID1A methylation"
arid1a_metrics_tumor$AUC <- auc_values_tumor[14]

tumor_lentele_atskiru <- merge(results_tumor,arid1a_metrics_tumor,
                               by = intersect(names(results_tumor),
                                              names(arid1a_metrics_tumor)), all = TRUE )
tumor_lentele_atskiruEN <- tumor_lentele_atskiru #for en version
tumor_lentele_atskiru <- tumor_lentele_atskiru %>%
  mutate(Predictor = factor(Predictor,
                            levels = c(raiska, "HOPX", "ALX4","CDX2", "ARID1A methylation"))) %>%
  arrange(Predictor)

#rename to lithuanian
colnames(tumor_lentele_atskiru) <- c("Biožymuo",  "plotas po kreive", "slenkstinė vertė",
                                     "tikslumas", "jautrumas", "specifiškumas")

tumor_lentele_atskiru$Biožymuo <- c(raiska, "HOPX metilinimas", "ALX4 metilinimas", 
                                    "CDX2 metilinimas", "ARID1A metilinimas")

#order according to AUC value
tumor_lentele_atskiru<- tumor_lentele_atskiru[order(-tumor_lentele_atskiru$`plotas po kreive`), ]
#nice formating of the Table metrics for ROC OC
gt_table_tumor <- tumor_lentele_atskiru %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai", 
    subtitle = "Gerybinių pakitimų atskyrimas nuo KV atvejų") %>%
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
gtsave(gt_table_tumor, filename = "met_exprs_roctable_OC_output20260130.png")

#save both the same size together
roc_image2   <- image_read("met_exprs_roc_OC_output20260130.png")
table_image2 <- image_read("met_exprs_roctable_OC_output20260130.png")

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
  "met_exprs_ROC_TABLE_OC_outputs20260130.png"
)

##roc figure OVCa######################################
roc_plotEN <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_tumor[["NOTCH1"]], print.auc = F, col = "#dcbeff",
           cex.main=1,
           main ="Separation of benign ovarian tumors form ovarian cancer", 
           xlab = "1 - Specificity",   # Custom x-axis label (e.g., in Lithuanian)
           ylab = "Sensitivity", 
           legacy.axes = T) #7
  lines(roc_results_tumor[["NOTCH2"]], col = "#911eb4", lwd =2) #6
  lines(roc_results_tumor[["NOTCH3"]], col ="#ffd8b1", lwd =2) #8
  lines(roc_results_tumor[["NOTCH4"]], col = "#42d4f4", lwd =2) #4-5
  lines(roc_results_tumor[["ARID1A"]], col = "#fabed4", lwd =2) #9
  lines(roc_results_tumor[["CTNNB1"]], col = "#f032e6", lwd =3.5) #1
  lines(roc_results_tumor[["FBXW7"]], col = "#f58231", lwd =2) #2-3
  lines(roc_results_tumor[["JAG2"]], col = "#a9a9a9", lwd =2) #2-3
  lines(roc_results_tumor[["DLL1"]], col = "#469990", lwd =2) #4-5
  lines(roc_results_tumor[["HES1"]], col = "#808000", lwd =2) #2-3
  lines(roc_results_tumor[["HOPX"]], col = "black", lwd =2) #2-3
  lines(roc_results_tumor[["ALX4"]], col = "red", lwd =2) #2-3
  lines(roc_results_tumor[["CDX2"]], col = "darkred", lwd =2) #2-3
  lines(roc_results_tumor[["ARID1A_met"]], col = "deeppink", lwd =2) #2-3
  # Add legend
  legend("bottomright", legend = c( expression(italic("CTNNB1")), #f032e6
                                    expression(italic("FBXW7")), #f58231
                                    expression(italic("HES1")), #808000
                                    expression(italic("DLL1")), #469990
                                    expression(italic("NOTCH4")),#42d4f4
                                    expression(italic("NOTCH2")),#911eb4
                                    expression(italic("NOTCH1")), #dcbeff
                                    expression(italic("NOTCH3")),#ffd8b1
                                    expression(italic("ARID1A")),#fabed4
                                    expression(italic("JAG2")),
                                    #a9a9a9
                                    expression(italic("HOPX methylation")),
                                    expression(italic("ALX4 methylation")),
                                    expression(italic("CDX2 methylation")),
                                    expression(italic("ARID1A methylation"))
                                    
  ),
  col = c("#f032e6", "#f58231","#808000", "#469990", "#42d4f4",
          "#911eb4", "#dcbeff", "#ffd8b1", "#fabed4", "#a9a9a9",
          "black", "red", "darkred", "deeppink"), lty = 1, 
  cex = 0.7, lwd =3)
}
roc_plotEN()

png("met_exprs_roc_OC_outputEN20260121.png",
    width = 15, height = 15, res = 510, units = "cm")
roc_plotEN()
dev.off()


#rename  methylation
tumor_lentele_atskiruEN$Predictor <- c(raiska, "HOPX methylation", "ALX4 methylation", 
                                       "CDX2 methylation", "ARID1A methylation")

#order according to AUC value
tumor_lentele_atskiruEN<- tumor_lentele_atskiruEN[order(-tumor_lentele_atskiruEN$AUC), ]
#nice formating of the Table metrics for ROC OC
gt_table_tumorEN <- tumor_lentele_atskiruEN %>%
  gt() %>%
  tab_header(
    title = "ROC metrics", 
    subtitle = "Separation of benign ovarian tumors form ovarian cancer") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Predictor))
  )
#show
gt_table_tumorEN

#there is no other convenient way to save gt outputs
gtsave(gt_table_tumorEN,
       filename = "met_exprs_roctable_OC_outputEN20260121.png")

#Combine the images 
roc_image2   <- image_read("met_exprs_roc_OC_outputEN20260121.png")
table_image2 <- image_read("met_exprs_roctable_OC_outputEN20260121.png")

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
  "met_exprs_ROC_TABLE_OC_outputsEN20260130.png"
)


#ROC HGSOC vs BENIGN###########################################
# HGSOC vs benign df
KN_BENIGN_HGSOC <- KN_data[KN_data$Grupė_Ieva != "Other", ] 
KN_BENIGN_HGSOC$Grupė_Ieva <- droplevels(KN_BENIGN_HGSOC$Grupė_Ieva)
table(KN_BENIGN_HGSOC$Grupė_Ieva) #51 left 42 vs 9

#tumor must be a factor with reference: in this command first  level is a reference
KN_BENIGN_HGSOC$Grupė_Ieva <- factor(KN_BENIGN_HGSOC$Grupė_Ieva, levels = c("Benign", "HGSOC"))

roc_results_tumor_bh<- lapply(biomarkers, function(col) {
  roc(response = KN_BENIGN_HGSOC$Grupė_Ieva, predictor = KN_BENIGN_HGSOC[[col]])
})
names(roc_results_tumor_bh) <- biomarkers
roc_results_tumor_bh
#extract the aucs
auc_values_tumor_bh <- sapply(roc_results_tumor_bh, function(roc_obj) {auc(roc_obj)})
auc_values_tumor_bh #extracted aucs

##roc figure HGSOC vs benign ##########################
roc_plot2 <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_tumor_bh[["NOTCH1"]], print.auc = F, col = "#dcbeff",
           cex.main=1, main ="Gerybinių pakitimų atskyrimas nuo HGSOC atvejų", 
           xlab = "1 - Specifiškumas",   # Custom x-axis label (e.g., in Lithuanian)
           ylab = "Jautrumas", 
           legacy.axes = T) #7
  lines(roc_results_tumor_bh[["NOTCH2"]], col = "#911eb4", lwd =2) #6
  lines(roc_results_tumor_bh[["NOTCH3"]], col ="#ffd8b1", lwd =2) #8
  lines(roc_results_tumor_bh[["NOTCH4"]], col = "#42d4f4", lwd =2) #4-5
  lines(roc_results_tumor_bh[["ARID1A"]], col = "#fabed4", lwd =2) #9
  lines(roc_results_tumor_bh[["CTNNB1"]], col = "#f032e6", lwd =3.5) #1
  lines(roc_results_tumor_bh[["FBXW7"]], col = "#f58231", lwd =2) #2-3
  lines(roc_results_tumor_bh[["JAG2"]], col = "#a9a9a9", lwd =2) #2-3
  lines(roc_results_tumor_bh[["DLL1"]], col = "#469990", lwd =2) #4-5
  lines(roc_results_tumor_bh[["HES1"]], col = "#808000", lwd =2) #2-3
  lines(roc_results_tumor_bh[["HOPX"]], col = "black", lwd =2) #2-3
  lines(roc_results_tumor_bh[["ALX4"]], col = "red", lwd =2) #2-3
  lines(roc_results_tumor_bh[["CDX2"]], col = "darkred", lwd =2, lty = 2) #2-3
  lines(roc_results_tumor_bh[["ARID1A_met"]], col = "deeppink", lwd =2) #2-3
  # Add legend
  legend("bottomright", legend = c( expression(italic("CTNNB1")), #f032e6
                                    expression(italic("FBXW7")), #f58231
                                    expression(italic("HES1")), #808000
                                    expression(italic("DLL1")), #469990
                                    expression(italic("NOTCH4")),#42d4f4
                                    expression(italic("NOTCH2")),#911eb4
                                    expression(italic("NOTCH1")), #dcbeff
                                    expression(italic("NOTCH3")),#ffd8b1
                                    expression(italic("ARID1A")),#fabed4
                                    expression(italic("JAG2")),
                                    #a9a9a9
                                    expression(italic("HOPX metilinimas")),
                                    expression(italic("ALX4 metilinimas")),
                                    expression(italic("CDX2 metilinimas")),
                                    expression(italic("ARID1A metilinimas"))
                                    
  ),
  col = c("#f032e6", "#f58231","#808000", "#469990", "#42d4f4",
          "#911eb4", "#dcbeff", "#ffd8b1", "#fabed4", "#a9a9a9",
          "black", "red", "darkred", "deeppink"), lty = 1, 
  cex = 0.75, lwd =3)
}
roc_plot2()
# Save the plot as a PNG file
png("met_exrs_roc_HGSOC_output20260121.png",
    width = 15, height = 15, res = 510, units = "cm")
roc_plot2()
dev.off()

##Table for ROC metrics, HGSOC vs benign#################################
coords_results_tumor_bh <- lapply(roc_results_tumor_bh, function(roc_obj) {
  coords(roc_obj, "best", ret = c("threshold", "accuracy", "sensitivity", "specificity"), transpose = FALSE)
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

#fix ARID1A
arid1a_metrics_tumor_bh <- coords(roc_results_tumor_bh[["ARID1A_met"]],
                                  best.method="closest.toplef", "best",
                                  ret=c("threshold", "accuracy", "sensitivity",
                                        "specificity"), 
                                  transpose = FALSE) #jei INF
arid1a_metrics_tumor_bh
arid1a_metrics_tumor_bh$Predictor <- "ARID1A methylation"
arid1a_metrics_tumor_bh$AUC <- auc_values_tumor_bh[14]

tumor_lentele_atskiru_bh <- merge(results_tumor_bh,arid1a_metrics_tumor_bh,
                                  by = intersect(names(results_tumor_bh),
                                                 names(arid1a_metrics_tumor_bh)), all = TRUE )
tumor_lentele_atskiru_bhEN <- tumor_lentele_atskiru_bh #for english version
tumor_lentele_atskiru_bh <- tumor_lentele_atskiru_bh %>%
  mutate(Predictor = factor(Predictor, levels = c(raiska, "HOPX", "ALX4","CDX2", "ARID1A methylation"))) %>%
  arrange(Predictor)
#rename to lithuanian variables
colnames(tumor_lentele_atskiru_bh) <- c("Biožymuo", "plotas po kreive", "slenkstinė vertė",
                                        "tikslumas", "jautrumas", "specifiškumas" 
                                       )

tumor_lentele_atskiru_bh$Biožymuo <- c(raiska, "HOPX metilinimas", "ALX4 metilinimas", 
                                       "CDX2 metilinimas", "ARID1A metilinimas")
#order according to AUC value
tumor_lentele_atskiru_bh<- tumor_lentele_atskiru_bh[order(-tumor_lentele_atskiru_bh$`plotas po kreive`), ]

#gt table
gt_table_tumor_bh <- tumor_lentele_atskiru_bh %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai", 
    subtitle = "Gerybinių pakitimų atskyrimas nuo HGSOC atvejų") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biožymuo))
  )
#show
gt_table_tumor_bh

#there is no other convenient way to save gt outputs
gtsave(gt_table_tumor_bh, filename = "met_exprs_roctable_HGSOC_output20260121.png")

#Combine the images 
roc_image2   <- image_read("met_exrs_roc_HGSOC_output20260121.png")
table_image2 <- image_read("met_exprs_roctable_HGSOC_output20260121.png")

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
  "met_exrs_ROCTABLE_HGSOC_output20260130.png"
)

## HGSOC vs benign ROC EN version #########################
roc_plot2EN <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_tumor_bh[["NOTCH1"]], print.auc = F, col = "#dcbeff",
           cex.main=1, main ="Separation of benign ovarian tumors form HGSOC", 
           xlab = "1 - Specificity",   # Custom x-axis label (e.g., in Lithuanian)
           ylab = "Sensitivity", 
           legacy.axes = T) #7
  lines(roc_results_tumor_bh[["NOTCH2"]], col = "#911eb4", lwd =2) #6
  lines(roc_results_tumor_bh[["NOTCH3"]], col ="#ffd8b1", lwd =2) #8
  lines(roc_results_tumor_bh[["NOTCH4"]], col = "#42d4f4", lwd =2) #4-5
  lines(roc_results_tumor_bh[["ARID1A"]], col = "#fabed4", lwd =2) #9
  lines(roc_results_tumor_bh[["CTNNB1"]], col = "#f032e6", lwd =3.5) #1
  lines(roc_results_tumor_bh[["FBXW7"]], col = "#f58231", lwd =2) #2-3
  lines(roc_results_tumor_bh[["JAG2"]], col = "#a9a9a9", lwd =2) #2-3
  lines(roc_results_tumor_bh[["DLL1"]], col = "#469990", lwd =2) #4-5
  lines(roc_results_tumor_bh[["HES1"]], col = "#808000", lwd =2) #2-3
  lines(roc_results_tumor_bh[["HOPX"]], col = "black", lwd =2) #2-3
  lines(roc_results_tumor_bh[["ALX4"]], col = "red", lwd =2) #2-3
  lines(roc_results_tumor_bh[["CDX2"]], col = "darkred", lwd =2, lty = 2) #2-3
  lines(roc_results_tumor_bh[["ARID1A_met"]], col = "deeppink", lwd =2) #2-3
  # Add legend
  legend("bottomright", legend = c( expression(italic("CTNNB1")), #f032e6
                                    expression(italic("FBXW7")), #f58231
                                    expression(italic("HES1")), #808000
                                    expression(italic("DLL1")), #469990
                                    expression(italic("NOTCH4")),#42d4f4
                                    expression(italic("NOTCH2")),#911eb4
                                    expression(italic("NOTCH1")), #dcbeff
                                    expression(italic("NOTCH3")),#ffd8b1
                                    expression(italic("ARID1A")),#fabed4
                                    expression(italic("JAG2")),
                                    #a9a9a9
                                    expression(italic("HOPX methylation")),
                                    expression(italic("ALX4 methylation")),
                                    expression(italic("CDX2 methylation")),
                                    expression(italic("ARID1A methylation"))
                                    
  ),
  col = c("#f032e6", "#f58231","#808000", "#469990", "#42d4f4",
          "#911eb4", "#dcbeff", "#ffd8b1", "#fabed4", "#a9a9a9",
          "black", "red", "darkred", "deeppink"), lty = 1, 
  cex = 0.75, lwd =3)
}
roc_plot2EN()

png("met_exprs_roc_HGSOC_outputEN20260121.png",
    width = 15, height = 15, res = 510, units = "cm")
roc_plot2EN()
dev.off()

#tumor_lentele_atskiru_bhEN$Predictor <- c(raiska, "HOPX methylation", "ALX4 methylation", 
#                                          "CDX2 methylation", "ARID1A methylation")
#order according to AUC value
tumor_lentele_atskiru_bhEN<- tumor_lentele_atskiru_bhEN[order(-tumor_lentele_atskiru_bhEN$AUC), ]

#gt table
gt_table_tumor_bhEN <- tumor_lentele_atskiru_bhEN %>%
  gt() %>%
  tab_header(
    title = "ROC metrics", 
    subtitle = "Separation of benign ovarian tumors form HGSOC") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Predictor))
  )
#show
gt_table_tumor_bhEN

#there is no other convenient way to save gt outputs
gtsave(gt_table_tumor_bhEN, filename = "met_exprs_roctable_HGSOC_outputEN20260121.png")

#combine images
roc_image2   <- image_read("met_exprs_roc_HGSOC_outputEN20260121.png")
table_image2 <- image_read("met_exprs_roctable_HGSOC_outputEN20260121.png")

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
  "met_exrs_ROCTABLE_HGSOC_outputEN20260130.png"
)
#ROC HGSOC vs other###########################################
# HGSOC vs other df
KN_OTHER_HGSOC <- KN_data[KN_data$Grupė_Ieva != "Benign", ] 
KN_OTHER_HGSOC$Grupė_Ieva <- droplevels(KN_OTHER_HGSOC$Grupė_Ieva)
table(KN_OTHER_HGSOC$Grupė_Ieva) #51 left 42 vs 9

#tumor must be a factor with reference: in this command first  level is a reference
KN_OTHER_HGSOC$Grupė_Ieva <- factor(KN_OTHER_HGSOC$Grupė_Ieva, levels = c("Other", "HGSOC"))

roc_results_tumor_oh<- lapply(biomarkers, function(col) {
  roc(response = KN_OTHER_HGSOC$Grupė_Ieva, predictor = KN_OTHER_HGSOC[[col]])
})
names(roc_results_tumor_oh) <- biomarkers
roc_results_tumor_oh
#extract the aucs
auc_values_tumor_oh <- sapply(roc_results_tumor_oh, function(roc_obj) {auc(roc_obj)})
auc_values_tumor_oh #extracted aucs

##roc figure HGSOC vs others ###################################
roc_plot3 <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_tumor_oh[["NOTCH1"]], print.auc = F, col = "#dcbeff",
           cex.main=1, main ="HGSOC navikų atskyrimas nuo kitų KV atvejų", 
           xlab = "1 - Specifiškumas",   # Custom x-axis label (e.g., in Lithuanian)
           ylab = "Jautrumas", 
           legacy.axes = T) #7
  lines(roc_results_tumor_oh[["NOTCH2"]], col = "#911eb4", lwd =2) #6
  lines(roc_results_tumor_oh[["NOTCH3"]], col ="#ffd8b1", lwd =2) #8
  lines(roc_results_tumor_oh[["NOTCH4"]], col = "#42d4f4", lwd =2) #4-5
  lines(roc_results_tumor_oh[["ARID1A"]], col = "#fabed4", lwd =2) #9
  lines(roc_results_tumor_oh[["CTNNB1"]], col = "#f032e6", lwd =3.5) #1
  lines(roc_results_tumor_oh[["FBXW7"]], col = "#f58231", lwd =2) #2-3
  lines(roc_results_tumor_oh[["JAG2"]], col = "#a9a9a9", lwd =2) #2-3
  lines(roc_results_tumor_oh[["DLL1"]], col = "#469990", lwd =2) #4-5
  lines(roc_results_tumor_oh[["HES1"]], col = "#808000", lwd =2) #2-3
  lines(roc_results_tumor_oh[["HOPX"]], col = "black", lwd =2) #2-3
  lines(roc_results_tumor_oh[["ALX4"]], col = "red", lwd =2) #2-3
  lines(roc_results_tumor_oh[["CDX2"]], col = "darkred", lwd =2, lty = 2) #2-3
  lines(roc_results_tumor_oh[["ARID1A_met"]], col = "deeppink", lwd =2) #2-3
  # Add legend
  legend("bottomright", legend = c( expression(italic("CTNNB1")), #f032e6
                                    expression(italic("FBXW7")), #f58231
                                    expression(italic("HES1")), #808000
                                    expression(italic("DLL1")), #469990
                                    expression(italic("NOTCH4")),#42d4f4
                                    expression(italic("NOTCH2")),#911eb4
                                    expression(italic("NOTCH1")), #dcbeff
                                    expression(italic("NOTCH3")),#ffd8b1
                                    expression(italic("ARID1A")),#fabed4
                                    expression(italic("JAG2")),
                                    #a9a9a9
                                    expression(italic("HOPX metilinimas")),
                                    expression(italic("ALX4 metilinimas")),
                                    expression(italic("CDX2 metilinimas")),
                                    expression(italic("ARID1A metilinimas"))
                                    
  ),
  col = c("#f032e6", "#f58231","#808000", "#469990", "#42d4f4",
          "#911eb4", "#dcbeff", "#ffd8b1", "#fabed4", "#a9a9a9",
          "black", "red", "darkred", "deeppink"), lty = 1, 
  cex = 0.7, lwd =3)
}

roc_plot3()
# Save the plot as a PNG file
png("met_exprs_ROC_HGSOC_others_output20260121.png",
    width = 15, height = 15, res = 510, units = "cm")
roc_plot3()
dev.off()

##Table for ROC metrics, HGSOC vs BENIGN###################
coords_results_tumor_oh <- lapply(roc_results_tumor_oh, function(roc_obj) {
  coords(roc_obj, "best", best.method="closest.toplef",
         ret = c("threshold", "accuracy", "sensitivity","specificity"),
         transpose = FALSE)
})
coords_results_tumor_oh

#fix NOTCH3
notch3_metrics_tumor_oh <- coords(roc_results_tumor_oh[["NOTCH3"]],
                                  best.method=c("youden"), "best",
                                  ret=c("threshold", "accuracy", "sensitivity",
                                        "specificity"), 
                                  transpose = FALSE) #jei INF
notch3_metrics_tumor_oh
notch3_metrics_tumor_oh$AUC <- auc_values_tumor_bh[3]
notch3_metrics_tumor_oh$Predictor <- "NOTCH3"

# Create a dataframe combining AUC values and coordinates results
results_tumor_oh<- data.frame(
  Predictor = c("NOTCH1", "NOTCH2", "NOTCH4", "ARID1A" ,"CTNNB1",
                "FBXW7" , "JAG2" ,  "DLL1" ,  "HES1"  , "HOPX", "ALX4","CDX2", "ARID1A_met"),
  AUC = auc_values_tumor_oh[c(1:2, 4:14)],
  do.call(rbind, coords_results_tumor_oh[c(1:2, 4:14)]) 
)
rownames(results_tumor_oh) <- NULL


# Reorder df2 to match df1
notch3_metrics_tumor_oh <- notch3_metrics_tumor_oh[, names(results_tumor_oh)]
#MERGE
tumor_lentele_atskiru_oh <- rbind(results_tumor_oh,notch3_metrics_tumor_oh)

results_tumor_ohEN <- tumor_lentele_atskiru_oh # for en version
colnames(tumor_lentele_atskiru_oh) <- c("Biožymuo", "plotas po kreive", "slenkstinė vertė",
                                        "tikslumas", "jautrumas", "specifiškumas" 
                                       )
#reorder
tumor_lentele_atskiru_oh <- tumor_lentele_atskiru_oh[match(c(raiska, metilinimas), tumor_lentele_atskiru_oh$Biožymuo), ]
tumor_lentele_atskiru_oh$Biožymuo <- c(raiska, "HOPX metilinimas", "ALX4 metilinimas", 
                                       "CDX2 metilinimas", "ARID1A metilinimas")
#order according to AUC value
tumor_lentele_atskiru_oh<- tumor_lentele_atskiru_oh[order(-tumor_lentele_atskiru_oh$`plotas po kreive`), ]

#gt table
gt_table_tumor_oh <- tumor_lentele_atskiru_oh %>%
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
gt_table_tumor_oh

#there is no other convenient way to save gt outputs
gtsave(gt_table_tumor_oh, filename = "met_exprs_table_HGSOC_others_output20260121.png")

#Combine the images
roc_image2   <- image_read("met_exprs_ROC_HGSOC_others_output20260121.png")
table_image2 <- image_read("met_exprs_table_HGSOC_others_output20260121.png")

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
  "met_exprs_ROCTABLE_HGSOC_others_output202602130.png"
)

##roc figure HGSOC vs others ###################################
roc_plot3EN <- function() {
  par(pty = "s") #sets square
  plot.roc(roc_results_tumor_oh[["NOTCH1"]], print.auc = F, col = "#dcbeff",
           cex.main=1, main ="Separation of HGSOC form other ovarian cancer cases", 
           xlab = "1 - Specificity",   # Custom x-axis label (e.g., in Lithuanian)
           ylab = "Sensitivity", 
           legacy.axes = T) #7
  lines(roc_results_tumor_oh[["NOTCH2"]], col = "#911eb4", lwd =2) #6
  lines(roc_results_tumor_oh[["NOTCH3"]], col ="#ffd8b1", lwd =2) #8
  lines(roc_results_tumor_oh[["NOTCH4"]], col = "#42d4f4", lwd =2) #4-5
  lines(roc_results_tumor_oh[["ARID1A"]], col = "#fabed4", lwd =2) #9
  lines(roc_results_tumor_oh[["CTNNB1"]], col = "#f032e6", lwd =3.5) #1
  lines(roc_results_tumor_oh[["FBXW7"]], col = "#f58231", lwd =2) #2-3
  lines(roc_results_tumor_oh[["JAG2"]], col = "#a9a9a9", lwd =2) #2-3
  lines(roc_results_tumor_oh[["DLL1"]], col = "#469990", lwd =2) #4-5
  lines(roc_results_tumor_oh[["HES1"]], col = "#808000", lwd =2) #2-3
  lines(roc_results_tumor_oh[["HOPX"]], col = "black", lwd =2) #2-3
  lines(roc_results_tumor_oh[["ALX4"]], col = "red", lwd =2) #2-3
  lines(roc_results_tumor_oh[["CDX2"]], col = "darkred", lwd =2, lty = 2) #2-3
  lines(roc_results_tumor_oh[["ARID1A_met"]], col = "deeppink", lwd =2) #2-3
  # Add legend
  legend("bottomright", legend = c( expression(italic("CTNNB1")), #f032e6
                                    expression(italic("FBXW7")), #f58231
                                    expression(italic("HES1")), #808000
                                    expression(italic("DLL1")), #469990
                                    expression(italic("NOTCH4")),#42d4f4
                                    expression(italic("NOTCH2")),#911eb4
                                    expression(italic("NOTCH1")), #dcbeff
                                    expression(italic("NOTCH3")),#ffd8b1
                                    expression(italic("ARID1A")),#fabed4
                                    expression(italic("JAG2")),
                                    #a9a9a9
                                    expression(italic("HOPX methylation")),
                                    expression(italic("ALX4 methylation")),
                                    expression(italic("CDX2 methylation")),
                                    expression(italic("ARID1A methylation"))
                                    
  ),
  col = c("#f032e6", "#f58231","#808000", "#469990", "#42d4f4",
          "#911eb4", "#dcbeff", "#ffd8b1", "#fabed4", "#a9a9a9",
          "black", "red", "darkred", "deeppink"), lty = 1, 
  cex = 0.7, lwd =3)
}

roc_plot3EN()

png("met_exprs_roc_HGSOC_Others_outputEN20260121.png", 
    width = 15, height = 15, res = 510, units = "cm")
roc_plot3EN()
dev.off()

results_tumor_ohEN$Predictor <- c("NOTCH1", "NOTCH2",  "NOTCH4", "ARID1A", 
                                  "CTNNB1", "FBXW7",  "JAG2" ,  "DLL1",
                                  "HES1", "HOPX", "ALX4", 
                                  "CDX2", "ARID1A methylation", "NOTCH3")
#order according to AUC value
results_tumor_ohEN<- results_tumor_ohEN[order(-results_tumor_ohEN$AUC), ]

#gt table
gt_table_tumor_ohEN <- results_tumor_ohEN %>%
  gt() %>%
  tab_header(
    title = "ROC metrics", 
    subtitle = "Separation of HGSOC form other ovarian cancer cases") %>%
  fmt_number(
    columns = everything(),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Predictor))
  )
#show
gt_table_tumor_ohEN

#there is no other convenient way to save gt outputs
gtsave(gt_table_tumor_ohEN,
       filename = "met_exprs_table_HGSOC_others_outputEN20260121.png")

#Combine the images
roc_image2   <- image_read("met_exprs_roc_HGSOC_Others_outputEN20260121.png")
table_image2 <- image_read("met_exprs_table_HGSOC_others_outputEN20260121.png")

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
  "met_exprs_ROCTABLE_HGSOC_others_outputEN202602130.png"
)

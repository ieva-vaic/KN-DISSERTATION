#KN-DISSERTATION project. MET_EXPRS - methylation, expression, mutation data
#ROC analisies, combined biomarkers 
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
#ROC COMBINATIONS OVCa ###########################
#this will have the problem of the one missing expression
#remove the missing case
KN_expr_full <- KN_data %>% filter(patient_id_aud != "KN-100")
#MODELS OVCa#########################################
##only expression###################################
expr_tumor <- KN_expr_full[colnames(KN_expr_full) %in% raiska]
group1 <- factor(KN_expr_full$tumor, levels = c("Benign", "OC") )
logistic.model_1<- glm(group1 ~ ., data = expr_tumor, family = "binomial") #will not converge
#try another model to see if converges
brglm.model_1 <- glm(group1 ~ ., data = expr_tumor, family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_1 <- predict.glm(brglm.model_1, type='response') 
roc_curve1 <- roc(group1, predicted_probs_1)
auc(roc_curve1)
plot(roc_curve1) #auc =  1
coords1 <- coords(roc_curve1, "best",
                  ret=c("threshold", "accuracy", "sensitivity", "specificity",
                        "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords1

##only methylation##################################
#keeping the same dataset for consistency
expr_hb1.2 <- KN_expr_full[colnames(KN_expr_full) %in% metilinimas]

logistic.model_1.2<- glm(group1 ~ ., data = expr_hb1.2, family = "binomial") #will not converge

predicted_probs_1.2 <- predict.glm(logistic.model_1.2, type='response') 
roc_curve1.2 <- roc(group1, predicted_probs_1.2)
auc(roc_curve1.2)
plot(roc_curve1.2) #auc =  0,88
coords1.2 <- coords(roc_curve1.2, "best", best.method= "closest.topleft",
                    ret=c("threshold", "accuracy", "sensitivity", "specificity",
                          "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords1.2

##methylation + gene expression##################################
# keeping the same dataset for consistency
expr_hb1.3 <- KN_expr_full[colnames(KN_expr_full) %in% biomarkers]

logistic.model_1.3<- glm(group1 ~ ., data = expr_hb1.3, family = "binomial") #will not converge
#to make it converge
brglm.model_1.3 <- glm(group1 ~ ., data = expr_hb1.3, family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_1.3 <- predict.glm(brglm.model_1.3, type='response') 
roc_curve1.3<- roc(group1, predicted_probs_1.3)
auc(roc_curve1.3)
plot(roc_curve1.3) #auc =  0.998
coords1.3 <- coords(roc_curve1.3, "best", 
                    ret=c("threshold", "accuracy", "sensitivity", "specificity",
                          "precision", "npv","tpr", "fpr"), transpose = FALSE)
coords1.3

##only NOTCH genes##################################
NOTCH_genes <- c("NOTCH1", "NOTCH2", "NOTCH3" ,"NOTCH4", "JAG2","DLL1","HES1" )
# keeping the same dataset for consistency
expr_hb1.4 <- KN_expr_full[colnames(KN_expr_full) %in% NOTCH_genes]

logistic.model_1.4<- glm(group1 ~ ., data = expr_hb1.4, family = "binomial") #will not converge
#to make it converge
brglm.model_1.4 <- glm(group1 ~ ., data = expr_hb1.4,
                       family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_1.4 <- predict.glm(brglm.model_1.4, type='response') 
roc_curve1.4<- roc(group1, predicted_probs_1.4)
auc(roc_curve1.4)
plot(roc_curve1.4) #auc =  0.982
coords1.4 <- coords(roc_curve1.4, best.method= "closest.topleft", "best", 
                    ret=c("threshold", "accuracy", "sensitivity","specificity", 
                          "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords1.4

##only HOX methylation###############################################
HOX_genes <- c("HOPX","ALX4","CDX2"  )
# keeping the same dataset for consistency
expr_hb1.5 <- KN_expr_full[colnames(KN_expr_full) %in% HOX_genes]

logistic.model_1.5<- glm(group1 ~ ., data = expr_hb1.5, family = "binomial") #will not converge
predicted_probs_1.5 <- predict.glm(logistic.model_1.5, type='response') 
roc_curve1.5<- roc(group1, predicted_probs_1.5)
auc(roc_curve1.5)
plot(roc_curve1.5) #auc =  0.982
coords1.5<- coords(roc_curve1.5, "best", ret=c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv",
                                               "tpr", "fpr"), transpose = FALSE)
coords1.5

##CA125 OVCa for comparisons#################################
# keeping the same dataset for consistency
KN_CA <- KN_expr_full[!is.na(KN_expr_full$CA125_f), ] #remove empty
table(KN_CA$tumor)
KN_CA$CA125_fN <- as.numeric(factor(KN_CA$CA125_f))- 1
roc_curve_CA <- roc(KN_CA$tumor, KN_CA$CA125_fN, direction = ">")
plot(roc_curve_CA) #auc = 0.765
auc(roc_curve_CA)
coords_ca <- coords(roc_curve_CA, "best", ret=c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv",
                                                "tpr", "fpr"), transpose = FALSE)
coords_ca

##PLOT ROC COMBINATIONS OVCa  ######################################
roc_plot_3 <- function() {
par(pty = "s") #sets square
plot.roc(roc_curve1, print.auc = F, col = "#911eb4",
         cex.main=0.8, main ="Gerybinių pokyčių atskyrimas nuo KV atvejų",
         xlab = "Specifiškumas",   # Custom x-axis label (e.g., in Lithuanian)
         ylab = "Jautrumas") #7
lines(roc_curve1.2, col = "#dcbeff", lwd =2, ) #6
lines(roc_curve1.3, col ="#fabed4", lwd =2, lty = 2) #8
lines(roc_curve1.4, col ="darkred", lwd =2) 
lines(roc_curve1.5, col ="deeppink", lwd =2) 
lines(roc_curve_CA, col = "grey", lwd = 2)

# Add legend
legend("bottomright", legend = c( expression(italic("Genų raiškos žymenų kombinacija")), 
                                  expression(italic("Promotorių metilinimo statuso kombinacija")),
                                  expression(italic("Genų raiškos žymenų + promotorių metilinimo statuso kombinacija ")), 
                                  expression(italic("NOTCH genų raiškos žymenų kombinacija")),
                                  expression(italic("HOX promotorių metilinimo statuso kombinacija")),
                                  expression(italic("Serumo CA125 biožymens statusas")))
       ,
       col = c("#911eb4","#dcbeff", "#fabed4", "darkred",
               "deeppink", "grey" ), lty = 1, 
       cex = 0.8, lwd =3)
}
# Save the plot to a variable
saved_plot <- recordPlot(roc_plot_3())

# Save the plot as a PNG file
png("metexprs_roc_combinations_OVCa_output20250625.png", width = 1000, height = 1000, res = 150)
roc_plot_3()
dev.off()

##TABLE ROC COMBINATIONS OVCa ###########################  
results_roc1<- data.frame(
  Biožymenys = c("Genų raiškos žymenų kombinacija", 
                "Promotorių metilinimo statuso kombinacija",
                "Genų raiškos žymenų + promotorių metilinimo statuso kombinacija", 
                "NOTCH genų raiškos žymenų kombinacija",
                "HOX promotorių metilinimo statuso kombinacija",
                "Serumo CA125 biožymens statusas"),
  `plotas po kreive` = c(roc_curve1$auc, roc_curve1.2$auc, roc_curve1.3$auc, roc_curve1.4$auc, roc_curve1.5$auc, roc_curve_CA$auc), 
  `slenkstinė vertė` =c(coords1$threshold, coords1.2$threshold , coords1.3$threshold,
                        coords1.4$threshold, coords1.5$threshold, coords_ca$threshold ),
  tikslumas = c(coords1$accuracy, coords1.2$accuracy , coords1.3$accuracy,
               coords1.4$accuracy, coords1.5$accuracy, coords_ca$accuracy ),
  jautrumas = c(coords1$sensitivity, coords1.2$sensitivity, coords1.3$sensitivity,
                  coords1.4$sensitivity, coords1.5$sensitivity, coords_ca$sensitivity),
  specifiškumas = c(coords1$specificity, coords1.2$specificity, coords1.3$specificity,
                  coords1.4$specificity, coords1.5$specificity, coords_ca$specificity),
  ppv  = c(coords1$precision, coords1.2$precision, coords1.3$precision,
                 coords1.4$precision,coords1.5$precision,coords_ca$precision ),
  npv  = c(coords1$npv, coords1.2$npv, coords1.3$npv,
           coords1.4$npv, coords1.5$npv, coords_ca$npv),
  tpr  = c(coords1$tpr, coords1.2$tpr, coords1.3$tpr,
           coords1.4$tpr,coords1.5$tpr,coords_ca$tpr),
  fpr  = c(coords1$fpr, coords1.2$fpr, coords1.3$fpr,
           coords1.4$fpr, coords1.5$fpr,coords_ca$fpr),
  check.names = FALSE
)
rownames(results_roc1) <- NULL
results_roc1

gt_table <- results_roc1 %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai",
    subtitle = "Gerybinių pokyčių atskyrimas nuo KV atvejų"
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
gt_table

#there is no other convieneat way to save gt outputs
gtsave(gt_table, filename = "metexprs_table_combinations_OVCa_output20250625.png", vheight = 800)

#Combine the images
roc_image2<- image_read("metexprs_roc_combinations_OVCa_output20250625.png")
table_image2 <- image_read("metexprs_table_combinations_OVCa_output20250625.png")

combined_image2 <- image_append(c(roc_image2, table_image2), stack = F)

# Save the combined image
image_write(combined_image2, 
            "metexprs_roctable_combinations_OVCa_output20250625.png")
# center
# Find the max width to align both
roc_info <- image_info(roc_image2)
table_info <- image_info(table_image2)
max_width <- max(roc_info$width, table_info$width)

# Pad each image to the max width
roc_image2_padded <- image_extent(roc_image2, geometry = geometry_area(max_width, roc_info$height), gravity = "center", color = "white")
table_image2_padded <- image_extent(table_image2, geometry = geometry_area(max_width, table_info$height), gravity = "center", color = "white")

# Now append vertically
combined_image2 <- image_append(c(roc_image2_padded, table_image2_padded), stack = T)

# Save the combined image
image_write(combined_image2, 
            "metexprs_roc_combinations_OVCa_v2_output20250909.png")

#non staked version
combined_image2.2 <- image_append(c(roc_image2, table_image2), stack = F)

# Save the combined image
image_write(combined_image2.2, 
            "metexprs_roc_combinations_OVCa_v1_output20250909.png")
##delong tests with CA125 OVCa#################################
roc.test(roc_curve_CA, roc_curve1)#0.03112 10 genes
roc.test(roc_curve_CA, roc_curve1.2)#0.3226 4 methyl
roc.test(roc_curve_CA, roc_curve1.3)#0.03112 14 biomakers
roc.test(roc_curve_CA, roc_curve1.4)#0.04916 notch
roc.test(roc_curve_CA, roc_curve1.5)#0.4405 #methyl

#ROC COMBINATIONS HGSOC###################
# HGSOC vs benign df
KN_BENIGN_HGSOC <- KN_data[KN_data$Grupė_Ieva != "Other", ] 
KN_BENIGN_HGSOC$Grupė_Ieva <- droplevels(KN_BENIGN_HGSOC$Grupė_Ieva)
table(KN_BENIGN_HGSOC$Grupė_Ieva) #51 left 42 vs 9
#tumor must be a factor with reference: in this command first level is a reference
KN_BENIGN_HGSOC$Grupė_Ieva <- factor(KN_BENIGN_HGSOC$Grupė_Ieva, levels = c("Benign", "HGSOC"))
#MODELS:
##only expression########################################
expr_HB <- KN_BENIGN_HGSOC[colnames(KN_BENIGN_HGSOC) %in% raiska]
group_HB<- factor(KN_BENIGN_HGSOC$Grupė_Ieva, levels = c("Benign", "HGSOC") )
logistic.model_2<- glm(group_HB ~ ., data = expr_HB, family = "binomial") #will not converge
#to make it converge
brglm.model_2 <- glm(group_HB ~ ., data = expr_HB, family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_2 <- predict.glm(brglm.model_2, type='response') 
roc_curve2 <- roc(group_HB, predicted_probs_2)
auc(roc_curve2)
plot(roc_curve2) #auc =  1
coords2 <- coords(roc_curve2, "best", 
                  ret=c("threshold", "accuracy", "sensitivity", "specificity",
                        "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords2

##only methylation#####################################
met_HB <- KN_BENIGN_HGSOC[colnames(KN_BENIGN_HGSOC) %in% metilinimas]
logistic.model_2.2<- glm(group_HB ~ ., data = met_HB, family = "binomial") #will not converge
predicted_probs_2.2 <- predict.glm(logistic.model_2.2, type='response') 
roc_curve2.2 <- roc(group_HB, predicted_probs_2.2)
auc(roc_curve2.2)
plot(roc_curve2.2) #auc =  0.87
coords2.2 <- coords(roc_curve2.2, "best", 
                    ret=c("threshold", "accuracy", "sensitivity", "specificity",
                          "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords2.2

##gene expression + methylation##################################
expr_hb2.3 <- KN_BENIGN_HGSOC[colnames(KN_BENIGN_HGSOC) %in% biomarkers]

logistic.model_2.3<- glm(group_HB ~ ., data = expr_hb2.3, family = "binomial") #will not converge
#to make it converge
brglm.model_2.3 <- glm(group_HB ~ ., data = expr_hb2.3, family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_2.3 <- predict.glm(brglm.model_2.3, type='response') 
roc_curve2.3<- roc(group_HB, predicted_probs_2.3)
auc(roc_curve2.3)
plot(roc_curve2.3) #auc =  1
coords2.3<- coords(roc_curve2.3, "best", 
                   ret=c("threshold", "accuracy", "sensitivity", "specificity",
                         "precision", "npv","tpr", "fpr"), transpose = FALSE)
coords2.3

##only NOTCH genes#############################
expr_hb2.4 <- KN_BENIGN_HGSOC[colnames(KN_BENIGN_HGSOC) %in% NOTCH_genes]

logistic.model_2.4<- glm(group_HB ~ ., data = expr_hb2.4, family = "binomial") #will not converge
#to make it converge
brglm.model_2.4 <- glm(group_HB ~ ., data = expr_hb2.4, family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_2.4 <- predict.glm(brglm.model_2.4, type='response') 
roc_curve2.4<- roc(group_HB, predicted_probs_2.4)
auc(roc_curve2.4)
plot(roc_curve2.4) #auc =  1
coords2.4<- coords(roc_curve2.4, "best", 
                   ret=c("threshold", "accuracy", "sensitivity", "specificity",
                         "precision", "npv","tpr", "fpr"), transpose = FALSE)
coords2.4

##only HOX methylation #####################################
expr_hb2.5 <- KN_BENIGN_HGSOC[colnames(KN_BENIGN_HGSOC) %in% HOX_genes]

logistic.model_2.5 <- glm(group_HB ~ ., data = expr_hb2.5, family = "binomial") #will not converge
predicted_probs_2.5 <- predict.glm(logistic.model_2.5, type='response') 
roc_curve2.5 <- roc(group_HB, predicted_probs_2.5)
auc(roc_curve2.5)
plot(roc_curve2.5) #auc =  0,85
coords2.5<- coords(roc_curve2.5, "best", 
                   ret=c("threshold", "accuracy", "sensitivity", "specificity",
                         "precision", "npv","tpr", "fpr"), transpose = FALSE)
coords2.5

##CA125 HGSOC for comparisons###########################
KN_CA2 <- KN_BENIGN_HGSOC[!is.na(KN_BENIGN_HGSOC$CA125_f), ] #remove empty
table(KN_CA2$Grupė_Ieva)
KN_CA2$CA125_fN <- as.numeric(factor(KN_CA2$CA125_f))- 1
roc_curve_CA2 <- roc(KN_CA2$Grupė_Ieva, KN_CA2$CA125_fN , direction = ">")
plot(roc_curve_CA2) #auc = 0.765
auc(roc_curve_CA2)
coords_ca2 <- coords(roc_curve_CA2, "best", 
                     ret=c("threshold", "accuracy", "sensitivity", "specificity",
                           "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coords_ca2

##PLOT ROC COMBINATIONS HGSOC ##########################
#HGSOC vs benign models
roc_plot_4 <- function() {
par(pty = "s") #sets square
plot.roc(roc_curve2, print.auc = F, col = "#911eb4", lty = 2,
         cex.main=0.8, main ="Gerybinių pokyčių atskyrimas nuo HGSOC atvejų",
         xlab = "Specifiškumas",   # Custom x-axis label (e.g., in Lithuanian)
         ylab = "Jautrumas") #7
lines(roc_curve2.2, col = "#dcbeff", lwd =2 ) #6
lines(roc_curve2.3, col ="#fabed4", lwd =2, lty = 4) #8
lines(roc_curve2.4, col ="darkred", lwd =2, lty = 3) 
lines(roc_curve2.5, col ="deeppink", lwd =2) 
lines(roc_curve_CA2, col = "grey", lwd = 2)

# Add legend
legend("bottomright", legend = c( expression(italic("Genų raiškos žymenų kombinacija")), 
                                  expression(italic("Promotorių metilinimo statuso kombinacija")),
                                  expression(italic("Genų raiškos žymenų + promotorių metilinimo statuso kombinacija ")), 
                                  expression(italic("NOTCH genų raiškos žymenų kombinacija")),
                                  expression(italic("HOX promotorių metilinimo statuso kombinacija")),
                                  expression(italic("Serumo CA125 biožymens statusas")))
       ,
       col = c("#911eb4","#dcbeff", "#fabed4", "darkred",
               "deeppink", "grey" ), lty = 1, 
       cex = 0.8, lwd =3)
}
# Save the plot as a PNG file
png("metexprs_roc_HGSOC_output20250623.png", width = 1000, height = 1000, res = 150)
roc_plot_4()
dev.off()

##TABLE ROC COMBINATIONS HGSOC #####################
results_roc2<- data.frame(
  Biožymenys = c("Genų raiškos žymenų kombinacija", 
                 "Promotorių metilinimo statuso kombinacija",
                 "Genų raiškos žymenų + promotorių metilinimo statuso kombinacija", 
                 "NOTCH genų raiškos žymenų kombinacija",
                 "HOX promotorių metilinimo statuso kombinacija",
                 "Serumo CA125 biožymens statusas"),
  `plotas po kreive` = c(roc_curve2$auc, roc_curve2.2$auc, roc_curve2.3$auc, roc_curve2.4$auc, roc_curve2.5$auc, roc_curve_CA2$auc), 
  `slenkstinė vertė` = c(coords2$threshold, coords2.2$threshold , coords2.3$threshold,
                         coords2.4$threshold, coords2.5$threshold, coords_ca2$threshold ),
  tikslumas = c(coords2$accuracy, coords2.2$accuracy , coords2.3$accuracy,
               coords2.4$accuracy, coords2.5$accuracy, coords_ca2$accuracy ),
  jautrumas = c(coords2$sensitivity, coords2.2$sensitivity, coords2.3$sensitivity,
                  coords2.4$sensitivity, coords2.5$sensitivity, coords_ca2$sensitivity),
  specifiškumas = c(coords2$specificity, coords2.2$specificity, coords2.3$specificity,
                  coords2.4$specificity, coords2.5$specificity, coords_ca2$specificity),
  ppv  = c(coords2$precision, coords2.2$precision, coords2.3$precision,
                 coords2.4$precision,coords2.5$precision,coords_ca2$precision ),
  npv  = c(coords2$npv, coords2.2$npv, coords2.3$npv,
           coords2.4$npv, coords2.5$npv, coords_ca2$npv),
  tpr  = c(coords2$tpr, coords2.2$tpr, coords2.3$tpr,
           coords2.4$tpr,coords2.5$tpr,coords_ca2$tpr),
  fpr  = c(coords2$fpr, coords2.2$fpr, coords2.3$fpr,
           coords2.4$fpr, coords2.5$fpr,coords_ca2$fpr),
  check.names = FALSE
)
rownames(results_roc2) <- NULL
results_roc2

gt_table2 <- results_roc2 %>%
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
gt_table2

#there is no other convenient way to save gt outputs
gtsave(gt_table2,vwidth = 800,
       filename = "metexprs_table_HGSOC_output20250623.png")

#Combine the images
roc_image2<- image_read("metexprs_roc_HGSOC_output20250623.png")
table_image2 <- image_read("metexprs_table_HGSOC_output20250623.png")


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
            "metexprs_Roctable_HGSOC_output20250909.png")

##delong tests with CA125 HGSOC##########################
roc.test(roc_curve_CA2, roc_curve2)#0.0397 10 genes
roc.test(roc_curve_CA2, roc_curve2.2)#0.4493 4 methyl
roc.test(roc_curve_CA2, roc_curve2.3)#0.0397 14 biomakers
roc.test(roc_curve_CA2, roc_curve2.4)#0.0397 notch
roc.test(roc_curve_CA2, roc_curve2.5)#0.5977 #methyl HOX
#compare together
roc.test(roc_curve2.4, roc_curve2.5, method = "delong")#0.004806 #methyl vs notch
roc.test(roc_curve2.4, roc_curve2.2, method = "delong")#0.004806 #methyl 4 vs notch

#ROC COMBINATIONS hgsoc vs other #####################
# HGSOC vs other df
KN_OTHER_HGSOC <- KN_expr_full[KN_expr_full$Grupė_Ieva != "Benign", ] 
KN_OTHER_HGSOC$Grupė_Ieva <- droplevels(KN_OTHER_HGSOC$Grupė_Ieva)
table(KN_OTHER_HGSOC$Grupė_Ieva) #51 left 42 vs 9
#tumor must be a factor with reference: in this command first  level is a reference
KN_OTHER_HGSOC$Grupė_Ieva <- factor(KN_OTHER_HGSOC$Grupė_Ieva, levels = c("Other", "HGSOC"))
#13 other, 42 HGSOC
#MODELS:
##only expression#######################################
expr_HO <- KN_OTHER_HGSOC[colnames(KN_OTHER_HGSOC) %in% raiska]
group_OB<- factor(KN_OTHER_HGSOC$Grupė_Ieva, levels = c("Other", "HGSOC") )
logistic.model_HO_1<- glm(group_OB ~ ., data = expr_HO, family = "binomial") #will not converge
#to make it converge
brglm.model_x <- glm(group_OB ~ ., data = expr_HO, family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_x <- predict.glm(brglm.model_x, type='response') 
roc_curvex <- roc(group_OB, predicted_probs_x)
auc(roc_curvex)
plot(roc_curvex) #auc =  1
coordsx <- coords(roc_curvex, "best", 
                  ret=c("threshold", "accuracy", "sensitivity", "specificity",
                        "precision", "npv","tpr", "fpr"), transpose = FALSE)
coordsx

##only methylation#############################################
met_OB <- KN_OTHER_HGSOC[colnames(KN_OTHER_HGSOC) %in% metilinimas]
logistic.model_x1<- glm(group_OB ~ ., data = met_OB, family = "binomial") #will not converge
predicted_probs_x.1 <- predict.glm(logistic.model_x1, type='response') 
roc_curvex.1 <- roc(group_OB, predicted_probs_x.1)
auc(roc_curvex.1)
plot(roc_curvex.1) #auc =  0.87
coordsx.1 <- coords(roc_curvex.1, "best",
                    ret=c("threshold", "accuracy", "sensitivity", "specificity",
                          "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coordsx.1

##gene expression + methylation#########################################
expr_obx.3 <- KN_OTHER_HGSOC[colnames(KN_OTHER_HGSOC) %in% biomarkers]

#logistic.model_x.3<- glm(group_OB ~ ., data = expr_obx.3, family = "binomial") #will not converge
#to make it converge
brglm.model_x.3 <- glm(group_OB ~ ., data = expr_obx.3, family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_x.3 <- predict.glm(brglm.model_x.3, type='response') 
roc_curvex.3<- roc(group_OB, predicted_probs_x.3)
auc(roc_curvex.3)
plot(roc_curvex.3) #auc =  1
coordsx.3<- coords(roc_curvex.3, "best", 
                   ret=c("threshold", "accuracy", "sensitivity", "specificity",
                         "precision", "npv","tpr", "fpr"), transpose = FALSE)
coordsx.3

##only NOTCH genes######################################
expr_ob2.4 <- KN_OTHER_HGSOC[colnames(KN_OTHER_HGSOC) %in% NOTCH_genes]

logistic.model_x.4<- glm(group_OB ~ ., data = expr_ob2.4, family = "binomial") #will not converge
#to make it converge
#brglm.model_x.4 <- glm(group_OB ~ ., data = expr_ob2.4, family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_x.4 <- predict.glm(logistic.model_x.4, type='response') 
roc_curvex.4<- roc(group_OB, predicted_probs_x.4)
auc(roc_curvex.4)
plot(roc_curvex.4) #auc =  1
coordsx.4<- coords(roc_curvex.4, "best",
                   ret=c("threshold", "accuracy", "sensitivity", "specificity",
                         "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coordsx.4

##only HOX methylation ###########################################
expr_ohx.5 <- KN_OTHER_HGSOC[colnames(KN_OTHER_HGSOC) %in% HOX_genes]

logistic.model_x.5 <- glm(group_OB ~ ., data = expr_ohx.5, family = "binomial") #will not converge
predicted_probs_X.5 <- predict.glm(logistic.model_x.5, type='response') 
roc_curveX.5 <- roc(group_OB, predicted_probs_X.5)
auc(roc_curveX.5)
plot(roc_curveX.5) #auc =  0,85
coordsX.5<- coords(roc_curveX.5, "best",
                   ret=c("threshold", "accuracy", "sensitivity", "specificity",
                         "precision", "npv", "tpr", "fpr"), transpose = FALSE)
coordsX.5

## CA125 HGSOC vs OTHERS for comparisons##################################
KN_CA2X <- KN_OTHER_HGSOC[!is.na(KN_OTHER_HGSOC$CA125_f), ] #remove empty
table(KN_CA2X$Grupė_Ieva)
KN_CA2X$CA125_fN <- as.numeric(factor(KN_CA2X$CA125_f))- 1
roc_curve_CA2X <- roc(KN_CA2X$Grupė_Ieva, KN_CA2X$CA125_fN , direction = ">")
plot(roc_curve_CA2X) #auc = 0.765
auc(roc_curve_CA2X)
coords_ca2X <- coords(roc_curve_CA2X, "best", ret=c("threshold", "accuracy", "sensitivity", "specificity", "precision", "npv",
                                                  "tpr", "fpr"), transpose = FALSE)
coords_ca2X

## PLOT ROC COMBINATIONS HGSOC vs OTHERS##########################
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
  lines(roc_curve_CA2X, col = "grey", lwd = 2)
  
  # Add legend
  legend("bottomright", legend = c( expression(italic("Genų raiškos žymenų kombinacija")), 
                                    expression(italic("Promotorių metilinimo statuso kombinacija")),
                                    expression(italic("Genų raiškos žymenų + promotorių metilinimo statuso kombinacija ")), 
                                    expression(italic("NOTCH genų raiškos žymenų kombinacija")),
                                    expression(italic("HOX promotorių metilinimo statuso kombinacija")),
                                    expression(italic("Serumo CA125 biožymens statusas")))
         ,
         col = c("#911eb4","#dcbeff", "#fabed4", "darkred",
                 "deeppink", "grey" ), lty = 1, 
         cex = 0.73, lwd =3)
}

# Save the plot as a PNG file
png("metexprs_roc_HGSOC_OTHERS_MODELS_output20250610.png", width = 1000, height = 1000, res = 150)
roc_plot_5()
dev.off()

##TABLE ROC COMBINATIONS HGSOC vs OTHERS#####################
results_rocx<- data.frame(
  Biožymenys = c("Genų raiškos žymenų kombinacija", 
                 "Promotorių metilinimo statuso kombinacija",
                 "Genų raiškos žymenų + promotorių metilinimo statuso kombinacija", 
                 "NOTCH genų raiškos žymenų kombinacija",
                 "HOX promotorių metilinimo statuso kombinacija",
                 "Serumo CA125 biožymens statusas"),
  `plotas po kreive` = c(roc_curvex$auc, roc_curvex.1$auc, roc_curvex.3$auc, 
                         roc_curvex.4$auc, roc_curveX.5$auc, roc_curve_CA2X$auc), 
  `slenkstinė vertė` = c(coordsx$threshold, coordsx.1$threshold , coordsx.3$threshold,
                         coordsx.4$threshold, coordsX.5$threshold, coords_ca2X$threshold ),
  tikslumas = c(coordsx$accuracy, coordsx.1$accuracy , coordsx.3$accuracy,
                coordsx.4$accuracy, coordsX.5$accuracy, coords_ca2X$accuracy ),
  jautrumas = c(coordsx$sensitivity, coordsx.1$sensitivity, coordsx.3$sensitivity,
                coordsx.4$sensitivity, coordsX.5$sensitivity, coords_ca2X$sensitivity),
  specifiškumas = c(coordsx$specificity, coordsx.1$specificity, coordsx.3$specificity,
                    coordsx.4$specificity, coordsX.5$specificity, coords_ca2X$specificity),
  ppv  = c(coordsx$precision, coordsx.1$precision, coordsx.3$precision,
           coordsx.4$precision,coordsX.5$precision,coords_ca2X$precision ),
  npv  = c(coordsx$npv, coordsx.1$npv, coordsx.3$npv,
           coordsx.4$npv, coordsX.5$npv, coords_ca2X$npv),
  tpr  = c(coordsx$tpr, coordsx.1$tpr, coordsx.3$tpr,
           coordsx.4$tpr,coordsX.5$tpr,coords_ca2X$tpr),
  fpr  = c(coordsx$fpr, coordsx.1$fpr, coordsx.3$fpr,
           coordsx.4$fpr, coordsX.5$fpr,coords_ca2X$fpr),
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
       filename = "metexprs_table_HGSOC_OTHERS_MODELS_output20250610.png")

#Combine the images
roc_image2<- image_read("metexprs_roc_HGSOC_OTHERS_MODELS_output20250610.png")
table_image2 <- image_read("metexprs_table_HGSOC_OTHERS_MODELS_output20250610.png")

# Find the max width to align both
roc_info <- image_info(roc_image2)
table_info <- image_info(table_image2)
max_width <- max(roc_info$width, table_info$width)

# Pad each image to the max width
roc_image2_padded <- image_extent(roc_image2, geometry = geometry_area(max_width, roc_info$height), gravity = "center", color = "white")
table_image2_padded <- image_extent(table_image2, geometry = geometry_area(max_width, table_info$height), gravity = "center", color = "white")

# Now append vertically
combined_image2 <- image_append(c(roc_image2, table_image2), stack = F)

# Save the combined image
image_write(combined_image2, 
            "metexprs_tableroc_HGSOC_OTHERS_MODELS_output20250610.png")
##delong tests with CA125 for comparisons######################################
roc.test(roc_curve_CA2X, roc_curvex)#9.907e-12 10 genes
roc.test(roc_curve_CA2X, roc_curvex.1)#0.3469 4 methyl
roc.test(roc_curve_CA2X, roc_curvex.3)#8.825e-14 14 biomakers
roc.test(roc_curve_CA2X, roc_curvex.4)#2.548e-10 notch
roc.test(roc_curve_CA2X, roc_curveX.5)#0.5977 #methyl
roc.test(roc_curvex.3, roc_curvex, method = "delong")#0.4405 #methyl

#save objects####################################################
# Save the list
save( roc_curve1,
      roc_curve1.2,
      roc_curve1.3,
      roc_curve1.4,
      roc_curve1.5,
      roc_curve_CA,
      roc_curve2,
      roc_curve2.2,
      roc_curve2.3,
      roc_curve2.4,
      roc_curve2.5,
      roc_curve_CA2,
      roc_curvex,
      roc_curvex.1,
      roc_curvex.3,
      roc_curvex.4,
      roc_curveX.5,
      roc_curve_CA2X,
 file = "roc_list20250826.RData")

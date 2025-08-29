#KN-DISSERTATION project. 10 gene panel - genes selected via statistical analysis of TCGA and GTEx data
#KN 10 gene analysis data (cv 2 %, one threshold, 65 cases!, restored KN-59) 2025-03-05
#ROCs: Iterate logistic regression between each variable 
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
#set wd for outputs
setwd("C:/Users/Ieva/rprojects/outputs_all/DISS/")
#data###################################################################
OC_full <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES//OC_10_genes_clean_2025_02_14.RDS")
OC_full <- OC_full[c(OC_full$KN != "KN-100"), ]
#expression genes list
expression <- c("EXO1", "RAD50","PPT2", "LUC7L2","PKP3", "CDCA5","ZFPL1","VPS33B", "GRB7","TCEAL4")

#HGSOC vs BENIGN DF
OC_HGSOC_BENIGN<- OC_full[c(OC_full$type != "Other"),] #51 cases left
OC_HGSOC_BENIGN$tumor <- relevel(factor(OC_HGSOC_BENIGN$type), ref = "Benign")

#HGSOC VS OTHERS DF
OC_HGSOC_OTHERS<- OC_full[c(OC_full$type != "Benign"),] #56 cases left
OC_HGSOC_OTHERS$tumor <- relevel(factor(OC_HGSOC_OTHERS$type), ref = "Other")

#manually make combination of the best two HGSOC vs benign#####
genes2 <- c( "GRB7", "TCEAL4")
expr_tumor2 <- OC_HGSOC_BENIGN[colnames(OC_HGSOC_BENIGN) %in% genes2]
brglm.model_2 <- glm(OC_HGSOC_BENIGN$tumor ~ ., data = expr_tumor2, family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_2 <- predict.glm(brglm.model_2, type='response') 
pred_data2 <- OC_HGSOC_BENIGN[(rownames(OC_HGSOC_BENIGN) #remove incoplete rows
                               %in% names(predicted_probs_2)), ]
roc_curve2 <- roc(pred_data2$tumor, predicted_probs_2)
AUC2 <- auc(roc_curve2)
coords2 <- coords(roc_curve2,
                  "best", ret=c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
coords2$AUC <- AUC2
coords2

#HGSOC vs Benign to loop over pairs################
# Generate all possible pairs of column names
var_pairs <- combn(expression, 2, simplify = FALSE)
# Create an empty dataframe to store results
results_df <- data.frame(Gene1 = character(),
                         Gene2 = character(),
                         AUC = numeric(),
                         Accuracy = numeric(),
                         Sensitivity = numeric(),
                         Specificity = numeric(),
                         benign_n = integer(),
                         OC_n = integer(),
                         stringsAsFactors = FALSE)

# Loop over each gene pair
for (pair in var_pairs) {
  genes2 <- pair  # Extract gene names from pair
  
  # Select expression data for the current gene pair
  expr_tumor2 <- OC_HGSOC_BENIGN[colnames(OC_HGSOC_BENIGN) %in% genes2]
  
  # Run logistic regression
  brglm.model_2 <- glm(OC_HGSOC_BENIGN$tumor ~ ., data = expr_tumor2, 
                       family = binomial("logit"), method = "brglm_fit")
  
  # Get predicted probabilities
  predicted_probs_2 <- predict.glm(brglm.model_2, type = 'response')
  
  # Filter dataset to remove incomplete rows
  pred_data2 <- OC_HGSOC_BENIGN[rownames(OC_HGSOC_BENIGN) %in% names(predicted_probs_2), ]
  
  # Compute tumor group counts
  tumor_counts <- table(pred_data2$tumor)
  group_benign_n <- ifelse("Benign" %in% names(tumor_counts), tumor_counts["Benign"], 0)
  group_OC_n <- ifelse("HGSOC" %in% names(tumor_counts), tumor_counts["HGSOC"], 0)
  
  # Print table for each pair
  print(paste("Gene Pair:", genes2[1], "-", genes2[2]))
  print(tumor_counts)  # This prints the table for each gene pair
  
  # Compute ROC curve and AUC
  roc_curve2 <- roc(pred_data2$tumor, predicted_probs_2)
  AUC2 <- auc(roc_curve2)
  
  # Get best accuracy, sensitivity, specificity
  coords2 <- coords(roc_curve2, "best", ret = c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
  
  # Append results to dataframe
  results_df <- rbind(results_df, data.frame(Gene1 = genes2[1], 
                                             Gene2 = genes2[2], 
                                             AUC = AUC2, 
                                             Accuracy = coords2$accuracy, 
                                             Sensitivity = coords2$sensitivity, 
                                             Specificity = coords2$specificity,
                                             benign_n = group_benign_n,
                                             OC_n = group_OC_n))
}

# Sort results by AUC in descending order
results_df <- results_df %>%
  arrange(desc(AUC))

# Print final results dataframe
print(results_df)

#write.csv(results_df, "hgsc_benign_pair_combinations20250307.csv")

#HGSOC vs Benign to loop over trios#########
var_pairs3 <- combn(expression, 3, simplify = FALSE)
# Create an empty dataframe to store results
results_df3 <- data.frame(Gene1 = character(),
                         Gene2 = character(),
                         Gene3 = character(),
                         AUC = numeric(),
                         Accuracy = numeric(),
                         Sensitivity = numeric(),
                         Specificity = numeric(),
                         benign_n = integer(),
                         OC_n = integer(),
                         stringsAsFactors = FALSE)

# Loop over each gene TRIO
for (pair in var_pairs3) {
  genes3 <- pair  # Extract gene names from pair
  
  # Select expression data for the current gene pair
  expr_tumor2 <- OC_HGSOC_BENIGN[colnames(OC_HGSOC_BENIGN) %in% genes3]
  
  # Run logistic regression
  brglm.model_2 <- glm(OC_HGSOC_BENIGN$tumor ~ ., data = expr_tumor2, 
                       family = binomial("logit"), method = "brglm_fit")
  
  # Get predicted probabilities
  predicted_probs_2 <- predict.glm(brglm.model_2, type = 'response')
  
  # Filter dataset to remove incomplete rows
  pred_data2 <- OC_HGSOC_BENIGN[rownames(OC_HGSOC_BENIGN) %in% names(predicted_probs_2), ]
  
  # Compute tumor group counts
  tumor_counts <- table(pred_data2$tumor)
  group_benign_n <- ifelse("Benign" %in% names(tumor_counts), tumor_counts["Benign"], 0)
  group_OC_n <- ifelse("HGSOC" %in% names(tumor_counts), tumor_counts["HGSOC"], 0)
  
  # Print table for each pair
  print(paste("Gene Pair:", genes3[1], "-", genes3[2], "-", genes3[3]))
  print(tumor_counts)  # This prints the table for each gene pair
  
  # Compute ROC curve and AUC
  roc_curve2 <- roc(pred_data2$tumor, predicted_probs_2)
  AUC2 <- auc(roc_curve2)
  
  # Get best accuracy, sensitivity, specificity
  coords2 <- coords(roc_curve2, "best", best.method=c("closest.topleft"),
                      ret = c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
  
  # Append results to dataframe
  results_df3 <- rbind(results_df3, data.frame(Gene1 = genes3[1], 
                                             Gene2 = genes3[2], 
                                             Gene3 = genes3[3],
                                             AUC = AUC2, 
                                             Accuracy = coords2$accuracy, 
                                             Sensitivity = coords2$sensitivity, 
                                             Specificity = coords2$specificity,
                                             benign_n = group_benign_n,
                                             OC_n = group_OC_n))
}

# Sort results by AUC in descending order
results_df3 <- results_df3 %>%
  arrange(desc(AUC))

# Print final results dataframe
print(results_df3)

#write.csv(results_df3, "hgsc_benign_trio_combinations20250307.csv")

#HGSOC vs Benign to loop over fours #########
var_pairs4 <- combn(expression, 4, simplify = FALSE)
# Create an empty dataframe to store results
results_df4 <- data.frame(Gene1 = character(),
                          Gene2 = character(),
                          Gene3 = character(),
                          Gene4 = character(),
                          AUC = numeric(),
                          Accuracy = numeric(),
                          Sensitivity = numeric(),
                          Specificity = numeric(),
                          benign_n = integer(),
                          OC_n = integer(),
                          stringsAsFactors = FALSE)

# Loop over each gene pair
for (pair in var_pairs4) {
  genes3 <- pair  # Extract gene names from pair
  
  # Select expression data for the current gene pair
  expr_tumor2 <- OC_HGSOC_BENIGN[colnames(OC_HGSOC_BENIGN) %in% genes3]
  
  # Run logistic regression
  brglm.model_2 <- glm(OC_HGSOC_BENIGN$tumor ~ ., data = expr_tumor2, 
                       family = binomial("logit"), method = "brglm_fit")
  
  # Get predicted probabilities
  predicted_probs_2 <- predict.glm(brglm.model_2, type = 'response')
  
  # Filter dataset to remove incomplete rows
  pred_data2 <- OC_HGSOC_BENIGN[rownames(OC_HGSOC_BENIGN) %in% names(predicted_probs_2), ]
  
  # Compute tumor group counts
  tumor_counts <- table(pred_data2$tumor)
  group_benign_n <- ifelse("Benign" %in% names(tumor_counts), tumor_counts["Benign"], 0)
  group_OC_n <- ifelse("HGSOC" %in% names(tumor_counts), tumor_counts["HGSOC"], 0)
  
  # Print table for each pair
  print(paste("Gene Pair:", genes3[1], "-", genes3[2], "-", genes3[3], "-", genes3[4]))
  print(tumor_counts)  # This prints the table for each gene pair
  
  # Compute ROC curve and AUC
  roc_curve2 <- roc(pred_data2$tumor, predicted_probs_2)
  AUC2 <- auc(roc_curve2)
  
  # Get best accuracy, sensitivity, specificity
  coords2 <- coords(roc_curve2, "best", best.method=c("closest.topleft"),
                    ret = c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
  
  # Append results to dataframe
  results_df4 <- rbind(results_df4, data.frame(Gene1 = genes3[1], 
                                               Gene2 = genes3[2], 
                                               Gene3 = genes3[3],
                                               Gene4 = genes3[3],
                                               AUC = AUC2, 
                                               Accuracy = coords2$accuracy, 
                                               Sensitivity = coords2$sensitivity, 
                                               Specificity = coords2$specificity,
                                               benign_n = group_benign_n,
                                               OC_n = group_OC_n))
}

# Sort results by AUC in descending order
results_df4 <- results_df4 %>%
  arrange(desc(AUC))

# Print final results dataframe
print(results_df4)

#write.csv(results_df4, "hgsoc_benign_four_combinations20250307.csv")

#OC vs Benign to loop over pairs################
# Generate all possible pairs of column names
var_pairs <- combn(expression, 2, simplify = FALSE)
# Create an empty dataframe to store results
results_df_oc <- data.frame(Gene1 = character(),
                         Gene2 = character(),
                         AUC = numeric(),
                         Accuracy = numeric(),
                         Sensitivity = numeric(),
                         Specificity = numeric(),
                         benign_n = integer(),
                         OC_n = integer(),
                         stringsAsFactors = FALSE)

# Loop over each gene pair
for (pair in var_pairs) {
  genes2 <- pair  # Extract gene names from pair
  
  # Select expression data for the current gene pair
  expr_tumor2 <- OC_full[colnames(OC_full) %in% genes2]
  
  # Run logistic regression
  brglm.model_2 <- glm(OC_full$tumor ~ ., data = expr_tumor2, 
                       family = binomial("logit"), method = "brglm_fit")
  
  # Get predicted probabilities
  predicted_probs_2 <- predict.glm(brglm.model_2, type = 'response')
  
  # Filter dataset to remove incomplete rows
  pred_data2 <- OC_full[rownames(OC_full) %in% names(predicted_probs_2), ]
  
  # Compute tumor group counts
  tumor_counts <- table(pred_data2$tumor)
  group_benign_n <- ifelse("Benign" %in% names(tumor_counts), tumor_counts["Benign"], 0)
  group_OC_n <- ifelse("OC" %in% names(tumor_counts), tumor_counts["OC"], 0)
  
  # Print table for each pair
  print(paste("Gene Pair:", genes2[1], "-", genes2[2]))
  print(tumor_counts)  # This prints the table for each gene pair
  
  # Compute ROC curve and AUC
  roc_curve2 <- roc(pred_data2$tumor, predicted_probs_2)
  AUC2 <- auc(roc_curve2)
  
  # Get best accuracy, sensitivity, specificity
  coords2 <- coords(roc_curve2, "best", ret = c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
  
  # Append results to dataframe
  results_df_oc <- rbind(results_df_oc, data.frame(Gene1 = genes2[1], 
                                             Gene2 = genes2[2], 
                                             AUC = AUC2, 
                                             Accuracy = coords2$accuracy, 
                                             Sensitivity = coords2$sensitivity, 
                                             Specificity = coords2$specificity,
                                             benign_n = group_benign_n,
                                             OC_n = group_OC_n))
}

# Sort results by AUC in descending order
results_df_oc <- results_df_oc %>%
  arrange(desc(AUC))

# Print final results dataframe
print(results_df_oc)

#write.csv(results_df_oc, "oc_benign_pair_combinations20250307.csv")

#OC vs Benign to loop over trios################
# Generate all possible trios of column names
var_pairs3 <- combn(expression, 3, simplify = FALSE)
# Create an empty dataframe to store results
results_df_oc3 <- data.frame(Gene1 = character(),
                            Gene2 = character(),
                            Gene3 = character(),
                            AUC = numeric(),
                            Accuracy = numeric(),
                            Sensitivity = numeric(),
                            Specificity = numeric(),
                            benign_n = integer(),
                            OC_n = integer(),
                            stringsAsFactors = FALSE)

# Loop over each gene pair
for (pair in var_pairs3) {
  genes2 <- pair  # Extract gene names from pair
  
  # Select expression data for the current gene pair
  expr_tumor2 <- OC_full[colnames(OC_full) %in% genes2]
  
  # Run logistic regression
  brglm.model_2 <- glm(OC_full$tumor ~ ., data = expr_tumor2, 
                       family = binomial("logit"), method = "brglm_fit")
  
  # Get predicted probabilities
  predicted_probs_2 <- predict.glm(brglm.model_2, type = 'response')
  
  # Filter dataset to remove incomplete rows
  pred_data2 <- OC_full[rownames(OC_full) %in% names(predicted_probs_2), ]
  
  # Compute tumor group counts
  tumor_counts <- table(pred_data2$tumor)
  group_benign_n <- ifelse("Benign" %in% names(tumor_counts), tumor_counts["Benign"], 0)
  group_OC_n <- ifelse("OC" %in% names(tumor_counts), tumor_counts["OC"], 0)
  
  # Print table for each pair
  print(paste("Gene Pair:", genes2[1], "-", genes2[2], "-", genes2[3]))
  print(tumor_counts)  # This prints the table for each gene pair
  
  # Compute ROC curve and AUC
  roc_curve2 <- roc(pred_data2$tumor, predicted_probs_2)
  AUC2 <- auc(roc_curve2)
  
  # Get best accuracy, sensitivity, specificity
  coords2 <- coords(roc_curve2, "best", best.method=c("closest.topleft"),
  ret = c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
  
  # Append results to dataframe
  results_df_oc3 <- rbind(results_df_oc3, data.frame(Gene1 = genes2[1], 
                                                   Gene2 = genes2[2], 
                                                   Gene3 = genes2[3],
                                                   AUC = AUC2, 
                                                   Accuracy = coords2$accuracy, 
                                                   Sensitivity = coords2$sensitivity, 
                                                   Specificity = coords2$specificity,
                                                   benign_n = group_benign_n,
                                                   OC_n = group_OC_n))
}

# Sort results by AUC in descending order
results_df_oc3 <- results_df_oc3 %>%
  arrange(desc(AUC))

# Print final results dataframe
print(results_df_oc3)

#write.csv(results_df_oc3, "oc_benign_trio_combinations20250307.csv")

#OC vs Benign to loop over fours #########
var_pairs4 <- combn(expression, 4, simplify = FALSE)
# Create an empty dataframe to store results
results_df4oc <- data.frame(Gene1 = character(),
                          Gene2 = character(),
                          Gene3 = character(),
                          Gene4 = character(),
                          AUC = numeric(),
                          Accuracy = numeric(),
                          Sensitivity = numeric(),
                          Specificity = numeric(),
                          benign_n = integer(),
                          OC_n = integer(),
                          stringsAsFactors = FALSE)

# Loop over each gene pair
for (pair in var_pairs4) {
  genes3 <- pair  # Extract gene names from pair
  
  # Select expression data for the current gene pair
  expr_tumor2 <- OC_full[colnames(OC_full) %in% genes3]
  
  # Run logistic regression
  brglm.model_2 <- glm(OC_full$tumor ~ ., data = expr_tumor2, 
                       family = binomial("logit"), method = "brglm_fit")
  
  # Get predicted probabilities
  predicted_probs_2 <- predict.glm(brglm.model_2, type = 'response')
  
  # Filter dataset to remove incomplete rows
  pred_data2 <- OC_full[rownames(OC_full) %in% names(predicted_probs_2), ]
  
  # Compute tumor group counts
  tumor_counts <- table(pred_data2$tumor)
  group_benign_n <- ifelse("Benign" %in% names(tumor_counts), tumor_counts["Benign"], 0)
  group_OC_n <- ifelse("OC" %in% names(tumor_counts), tumor_counts["OC"], 0)
  
  # Print table for each pair
  print(paste("Gene Pair:", genes3[1], "-", genes3[2], "-", genes3[3], "-", genes3[4]))
  print(tumor_counts)  # This prints the table for each gene pair
  
  # Compute ROC curve and AUC
  roc_curve2 <- roc(pred_data2$tumor, predicted_probs_2)
  AUC2 <- auc(roc_curve2)
  
  # Get best accuracy, sensitivity, specificity
  coords2 <- coords(roc_curve2, "best", best.method=c("closest.topleft"),
                    ret = c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
  
  # Append results to dataframe
  results_df4oc <- rbind(results_df4oc, data.frame(Gene1 = genes3[1], 
                                               Gene2 = genes3[2], 
                                               Gene3 = genes3[3],
                                               Gene4 = genes3[4],
                                               AUC = AUC2, 
                                               Accuracy = coords2$accuracy, 
                                               Sensitivity = coords2$sensitivity, 
                                               Specificity = coords2$specificity,
                                               benign_n = group_benign_n,
                                               OC_n = group_OC_n))
}

# Sort results by AUC in descending order
results_df4oc <- results_df4oc %>%
  arrange(desc(AUC))

# Print final results dataframe
print(results_df4oc)

#write.csv(results_df4oc, "OC_benign_four_combinations20250305.csv")

#HGSOC vs OTHERS to loop over pairs################
# Generate all possible pairs of column names
var_pairs <- combn(expression, 2, simplify = FALSE)
# Create an empty dataframe to store results
results_df_oc <- data.frame(Gene1 = character(),
                            Gene2 = character(),
                            AUC = numeric(),
                            Accuracy = numeric(),
                            Sensitivity = numeric(),
                            Specificity = numeric(),
                            benign_n = integer(),
                            OC_n = integer(),
                            stringsAsFactors = FALSE)

# Loop over each gene pair
for (pair in var_pairs) {
  genes2 <- pair  # Extract gene names from pair
  
  # Select expression data for the current gene pair
  expr_tumor2 <- OC_HGSOC_OTHERS[colnames(OC_HGSOC_OTHERS) %in% genes2]
  
  # Run logistic regression
  brglm.model_2 <- glm(OC_HGSOC_OTHERS$tumor ~ ., data = expr_tumor2, 
                       family = binomial("logit"), method = "brglm_fit")
  
  # Get predicted probabilities
  predicted_probs_2 <- predict.glm(brglm.model_2, type = 'response')
  
  # Filter dataset to remove incomplete rows
  pred_data2 <- OC_HGSOC_OTHERS[rownames(OC_HGSOC_OTHERS) %in% names(predicted_probs_2), ]
  
  # Compute tumor group counts
  tumor_counts <- table(pred_data2$tumor)
  group_benign_n <- ifelse("Other" %in% names(tumor_counts), tumor_counts["Other"], 0)
  group_OC_n <- ifelse("HGSOC" %in% names(tumor_counts), tumor_counts["HGSOC"], 0)
  
  # Print table for each pair
  print(paste("Gene Pair:", genes2[1], "-", genes2[2]))
  print(tumor_counts)  # This prints the table for each gene pair
  
  # Compute ROC curve and AUC
  roc_curve2 <- roc(pred_data2$tumor, predicted_probs_2)
  AUC2 <- auc(roc_curve2)
  
  # Get best accuracy, sensitivity, specificity
  coords2 <- coords(roc_curve2, "best",best.method=c("closest.topleft"),
                    ret = c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
  
  # Append results to dataframe
  results_df_oc <- rbind(results_df_oc, data.frame(Gene1 = genes2[1], 
                                                   Gene2 = genes2[2], 
                                                   AUC = AUC2, 
                                                   Accuracy = coords2$accuracy, 
                                                   Sensitivity = coords2$sensitivity, 
                                                   Specificity = coords2$specificity,
                                                   benign_n = group_benign_n,
                                                   OC_n = group_OC_n))
}

# Sort results by AUC in descending order
results_df_oc <- results_df_oc %>%
  arrange(desc(AUC))

# Print final results dataframe
print(results_df_oc)

#write.csv(results_df_oc, "hgsoc_others_pair_combinations20250307.csv")

#HGSOC vs Others to loop over trios################
# Generate all possible trios of column names
var_pairs3 <- combn(expression, 3, simplify = FALSE)
# Create an empty dataframe to store results
results_df_other3 <- data.frame(Gene1 = character(),
                             Gene2 = character(),
                             Gene3 = character(),
                             AUC = numeric(),
                             Accuracy = numeric(),
                             Sensitivity = numeric(),
                             Specificity = numeric(),
                             benign_n = integer(),
                             OC_n = integer(),
                             stringsAsFactors = FALSE)

# Loop over each gene pair
for (pair in var_pairs3) {
  genes2 <- pair  # Extract gene names from pair
  
  # Select expression data for the current gene pair
  expr_tumor2 <- OC_HGSOC_OTHERS[colnames(OC_HGSOC_OTHERS) %in% genes2]
  
  # Run logistic regression
  brglm.model_2 <- glm(OC_HGSOC_OTHERS$tumor ~ ., data = expr_tumor2, 
                       family = binomial("logit"), method = "brglm_fit")
  
  # Get predicted probabilities
  predicted_probs_2 <- predict.glm(brglm.model_2, type = 'response')
  
  # Filter dataset to remove incomplete rows
  pred_data2 <- OC_HGSOC_OTHERS[rownames(OC_HGSOC_OTHERS) %in% names(predicted_probs_2), ]
  
  # Compute tumor group counts
  tumor_counts <- table(pred_data2$tumor)
  group_benign_n <- ifelse("Other" %in% names(tumor_counts), tumor_counts["Other"], 0)
  group_OC_n <- ifelse("HGSOC" %in% names(tumor_counts), tumor_counts["HGSOC"], 0)
  
  # Print table for each pair
  print(paste("Gene Pair:", genes2[1], "-", genes2[2], "-", genes2[3]))
  print(tumor_counts)  # This prints the table for each gene pair
  
  # Compute ROC curve and AUC
  roc_curve2 <- roc(pred_data2$tumor, predicted_probs_2)
  AUC2 <- auc(roc_curve2)
  
  # Get best accuracy, sensitivity, specificity
  coords2 <- coords(roc_curve2, "best", best.method=c("closest.topleft"),
                    ret = c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
  
  # Append results to dataframe
  results_df_other3 <- rbind(results_df_other3, data.frame(Gene1 = genes2[1], 
                                                     Gene2 = genes2[2], 
                                                     Gene3 = genes2[3],
                                                     AUC = AUC2, 
                                                     Accuracy = coords2$accuracy, 
                                                     Sensitivity = coords2$sensitivity, 
                                                     Specificity = coords2$specificity,
                                                     benign_n = group_benign_n,
                                                     OC_n = group_OC_n))
}

# Sort results by AUC in descending order
results_df_other3 <- results_df_other3 %>%
  arrange(desc(AUC))

# Print final results dataframe
print(results_df_other3)

#write.csv(results_df_other3, 
#"C:/Users/user/Documents/rprojects/ALL OUTPUTS/KN_10_genes/hgsoc_others_trio_combinations20250307.csv")

#HGSOC vs Others to loop over fours #########
var_pairs4 <- combn(expression, 4, simplify = FALSE)
# Create an empty dataframe to store results
results_df4other <- data.frame(Gene1 = character(),
                            Gene2 = character(),
                            Gene3 = character(),
                            Gene4 = character(),
                            AUC = numeric(),
                            Accuracy = numeric(),
                            Sensitivity = numeric(),
                            Specificity = numeric(),
                            benign_n = integer(),
                            OC_n = integer(),
                            stringsAsFactors = FALSE)

# Loop over each gene pair
for (pair in var_pairs4) {
  genes3 <- pair  # Extract gene names from pair
  
  # Select expression data for the current gene pair
  expr_tumor2 <- OC_HGSOC_OTHERS[colnames(OC_HGSOC_OTHERS) %in% genes3]
  
  # Run logistic regression
  brglm.model_2 <- glm(OC_HGSOC_OTHERS$tumor ~ ., data = expr_tumor2, 
                       family = binomial("logit"), method = "brglm_fit")
  
  # Get predicted probabilities
  predicted_probs_2 <- predict.glm(brglm.model_2, type = 'response')
  
  # Filter dataset to remove incomplete rows
  pred_data2 <- OC_HGSOC_OTHERS[rownames(OC_HGSOC_OTHERS) %in% names(predicted_probs_2), ]
  
  # Compute tumor group counts
  tumor_counts <- table(pred_data2$tumor)
  group_benign_n <- ifelse("Other" %in% names(tumor_counts), tumor_counts["Other"], 0)
  group_OC_n <- ifelse("HGSOC" %in% names(tumor_counts), tumor_counts["HGSOC"], 0)
  
  # Print table for each pair
  print(paste("Gene Pair:", genes3[1], "-", genes3[2], "-", genes3[3], "-", genes3[4]))
  print(tumor_counts)  # This prints the table for each gene pair
  
  # Compute ROC curve and AUC
  roc_curve2 <- roc(pred_data2$tumor, predicted_probs_2)
  AUC2 <- auc(roc_curve2)
  
  # Get best accuracy, sensitivity, specificity
  coords2 <- coords(roc_curve2, "best", best.method=c("closest.topleft"),
                    ret = c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
  
  #fix the two line pronlem
  if (nrow(coords2) > 1) {
    # If multiple best thresholds, prioritize Youden first, then Closest to (0,1)
    coords2 <- coords(roc_curve2, "best", best.method = "youden",
                      ret = c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
    
    #if (nrow(coords2) > 1) {
     # coords2 <- coords(roc_curve2, "best", best.method = "closest.topleft")}
    }
  
  
  # Append results to dataframe
  results_df4other <- rbind(results_df4other, data.frame(Gene1 = genes3[1], 
                                                   Gene2 = genes3[2], 
                                                   Gene3 = genes3[3],
                                                   Gene4 = genes3[4],
                                                   AUC = AUC2, 
                                                   Accuracy = coords2$accuracy, 
                                                   Sensitivity = coords2$sensitivity, 
                                                   Specificity = coords2$specificity,
                                                   benign_n = group_benign_n,
                                                   OC_n = group_OC_n))
}

# Sort results by AUC in descending order
results_df4other <- results_df4other %>%
  arrange(desc(AUC))

# Print final results dataframe
str(results_df4other)

#write.csv(results_df4other, "hgosc_other_four_combinations20250307.csv")

#HGSOC vs Others to loop over fives #########
var_pairs5 <- combn(expression, 5, simplify = FALSE)
# Create an empty dataframe to store results
results_df5other <- data.frame(Gene1 = character(),
                               Gene2 = character(),
                               Gene3 = character(),
                               Gene4 = character(),
                               Gene5 = character(),
                               AUC = numeric(),
                               Accuracy = numeric(),
                               Sensitivity = numeric(),
                               Specificity = numeric(),
                               benign_n = integer(),
                               OC_n = integer(),
                               stringsAsFactors = FALSE)

# Loop over each gene pair
for (pair in var_pairs5) {
  genes5 <- pair  # Extract gene names from pair
  
  # Select expression data for the current gene pair
  expr_tumor2 <- OC_HGSOC_OTHERS[colnames(OC_HGSOC_OTHERS) %in% genes5]
  
  # Run logistic regression
  brglm.model_2 <- glm(OC_HGSOC_OTHERS$tumor ~ ., data = expr_tumor2, 
                       family = binomial("logit"), method = "brglm_fit")
  
  # Get predicted probabilities
  predicted_probs_2 <- predict.glm(brglm.model_2, type = 'response')
  
  # Filter dataset to remove incomplete rows
  pred_data2 <- OC_HGSOC_OTHERS[rownames(OC_HGSOC_OTHERS) %in% names(predicted_probs_2), ]
  
  # Compute tumor group counts
  tumor_counts <- table(pred_data2$tumor)
  group_benign_n <- ifelse("Other" %in% names(tumor_counts), tumor_counts["Other"], 0)
  group_OC_n <- ifelse("HGSOC" %in% names(tumor_counts), tumor_counts["HGSOC"], 0)
  
  # Print table for each pair
  print(paste("Gene Pair:", genes5[1], "-", genes5[2], "-", genes5[3], "-", genes5[4], "-", genes5[5]))
  print(tumor_counts)  # This prints the table for each gene pair
  
  # Compute ROC curve and AUC
  roc_curve2 <- roc(pred_data2$tumor, predicted_probs_2)
  AUC2 <- auc(roc_curve2)
  
  # Get best accuracy, sensitivity, specificity
  coords2 <- coords(roc_curve2, "best", best.method=c("closest.topleft"),
                    ret = c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
  
  #fix the two line pronlem
  if (nrow(coords2) > 1) {
    # If multiple best thresholds, prioritize Youden first, then Closest to (0,1)
    coords2 <- coords(roc_curve2, "best", best.method = "youden",
                      ret = c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
    
    #if (nrow(coords2) > 1) {
    # coords2 <- coords(roc_curve2, "best", best.method = "closest.topleft")}
  }
  
  
  # Append results to dataframe
  results_df5other <- rbind(results_df5other, data.frame(Gene1 = genes5[1], 
                                                         Gene2 = genes5[2], 
                                                         Gene3 = genes5[3],
                                                         Gene4 = genes5[4],
                                                         Gene5 = genes5[5],
                                                         AUC = AUC2, 
                                                         Accuracy = coords2$accuracy, 
                                                         Sensitivity = coords2$sensitivity, 
                                                         Specificity = coords2$specificity,
                                                         benign_n = group_benign_n,
                                                         OC_n = group_OC_n))
}

# Sort results by AUC in descending order
results_df5other <- results_df5other %>%
  arrange(desc(AUC))

# Print final results dataframe
head(results_df5other)

#write.csv(results_df5other, "hgosc_other_five_combinations20250425.csv")

#HGSOC vs Others to loop over six #########
var_pairs6 <- combn(expression, 6, simplify = FALSE)
# Create an empty dataframe to store results
results_df6other <- data.frame(Gene1 = character(),
                               Gene2 = character(),
                               Gene3 = character(),
                               Gene4 = character(),
                               Gene5 = character(),
                               Gene6 = character(),
                               AUC = numeric(),
                               Accuracy = numeric(),
                               Sensitivity = numeric(),
                               Specificity = numeric(),
                               benign_n = integer(),
                               OC_n = integer(),
                               stringsAsFactors = FALSE)

# Loop over each gene pair
for (pair in var_pairs6) {
  genes6 <- pair  # Extract gene names from pair
  
  # Select expression data for the current gene pair
  expr_tumor2 <- OC_HGSOC_OTHERS[colnames(OC_HGSOC_OTHERS) %in% genes6]
  
  # Run logistic regression
  brglm.model_2 <- glm(OC_HGSOC_OTHERS$tumor ~ ., data = expr_tumor2, 
                       family = binomial("logit"), method = "brglm_fit")
  
  # Get predicted probabilities
  predicted_probs_2 <- predict.glm(brglm.model_2, type = 'response')
  
  # Filter dataset to remove incomplete rows
  pred_data2 <- OC_HGSOC_OTHERS[rownames(OC_HGSOC_OTHERS) %in% names(predicted_probs_2), ]
  
  # Compute tumor group counts
  tumor_counts <- table(pred_data2$tumor)
  group_benign_n <- ifelse("Other" %in% names(tumor_counts), tumor_counts["Other"], 0)
  group_OC_n <- ifelse("HGSOC" %in% names(tumor_counts), tumor_counts["HGSOC"], 0)
  
  # Print table for each pair
  print(paste("Gene Pair:", genes6[1], "-", genes6[2], "-", genes6[3],
              "-", genes6[4], "-", genes6[5], "-", genes6[6]))
  print(tumor_counts)  # This prints the table for each gene pair
  
  # Compute ROC curve and AUC
  roc_curve2 <- roc(pred_data2$tumor, predicted_probs_2)
  AUC2 <- auc(roc_curve2)
  
  # Get best accuracy, sensitivity, specificity
  coords2 <- coords(roc_curve2, "best", best.method=c("closest.topleft"),
                    ret = c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
  
  #fix the two line pronlem
  if (nrow(coords2) > 1) {
    # If multiple best thresholds, prioritize Youden first, then Closest to (0,1)
    coords2 <- coords(roc_curve2, "best", best.method = "youden",
                      ret = c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
    
    #if (nrow(coords2) > 1) {
    # coords2 <- coords(roc_curve2, "best", best.method = "closest.topleft")}
  }
  
  
  # Append results to dataframe
  results_df6other <- rbind(results_df6other, data.frame(Gene1 = genes6[1], 
                                                         Gene2 = genes6[2], 
                                                         Gene3 = genes6[3],
                                                         Gene4 = genes6[4],
                                                         Gene5 = genes6[5],
                                                         Gene6 = genes6[6],
                                                         AUC = AUC2, 
                                                         Accuracy = coords2$accuracy, 
                                                         Sensitivity = coords2$sensitivity, 
                                                         Specificity = coords2$specificity,
                                                         benign_n = group_benign_n,
                                                         OC_n = group_OC_n))
}

# Sort results by AUC in descending order
results_df6other <- results_df6other %>%
  arrange(desc(AUC))

# Print final results dataframe
head(results_df6other)

#write.csv(results_df6other, "hgosc_other_six_combinations20250425.csv")

#HGSOC vs Others to loop over sevens #########
var_pairs7 <- combn(expression, 7, simplify = FALSE)
# Create an empty dataframe to store results
results_df7other <- data.frame(Gene1 = character(),
                               Gene2 = character(),
                               Gene3 = character(),
                               Gene4 = character(),
                               Gene5 = character(),
                               Gene6 = character(),
                               Gene7 = character(),
                               AUC = numeric(),
                               Accuracy = numeric(),
                               Sensitivity = numeric(),
                               Specificity = numeric(),
                               benign_n = integer(),
                               OC_n = integer(),
                               stringsAsFactors = FALSE)

# Loop over each gene pair
for (pair in var_pairs7) {
  genes7 <- pair  # Extract gene names from pair
  
  # Select expression data for the current gene pair
  expr_tumor2 <- OC_HGSOC_OTHERS[colnames(OC_HGSOC_OTHERS) %in% genes7]
  
  # Run logistic regression
  brglm.model_2 <- glm(OC_HGSOC_OTHERS$tumor ~ ., data = expr_tumor2, 
                       family = binomial("logit"), method = "brglm_fit")
  
  # Get predicted probabilities
  predicted_probs_2 <- predict.glm(brglm.model_2, type = 'response')
  
  # Filter dataset to remove incomplete rows
  pred_data2 <- OC_HGSOC_OTHERS[rownames(OC_HGSOC_OTHERS) %in% names(predicted_probs_2), ]
  
  # Compute tumor group counts
  tumor_counts <- table(pred_data2$tumor)
  group_benign_n <- ifelse("Other" %in% names(tumor_counts), tumor_counts["Other"], 0)
  group_OC_n <- ifelse("HGSOC" %in% names(tumor_counts), tumor_counts["HGSOC"], 0)
  
  # Print table for each pair
  print(paste("Gene Pair:", genes7[1], "-", genes7[2], "-", genes7[3],
              "-", genes7[4], "-", genes7[5], "-", genes7[6], "-", genes7[7]))
  print(tumor_counts)  # This prints the table for each gene pair
  
  # Compute ROC curve and AUC
  roc_curve2 <- roc(pred_data2$tumor, predicted_probs_2)
  AUC2 <- auc(roc_curve2)
  
  # Get best accuracy, sensitivity, specificity
  coords2 <- coords(roc_curve2, "best", best.method=c("closest.topleft"),
                    ret = c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
  
  #fix the two line pronlem
  if (nrow(coords2) > 1) {
    # If multiple best thresholds, prioritize Youden first, then Closest to (0,1)
    coords2 <- coords(roc_curve2, "best", best.method = "youden",
                      ret = c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
    
    #if (nrow(coords2) > 1) {
    # coords2 <- coords(roc_curve2, "best", best.method = "closest.topleft")}
  }
  
  
  # Append results to dataframe
  results_df7other <- rbind(results_df7other, data.frame(Gene1 = genes7[1], 
                                                         Gene2 = genes7[2], 
                                                         Gene3 = genes7[3],
                                                         Gene4 = genes7[4],
                                                         Gene5 = genes7[5],
                                                         Gene6 = genes7[6],
                                                         Gene6 = genes7[7],
                                                         AUC = AUC2, 
                                                         Accuracy = coords2$accuracy, 
                                                         Sensitivity = coords2$sensitivity, 
                                                         Specificity = coords2$specificity,
                                                         benign_n = group_benign_n,
                                                         OC_n = group_OC_n))
}

# Sort results by AUC in descending order
results_df7other <- results_df7other %>%
  arrange(desc(AUC))

# Print final results dataframe
head(results_df7other)

#write.csv(results_df7other,  "hgosc_other_seven_combinations20250425.csv")

#HGSOC vs Others to loop over eights #########
var_pairs8 <- combn(expression, 8, simplify = FALSE)
# Create an empty dataframe to store results
results_df8other <- data.frame(Gene1 = character(),
                               Gene2 = character(),
                               Gene3 = character(),
                               Gene4 = character(),
                               Gene5 = character(),
                               Gene6 = character(),
                               Gene7 = character(),
                               Gene8 = character(),
                               AUC = numeric(),
                               Accuracy = numeric(),
                               Sensitivity = numeric(),
                               Specificity = numeric(),
                               benign_n = integer(),
                               OC_n = integer(),
                               stringsAsFactors = FALSE)

# Loop over each gene pair
for (pair in var_pairs8) {
  genes8 <- pair  # Extract gene names from pair
  
  # Select expression data for the current gene pair
  expr_tumor2 <- OC_HGSOC_OTHERS[colnames(OC_HGSOC_OTHERS) %in% genes8]
  
  # Run logistic regression
  brglm.model_2 <- glm(OC_HGSOC_OTHERS$tumor ~ ., data = expr_tumor2, 
                       family = binomial("logit"), method = "brglm_fit")
  
  # Get predicted probabilities
  predicted_probs_2 <- predict.glm(brglm.model_2, type = 'response')
  
  # Filter dataset to remove incomplete rows
  pred_data2 <- OC_HGSOC_OTHERS[rownames(OC_HGSOC_OTHERS) %in% names(predicted_probs_2), ]
  
  # Compute tumor group counts
  tumor_counts <- table(pred_data2$tumor)
  group_benign_n <- ifelse("Other" %in% names(tumor_counts), tumor_counts["Other"], 0)
  group_OC_n <- ifelse("HGSOC" %in% names(tumor_counts), tumor_counts["HGSOC"], 0)
  
  # Print table for each pair
  print(paste("Gene Pair:", genes8[1], "-", genes8[2], "-", genes8[3],
              "-", genes8[4], "-", genes8[5], "-", genes8[6], "-", genes8[7],
              "-", genes8[8]))
  print(tumor_counts)  # This prints the table for each gene pair
  
  # Compute ROC curve and AUC
  roc_curve2 <- roc(pred_data2$tumor, predicted_probs_2)
  AUC2 <- auc(roc_curve2)
  
  # Get best accuracy, sensitivity, specificity
  coords2 <- coords(roc_curve2, "best", best.method=c("closest.topleft"),
                    ret = c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
  
  #fix the two line pronlem
  # If multiple rows, pick the one with max accuracy, then sensitivity if needed
  if (nrow(coords2) > 1) {
    # Keep only rows with the max accuracy
    max_acc <- max(coords2$accuracy)
    coords2 <- coords2[coords2$accuracy == max_acc, ]
    
    # If still multiple, keep one with highest sensitivity
    if (nrow(coords2) > 1) {
      max_sens <- max(coords2$sensitivity)
      coords2 <- coords2[coords2$sensitivity == max_sens, ]
    }
    
    # If still multiple, just take the first one (or you could sort by threshold, etc.)
    if (nrow(coords2) > 1) {
      coords2 <- coords2[1, , drop = FALSE]
    }
  }
  
  
  # Append results to dataframe
  results_df8other <- rbind(results_df8other, data.frame(Gene1 = genes8[1], 
                                                         Gene2 = genes8[2], 
                                                         Gene3 = genes8[3],
                                                         Gene4 = genes8[4],
                                                         Gene5 = genes8[5],
                                                         Gene6 = genes8[6],
                                                         Gene7 = genes8[7],
                                                         Gene8 = genes8[8],
                                                         AUC = AUC2, 
                                                         Accuracy = coords2$accuracy, 
                                                         Sensitivity = coords2$sensitivity, 
                                                         Specificity = coords2$specificity,
                                                         benign_n = group_benign_n,
                                                         OC_n = group_OC_n))
}

# Sort results by AUC in descending order
results_df8other <- results_df8other %>%
  arrange(desc(AUC))

# Print final results dataframe
head(results_df8other)

#write.csv(results_df8other, "hgosc_other_eight_combinations20250425.csv")

#HGSOC vs Others to loop over nines #########
var_pairs9 <- combn(expression, 9, simplify = FALSE)
# Create an empty dataframe to store results
results_df9other <- data.frame(Gene1 = character(),
                               Gene2 = character(),
                               Gene3 = character(),
                               Gene4 = character(),
                               Gene5 = character(),
                               Gene6 = character(),
                               Gene7 = character(),
                               Gene8 = character(),
                               Gene9 = character(),
                               AUC = numeric(),
                               Accuracy = numeric(),
                               Sensitivity = numeric(),
                               Specificity = numeric(),
                               benign_n = integer(),
                               OC_n = integer(),
                               stringsAsFactors = FALSE)

# Loop over each gene pair
for (pair in var_pairs9) {
  genes9 <- pair  # Extract gene names from pair
  
  # Select expression data for the current gene pair
  expr_tumor2 <- OC_HGSOC_OTHERS[colnames(OC_HGSOC_OTHERS) %in% genes9]
  
  # Run logistic regression
  brglm.model_2 <- glm(OC_HGSOC_OTHERS$tumor ~ ., data = expr_tumor2, 
                       family = binomial("logit"), method = "brglm_fit")
  
  # Get predicted probabilities
  predicted_probs_2 <- predict.glm(brglm.model_2, type = 'response')
  
  # Filter dataset to remove incomplete rows
  pred_data2 <- OC_HGSOC_OTHERS[rownames(OC_HGSOC_OTHERS) %in% names(predicted_probs_2), ]
  
  # Compute tumor group counts
  tumor_counts <- table(pred_data2$tumor)
  group_benign_n <- ifelse("Other" %in% names(tumor_counts), tumor_counts["Other"], 0)
  group_OC_n <- ifelse("HGSOC" %in% names(tumor_counts), tumor_counts["HGSOC"], 0)
  
  # Print table for each pair
  print(paste("Gene Pair:", genes9[1], "-", genes9[2], "-", genes9[3],
              "-", genes9[4], "-", genes9[5], "-", genes9[6], "-", genes9[7],
              "-", genes9[8],"-", genes9[9]))
  print(tumor_counts)  # This prints the table for each gene pair
  
  # Compute ROC curve and AUC
  roc_curve2 <- roc(pred_data2$tumor, predicted_probs_2)
  AUC2 <- auc(roc_curve2)
  
  # Get best accuracy, sensitivity, specificity
  coords2 <- coords(roc_curve2, "best", best.method=c("closest.topleft"),
                    ret = c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
  
  #fix the two line pronlem
  # If multiple rows, pick the one with max accuracy, then sensitivity if needed
  if (nrow(coords2) > 1) {
    # Keep only rows with the max accuracy
    max_acc <- max(coords2$accuracy)
    coords2 <- coords2[coords2$accuracy == max_acc, ]
    
    # If still multiple, keep one with highest sensitivity
    if (nrow(coords2) > 1) {
      max_sens <- max(coords2$sensitivity)
      coords2 <- coords2[coords2$sensitivity == max_sens, ]
    }
    
    # If still multiple, just take the first one (or you could sort by threshold, etc.)
    if (nrow(coords2) > 1) {
      coords2 <- coords2[1, , drop = FALSE]
    }
  }
  
  
  # Append results to dataframe
  results_df9other <- rbind(results_df9other, data.frame(Gene1 = genes9[1], 
                                                         Gene2 = genes9[2], 
                                                         Gene3 = genes9[3],
                                                         Gene4 = genes9[4],
                                                         Gene5 = genes9[5],
                                                         Gene6 = genes9[6],
                                                         Gene7 = genes9[7],
                                                         Gene8 = genes9[8],
                                                         Gene9 = genes9[9],
                                                         AUC = AUC2, 
                                                         Accuracy = coords2$accuracy, 
                                                         Sensitivity = coords2$sensitivity, 
                                                         Specificity = coords2$specificity,
                                                         benign_n = group_benign_n,
                                                         OC_n = group_OC_n))
}

# Sort results by AUC in descending order
results_df9other <- results_df9other %>%
  arrange(desc(AUC))

# Print final results dataframe
head(results_df9other)

#write.csv(results_df9other,"hgosc_other_nine_combinations20250425.csv")

#all 10 HGSOC vs other#####
expr_tumor10 <- OC_HGSOC_OTHERS[colnames(OC_HGSOC_OTHERS) %in% expression]
brglm.model_10 <- glm(OC_HGSOC_OTHERS$tumor ~ ., data = expr_tumor10,
                      family = binomial("logit"), method = "brglm_fit") #no problems
predicted_probs_10<- predict.glm(brglm.model_10, type='response') 
pred_data10 <- OC_HGSOC_OTHERS[(rownames(OC_HGSOC_OTHERS) #remove incoplete rows
                               %in% names(predicted_probs_10)), ]
roc_curve10 <- roc(pred_data10$tumor, predicted_probs_10)
AUC10 <- auc(roc_curve10)
coords10 <- coords(roc_curve10,
 "best", ret=c("accuracy", "sensitivity", "specificity"), transpose = FALSE)
coords10$AUC <- AUC10
coords10

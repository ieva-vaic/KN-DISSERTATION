#KN-DISSERTATION project. Genes selection via statistical analysis of TCGA and GTEx data 
#(only scripts producing lithuanian figures, included in the dissertation)
##THIS IS TCGA-OV-RISK-GENES project script No. 13
#TEST the risk model
#ROCs, boxplots with TCGA, test data
# Load packages ##########################################
Sys.setenv(LANG = "en")
library(tidyverse)
library(data.table)
library(glmnet)
library(tidyverse)
library(gplots)
library(survminer)
library(survivalROC)
library(survival)
library(gridExtra)
library(grid)  
library(timeROC)
library(RColorBrewer) 
library(ggprism)
library(patchwork)
library(rstatix) 
library(magick)
library(gt)
library(multcomp)
library(FSA)
library(cowplot)
#set directory of the data
setwd("../TCGA-OV-RISK-PROJECT/Public data RDSs/")
# Load test data ###################################
gtex_counts_test <- readRDS("test_gtcga_normcounts_prot_2025.RDS")
#filter for lasso genes
gtex_genes <- readRDS("gtcga_elastic_2025.RDS")
gtex_filtered_counts_test <- gtex_counts_test[colnames(gtex_counts_test) %in% gtex_genes] 
#filter for TCGA cases
gtex_filtered_counts_test2 <- gtex_filtered_counts_test %>%
  dplyr::filter(grepl("^TCGA", rownames(gtex_filtered_counts_test)))
dim(gtex_filtered_counts_test2)
# Load clinical data ###################################
clin_df <- readRDS("joinedTCGA_XENA_clinical2025.RDS")
#fix survival
clin_df$deceased = clin_df$vital_status == "Dead"
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)
#Filter clinical data###################################
test_ids <- rownames(gtex_filtered_counts_test2)
clin_df <- clin_df[clin_df$barcode %in% test_ids, ]  #79 samples
rownames(clin_df) <- clin_df$barcode

#Join clinical (mostly survival) data with train data ########################
#add all genes to clin_df
colnames(clin_df)
gtex_filtered_counts_test2$barcode <- rownames(gtex_filtered_counts_test2)
clin_df$barcode == gtex_filtered_counts_test2$barcode 
clin_df_joined_test <- left_join(clin_df, gtex_filtered_counts_test2, by = "barcode")
rownames(clin_df_joined_test) <- clin_df_joined_test$barcode
colnames(clin_df_joined_test)

#GTEX vs TCGA, boxplot ###########################
#CREATE GROUPINGS ACCORDING TO DATA#
snames = rownames(gtex_counts_test)
group = substr(snames, 1, 4)
group = as.factor(group)
levels(group) <- c("GTEx", "TCGA-OV")
gtex_counts_test2 <- gtex_counts_test
gtex_counts_test2$group <- group
#get genes of interest
expression <- c( "EXO1",   "RAD50",  "PPT2",   "LUC7L2", "PKP3",
                 "CDCA5",  "ZFPL1" , "VPS33B", "GRB7",   "TCEAL4")
#get long df
gtcga_table_full <- reshape2::melt(gtex_counts_test2[, colnames(gtex_counts_test2) %in%
                                                        c("group", expression )],
                                   id.vars="group",
                                   measure.vars= expression)
#get t test
t.test_gtex <- gtcga_table_full %>%
  group_by(variable) %>%
  t_test(value ~ group,
         p.adjust.method = "BH", 
         var.equal = TRUE, #stjudents
         paired = FALSE, 
         detailed=TRUE 
  ) # Format p-values to remove scientific notation
t.test_gtex
#make a tibble of p values: ~group1, ~group2, ~p.adj,   ~y.position, ~variable
t.test_gtex_tibble <- t.test_gtex %>% 
  dplyr::select(group1, group2, p, variable) %>%
  mutate(
    y.position = c(6, 8, 8, 8.5, 10, 
                   8, 6, 6, 13, 11) #choose where to plot p values
  )
t.test_gtex_tibble$p_custom <- ifelse(t.test_gtex_tibble$p < 0.001, 
                                      "p < 0.001", 
                                      paste0("p = ", sprintf("%.3f",
                                                             each.vs.ref_sig$pj)))
#get colors 
custom_colors <- c("GTEx" = "darkgreen","TCGA-OV" = "purple") 
#plot
gtex_plot <- ggplot(gtcga_table_full, aes(x=group , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = group )) +
  geom_jitter(aes(color = group ), size=1, alpha=0.5) +
  ylab(label = expression("Genų raiška")) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(t.test_gtex_tibble, label = "p_custom") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  ggtitle("Genų raiška testavimo imtyje")+
  coord_cartesian(clip = "off")  

#show plot
gtex_plot
#save short
png("C:/Users/Ieva/rprojects/outputs_all/DISS/test_barplot_20260127.png",
    width = 15, height = 10, res = 300, units = "cm") 
gtex_plot +
  theme(
    strip.text.x = element_text(size = 8, face = "bold.italic"),
    axis.text = element_text(size = 7),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 9, hjust = 0.5),
    plot.margin = margin(5, 5, 5, 5)
  )
dev.off() 

#save gtex vs tca test plot
# png("C:/Users/Ieva/rprojects/outputs_all/DISS/test_barplot_20260121.png",
#     width = 18, height = 17, res = 500, units = "cm") # width and height in pixels, resolution in dpi
# gtex_plot #
# dev.off() # Close the PNG device

# Load cox model ###################################
cvfit <- readRDS("coxnet_cvfit_2025.RDS")
cox_fitx <- readRDS("coxnet_fit_2025.RDS")

#get coeficients - gene names
coef_x <- coef(cox_fitx, s = 0.088)
head(coef_x)
coef_x = coef_x[coef_x[,1] != 0,] 
res_coef_cox_names = names(coef_x) # get names of the (non-zero) variables.
res_coef_cox_names #10

#get FC
gtex_counts_test_my <- gtex_counts_test2[, colnames(gtex_counts_test2) %in%
                                           c("group", expression )]
#drop group 
gtex_counts_test_my$group
group <- gtex_counts_test_my$group
gtex_counts_test_my <- gtex_counts_test_my[, !(names(gtex_counts_test_my) %in% "group")] 
#get group means
group1_mean <- colMeans(gtex_counts_test_my[group == "GTEx", ])
group2_mean <- colMeans(gtex_counts_test_my[group == "TCGA-OV", ])
#get logFC  and FC
log2FC <- group2_mean - group1_mean
FC <- 2^log2FC
#show
fc_table <- data.frame(
  Gene = colnames(gtex_counts_test_my),
  log2FC = log2FC,
  FoldChange = FC
)
#show fc
print(fc_table)

#RISK SCORE ##################################
# - cv_model: Your trained cv.glmnet model
cv_model <- cvfit
# - gene_data: A data frame (or matrix) with the expression values for your genes 
# (rows = samples, columns = genes)
# It should have the same feature names as used in the model (genes).
gene_data_test <- clin_df_joined_test[, res_coef_cox_names]
# - gene_list: A vector with the list of genes you're interested in.
res_coef_cox_names
# coefs
# Extract coefficients at optimal lambda
coefs <- coef_x
coefs_df <- as.data.frame(as.matrix(coefs))
coefs_df <- rownames_to_column(coefs_df, var = "Feature")
colnames(coefs_df)[2] <- "Coefficient"
print(coefs_df)
# Calculate the risk score: linear combination of gene expressions and coefficients
# Risk score = sum( gene_expression * coefficient )
risk_scores_test <- rowSums(sweep(gene_data_test, 2, coefs_df$Coefficient, "*"))
# View the risk scores
print(risk_scores_test) # now I have some risk scores
#add risk scores to the clin_df_joined_test
clin_df_joined_test$RiskScore <- risk_scores_test[rownames(clin_df_joined_test)]
#create df wih survival data
surv_df_test <- clin_df_joined_test[, colnames(clin_df_joined_test) %in%
        c("deceased", "overall_survival", res_coef_cox_names, "RiskScore")]
surv_df_test <- surv_df_test %>%
  dplyr::rename(censor = deceased, surv_time = overall_survival) 
#TIME ROC, Risk score #######################
#create features for timeroc
nobs <- NROW(surv_df_test)
#make surf df but only of my genes!
time <- surv_df_test$surv_time
event <- surv_df_test$censor
#time roc
t_eval <- c(365, 1095, 1825)  # time points
roc_result <- timeROC(
  T = time,       # Survival time from df
  delta = event, # Event indicator from df
  marker = surv_df_test[, "RiskScore"], # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = t_eval,    # Time points for ROC
  iid = TRUE         # Compute confidence intervals
)
roc_result
#time rocs for separate biomarkers ######################################

coxnet.df <- surv_df_test[, (colnames(surv_df_test) %in% res_coef_cox_names)]
dim(coxnet.df)

rez_list <- apply(coxnet.df, 2, timeROC,
                   T = time,       # Survival time from df
                   delta = event, # Event indicator from df
                   #marker  # Predictor already in the df
                   cause = 1,         # Event of interest
                   times = t_eval,    # Time points for ROC
                   iid = TRUE )        # Compute confidence intervals)

auc_table <- map_dfr(names(rez_list), function(gene) {
  roc <- rez_list[[gene]]
  
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

auc_table
##plot at year 1###############
# Choose target time
target_time <- 365
time_index <- which(rez_list[[1]]$times == target_time)

# Set up base plot with gene 1
plot(
  rez_list[[1]]$FP[, time_index],
  rez_list[[1]]$TP[, time_index],
  type = "l",
  col = 1,
  lwd = 2,
  xlab = "1 - Specificity (FPR)",
  ylab = "Sensitivity (TPR)",
  main = paste("Time-dependent ROC Curves at", target_time, "days"),
  xlim = c(0, 1),
  ylim = c(0, 1),
  asp = 1
)

# Add ROC lines for all genes
for (i in 2:length(rez_list)) {
  lines(
    rez_list[[i]]$FP[, time_index],
    rez_list[[i]]$TP[, time_index],
    col = i,
    lwd = 2
  )
}

# Add risk score ROC line in bold black
lines(
  roc_result$FP[, time_index],
  roc_result$TP[, time_index],
  col = "maroon",
  lwd = 3,
  lty = 1
)

# Add diagonal reference line
abline(0, 1, lty = 2, col = "gray")

# Build legend names: italic gene names + "risk score"
legend_labels <- c(
  parse(text = paste0("italic('", names(rez_list), "')")),
  "Risk Score"
)

# Add legend (last color is black for risk score)
legend(
  "bottomright",
  legend = legend_labels,
  col = c(1:length(rez_list), "maroon"),
  lwd = c(rep(2, length(rez_list)), 3),
  cex = 0.6,
  bty = "n"
)
##plot at year 3###############
# Choose target time
target_time <- 1095    
time_index <- which(rez_list[[1]]$times == target_time)

# Set up base plot with gene 1
plot(
  rez_list[[1]]$FP[, time_index],
  rez_list[[1]]$TP[, time_index],
  type = "l",
  col = 1,
  lwd = 2,
  xlab = "1 - Specificity (FPR)",
  ylab = "Sensitivity (TPR)",
  main = paste("Time-dependent ROC Curves at", target_time, "days"),
  xlim = c(0, 1),
  ylim = c(0, 1),
  asp = 1
)

# Add ROC lines for all genes
for (i in 2:length(rez_list)) {
  lines(
    rez_list[[i]]$FP[, time_index],
    rez_list[[i]]$TP[, time_index],
    col = i,
    lwd = 2
  )
}

# Add risk score ROC line in bold black
lines(
  roc_result$FP[, time_index],
  roc_result$TP[, time_index],
  col = "maroon",
  lwd = 3,
  lty = 1
)

# Add diagonal reference line
abline(0, 1, lty = 2, col = "gray")

# Build legend names: italic gene names + "risk score"
legend_labels <- c(
  parse(text = paste0("italic('", names(rez_list), "')")),
  "Risk Score"
)

# Add legend (last color is black for risk score)
legend(
  "bottomright",
  legend = legend_labels,
  col = c(1:length(rez_list), "maroon"),
  lwd = c(rep(2, length(rez_list)), 3),
  cex = 0.6,
  bty = "n"
)

## plot at year 5 with saving###########################
# Choose target time
target_time <- 1825   # choose year 5
time_index <- which(rez_list[[1]]$times == target_time)

# plot with save
png("C:/Users/Ieva/rprojects/outputs_all/DISS/10_genų_timeROC_test20251124.png",
    width = 1000, height = 1000, res = 200)

par(pty="s")

# Base plot with first gene
plot(
  rez_list[[1]]$FP[, time_index],
  rez_list[[1]]$TP[, time_index],
  type = "l",
  col = 1,
  lwd = 2,
  xlab = "1 - Specifiškumas",
  ylab = "Jautrumas",
  main = paste("Nuo laiko priklausomos ROC kreivės,\n5 metai po diagnozės testavimo imtyje"),
  xlim = c(0, 1),
  ylim = c(0, 1),
  asp = 1
)

# Add ROC lines for all genes
for (i in 2:length(rez_list)) {
  lines(
    rez_list[[i]]$FP[, time_index],
    rez_list[[i]]$TP[, time_index],
    col = i,
    lwd = 2
  )
}

# Add risk score ROC line in bold
lines(
  roc_result$FP[, time_index],
  roc_result$TP[, time_index],
  col = "maroon",
  lwd = 3,
  lty = 1
)

# Add diagonal reference line
abline(0, 1, lty = 2, col = "gray")

# Build legend labels without AUC
legend_labels <- c(paste0("italic('", names(rez_list), "')"),
                   "'Rizikos balas'")

# Add legend
legend(
  "bottomright",
  legend = parse(text = legend_labels),
  col = c(1:length(rez_list), "maroon"),
  lwd = c(rep(2, length(rez_list)), 3),
  cex = 0.6,
  bty = "n"
)

# Add panel label
mtext("D", side = 3, line = 2.5, adj = -0.2, font = 2, cex = 1.5)

dev.off()

#COORDS table #####################################################

# extract best sens/spec at a given time
extract_coords <- function(roc, gene, timepoint = 1825) {
  idx <- which(roc$times == timepoint)
  
  # sensitivity & specificity across thresholds
  sens <- as.vector(roc$TP[, idx])
  spec <- as.vector(1 - roc$FP[, idx])
  
  # Youden index to pick best threshold
  youden <- sens + spec - 1
  best_idx <- which.max(youden)
  
  tibble(
    gene = gene,
    time = roc$times[idx],
    auc = roc$AUC[idx] ,
    sens = sens[best_idx],
    spec = spec[best_idx],
    cutoff = roc$cutoffs[best_idx]
  )
}

# Apply to each biomarker
sens_spec_auc_60 <- map_dfr(names(rez_list), ~ extract_coords(rez_list[[.x]], .x, 1825))

# Apply to risk score
risk_score_row <- extract_coords(roc_result, "10 Genų raiškos rizikos balas", 1825)

# Combine into one tibble
sens_spec_auc_60_all <- bind_rows(sens_spec_auc_60, risk_score_row)
# View
sens_spec_auc_60_all

#make prety table
# Prepare table: rename columns for display
roc_table_display <- sens_spec_auc_60_all %>%
  rename(
    Biožymuo = gene,
    Laikas_dienom = time,
    AUC = auc,
    Jautrumas = sens,
    Specifiškumas = spec
  )

# Create gt table
gt_table_roc_60 <- roc_table_display %>%
  gt() %>%
  tab_header(
    title = "ROC kriterijai",
    subtitle = "Prognostinis tikslumas 5 metų išgyvenimui"
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
       filename = "C:/Users/Ieva/rprojects/outputs_all/DISS/TEST_timeroc_table_20251002.png")

#Combine the images
roc_image <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/10_genų_timeROC_test20251124.png")
table_image <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/TEST_timeroc_table_20251002.png")

combined_image <- image_append(c(roc_image, table_image), stack = F)

# Save the combined image
image_write(combined_image, 
            "C:/Users/Ieva/rprojects/outputs_all/DISS/TESTROC_W_TABLE20251124.png")

#KM plot with RISK SCORE#################
#Kaplan-meier plot ##################################
# Calculate the median risk score
median_risk <- median(clin_df_joined_test$RiskScore, na.rm = TRUE) #-0.05914624
# Create a new factor column based on the median value
clin_df_joined_test$RiskGroup <- ifelse(clin_df_joined_test$RiskScore <= median_risk, "Low Risk", "High Risk")
#Create a survival object
surv_object <- Surv(time = clin_df_joined_test$overall_survival, event = clin_df_joined_test$deceased )
# Fit a Kaplan-Meier model
km_fit <- survfit(surv_object ~ RiskGroup, data = clin_df_joined_test)
# Compute log-rank test
surv_test <- survdiff(Surv(overall_survival, deceased) ~ RiskGroup, data = clin_df_joined_test)
pval_num <- 1 - pchisq(surv_test$chisq, length(surv_test$n) - 1)
# Format nicely
if (pval_num < 0.001) {
  pval_text <- "Log-rank p < 0.001"
} else {
  pval_text <- paste0("Log-rank p = ", signif(pval_num, 3))
}

# Plot
train_surv <- ggsurvplot(
  km_fit, 
  data = clin_df_joined_test, 
  pval = FALSE,  # disable built-in
  risk.table = TRUE,
  risk.table.title = "Pacienčių skaičius rizikos grupėje",
  title = "Didelės vs. mažos rizikos atvejai testavimo imtyje",
  xlab = "Bendras išgyvenamumo laikas",
  ylab = "Išgyvenamumo tikimybė",
  palette = c("turquoise", "deeppink"),
  legend.title = "Rizikos grupė", 
  legend.labs = c("Mažos rizikos balas", "Didelės rizikos balas")
)

# Add subtitle
train_surv$plot <- train_surv$plot + labs(subtitle = pval_text)
# Add a big "B" to the top-left corner
print(train_surv)
#save km plot
png("C:/Users/Ieva/rprojects/outputs_all/DISS/dis_lt_km_test20260123.png",
    width = 15, height = 13, res = 100, units = "cm") # width and height in pixels, resolution in dpi
train_surv #
dev.off() # Close the PNG device

# #Forest plot################
# colnames(clin_df_joined_test)
# clin_df_joined_test$neoplasmhistologicgrade
# #first, fix stage
# clin_df_joined_test <- clin_df_joined_test %>%
#   mutate(stage_early_late = case_when(
#     clinicalstage2 %in% c("Stage I", "Stage II") ~ "early stage",
#     clinicalstage2 %in% c("Stage III", "Stage IV") ~ "late stage",
#     TRUE ~ NA_character_
#   ))
# #first, fix grade
# clin_df_joined_test$neoplasmhistologicgrade <- recode(clin_df_joined_test$neoplasmhistologicgrade,
#                                                  "GX" = NA_character_)
# vars <- c("RiskScore", "stage_early_late",
#           "neoplasmhistologicgrade", "ageatinitialpathologicdiagnosis" )
# 
# #model
# models <- vars %>%       # begin with variables of interest
#   str_c("deceased ~ ", .) %>%   # combine each variable into formula ("outcome ~ variable of interest")
#   # iterate through each univariate formula
#   map(                               
#     .f = ~glm(                       # pass the formulas one-by-one to glm()
#       formula = as.formula(.x),      # within glm(), the string formula is .x
#       family = "binomial",           # specify type of glm (logistic)
#       data = clin_df_joined_test)) %>%          # dataset
#   
#   # tidy up each of the glm regression outputs from above
#   map(
#     .f = ~tidy(
#       .x, 
#       exponentiate = TRUE,           # exponentiate 
#       conf.int = TRUE)) %>%          # return confidence intervals
#   
#   # collapse the list of regression outputs in to one data frame
#   bind_rows() %>% 
#   
#   # round all numeric columns
#   mutate(across(where(is.numeric), round, digits = 2)) %>%
#   filter(term != "(Intercept)")
# 
# models #simple table models
# 
# # fix terms
# models$term <-  c("Risk Score", "FIGO Stage late vs early",
#                   "Neoplasm histologic grade G3 vs G2", "Age at initial pathologic diagnosis")
# term_order <-  c("Risk Score", "FIGO Stage late vs early",
#                  "Neoplasm histologic grade G3 vs G2", "Age at initial pathologic diagnosis")
# # Apply the order to the term variable in the models data frame
# models$term <- factor(models$term, levels = term_order)
# models
# 
# #forest plot
# p <- 
#   models |>
#   ggplot(aes(y = fct_rev(term))) + 
#   theme_classic()+
#   geom_point(aes(x=estimate), shape=15, size=3) +
#   geom_linerange(aes(xmin=conf.low, xmax=conf.high)) +
#   geom_vline(xintercept = 1, linetype="dashed") +
#   labs(x="Odds Ratio, 95% CI", y="")+
#   theme(
#     axis.text.y = element_text(face = "italic")  # Make y-axis labels italic
#   )+
#   scale_x_continuous(
#     breaks = c(1, seq(-3,  10, by = 1))  # Increase the number of x-axis ticks
#   )
# 
# 
# odds1 <- p +geom_text(aes(y = term, x = 2, label = sprintf("%0.2f", round(estimate, digits = 2))), size = 4) + ## decimal places
#   annotate("text", x = 8, y= 4.5,  label = "OR", size = 4, fontface = "bold", hjust = 0)+
#   geom_text(aes(y = term, x =6, label = paste(conf.low, "-", conf.high)), size = 4) + 
#   annotate("text", x =5, y= 4.5,  label = "95 % CI", size =4, fontface = "bold", hjust = 0)+
#   geom_text(aes(y = term, x =4, label = sprintf("%0.2f", round(p.value, digits = 2))), size = 4) + 
#   annotate("text", x = 7, y= 4.5,  label = "p value", size = 4, fontface = "bold", hjust = 0)+
#   labs(title = "Hazard Ratios for Mortality")
# odds1

#TCGA - stage, boxplot ###########################
#CREATE GROUPINGS ACCORDING TO DATA#
table(clin_df_joined_test$clinicalstage2, useNA ="a")
table(clin_df_joined_test$stage_early_late, useNA ="a")
#filter for data of interest
colnames(clin_df_joined_test)
stage_test <- clin_df_joined_test[, colnames(clin_df_joined_test) %in% c(expression,
                                                                         "clinicalstage2", "barcode" )]
stage_test <- stage_test[!is.na(stage_test$clinicalstage2), ]
##anova, stage################################
#levene test for variance, stage, tcga
levene_rez <- stage_test[, c(2:12)] %>%
  pivot_longer(cols = -clinicalstage2 , names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = car::leveneTest(value ~ clinicalstage2 )$`Pr(>F)`[1],
    .groups = "drop"
  )
levene_rez  #all normal
#anova

results <- lapply(names(stage_test)[c(3:12)], function(var) {
  formula <- as.formula(paste(var, "~ clinicalstage2"))
  fit <- aov(formula, data = stage_test)
  
  list(
    variable = var,
    summary = summary(fit),
    shapiro = shapiro.test(residuals(fit))
  )
})

names(results) <- names(stage_test)[c(3:12)]
#make df of results
results_df <- lapply(results, function(x) {
  data.frame(
    variable   = x$variable,
    anova_p    = x$summary[[1]][["Pr(>F)"]][1],
    shapiro_p  = x$shapiro$p.value
  )
}) %>%
  bind_rows()

results_df
#show not normal
non_normal_df <- results_df %>% filter(shapiro_p <= 0.05)
non_normal_df
#show significant
results_df %>% filter(anova_p <= 0.05) #GRB7 

##stage Kruskal test######################################################
kruskal_results <- lapply(non_normal_df$variable, function(var) {
  df <- stage_test[, c(var, "clinicalstage2")]
  # df <- na.omit(df)  # remove NAs
  test <- kruskal.test(as.formula(paste(var, "~ clinicalstage2")), data = df)
  
  data.frame(
    variable = var,
    kruskal_p = test$p.value
  )
}) %>%
  bind_rows()

kruskal_results# not significant at all
##plot stage #########
#get long df
stage_table_full <- reshape2::melt(stage_test[, colnames(stage_test) %in%
                                                c("clinicalstage2", expression )],
                                   id.vars="clinicalstage2",
                                   measure.vars= expression)
stage_table_full <- stage_table_full %>%
  mutate(clinicalstage2 = case_when(
    clinicalstage2 == "Stage I" ~ "I stadija",
    clinicalstage2 == "Stage II" ~ "II stadija",
    clinicalstage2 == "Stage III" ~ "III stadija",
    clinicalstage2 == "Stage IV" ~ "IV stadija",
    TRUE ~ NA_character_
  ))

#get colors 
custom_colors_stage <- c("II stadija" = "lightgreen","III stadija" = "green",
                         "IV stadija" = "darkgreen") 

#plot
stgae_plot <- ggplot(stage_table_full, aes(x=clinicalstage2 , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = clinicalstage2 )) +
  geom_jitter(aes(color = clinicalstage2 ), size=1, alpha=0.5) +
  ylab(label = expression("Genų raiška")) + 
  facet_wrap(.~ factor(variable), nrow = 2, scales = "free") +
  theme_minimal()+
  theme(
    strip.text.x = element_text( size = 12, face = "bold.italic"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values = custom_colors_stage) +
  scale_color_manual(values = custom_colors_stage) +
  ggtitle("Genų raiška pagal su stadiją testavimo imtyje")+
  coord_cartesian(clip = "off")  

stgae_plot

stage_plot2 <- ggdraw(stgae_plot) +
  draw_plot_label(label = "B", x = 0, y = 1, hjust = 0, vjust = 1, size = 20, fontface = "bold")
#save short
png("C:/Users/Ieva/rprojects/outputs_all/DISS/test_stage_barplot_20260129.png",
    width = 15, height = 11, res = 300, units = "cm") 
stage_plot2 
dev.off() 
#save stage boxplot 
# png("C:/Users/Ieva/rprojects/outputs_all/DISS/test_stage_barplot_20260121.png",
#     width = 32, height = 19, res = 500, units = "cm") # width and height in pixels, resolution in dpi
# stage_plot2 #
# dev.off() # Close the PNG device

#TCGA - grade, boxplot ###########################
#CREATE GROUPINGS ACCORDING TO DATA#
table(clin_df_joined_test$neoplasmhistologicgrade, useNA ="a")
clin_df_joined_test <- clin_df_joined_test %>%
  mutate(neoplasmhistologicgrade = case_when(
    neoplasmhistologicgrade %in% c( "GX") ~ NA_character_,
    TRUE ~ neoplasmhistologicgrade
  ))

#filter for data of interest
colnames(clin_df_joined_test)
grade_test <- clin_df_joined_test[, colnames(clin_df_joined_test) %in% c(expression,
                                                                         "neoplasmhistologicgrade", "barcode" )]
grade_test <- grade_test[!is.na(grade_test$neoplasmhistologicgrade), ]
#get long df
grade_table_full <- reshape2::melt(grade_test[, colnames(grade_test) %in%
                                                c("neoplasmhistologicgrade", expression )],
                                   id.vars="neoplasmhistologicgrade",
                                   measure.vars= expression)
#get t test
t.test_grade <- grade_table_full %>%
  group_by(variable) %>%
  t_test(value ~ neoplasmhistologicgrade,
         p.adjust.method = "BH", 
         var.equal = TRUE, #stjudents
         paired = FALSE, 
         detailed=TRUE 
  ) # Format p-values to remove scientific notation
t.test_grade
#make a tibble of p values: ~group1, ~group2, ~p.adj,   ~y.position, ~variable
t.test_grade_tibble <- t.test_grade %>% 
  dplyr::select(group1, group2, p, variable) %>%
  mutate(
    y.position = c(6, 8, 8, 8.5, 13, 
                   8, 6, 6, 8, 11) #choose where to plot p values
  )%>%
  filter(p < 0.05)
t.test_grade_tibble

#get colors 
custom_colors_stage <- c("G2" = "lightgreen", "G3"= "darkgreen") 
#plot
grade_plot <- ggplot(grade_table_full, aes(x=neoplasmhistologicgrade , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = neoplasmhistologicgrade )) +
  geom_jitter(aes(color = neoplasmhistologicgrade ), size=1, alpha=0.5) +
  ylab(label = expression("Genų raiška")) + 
  facet_wrap(.~ factor(variable), nrow = 2, scales = "free") +
  add_pvalue(t.test_grade_tibble, label = "p") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values = custom_colors_stage) +
  scale_color_manual(values = custom_colors_stage) +
  ggtitle("Genų raiška pagal su naviko diferenciacijos laipsnį testavimo imtyje")+
  coord_cartesian(clip = "off")  
grade_plot

grade_plot2 <- ggdraw(grade_plot) +
  draw_plot_label(label = "B", x = 0, y = 1, hjust = 0, vjust = 1, size = 20, fontface = "bold")

grade_plot2
#save short
png("C:/Users/Ieva/rprojects/outputs_all/DISS/test_grade_barplot_20260129.png",
    width = 15, height = 11, res = 300, units = "cm") 
grade_plot2
dev.off() 
#save grade boxplot
# png("C:/Users/Ieva/rprojects/outputs_all/DISS/test_grade_barplot_20260121.png",
#     width = 32, height = 19, res = 500, units = "cm") # width and height in pixels, resolution in dpi
# grade_plot2 #
# dev.off() # Close the PNG device

#Coreliation with age#################
table(clin_df_joined_test$age_at_diagnosis, useNA ="a")
age_test <- clin_df_joined_test[, colnames(clin_df_joined_test) %in% c(expression,
                                                                       "age_at_diagnosis")]

stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "bottom", size = 4)
plot_list <- lapply(expression, function(gene) {
  ggplot(age_test, aes_string(x = "age_at_diagnosis", y = gene)) +
    geom_point(color = "steelblue") +
    geom_smooth(method = "lm", se = FALSE, color = "darkred") +
    stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "bottom", size = 4) +
    labs(
      title = bquote("Amžius vs"~italic(.(gene))),
      x = "Amžius metais",
      y = bquote("Raiška"~italic(.(gene)))
    ) +
    theme_minimal()
})

# Show plots individually
for (p in plot_list) print(p)

combined_plot <- wrap_plots(plot_list, ncol = 2) +
  plot_annotation(title = "Ryšys tarp amžiaus ir genų raiškos testavimo imtyje")

# Save to PNG
ggsave(
  filename = "C:/Users/Ieva/rprojects/outputs_all/DISS/age_vs_expressiontest20250619.png",
  plot = combined_plot,
  width = 8,       # adjust width as needed
  height = 11,       # adjust height as needed
  dpi = 300         # high resolution
)

#TCGA - lymphatic invasion, boxplot ###########################
table(clin_df_joined_test$lymphaticinvasion, useNA ="a")
clin_df_joined_test$lymphaticinvasion <- factor(clin_df_joined_test$lymphaticinvasion)
levels(clin_df_joined_test$lymphaticinvasion) <- c("nėra invazijos", "invazija į limfmazgius")


#filter for data of interest
lymph_test <- clin_df_joined_test[, colnames(clin_df_joined_test) %in% c(expression,
                                                                "lymphaticinvasion" )]
lymph_test <- lymph_test[!is.na(lymph_test$lymphaticinvasion), ]
#get long df
lymph_table_full <- reshape2::melt(lymph_test[, colnames(lymph_test) %in%
                                                 c("lymphaticinvasion", expression )],
                                   id.vars="lymphaticinvasion",
                                   measure.vars= expression)
#get t test
lymph_test_lymph <- lymph_table_full %>%
  group_by(variable) %>%
  t_test(value ~ lymphaticinvasion,
         p.adjust.method = "BH", 
         var.equal = TRUE, #stjudents
         paired = FALSE, 
         detailed=TRUE 
  ) # Format p-values to remove scientific notation
lymph_test_lymph
#make a tibble of p values: ~group1, ~group2, ~p.adj,   ~y.position, ~variable
lymph_test_lymph_tibble <- lymph_test_lymph %>% 
  dplyr::select(group1, group2, p, variable) %>%
  mutate(
    y.position = c(5, 5, 5, 7, 13, 
                   8, 5, 5, 5, 11) #choose where to plot p values
  )%>%
  filter(p < 0.05)
lymph_test_lymph_tibble

#get colors 
custom_colors_lymph <- c("invazija į limfmazgius" = "darkviolet",
                         "nėra invazijos" = "violet") 
#plot
lymph_plot <- ggplot(lymph_table_full, aes(x=lymphaticinvasion , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = lymphaticinvasion )) +
  geom_jitter(aes(color = lymphaticinvasion ), size=1, alpha=0.5) +
  ylab(label = expression("Genų raiška")) + 
  facet_wrap(.~ factor(variable), nrow = 2, scales = "free") +
  add_pvalue(lymph_test_lymph_tibble, label = "p") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text( size = 12, face = "bold.italic"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values = custom_colors_lymph) +
  scale_color_manual(values = custom_colors_lymph) +
  ggtitle("Genų raiška pagal su invaziją į limfmazgius testavimo imtyje")+
  coord_cartesian(clip = "off")  

#show
lymph_plot

lymph_plot2 <- ggdraw(lymph_plot) +
  draw_plot_label(label = "B", x = 0, y = 1, hjust = 0, vjust = 1, size = 20, fontface = "bold")

lymph_plot2

#save short
png("C:/Users/Ieva/rprojects/outputs_all/DISS/test_lymph_barplot_20260129.png",
    width = 15, height = 11, res = 300, units = "cm") 
lymph_plot2 
dev.off() 
#save lymphatic invasion plot
# png("C:/Users/Ieva/rprojects/outputs_all/DISS/test_lymph_barplot_20260121.png",
#     width = 32, height = 19, res = 500, units = "cm") # width and height in pixels, resolution in dpi
# lymph_plot2 #
# dev.off() # Close the PNG device

#TCGA - residual disease, barplot #######################
table(clin_df_joined_test$tumorresidualdisease, useNA ="a")
table(clin_df_joined_test$residualtumor, useNA ="a")
clin_df_joined_test$tumorresidualdisease <- factor(clin_df_joined_test$tumorresidualdisease)
levels(clin_df_joined_test$tumorresidualdisease) <- c(">20 mm", "1-10 mm ","11-20 mm",
                                                 "Be likutinio naviko")


#filter for data of interest
residual_test <- clin_df_joined_test[, colnames(clin_df_joined_test) %in% c(expression,
                                                                   "tumorresidualdisease" )]
residual_test <- residual_test[!is.na(residual_test$tumorresidualdisease), ]

#make residual tumor comaprisons############################################
table(residual_test$tumorresidualdisease, useNA = "a")
#levene test for variance, stage, tcga
levene_rez2 <- residual_test %>%
  pivot_longer(cols = -tumorresidualdisease , names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = car::leveneTest(value ~ tumorresidualdisease )$`Pr(>F)`[1],
    .groups = "drop"
  )
levene_rez2  #all normal, exept VPS3B

#make anova for vps33b
oneway.test(VPS33B ~ tumorresidualdisease, data = residual_test, var.equal = FALSE) #p sifnificant
#make Games-Howell posthoc test
residual_test %>% games_howell_test(VPS33B ~ tumorresidualdisease) #sifificant 3 groups

##anova, residual tumor################################
results2 <- lapply(names(residual_test)[c(2:11)], function(var) {
  formula <- as.formula(paste(var, "~ tumorresidualdisease"))
  fit <- aov(formula, data = residual_test)
  
  list(
    variable = var,
    summary = summary(fit),
    shapiro = shapiro.test(residuals(fit))
  )
})

names(results2) <- names(residual_test)[c(2:11)]
#make df of results2
results2_df <- lapply(results2, function(x) {
  data.frame(
    variable   = x$variable,
    anova_p    = x$summary[[1]][["Pr(>F)"]][1],
    shapiro_p  = x$shapiro$p.value
  )
}) %>%
  bind_rows()

results2_df
#show not normal
non_normal_df2 <- results2_df %>% filter(shapiro_p <= 0.05) 
non_normal_df2 #GRB7 ; zfpl1, exo1
#show significant
results2_df %>% filter(anova_p <= 0.05) #ppt2

#tuckey posthoc for ppt2
TukeyHSD(aov(PPT2 ~ tumorresidualdisease, data = residual_test)) #significant

##stage Kruskal test######################################################
kruskal_results2 <- lapply(non_normal_df2$variable, function(var) {
  df <- residual_test[, c(var, "tumorresidualdisease")]
  # df <- na.omit(df)  # remove NAs
  test <- kruskal.test(as.formula(paste(var, "~ tumorresidualdisease")), data = df)
  
  data.frame(
    variable = var,
    kruskal_p = test$p.value
  )
}) %>%
  bind_rows()

kruskal_results2# not significant

#add manually from anovas
t.test_res_tibble <-  tibble::tribble(
  ~group1, ~group2, ~p,   ~y.position, ~variable,
  "Be likutinio naviko",   ">20 mm", 0.0377, 5, "PPT2", #tukey
  "11-20 mm",   ">20 mm", 0.02, 6, "VPS33B", #Games-Howell
  "1-10 mm ",   "11-20 mm", 0.037 , 5, "VPS33B", #Games-Howell
  "Be likutinio naviko",   "11-20 mm", 0.014 , 4, "VPS33B", #Games-Howell
)

##plot residual##################
#get long df
residual_table_full <- reshape2::melt(residual_test[, colnames(residual_test) %in%
                                                       c("tumorresidualdisease", expression )],
                                      id.vars="tumorresidualdisease",
                                      measure.vars= expression)

#get colors 
custom_colors_res <- c(">20 mm" = "#E6E6FA",
                       "1-10 mm " = "#9966CC",
                       "11-20 mm" = "#7851A9",
                       "Be likutinio naviko" ="#673AB7") 
#plot
res_plot <- ggplot(residual_table_full, aes(x=tumorresidualdisease , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = tumorresidualdisease )) +
  geom_jitter(aes(color = tumorresidualdisease ), size=1, alpha=0.5) +
  ylab(label = expression("Genų raiška")) + 
  facet_wrap(.~ factor(variable), nrow = 2, scales = "free") +
  add_pvalue(t.test_res_tibble, label = "p") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7) )+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values = custom_colors_res) +
  scale_color_manual(values = custom_colors_res) +
  ggtitle("Genų raiška pagal likutinį naviką testavimo imtyje")+
  coord_cartesian(clip = "off")  
#show
res_plot

res_plot2 <- ggdraw(res_plot) +
  draw_plot_label(label = "B", x = 0, y = 1, hjust = 0, vjust = 1, size = 20, fontface = "bold")

res_plot2
#save short
png("C:/Users/Ieva/rprojects/outputs_all/DISS/test_res_barplot_20260129.png",
    width = 15, height = 12, res = 300, units = "cm") 
res_plot2 
dev.off() 
#save residual disease plot
# png("C:/Users/Ieva/rprojects/outputs_all/DISS/test_res_barplot_20260121.png",
#     width = 32, height = 19, res = 500, units = "cm") # width and height in pixels, resolution in dpi
# res_plot2 #
# dev.off() # Close the PNG device
#survival, separate genes ###########################################
colnames(clin_df_joined_test)
genes <-  c("EXO1", "RAD50", "PPT2", "LUC7L2", 
            "PKP3", "CDCA5", "ZFPL1", "VPS33B", 
            "GRB7", "TCEAL4")
genes_f <-  paste0(genes, "_f")
# Calculate the median expressions
clin_df_joined_test2 <- clin_df_joined_test %>%
  mutate(
    across(all_of(genes),
           ~ factor(if_else(. > median(., na.rm = TRUE), "High", "Low")),
           .names = "{.col}_f")
  )
colnames(clin_df_joined_test2)

##COX REGRESION, UNIVARIATE, using continuous variables######################
univ_results <- lapply(genes, function(gene) {
  formula <- as.formula(paste("Surv(overall_survival, deceased) ~", gene))
  cox_model <- coxph(formula, data = clin_df_joined_test2)
  sum_cox <- summary(cox_model)
  
  # Extract hazard ratio, 95% CI, and p-value
  if (!is.null(sum_cox$conf.int)) {
    data.frame(
      HR = sum_cox$conf.int[,"exp(coef)"],
      lower95 = sum_cox$conf.int[,"lower .95"],
      upper95 = sum_cox$conf.int[,"upper .95"],
      pvalue = sum_cox$coefficients[,"Pr(>|z|)"]
    )
  } else {
    data.frame(HR = NA, lower95 = NA, upper95 = NA, pvalue = NA)
  }
})

# Combine results
univ_df <- do.call(rbind, univ_results)
univ_df$Gene <- genes_f
univ_df <- univ_df[, c("Gene", "HR", "lower95", "upper95", "pvalue")]
print(univ_df)

##plot KM, separate ###################
plots_met <- list()

for (gene in genes_f) {
  # Clean gene name: remove "_f" suffix if present
  gene_clean <- gsub("_f$", "", gene)
  
  # Fit KM survival curve
  fit <- survfit(as.formula(paste("Surv(overall_survival, deceased) ~", gene)), data = clin_df_joined_test2)
  
  ## --- Extract univariable Cox HR + CI from univ_df ---
  hr_row_uni <- univ_df[univ_df$Gene == gene, ]
  if (nrow(hr_row_uni) == 0 || any(is.na(hr_row_uni$HR))) {
    hr_text_uni <- bquote("HR = NA")
  } else {
    hr_text_uni <- bquote("HR = " * .(sprintf("%.2f (95%% CI: %.2f–%.2f)", 
                                              hr_row_uni$HR, hr_row_uni$lower95, hr_row_uni$upper95)))
  }
  
  ## --- KM log-rank p-value ---
  survdiff_res <- survdiff(as.formula(paste("Surv(overall_survival, deceased) ~", gene)), 
                           data = clin_df_joined_test2)
  pval_km <- 1 - pchisq(survdiff_res$chisq, length(survdiff_res$n) - 1)
  if (pval_km < 0.05) {
    pval_expr <- bquote(bold("Log-rank p = " ~ .(sprintf("%.3f", pval_km))))
  } else {
    pval_expr <- bquote("Log-rank p = " ~ .(sprintf("%.3f", pval_km)))
  }
  
  ## --- Combine HR and log-rank p-value in subtitle ---
  subtitle_expr <- bquote(.(hr_text_uni) ~ "\n" ~ .(pval_expr))
  
  ## --- Legend labels ---
  legend_labels <- c("Low", "High")
  
  ## --- KM plot ---
  p <- ggsurvplot(
    fit,
    data = clin_df_joined_test2,
    pval = FALSE,
    risk.table = TRUE,
    legend.title = "",
    palette = c("blue", "red")
  )$plot +
    scale_color_manual(
      values = c("blue", "red"),
      labels = legend_labels
    ) +
    labs(
      title = bquote(italic(.(gene_clean))),
      subtitle = subtitle_expr
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 10, hjust = 0.5, lineheight = 1.1)
    )
  
  plots_met[[gene_clean]] <- p  # save using cleaned name
}

# Combine plots with overall title
combined_plot <- wrap_plots(plots_met, ncol = 5) +
  plot_annotation(
    title = "Survival analysis of test data",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

# Save to file
ggsave(
  filename = "C:/Users/Ieva/rprojects/outputs_all/DISS/KM_combined_TEST_data.png",
  plot = combined_plot,
  width = 25,
  height = 7,
  dpi = 300
)

#ENGLISH plot at year 5 with saving###########################
# Choose target time
target_time <- 1825   # choose year 5
time_index <- which(rez_list[[1]]$times == target_time)

# plot with save
png("C:/Users/Ieva/rprojects/outputs_all/DISS/10_genų_timeROC_test202501215EN.png",
    width = 1000, height = 1000, res = 200)

par(pty="s")

# Base plot with first gene
plot(
  rez_list[[1]]$FP[, time_index],
  rez_list[[1]]$TP[, time_index],
  type = "l",
  col = 1,
  lwd = 2,
  xlab = "1 - Specificity",
  ylab = "Sensitivity",
  main = paste("ROC curves,\n5 years from diagnosis, Test cohort"),
  xlim = c(0, 1),
  ylim = c(0, 1),
  asp = 1
)

# Add ROC lines for all genes
for (i in 2:length(rez_list)) {
  lines(
    rez_list[[i]]$FP[, time_index],
    rez_list[[i]]$TP[, time_index],
    col = i,
    lwd = 2
  )
}

# Add risk score ROC line in bold
lines(
  roc_result$FP[, time_index],
  roc_result$TP[, time_index],
  col = "maroon",
  lwd = 3,
  lty = 1
)

# Add diagonal reference line
abline(0, 1, lty = 2, col = "gray")

# Build legend labels without AUC
legend_labels <- c(paste0("italic('", names(rez_list), "')"),
                   "'Risk score'")

# Add legend
legend(
  "bottomright",
  legend = parse(text = legend_labels),
  col = c(1:length(rez_list), "maroon"),
  lwd = c(rep(2, length(rez_list)), 3),
  cex = 0.6,
  bty = "n"
)

# Add panel label
#mtext("B", side = 3, line = 2.5, adj = -0.2, font = 2, cex = 1.5)

dev.off()

#COORDS table #####################################################

# extract best sens/spec at a given time
extract_coords <- function(roc, gene, timepoint = 1825) {
  idx <- which(roc$times == timepoint)
  
  # sensitivity & specificity across thresholds
  sens <- as.vector(roc$TP[, idx])
  spec <- as.vector(1 - roc$FP[, idx])
  
  # Youden index to pick best threshold
  youden <- sens + spec - 1
  best_idx <- which.max(youden)
  
  tibble(
    gene = gene,
    time = roc$times[idx],
    auc = roc$AUC[idx] * 100,
    sens = sens[best_idx],
    spec = spec[best_idx],
    cutoff = roc$cutoffs[best_idx]
  )
}

# Apply to each biomarker
sens_spec_auc_60 <- map_dfr(names(rez_list), ~ extract_coords(rez_list[[.x]], .x, 1825))

# Apply to risk score
risk_score_row <- extract_coords(roc_result, "Risk score", 1825)

# Combine into one tibble
sens_spec_auc_60_all <- bind_rows(sens_spec_auc_60, risk_score_row)
# View
sens_spec_auc_60_all

#make prety table
# Prepare table: rename columns for display
roc_table_display <- sens_spec_auc_60_all %>%
  rename(
    Biomarker = gene,
    Time_Days = time,
    AUC = auc,
    Sensitivity = sens,
    Specificity = spec
  )

# Create gt table
gt_table_roc_60 <- roc_table_display %>%
  gt() %>%
  tab_header(
    title = "ROC metrics",
    subtitle = "Prognostic criteria for 5 year survival"
  ) %>%
  fmt_number(
    columns = vars(AUC, Sensitivity, Specificity),
    decimals = 3
  ) %>%
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = vars(Biomarker))
  )

# Show table
gt_table_roc_60

#there is no other convenient way to save gt outputs
gtsave(gt_table_roc_60,
       filename = "C:/Users/Ieva/rprojects/outputs_all/DISS/TEST_timeroc_table_20251023EN.png")

#Combine the images
roc_image <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/10_genų_timeROC_test202501215EN.png")
table_image <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/TEST_timeroc_table_20251023EN.png")

combined_image <- image_append(c(roc_image, table_image), stack = F)

# Save the combined image
image_write(combined_image, 
            "C:/Users/Ieva/rprojects/outputs_all/DISS/TESTROC_W_TABLE20251216EN.png")

#EN plot TCGA vs GTEx #######################
gtex_plotEN <- ggplot(gtcga_table_full, aes(x=group , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = group )) +
  geom_jitter(aes(color = group ), size=1, alpha=0.5) +
  ylab(label = expression("Gene expression")) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(t.test_gtex_tibble, label = "p_custom") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  ggtitle("Gene expression in test dataset")+
  coord_cartesian(clip = "off")  
#show plot
gtex_plotEN
#save gtex vs tca test plot
png("C:/Users/Ieva/rprojects/outputs_all/DISS/test_barplot_EN20260130.png",
    width = 15, height = 10, res = 300, units = "cm") # width and height in pixels, resolution in dpi
gtex_plotEN +
  theme(
    strip.text.x = element_text(size = 8, face = "bold.italic"),
    axis.text = element_text(size = 7),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 9, hjust = 0.5),
    plot.margin = margin(5, 5, 5, 5)
  )#
dev.off() # Close the PNG device


#EN KM PLOT TEST ####################
# Plot
train_survEN <- ggsurvplot(
  km_fit, 
  data = clin_df_joined_test, 
  pval = T,  # disable built-in
  risk.table = TRUE,
  #risk.table.title = "Pacienčių skaičius rizikos grupėje",
  title = "Survival analysis in test cohort",
  xlab = "Overall survival, days",
  #ylab = "Išgyvenamumo tikimybė",
  palette = c("turquoise", "deeppink"),
 # legend.title = "Rizikos grupė", 
  legend.labs = c("Low risk score", "High risk score")
)

#save km plot
png("C:/Users/Ieva/rprojects/outputs_all/DISS/dis_lt_km_test20251215EN.png",
    width = 800, height = 600, res = 120) # width and height in pixels, resolution in dpi
train_survEN #
dev.off() # Close the PNG device


#save boxplots together with the other cohort##############################
#do this after train plots are already generated 
img_top <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/train_barplot_20260121.png")
img_bottom <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/test_barplot_20260121.png")

combined_vertical <- image_append(
  c(img_top, img_bottom),
  stack = TRUE
)

image_write(
  combined_vertical,
  "C:/Users/Ieva/rprojects/outputs_all/DISS/TRAIN_TEST_boxplot_LT20260121.png"
)


#EN version 
#do this after train plots are already generated 
img_top <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/train_barplot_EN20260121.png")
img_bottom <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/test_barplot_EN20260121.png")

combined_vertical <- image_append(
  c(img_top, img_bottom),
  stack = TRUE
)

image_write(
  combined_vertical,
  "C:/Users/Ieva/rprojects/outputs_all/DISS/TRAIN_TEST_boxplot_EN20260121.png"
)

#SAVE FIGURES TOGETHER WITH TRAIN COHORTS#####################
#stage
#do this after train plots are already generated 
img_top <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/train_stage_barplot_20260129.png")
img_bottom <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/test_stage_barplot_20260129.png")

combined_vertical <- image_append(
  c(img_top, img_bottom),
  stack = TRUE
)

image_write(
  combined_vertical,
  "C:/Users/Ieva/rprojects/outputs_all/DISS/TRAIN_TEST_boxplot_STAGE20260129.png"
)
#grade
#do this after train plots are already generated 
img_top <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/train_grade_barplot_20260129.png")
img_bottom <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/test_grade_barplot_20260129.png")

combined_vertical <- image_append(
  c(img_top, img_bottom),
  stack = TRUE
)

image_write(
  combined_vertical,
  "C:/Users/Ieva/rprojects/outputs_all/DISS/TRAIN_TEST_boxplot_GRADE20260129.png"
)
#lymph nodes
#do this after train plots are already generated 
img_top <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/train_lymph_barplot_20260129.png")
img_bottom <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/test_lymph_barplot_20260129.png")

combined_vertical <- image_append(
  c(img_top, img_bottom),
  stack = TRUE
)

image_write(
  combined_vertical,
  "C:/Users/Ieva/rprojects/outputs_all/DISS/TRAIN_TEST_boxplot_LYMPH20260129.png"
)
#residual tumor
#do this after train plots are already generated 
img_top <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/train_res_barplot_20260129.png")
img_bottom <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/test_res_barplot_20260129.png")

combined_vertical <- image_append(
  c(img_top, img_bottom),
  stack = TRUE
)

image_write(
  combined_vertical,
  "C:/Users/Ieva/rprojects/outputs_all/DISS/TRAIN_TEST_boxplot_RESIDUAL_TUMOR20260129.png"
)

# save togerther gtex vs tcga###########

# Import images
img1 <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/train_barplot_20260127.png")
img2 <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/test_barplot_20260127.png")

# Combine vertically and save
image_append(c(img1, img2), stack = TRUE) %>%
  image_write("combined_train_test_barplot_20260130.png")

#EN
# Import images
img1 <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/train_barplot_EN20260129.png")
img2 <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/test_barplot_EN20260130.png")

# Combine vertically and save
image_append(c(img1, img2), stack = TRUE) %>%
  image_write("C:/Users/Ieva/rprojects/outputs_all/DISS/combined_train_test_barplot_EN20260130.png")

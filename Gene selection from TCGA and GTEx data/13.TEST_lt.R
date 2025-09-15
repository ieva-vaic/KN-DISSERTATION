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
  select(group1, group2, p, variable) %>%
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
  ggtitle("Genų raiška testavimo imtyje")
#show plot
gtex_plot
#save gtex vs tca test plot
png("C:/Users/Ieva/rprojects/outputs_all/DISS/test_barplot_20250618.png",
    width = 1500, height = 1200, res = 200) # width and height in pixels, resolution in dpi
gtex_plot #
dev.off() # Close the PNG device

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

##plot at year 5  with plotting###############
# Choose target time
target_time <- 1825        
time_index <- which(rez_list[[1]]$times == target_time)

png("C:/Users/Ieva/rprojects/outputs_all/DISS/10_genų_timeROC_test20250618.png",
    width = 1500, height = 1200, res = 200) # width and height in pixels, resolution in dpi
# Set up base plot with gene 1
par(pty="s")
# Set up base plot with gene 1
plot(
  rez_list[[1]]$FP[, time_index],
  rez_list[[1]]$TP[, time_index],
  type = "l",
  col = 1,
  lwd = 2,
  xlab = "Specifiškumas",
  ylab = "Jautrumas",
  main = paste("Nuo laiko priklausomos ROC kreivės,
5 metai po diagnozės testavimo imtyje"),
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

# Get AUCs for each gene at time_index
auc_list <- sapply(rez_list, function(x) x$AUC[time_index])
auc_risk <- roc_result$AUC[time_index]

# Build gene labels with italic names and AUCs
legend_labels <- mapply(function(name, auc) {
  paste0("italic('", name, "')~'(AUC = ", sprintf("%.3f", auc), ")'")
}, names(rez_list), auc_list)

# Add risk score with AUC
legend_labels <- c(legend_labels,
                   paste0("'Risk Score (AUC = ", sprintf("%.3f", auc_risk), ")'"))

# Add legend
legend(
  "bottomright",
  legend = parse(text = legend_labels),
  col = c(1:length(rez_list), "maroon"),
  lwd = c(rep(2, length(rez_list)), 3),
  cex = 0.6,
  bty = "n"
)
#add A label during png
mtext("D", side = 3, line = 2.5, adj = -0.2, font = 2, cex = 1.5)

#run plot
dev.off() # Close the PNG device
#KM plot with RISK SCORE#################
# Calculate the median risk score
median_risk <- median(clin_df_joined_test$RiskScore, na.rm = TRUE) #--0.07747393
# Create a new factor column based on the median value
clin_df_joined_test$RiskGroup <- ifelse(clin_df_joined_test$RiskScore <= median_risk,
                                   "Low Risk", "High Risk")
#Create a survival object
surv_object <- Surv(time = clin_df_joined_test$overall_survival,
                    event = clin_df_joined_test$deceased )
# Fit a Kaplan-Meier model
km_fit <- survfit(surv_object ~ RiskGroup, data = clin_df_joined_test)
# Plot the Kaplan-Meier curve using ggsurvplot
test_survplot <- ggsurvplot(km_fit, data = clin_df_joined_test, 
           pval = TRUE,  # Show p-value of the log-rank test
           risk.table = TRUE,  # Add risk table below the plot
           title = "Kaplan-Meier kreivė: Didelės vs. mažos rizikos atvejai testavimo imtyje",
           xlab = "Bendras išgyvenamumo laikas",
           ylab = "Išgyvenamumo tikimybė",
           palette = c("turquoise", "deeppink"),  # Color palette for groups
           legend.title = "Rizikos grupė", 
           legend.labs = c("Žema rizika", "Didelė rizika"))

png("C:/Users/Ieva/rprojects/outputs_all/DISS/dis_lt_km_test.png",
    width = 800, height = 600, res = 100) # width and height in pixels, resolution in dpi
test_survplot #
dev.off() # Close the PNG device

#Forest plot################
colnames(clin_df_joined_test)
clin_df_joined_test$neoplasmhistologicgrade
#first, fix stage
clin_df_joined_test <- clin_df_joined_test %>%
  mutate(stage_early_late = case_when(
    clinicalstage2 %in% c("Stage I", "Stage II") ~ "early stage",
    clinicalstage2 %in% c("Stage III", "Stage IV") ~ "late stage",
    TRUE ~ NA_character_
  ))
#first, fix grade
clin_df_joined_test$neoplasmhistologicgrade <- recode(clin_df_joined_test$neoplasmhistologicgrade,
                                                 "GX" = NA_character_)
vars <- c("RiskScore", "stage_early_late",
          "neoplasmhistologicgrade", "ageatinitialpathologicdiagnosis" )

#model
models <- vars %>%       # begin with variables of interest
  str_c("deceased ~ ", .) %>%   # combine each variable into formula ("outcome ~ variable of interest")
  # iterate through each univariate formula
  map(                               
    .f = ~glm(                       # pass the formulas one-by-one to glm()
      formula = as.formula(.x),      # within glm(), the string formula is .x
      family = "binomial",           # specify type of glm (logistic)
      data = clin_df_joined_test)) %>%          # dataset
  
  # tidy up each of the glm regression outputs from above
  map(
    .f = ~tidy(
      .x, 
      exponentiate = TRUE,           # exponentiate 
      conf.int = TRUE)) %>%          # return confidence intervals
  
  # collapse the list of regression outputs in to one data frame
  bind_rows() %>% 
  
  # round all numeric columns
  mutate(across(where(is.numeric), round, digits = 2)) %>%
  filter(term != "(Intercept)")

models #simple table models

# fix terms
models$term <-  c("Risk Score", "FIGO Stage late vs early",
                  "Neoplasm histologic grade G3 vs G2", "Age at initial pathologic diagnosis")
term_order <-  c("Risk Score", "FIGO Stage late vs early",
                 "Neoplasm histologic grade G3 vs G2", "Age at initial pathologic diagnosis")
# Apply the order to the term variable in the models data frame
models$term <- factor(models$term, levels = term_order)
models

#forest plot
p <- 
  models |>
  ggplot(aes(y = fct_rev(term))) + 
  theme_classic()+
  geom_point(aes(x=estimate), shape=15, size=3) +
  geom_linerange(aes(xmin=conf.low, xmax=conf.high)) +
  geom_vline(xintercept = 1, linetype="dashed") +
  labs(x="Odds Ratio, 95% CI", y="")+
  theme(
    axis.text.y = element_text(face = "italic")  # Make y-axis labels italic
  )+
  scale_x_continuous(
    breaks = c(1, seq(-3,  10, by = 1))  # Increase the number of x-axis ticks
  )


odds1 <- p +geom_text(aes(y = term, x = 2, label = sprintf("%0.2f", round(estimate, digits = 2))), size = 4) + ## decimal places
  annotate("text", x = 8, y= 4.5,  label = "OR", size = 4, fontface = "bold", hjust = 0)+
  geom_text(aes(y = term, x =6, label = paste(conf.low, "-", conf.high)), size = 4) + 
  annotate("text", x =5, y= 4.5,  label = "95 % CI", size =4, fontface = "bold", hjust = 0)+
  geom_text(aes(y = term, x =4, label = sprintf("%0.2f", round(p.value, digits = 2))), size = 4) + 
  annotate("text", x = 7, y= 4.5,  label = "p value", size = 4, fontface = "bold", hjust = 0)+
  labs(title = "Hazard Ratios for Mortality")
odds1

#TCGA - stage, boxplot ###########################
#CREATE GROUPINGS ACCORDING TO DATA#
table(clin_df_joined_test$clinicalstage2, useNA ="a")
table(clin_df_joined_test$stage_early_late, useNA ="a")
#filter for data of interest
colnames(clin_df_joined_test)
stage_test <- clin_df_joined_test[, colnames(clin_df_joined_test) %in% c(expression,
                                                                         "clinicalstage2", "barcode" )]
stage_test <- stage_test[!is.na(stage_test$clinicalstage2), ]
#get long df
stage_table_full <- reshape2::melt(stage_test[, colnames(stage_test) %in%
                                                c("clinicalstage2", expression )],
                                   id.vars="clinicalstage2",
                                   measure.vars= expression)
#get t test
t.test_stage <- stage_table_full %>%
  group_by(variable) %>%
  t_test(value ~ clinicalstage2,
         p.adjust.method = "BH", 
         var.equal = TRUE, #stjudents
         paired = FALSE, 
         detailed=TRUE 
  ) # Format p-values to remove scientific notation
t.test_stage
#make a tibble of p values: ~group1, ~group2, ~p.adj,   ~y.position, ~variable
t.test_stage_tibble <- t.test_stage %>% 
  select(group1, group2, p, variable) %>%
  mutate(
    y.position = c(6, 8, 8, 8.5, 13, 
                   8, 6, 6, 8, 11,
                   6, 8, 8, 8.5, 13, 
                   8, 6, 6, 8, 11,
                   6, 8, 8, 8.5, 13, 
                   8, 6, 6, 8, 11) #choose where to plot p values
  )%>%
  filter(p < 0.05)
t.test_stage_tibble

#get colors 
custom_colors_stage <- c("Stage II" = "lightgreen","Stage III" = "green", "Stage IV"= "darkgreen") 
#plot
stgae_plot <- ggplot(stage_table_full, aes(x=clinicalstage2 , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = clinicalstage2 )) +
  geom_jitter(aes(color = clinicalstage2 ), size=1, alpha=0.5) +
  ylab(label = expression("Genų raiška")) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(t.test_stage_tibble, label = "p") + #pvalue
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
  ggtitle("Genų raiškos koreliacija su stadija testavimo imtyje")

stgae_plot
#save stage boxplot 
png("C:/Users/Ieva/rprojects/outputs_all/DISS/test_stage_barplot_20250618.png",
    width = 2200, height = 1300, res = 200) # width and height in pixels, resolution in dpi
stgae_plot #
dev.off() # Close the PNG device

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
  select(group1, group2, p, variable) %>%
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
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
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
  ggtitle("Genų raiškos koreliacija su naviko diferenciacijos laipsniu testavimo imtyje")

grade_plot
#save grade boxplot
png("C:/Users/Ieva/rprojects/outputs_all/DISS/test_grade_barplot_20250618.png",
    width = 2200, height = 1300, res = 200) # width and height in pixels, resolution in dpi
grade_plot #
dev.off() # Close the PNG device

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
  select(group1, group2, p, variable) %>%
  mutate(
    y.position = c(6, 6, 6, 7, 13, 
                   8, 6, 6, 6, 11) #choose where to plot p values
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
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(lymph_test_lymph_tibble, label = "p") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values = custom_colors_lymph) +
  scale_color_manual(values = custom_colors_lymph) +
  ggtitle("Genų raiškos koreliacija su invazija į limfmazgius testavimo imtyje")
#show
lymph_plot
#save lymphatic invasion plot
png("C:/Users/Ieva/rprojects/outputs_all/DISS/test_lymph_barplot_20250619.png",
    width = 2450, height = 1300, res = 200) # width and height in pixels, resolution in dpi
lymph_plot #
dev.off() # Close the PNG device

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
#get long df
residual_table_full <- reshape2::melt(residual_test[, colnames(residual_test) %in%
                                                       c("tumorresidualdisease", expression )],
                                      id.vars="tumorresidualdisease",
                                      measure.vars= expression)
#get t test
t.test_res <- residual_table_full %>%
  group_by(variable) %>%
  t_test(value ~ tumorresidualdisease,
         p.adjust.method = "BH", 
         var.equal = TRUE, #stjudents
         paired = FALSE, 
         detailed=TRUE 
  ) # Format p-values to remove scientific notation
t.test_res
#make a tibble of p values: ~group1, ~group2, ~p.adj,   ~y.position, ~variable
t.test_res_tibble <- t.test_res %>% 
  select(group1, group2, p, variable) %>%
  mutate(
    y.position = c(9, 7, 7, 4, 5, 
                   8, 8, 7, 9, 8,
                   
                   9, 7, 7, 4, 8, 
                   8, 8, 7, 9, 8,
                   
                   9, 7, 7, 4, 8, 
                   8, 8, 7, 9, 8,
                   
                   9, 7, 7, 4, 8, 
                   8, 8, 7, 9, 8,
                   
                   9, 7, 7, 4, 8, 
                   8, 8, 7, 9, 8,
                   
                   9, 7, 7, 4, 8, 
                   8, 8, 7, 9, 8) #choose where to plot p values
  )%>%
  filter(p < 0.05)
t.test_res_tibble

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
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(t.test_res_tibble, label = "p") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1) )+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values = custom_colors_res) +
  scale_color_manual(values = custom_colors_res) +
  ggtitle("Genų raiškos koreliacija su invazija į limfmazgius testavimo imtyje")
#show
res_plot
#save residual disease plot
png("C:/Users/Ieva/rprojects/outputs_all/DISS/test_res_barplot_20250619.png",
    width = 2450, height = 1300, res = 170) # width and height in pixels, resolution in dpi
res_plot #
dev.off() # Close the PNG device
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


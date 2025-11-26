#KN-DISSERTATION project. Genes selection via statistical analysis of TCGA and GTEx data 
#(only scripts producing lithuanian figures, included in the dissertation)
#THIS IS TCGA-OV-RISK-GENES project script No. 12
#Make risk score, compare to clinical features
#ROCs, boxplots with TCGA, train data
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
library(rstatix) 
library(patchwork)
library(magick)
library(gt)
library(multcomp)
library(FSA)
library(cowplot)
#set directory of the data
setwd("../TCGA-OV-RISK-PROJECT/Public data RDSs/")
# Load train data ###################################
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot_2025.RDS")
#filter for lasso genes
gtex_genes <- readRDS("gtcga_elastic_2025.RDS")
gtex_filtered_counts_train <- gtex_counts_train[colnames(gtex_counts_train) %in% gtex_genes] 
#filter for TCGA cases
gtex_filtered_counts_train2 <- gtex_filtered_counts_train %>%
  dplyr::filter(grepl("^TCGA", rownames(gtex_filtered_counts_train)))
dim(gtex_filtered_counts_train2)
#gtex_filtered_counts_train2 <- as.data.frame(t(gtex_filtered_counts_train2))
#get genes of interest
expression <- c( "EXO1",   "RAD50",  "PPT2",   "LUC7L2", "PKP3",
                 "CDCA5",  "ZFPL1" , "VPS33B", "GRB7",   "TCEAL4")
# Load clinical data ###################################
clin_df <- readRDS("joinedTCGA_XENA_clinical2025.RDS")
#fix survival
clin_df$deceased = clin_df$vital_status == "Dead"
clin_df$overall_survival = ifelse(clin_df$deceased,
                                  clin_df$days_to_death,
                                  clin_df$days_to_last_follow_up)
#Filter clinical data###################################
train_ids <- rownames(gtex_filtered_counts_train2)
clin_df <- clin_df[clin_df$barcode %in% train_ids, ]  #336 samples
rownames(clin_df) <- clin_df$barcode
#Join clinical (mostly survival) with train data ########################
#add all genes to clin_df
colnames(clin_df)
gtex_filtered_counts_train2$barcode <- rownames(gtex_filtered_counts_train2)
clin_df$barcode == gtex_filtered_counts_train2$barcode 
clin_df_joined <- left_join(clin_df, gtex_filtered_counts_train2, by = "barcode")
rownames(clin_df_joined) <- clin_df_joined$barcode
colnames(clin_df_joined)
#saveRDS(clin_df_joined, "clin_df_joined_2025.RDS")
# Load cox model ###################################
cvfit <- readRDS("coxnet_cvfit_2025.RDS")
cox_fitx <- readRDS("coxnet_fit_2025.RDS")

#get coeficients - gene names
coef_x <- coef(cox_fitx, s = 0.088)
head(coef_x)
coef_x = coef_x[coef_x[,1] != 0,] 
res_coef_cox_names = names(coef_x) # get names of the (non-zero) variables.
res_coef_cox_names #10

#Value of coeficients plot ##################################
cv_model <- cvfit
# Extract coefficients at optimal lambda
coefs <- coef_x
coefs_df <- as.data.frame(as.matrix(coefs))
coefs_df <- rownames_to_column(coefs_df, var = "Feature")
colnames(coefs_df)[2] <- "Coefficient"
print(coefs_df)
# Plot coeficients ###########################
# get colors
n_features <- nrow(coefs_df)
colors <- colorRampPalette(brewer.pal(9, "Set3"))(n_features)  # Light pastel rainbow colors
# Pretty ggplot2 
plot_coef <- ggplot(coefs_df, aes(x = reorder(Feature, Coefficient), y = Coefficient, fill = Feature)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.8) +  # Colorful bars
  scale_fill_manual(values = colors) +  # Apply colors to bars
  coord_flip() +
  theme_minimal(base_size = 10) +
  labs( x = "Genai", y = "Koeficientų reikšmės") +
  theme(
    axis.text.y = element_text(face = "italic"),  # Italic feature names
    panel.grid.major.y = element_blank(),  # Remove gridlines for cleaner look
    panel.grid.minor = element_blank(),
    legend.position = "none"  # Hide legend for a clean visual
  )
#save value of coeficients plot
png("C:/Users/Ieva/rprojects/outputs_all/DISS/coeficient_plot_bar20250618.png",
    width = 900, height = 800, res = 200) 
plot_coef 
dev.off() 

#RISK SCORE ##################################
# - cv_model: Your trained cv.glmnet model
cv_model
# - gene_data: A data frame (or matrix) with the expression values for your genes 
# (rows = samples, columns = genes)
# It should have the same feature names as used in the model (genes).
gene_data <- clin_df_joined[, res_coef_cox_names]
# - gene_list: A vector with the list of genes you're interested in.
res_coef_cox_names
# coefs
coefs_df
# Calculate the risk score: linear combination of gene expressions and coefficients
# Risk score = sum( gene_expression * coefficient )
risk_scores <- rowSums(sweep(gene_data, 2, coefs_df$Coefficient, "*"))
# View the risk scores
print(risk_scores) # now I have some risk scores

#add risk scores to the clin_df_joined
clin_df_joined$RiskScore <- risk_scores[rownames(clin_df_joined)]

#survival time vs coeficient plot ##################################
clin_df_joined$overall_survival
clin_df_joined$RiskScore
clin_df_joined$vital_status

#Kaplan-meier plot RISK SCORE ##################################
# Calculate the median risk score
median_risk <- median(clin_df_joined$RiskScore, na.rm = TRUE) #-0.05914624
# Create a new factor column based on the median value
clin_df_joined$RiskGroup <- ifelse(clin_df_joined$RiskScore <= median_risk, "Low Risk", "High Risk")
#Create a survival object
surv_object <- Surv(time = clin_df_joined$overall_survival, event = clin_df_joined$deceased )
# Fit a Kaplan-Meier model
km_fit <- survfit(surv_object ~ RiskGroup, data = clin_df_joined)
# Compute log-rank test
surv_test <- survdiff(Surv(overall_survival, deceased) ~ RiskGroup, data = clin_df_joined)
pval_num <- 1 - pchisq(surv_test$chisq, length(surv_test$n) - 1)
# Format nicely
if (pval_num < 0.001) {
  pval_text <- "Long-rank p < 0.001"
} else {
  pval_text <- paste0("Long-rank p = ", signif(pval_num, 3))
}

# Plot
train_surv <- ggsurvplot(
  km_fit, 
  data = clin_df_joined, 
  pval = FALSE,  # disable built-in
  risk.table = TRUE,
  risk.table.title = "Pacienčių skaičius rizikos grupėje",
  title = "Didelės vs. mažos rizikos atvejai mokymosi imtyje",
  xlab = "Bendras išgyvenamumo laikas",
  ylab = "Išgyvenamumo tikimybė",
  palette = c("turquoise", "deeppink"),
  legend.title = "Rizikos grupė", 
  legend.labs = c("Mažas rizikos balas", "Didelis rizikos balas")
)

# Add subtitle
train_surv$plot <- train_surv$plot + labs(subtitle = pval_text)
print(train_surv)
#save km plot
png("C:/Users/Ieva/rprojects/outputs_all/DISS/dis_lt_km_train20250925.png",
    width = 800, height = 600, res = 100) # width and height in pixels, resolution in dpi
train_surv #
dev.off() # Close the PNG device

#SURVIVAL ROCs for separate markers ##################################
surv_df <- clin_df_joined[, colnames(clin_df_joined) %in%
                          c("deceased", "overall_survival", res_coef_cox_names)]
surv_df <- surv_df %>% dplyr::rename(censor = deceased, surv_time = overall_survival) 
#need these for survrock
nobs <- NROW(surv_df)
cutoff <- 365
#make surf df but only of my genes!
coxnet.df <- surv_df[, (colnames(surv_df) %in% res_coef_cox_names)]
dim(coxnet.df) #334  10 #kazkur pametu 2 zmones, kol kas neieskosiu - turbut tie kurie be deceased
time <- surv_df$surv_time
event <- surv_df$censor
rez_list <- apply(coxnet.df, 2, survivalROC, Stime = time, status = event, predict.time = cutoff, method="KM")
plot_list <- list()
# Loop through rez_list and create individual plots
for (i in seq_along(rez_list)) {
  p <- ggplot() +
    geom_line(aes(x = rez_list[[i]]$FP, y = rez_list[[i]]$TP)) +
    xlim(0, 1) + ylim(0, 1) +
    xlab(paste("FP", "\n", "AUC = ", round(rez_list[[i]]$AUC, 3))) +
    ylab("TP") +
    ggtitle(paste(names(rez_list)[i], ", Method = KM, Year = 1")) +
    theme_minimal() +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted")  # Add diagonal line
  
  # Store the plot in the plot_list
  plot_list[[i]] <- p
}
# Now, arrange all the plots in a grid and save the combined plot
combined_plot <- grid.arrange(grobs = plot_list, ncol = 2) 
#save separate plots
ggsave("C:/Users/Ieva/rprojects/outputs_all/DISS/SURVROC_SEPARATE20250410.png",
       combined_plot, width = 10, height = 15, dpi = 300)
#SURVIVAL ROC for my risk score ##################################
surv_df2 <- clin_df_joined[, colnames(clin_df_joined) %in%
                            c("deceased", "overall_survival", "RiskScore")]
surv_df2 <- surv_df2 %>% dplyr::rename(censor = deceased, surv_time = overall_survival) 
#need these for survrock
nobs <- NROW(surv_df2)
cutoff <- 365
#make surf df but only of my genes!
time2 <- surv_df2$surv_time
event2 <- surv_df2$censor
roc_result <- survivalROC(Stime = time2, 
                          status = event2, 
                          marker = surv_df2[, "RiskScore"], 
                          predict.time = cutoff, 
                          method = "KM")
# Plot the ROC curve for this single marker
p <- ggplot() +
  geom_line(aes(x = roc_result$FP, y = roc_result$TP)) +
  xlim(0, 1) + ylim(0, 1) +
  xlab(paste("FP", "\n", "AUC = ", round(roc_result$AUC, 3))) +
  ylab("TP") +
  ggtitle(paste("Marker: ", "RiskScore", ", Method = KM, Year = 1")) +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted")  # Add diagonal line

# Print the plot
print(p)

#TIME ROC for the risk score ###############################################
t_eval <- c(365, 1095, 1825)  # time points
roc_result <- timeROC(
  T = time2,       # Survival time from df
  delta = event2, # Event indicator from df
  marker = surv_df2[, "RiskScore"], # Predictor or risk score from df
  cause = 1,         # Event of interest
  times = t_eval,    # Time points for ROC
  iid = TRUE         # Compute confidence intervals
)

roc_result #time roc does basically the same as survival roc, but give se 
plot(roc_result, time = 1825) 
plot(roc_result, time = 1095) 
plot(roc_result, time = 365) 

#make prettier plot 
# Generate base ROC curve without default labels
plot(roc_result$FP[, "t=1825"], roc_result$TP[, "t=1825"], 
     type = "l", col = "blue", lwd = 2,
     xlab = "Specifiškumas", ylab = "Jautrumas",
     main = "10 genų raiškos modelio sąsaja su išgyvenamumu po 5 metų",
     xlim = c(0, 1), ylim = c(0, 1))

# Add a diagonal reference line
abline(a = 0, b = 1, lty = 2, col = "gray")

# Add AUC text
auc_val <- roc_result$AUC["t=1825"]
text(x = 0.6, y = 0.2, labels = paste0("AUC = ", round(auc_val, 3)), cex = 1.2)

#TIME ROC PLOTS for separate biomarkers ######################################
rez_list2 <- apply(coxnet.df, 2, timeROC,
                   T = time2,       # Survival time from df
                   delta = event2, # Event indicator from df
                   #marker  # Predictor already in the df
                   cause = 1,         # Event of interest
                   times = t_eval,    # Time points for ROC
                   iid = TRUE )        # Compute confidence intervals)

auc_table <- map_dfr(names(rez_list2), function(gene) {
  roc <- rez_list2[[gene]]
  
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
time_index <- which(rez_list2[[1]]$times == target_time)

# Set up base plot with gene 1
plot(
  rez_list2[[1]]$FP[, time_index],
  rez_list2[[1]]$TP[, time_index],
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
for (i in 2:length(rez_list2)) {
  lines(
    rez_list2[[i]]$FP[, time_index],
    rez_list2[[i]]$TP[, time_index],
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
  parse(text = paste0("italic('", names(rez_list2), "')")),
  "Risk Score"
)

# Add legend (last color is black for risk score)
legend(
  "bottomright",
  legend = legend_labels,
  col = c(1:length(rez_list2), "maroon"),
  lwd = c(rep(2, length(rez_list2)), 3),
  cex = 0.6,
  bty = "n"
)
##plot at year 3 ##########################
# Choose target time
target_time <- 1095      #choose year 3
time_index <- which(rez_list2[[1]]$times == target_time)

# Set up base plot with gene 1
plot(
  rez_list2[[1]]$FP[, time_index],
  rez_list2[[1]]$TP[, time_index],
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
for (i in 2:length(rez_list2)) {
  lines(
    rez_list2[[i]]$FP[, time_index],
    rez_list2[[i]]$TP[, time_index],
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

# Build legend names: italic gene names + "risk score"
legend_labels <- c(
  parse(text = paste0("italic('", names(rez_list2), "')")),
  "Risk Score"
)

# Add legend (last color for risk score)
legend(
  "bottomright",
  legend = legend_labels,
  col = c(1:length(rez_list2), "maroon"),
  lwd = c(rep(2, length(rez_list2)), 3),
  cex = 0.6,
  bty = "n"
)

## plot at year 5 with saving###########################
# Choose target time
target_time <- 1825   # choose year 5
time_index <- which(rez_list2[[1]]$times == target_time)

# plot with save
png("C:/Users/Ieva/rprojects/outputs_all/DISS/10_genų_timeROC_train20251124.png",
    width = 1000, height = 1000, res = 200)

par(pty="s")

# Base plot with first gene
plot(
  rez_list2[[1]]$FP[, time_index],
  rez_list2[[1]]$TP[, time_index],
  type = "l",
  col = 1,
  lwd = 2,
  xlab = "1 - Specifiškumas",
  ylab = "Jautrumas",
  main = paste("Nuo laiko priklausomos ROC kreivės,\n5 metai po diagnozės mokymosi imtyje"),
  xlim = c(0, 1),
  ylim = c(0, 1),
  asp = 1
)

# Add ROC lines for all genes
for (i in 2:length(rez_list2)) {
  lines(
    rez_list2[[i]]$FP[, time_index],
    rez_list2[[i]]$TP[, time_index],
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
legend_labels <- c(paste0("italic('", names(rez_list2), "')"),
                   "'Rizikos balas'")

# Add legend
legend(
  "bottomright",
  legend = parse(text = legend_labels),
  col = c(1:length(rez_list2), "maroon"),
  lwd = c(rep(2, length(rez_list2)), 3),
  cex = 0.6,
  bty = "n"
)

# Add panel label
mtext("C", side = 3, line = 2.5, adj = -0.2, font = 2, cex = 1.5)

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
    auc = roc$AUC[idx],
    sens = sens[best_idx],
    spec = spec[best_idx],
    cutoff = roc$cutoffs[best_idx]
  )
}

# Apply to each biomarker
sens_spec_auc_60 <- map_dfr(names(rez_list2), ~ extract_coords(rez_list2[[.x]], .x, 1825))

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
       filename = "C:/Users/Ieva/rprojects/outputs_all/DISS/TRAIN_timeroc_table_202501002.png")

#Combine the images
roc_image <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/10_genų_timeROC_train20251124.png")
table_image <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/TRAIN_timeroc_table_202501002.png")

combined_image <- image_append(c(roc_image, table_image), stack = F)

# Save the combined image
image_write(combined_image, 
            "C:/Users/Ieva/rprojects/outputs_all/DISS/TRAIN_ROC_W_TABLE20251124.png")

#Kaplan-meier plot, sepratae genes##################################
colnames(clin_df_joined) 
genes <-  c("EXO1", "RAD50", "PPT2", "LUC7L2", 
            "PKP3", "CDCA5", "ZFPL1", "VPS33B", 
            "GRB7", "TCEAL4")
genes_f <-  paste0(genes, "_f")
# Calculate the median expressions
clin_df_joined2 <- clin_df_joined %>%
  mutate(
    across(all_of(genes),
           ~ factor(if_else(. > median(., na.rm = TRUE), "High", "Low")),
           .names = "{.col}_f")
  )
colnames(clin_df_joined2)

##COX REGRESION, UNIVARIATE, using continuous variables######################
univ_results <- lapply(genes, function(gene) {
  formula <- as.formula(paste("Surv(overall_survival, deceased) ~", gene))
  cox_model <- coxph(formula, data = clin_df_joined2)
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
  fit <- survfit(as.formula(paste("Surv(overall_survival, deceased) ~", gene)), data = clin_df_joined2)
  
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
                           data = clin_df_joined2)
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
    data = clin_df_joined2,
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
    title = "Survival analysis of train data",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

# Save to file
ggsave(
  filename = "C:/Users/Ieva/rprojects/outputs_all/DISS/KM_combined_TRAIN_data.png",
  plot = combined_plot,
  width = 25,
  height = 7,
  dpi = 300
)

#With stage risk model ###########################
#choose variables
colnames(clin_df_joined)
table(clin_df_joined$clinicalstage2, useNA = "a")
clin_df_joined$RiskScore

#get anova
kruskal_test <- kruskal.test(clin_df_joined$RiskScore,
                             clin_df_joined$clinicalstage2, 
                             var.equal = F,
                             #alternative = "two.sided", 
                             na.rm = TRUE)
kruskal.testp_value <- kruskal_test$p.value
kruskal.testp_value
#get colors 
custom_colors <- c("Stage I" = "pink2", "Stage II" = "lightpink","Stage III" = "deeppink",
                   "Stage IV" = "darkviolet") 
#plot
stage_plot <- ggplot(clin_df_joined, aes(x=clinicalstage2 , y=RiskScore, fill = clinicalstage2)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = clinicalstage2 )) +
  geom_jitter(aes(color = clinicalstage2 ), size=1, alpha=0.5) +
  ylab(label = expression("Risk score")) + 
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
  ggtitle("Risk score correlation with stage ")+
  annotate("text", x = 2, y = max(clin_df_joined$RiskScore) + 0.11, 
           label = paste("p = ", format(kruskal.testp_value, digits = 3)), 
           size = 5, color = "black")

stage_plot

#With grade risk model###########################
table(clin_df_joined$neoplasmhistologicgrade, useNA = "a")
clin_df_joined$neoplasmhistologicgrade <- recode(clin_df_joined$neoplasmhistologicgrade,
                                                 "GB" = NA_character_,
                                                 "GX" = NA_character_,
                                                 "G4" = NA_character_)
# Remove NA values from both columns using complete.cases()
clin_df_joined_gr <- clin_df_joined[complete.cases(clin_df_joined$RiskScore,
                                                   clin_df_joined$neoplasmhistologicgrade), ]
clin_df_joined_gr$neoplasmhistologicgrade <- as.factor(clin_df_joined_gr$neoplasmhistologicgrade)
#fix stage while at it
clin_df_joined_gr <- clin_df_joined_gr %>%
  mutate(stage_early_late = case_when(
    clinicalstage2 %in% c("Stage I", "Stage II") ~ "early stage",
    clinicalstage2 %in% c("Stage III", "Stage IV") ~ "late stage",
    TRUE ~ NA_character_
  ))

ttest_grade <- t.test(clin_df_joined_gr$RiskScore ~ clin_df_joined_gr$neoplasmhistologicgrade ,  
                      var.equal = F,
                      alternative = "two.sided")
ttest_grade_p <- ttest_grade$p.value

#get colors 
custom_colors <- c("G2" = "lightpink","G3" = "deeppink") 
#plot
grade_plot <- ggplot(clin_df_joined, aes(x=neoplasmhistologicgrade ,
                                         y=RiskScore, fill = neoplasmhistologicgrade)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = neoplasmhistologicgrade )) +
  geom_jitter(aes(color = neoplasmhistologicgrade ), size=1, alpha=0.5) +
  ylab(label = expression("Risk score")) + 
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
  ggtitle("Risk score correlation with grade ")+
  annotate("text", x = 2, y = max(clin_df_joined$RiskScore) + 0.11, 
           label = paste("p = ", format(ttest_grade_p, digits = 3)), 
           size = 5, color = "black")

grade_plot

#GTEX vs TCGA, boxplot ###########################
#CREATE GROUPINGS ACCORDING TO DATA#
snames = rownames(gtex_counts_train)
group = substr(snames, 1, 4)
group = as.factor(group)
levels(group) <- c("GTEx", "TCGA-OV")
gtex_counts_train2 <- gtex_counts_train
gtex_counts_train2$group <- group
expression <- res_coef_cox_names
#get long df
gtcga_table_full <- reshape2::melt(gtex_counts_train2[, colnames(gtex_counts_train2) %in%
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
    y.position = c(6, 8, 8, 9, 10, 
                   8, 6, 6, 12, 12) #choose where to plot p values
  )
t.test_gtex_tibble$p_custom <- ifelse(t.test_gtex_tibble$p < 0.001, 
                                      "p < 0.001", 
                                      paste0("p = ", sprintf("%.3f",
                                                             each.vs.ref_sig$pj)))
#get colors 
custom_colors <- c("GTEx" = "turquoise","TCGA-OV" = "deeppink") 
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
  ggtitle("Genų raiška mokymosi imtyje")

gtex_plot
#save
png("C:/Users/Ieva/rprojects/outputs_all/DISS/train_barplot_20250618.png",
    width = 1500, height = 1200, res = 200) # width and height in pixels, resolution in dpi
gtex_plot #
dev.off() # Close the PNG device

#GET FC
gtex_counts_train_my <- gtex_counts_train2[, colnames(gtex_counts_train2) %in%
                                             c("group", expression )]
#drop group as comumn
gtex_counts_train_my$group
group <- gtex_counts_train_my$group
gtex_counts_train_my <- gtex_counts_train_my[, !(names(gtex_counts_train_my) %in% "group")] 
#get group means
group1_mean <- colMeans(gtex_counts_train_my[group == "GTEx", ])
group2_mean <- colMeans(gtex_counts_train_my[group == "TCGA-OV", ])
#get logFC  and FC
log2FC <- group2_mean - group1_mean
FC <- 2^log2FC
#show

fc_table_train <- data.frame(
  Gene = colnames(gtex_counts_train_my),
  log2FC = log2FC,
  FoldChange = FC
)

print(fc_table_train)
#TCGA - stage, boxplot #############################
#CREATE GROUPINGS ACCORDING TO DATA#
table(clin_df_joined$clinicalstage2, useNA ="a")

#filter for data of interest
colnames(clin_df_joined)
stage_train <- clin_df_joined[, colnames(clin_df_joined) %in% c(expression,
                                            "clinicalstage2", "barcode" )]
stage_train <- stage_train[!is.na(stage_train$clinicalstage2), ]

#make stage comaprisons############################################
#levene test for variance, stage, tcga
levene_rez <- stage_train[, c(2:12)] %>%
  pivot_longer(cols = -clinicalstage2 , names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = car::leveneTest(value ~ clinicalstage2 )$`Pr(>F)`[1],
    .groups = "drop"
  )
levene_rez  #all normal

##anova, stage################################
results <- lapply(names(stage_train)[c(3:12)], function(var) {
  formula <- as.formula(paste(var, "~ clinicalstage2"))
  fit <- aov(formula, data = stage_train)
  
  list(
    variable = var,
    summary = summary(fit),
    shapiro = shapiro.test(residuals(fit))
  )
})

names(results) <- names(stage_train)[c(3:12)]
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
non_normal_df <- results_df %>% filter(shapiro_p <= 0.05) #TCEAL4, LUC7L2, PKP3
non_normal_df
#show significant
results_df %>% filter(anova_p <= 0.05) #EXO1 

##stage Kruskal test######################################################
kruskal_results <- lapply(non_normal_df$variable, function(var) {
  df <- stage_train[, c(var, "clinicalstage2")]
  # df <- na.omit(df)  # remove NAs
  test <- kruskal.test(as.formula(paste(var, "~ clinicalstage2")), data = df)
  
  data.frame(
    variable = var,
    kruskal_p = test$p.value
  )
}) %>%
  bind_rows()

kruskal_results# EXO1 is only significant

##dunns tests (posthoc for kurskal)###################################
# Run Dunn post-hoc for all non-normal variables
dunn_results <- lapply(non_normal_df$variable, function(var) {
  df <- stage_train[, c(var, "clinicalstage2")]
  df <- na.omit(df)
  df$clinicalstage2 <- as.factor(df$clinicalstage2)
  
  test <- dunnTest(as.formula(paste(var, "~ clinicalstage2")), data = df, method = "holm")
  
  data.frame(
    variable   = var,
    comparison = test$res$Comparison,  # <- use this column instead of rownames
    Z          = test$res$Z,
    p_value    = test$res$P.adj
  )
}) %>%
  bind_rows()

# View results
dunn_results %>%
  filter(p_value < 0.1) #exo1 II vs IV 0.07

#get long df for plot
stage_table_train_full <- reshape2::melt(stage_train[, colnames(stage_train) %in%
                                                c("clinicalstage2", expression )],
                                   id.vars="clinicalstage2",
                                   measure.vars= expression)

stage_table_train_full <- stage_table_train_full %>%
  mutate(clinicalstage2 = case_when(
    clinicalstage2 == "Stage I" ~ "I stadija",
    clinicalstage2 == "Stage II" ~ "II stadija",
    clinicalstage2 == "Stage III" ~ "III stadija",
    clinicalstage2 == "Stage IV" ~ "IV stadija",
    TRUE ~ NA_character_
  ))

#add manually from anovas
t.train_stage_tibble <-  tibble::tribble(
    ~group1, ~group2, ~p,   ~y.position, ~variable,
    "II stadija",   "IV stadija", 0.0661, 6, "EXO1", #dunn
  )

#get colors 
custom_colors_stage <- c("I stadija" = "turquoise","II stadija" = "lightblue",
                         "III stadija" = "blue", "IV stadija"= "darkblue") 
##plot stage ###########
stage_plot <- ggplot(stage_table_train_full, aes(x=clinicalstage2 , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = clinicalstage2 )) +
  geom_jitter(aes(color = clinicalstage2 ), size=1, alpha=0.5) +
  ylab(label = expression("Genų raiška")) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(t.train_stage_tibble, label = "p") + #pvalue
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
  ggtitle("Genų raiškos koreliacija su stadija mokymosi imtyje")

stage_plot2 <- ggdraw(stage_plot) +
  draw_plot_label(label = "A", x = 0, y = 1, hjust = 0, vjust = 1, size = 20, fontface = "bold")

stage_plot2
#save stage plot
png("C:/Users/Ieva/rprojects/outputs_all/DISS/train_stage_barplot_20250925.png",
    width = 2200, height = 1300, res = 180) # width and height in pixels, resolution in dpi
stage_plot2 #
dev.off() # Close the PNG device

#TCGA - grade, boxplot ###########################
#CREATE GROUPINGS ACCORDING TO DATA#
table(clin_df_joined$neoplasmhistologicgrade, useNA ="a")

clin_df_joined$grade <- clin_df_joined$neoplasmhistologicgrade 
clin_df_joined <- clin_df_joined %>%
  mutate(grade = case_when(
    grade %in% c("GB", "GX", "G4") ~ NA_character_,
    TRUE ~ grade
  ))

#filter for data of interest
colnames(clin_df_joined)
grade_train <- clin_df_joined[, colnames(clin_df_joined) %in% c(expression,
                                "grade" )]
grade_train <- grade_train[!is.na(grade_train$grade), ]
#get long df
grade_table_full <- reshape2::melt(grade_train[, colnames(grade_train) %in%
                                                c("grade", expression )],
                                   id.vars="grade",
                                   measure.vars= expression)
#get t test
t.train_grade <- grade_table_full %>%
  group_by(variable) %>%
  t_test(value ~ grade,
         p.adjust.method = "BH", 
         var.equal = TRUE, #stjudents
         paired = FALSE, 
         detailed=TRUE 
  ) # Format p-values to remove scientific notation
t.train_grade
#make a tibble of p values: ~group1, ~group2, ~p.adj,   ~y.position, ~variable
t.train_grade_tibble <- t.train_grade %>% 
  dplyr::select(group1, group2, p, variable) %>%
  mutate(
    y.position = c(4, 4, 4, 8.5, 13, 
                   8, 8, 4,8, 8) #choose where to plot p values
  )%>%
  filter(p < 0.05)
t.train_grade_tibble

#get colors 
custom_colors_stage <- c("G1" = "lightblue", "G2" = "blue", "G3"= "darkblue") 
#plot
grade_plot <- ggplot(grade_table_full, aes(x=grade , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = grade )) +
  geom_jitter(aes(color = grade ), size=1, alpha=0.5) +
  ylab(label = expression("Genų raiška")) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(t.train_grade_tibble, label = "p") + #pvalue
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
  ggtitle("Genų raiškos koreliacija su naviko diferenciacijos laipsniu mokymosi imtyje")

grade_plot

grade_plot2 <- ggdraw(grade_plot) +
  draw_plot_label(label = "B", x = 0, y = 1, hjust = 0, vjust = 1, size = 20, fontface = "bold")

#save grade plot
png("C:/Users/Ieva/rprojects/outputs_all/DISS/train_grade_barplot_20251001.png",
    width = 2200, height = 1300, res = 200) 
grade_plot2 
dev.off() 

#Coreliation with age#################
table(clin_df_joined$age_at_diagnosis, useNA ="a")

age_train <- clin_df_joined[, colnames(clin_df_joined) %in% c(expression,
                                                              "age_at_diagnosis")]

stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "bottom", size = 4)
plot_list <- lapply(expression, function(gene) {
  ggplot(age_train, aes_string(x = "age_at_diagnosis", y = gene)) +
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

# Show them individually
for (p in plot_list) print(p)

combined_plot <- wrap_plots(plot_list, ncol = 2) +
  plot_annotation(title = "Ryšys tarp amžiaus ir genų raiškos mokymosi imtyje")

# Save to PNG
ggsave(
  filename = "C:/Users/Ieva/rprojects/outputs_all/DISS/age_vs_expressiontrain20250619.png",
  plot = combined_plot,
  width = 8,       # adjust width as needed
  height = 11,       # adjust height as needed
  dpi = 300         # high resolution
)

#TCGA - lymphovascular invasion, boxplot ##############

table(clin_df_joined$lymphaticinvasion, useNA ="a")
clin_df_joined$lymphaticinvasion <- factor(clin_df_joined$lymphaticinvasion)
levels(clin_df_joined$lymphaticinvasion) <- c("nėra invazijos", "invazija į limfmazgius")


#filter for data of interest
lymph_train <- clin_df_joined[, colnames(clin_df_joined) %in% c(expression,
                                                                "lymphaticinvasion" )]
lymph_train <- lymph_train[!is.na(lymph_train$lymphaticinvasion), ]
#get long df
lymph_table_full <- reshape2::melt(lymph_train[, colnames(lymph_train) %in%
                                                 c("lymphaticinvasion", expression )],
                                   id.vars="lymphaticinvasion",
                                   measure.vars= expression)
#get t test
t.train_lymph <- lymph_table_full %>%
  group_by(variable) %>%
  t_test(value ~ lymphaticinvasion,
         p.adjust.method = "BH", 
         var.equal = TRUE, #stjudents
         paired = FALSE, 
         detailed=TRUE 
  ) # Format p-values to remove scientific notation
t.train_lymph
#make a tibble of p values: ~group1, ~group2, ~p.adj,   ~y.position, ~variable
t.train_lymph_tibble <- t.train_lymph %>% 
  dplyr::select(group1, group2, p, variable) %>%
  mutate(
    y.position = c(6, 6, 6, 7, 13, 
                   8, 6, 6, 6, 11) #choose where to plot p values
  )%>%
  filter(p < 0.05)
t.train_lymph_tibble

#get colors 
custom_colors_lymph <- c("invazija į limfmazgius" = "hotpink",
                         "nėra invazijos" = "lightpink") 
#plot
lymph_plot <- ggplot(lymph_table_full, aes(x=lymphaticinvasion , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = lymphaticinvasion )) +
  geom_jitter(aes(color = lymphaticinvasion ), size=1, alpha=0.5) +
  ylab(label = expression("Genų raiška")) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(t.train_lymph_tibble, label = "p") + #pvalue
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
  ggtitle("Genų raiškos koreliacija su invazija į limfmazgius mokymosi imtyje")

lymph_plot

lymph_plot2 <- ggdraw(lymph_plot) +
  draw_plot_label(label = "A", x = 0, y = 1, hjust = 0, vjust = 1, size = 20, fontface = "bold")

lymph_plot2
#save lymphovascular invasion plot
png("C:/Users/Ieva/rprojects/outputs_all/DISS/train_lymph_barplot_20251001.png",
    width = 2450, height = 1300, res = 200) # width and height in pixels, resolution in dpi
lymph_plot2 #
dev.off() # Close the PNG device

#TCGA - resudual tumor, boxplot ##############
table(clin_df_joined$tumorresidualdisease, useNA ="a")
table(clin_df_joined$residualtumor, useNA ="a")
clin_df_joined$tumorresidualdisease <- factor(clin_df_joined$tumorresidualdisease)
levels(clin_df_joined$tumorresidualdisease) <- c(">20 mm", "1-10 mm ","11-20 mm",
                                                 "Be likutinio naviko")


#filter for data of interest
residual_train <- clin_df_joined[, colnames(clin_df_joined) %in% c(expression,
                                                                "tumorresidualdisease" )]
residual_train <- residual_train[!is.na(residual_train$tumorresidualdisease), ]

#make residual tumor comaprisons############################################
table(residual_train$tumorresidualdisease, useNA = "a")
#levene test for variance, stage, tcga
levene_rez2 <- residual_train %>%
  pivot_longer(cols = -tumorresidualdisease , names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = car::leveneTest(value ~ tumorresidualdisease )$`Pr(>F)`[1],
    .groups = "drop"
  )
levene_rez2  #all normal


##anova, residual tumor################################
results2 <- lapply(names(residual_train)[c(2:11)], function(var) {
  formula <- as.formula(paste(var, "~ tumorresidualdisease"))
  fit <- aov(formula, data = residual_train)
  
  list(
    variable = var,
    summary = summary(fit),
    shapiro = shapiro.test(residuals(fit))
  )
})

names(results2) <- names(residual_train)[c(2:11)]
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
non_normal_df2
#show significant
results2_df %>% filter(anova_p <= 0.05) #GRB7 

##stage Kruskal test######################################################
kruskal_results2 <- lapply(non_normal_df2$variable, function(var) {
  df <- residual_train[, c(var, "tumorresidualdisease")]
  # df <- na.omit(df)  # remove NAs
  test <- kruskal.test(as.formula(paste(var, "~ tumorresidualdisease")), data = df)
  
  data.frame(
    variable = var,
    kruskal_p = test$p.value
  )
}) %>%
  bind_rows()

kruskal_results2# GRB7  is only significant

##dunns tests (posthoc for kurskal)###################################
# Run Dunn post-hoc for all non-normal variables
dunn_results2 <- lapply(non_normal_df2$variable, function(var) {
  df <- residual_train[, c(var, "tumorresidualdisease")]
  df <- na.omit(df)
  df$tumorresidualdisease <- as.factor(df$tumorresidualdisease)
  
  test <- dunnTest(as.formula(paste(var, "~ tumorresidualdisease")), data = df, method = "holm")
  
  data.frame(
    variable   = var,
    comparison = test$res$Comparison,  # <- use this column instead of rownames
    Z          = test$res$Z,
    p_value    = test$res$P.adj
  )
}) %>%
  bind_rows()

# View results
dunn_results2 %>%
  filter(p_value < 0.1) #grb7 

#add manually from anovas
t.train_res_tibble <-  tibble::tribble(
  ~group1, ~group2, ~p,   ~y.position, ~variable,
  "Be likutinio naviko",   ">20 mm", 0.0091, 9, "GRB7", #dunn
)

#get long df
residual_table_full <- reshape2::melt(residual_train[, colnames(residual_train) %in%
                                                 c("tumorresidualdisease", expression )],
                                   id.vars="tumorresidualdisease",
                                   measure.vars= expression)
#get colors 
custom_colors_res <- c(">20 mm" = "#D8A7B1",
                       "1-10 mm " = "#FA8072",
                       "11-20 mm" = "#FF66CC",
                       "Be likutinio naviko" ="#F6CFCB") 
##plot residual tumor ###########
res_plot <- ggplot(residual_table_full, aes(x=tumorresidualdisease , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = tumorresidualdisease )) +
  geom_jitter(aes(color = tumorresidualdisease ), size=1, alpha=0.5) +
  ylab(label = expression("Genų raiška")) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(t.train_res_tibble, label = "p") + #pvalue
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
  ggtitle("Genų raiškos koreliacija su likutiniu naviku mokymosi imtyje")

res_plot

res_plot2 <- ggdraw(res_plot) +
  draw_plot_label(label = "B", x = 0, y = 1, hjust = 0, vjust = 1, size = 20, fontface = "bold")

res_plot2
#save residual tumor plot
png("C:/Users/Ieva/rprojects/outputs_all/DISS/train_res_barplot_20251001.png",
    width = 2450, height = 1300, res = 170) # width and height in pixels, resolution in dpi
res_plot2 #
dev.off() # Close the PNG device


#ENGLISH plot at year 5 with saving###########################
# Choose target time
target_time <- 1825   # choose year 5
time_index <- which(rez_list2[[1]]$times == target_time)

# plot with save
png("C:/Users/Ieva/rprojects/outputs_all/DISS/10_genų_timeROC_train20251023EN.png",
    width = 1000, height = 1000, res = 200)

par(pty="s")

# Base plot with first gene
plot(
  rez_list2[[1]]$FP[, time_index],
  rez_list2[[1]]$TP[, time_index],
  type = "l",
  col = 1,
  lwd = 2,
  xlab = "Specificity",
  ylab = "Sensitivity",
  main = paste("ROC curves,\n 5 years from diagnosis, Train cohort"),
  xlim = c(0, 1),
  ylim = c(0, 1),
  asp = 1
)

# Add ROC lines for all genes
for (i in 2:length(rez_list2)) {
  lines(
    rez_list2[[i]]$FP[, time_index],
    rez_list2[[i]]$TP[, time_index],
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
legend_labels <- c(paste0("italic('", names(rez_list2), "')"),
                   "'Risk score'")

# Add legend
legend(
  "bottomright",
  legend = parse(text = legend_labels),
  col = c(1:length(rez_list2), "maroon"),
  lwd = c(rep(2, length(rez_list2)), 3),
  cex = 0.6,
  bty = "n"
)

# Add panel label
mtext("A", side = 3, line = 2.5, adj = -0.2, font = 2, cex = 1.5)

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
sens_spec_auc_60 <- map_dfr(names(rez_list2), ~ extract_coords(rez_list2[[.x]], .x, 1825))

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
    Time_days = time,
    AUC = auc,
    Sensitivity = sens,
    Specificity = spec
  )

# Create gt table
gt_table_roc_60 <- roc_table_display %>%
  gt() %>%
  tab_header(
    title = "ROC criteria",
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
       filename = "C:/Users/Ieva/rprojects/outputs_all/DISS/TRAIN_timeroc_table_20250925EN.png")

#Combine the images
roc_image <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/10_genų_timeROC_train20251023EN.png")
table_image <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/TRAIN_timeroc_table_20250925EN.png")

combined_image <- image_append(c(roc_image, table_image), stack = F)

# Save the combined image
image_write(combined_image, 
            "C:/Users/Ieva/rprojects/outputs_all/DISS/TRAIN_ROC_W_TABLE20251023EN.png")

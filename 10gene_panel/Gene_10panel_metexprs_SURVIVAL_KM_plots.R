#KN-DISSERTATION project. 10 gene panel - genes selected via statistical analysis of TCGA and GTEx data
#KN 10 gene analysis data (cv 2 %, one threshold, 65 cases!, restored KN-59) 2025-02-17
#KN-DISSERTATION project. MET_EXPRS - methylation, expression, mutation data
#survival analysis
Sys.setenv(LANG = "en")
#libraries
library(survival)
library(survminer)
library(tidyverse)
library(grid)
library(gridExtra)
library(patchwork)
library(magick)
library(survivalROC)
library(purrr)
library(broom)
library(dplyr)
library(timeROC)
library(cowplot)
#set wd for plots
setwd("C:/Users/Ieva/rprojects/outputs_all/DISS/")
#get one file for all data
ALL_EXPRESSION_DF <- readRDS("../../OTHER DATA/KN-DISSERTATION FILES/ALL_expresssion_Df20250709.RDS")
#get survival data as of 2025-09-11
SURVIVAL_KN <- openxlsx::read.xlsx("../../OTHER DATA/KN-DISSERTATION FILES/KN_MIRTIES_FAILAS_20250911.xlsx")
#make only surv df
SURV <- SURVIVAL_KN[, c(2, 3,20, 21)]
head(SURV)
#join with main data
ALL_SURV_EXPRESSION <- left_join(ALL_EXPRESSION_DF, SURV, by = "patient_id_aud")
#remove random columns or rename
ALL_SURV_EXPRESSION <-  ALL_SURV_EXPRESSION[, -c(1, 15, 32)]
ALL_SURV_EXPRESSION <- ALL_SURV_EXPRESSION %>%
  rename(tumor = tumor.y ,
         KN = KN.y )
rownames(ALL_SURV_EXPRESSION) <- ALL_SURV_EXPRESSION$patient_id_aud
#remove NA from new surv data
ALL_SURV_EXPRESSION <- ALL_SURV_EXPRESSION%>%
  filter(!is.na(OS), !is.na(STATUS)) #58 observations
#make levels of my gene expression data
methylation <- c("ARID1A_met", "CDX2", "ALX4", "HOPX")
genes <- c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4",
           "ARID1A", "CTNNB1", "FBXW7",
           "JAG2", "DLL1", "HES1",
           "EXO1", "RAD50", "PPT2", "LUC7L2", 
           "PKP3", "CDCA5", "ZFPL1", "VPS33B", 
           "GRB7", "TCEAL4")
genes_f <- paste0(genes, "_f")
genes10 <- c("EXO1", "RAD50", "PPT2", "LUC7L2", 
             "PKP3", "CDCA5", "ZFPL1", "VPS33B", 
             "GRB7", "TCEAL4")
genes_notch <- genes[!genes %in% genes10]  #notch genes
#in plots i use names with _f so add it
genes10_f <- paste0(genes10, "_f")
genes_notch_f <- paste0(genes_notch, "_f")

#make levels of gene expresssion data
ALL_SURV_EXPRESSION <- ALL_SURV_EXPRESSION %>%
  mutate(
    across(all_of(genes),
           ~ factor(if_else(. > median(., na.rm = TRUE), "High", "Low")),
           .names = "{.col}_f")
  )

#rename AGE (comes up later in multivariable cox)
ALL_SURV_EXPRESSION <- ALL_SURV_EXPRESSION %>%
  rename(Age = Amžius )
colnames(ALL_SURV_EXPRESSION)

#make factors
ALL_SURV_EXPRESSION$Stage2 <- factor(ALL_SURV_EXPRESSION$Stage2)
ALL_EXPRESSION_DF$Stage4 <- factor(ALL_EXPRESSION_DF$Stage4)
ALL_EXPRESSION_DF$Grade2 <- factor(ALL_EXPRESSION_DF$Grade2)
ALL_EXPRESSION_DF$CA125_f <- factor(ALL_EXPRESSION_DF$CA125_f)

#remove non-OC (benign) cases#############################
OC_SURV_EXPRESSION <- ALL_SURV_EXPRESSION %>%
  filter(Grupė_Ieva != "Benign") #55 left#

#HGSOC ONLY##########################
HGSOC_SURV_EXPRESSION <- ALL_SURV_EXPRESSION %>%
  filter(Grupė_Ieva == "HGSOC") #41 left#
table(HGSOC_SURV_EXPRESSION$STATUS, useNA = "a") #17 dead

#SURVIVAL DATA##############################################
#get descriptive statistics on STATUS
table(ALL_SURV_EXPRESSION$STATUS, useNA = "a") #19 was dead
table(ALL_SURV_EXPRESSION$tumor
      , ALL_SURV_EXPRESSION$STATUS, useNA = "a")
fisher.test(table(ALL_SURV_EXPRESSION$tumor
                  , ALL_SURV_EXPRESSION$STATUS, useNA = "a"))

#get descriptive statistics on OS
shapiro.test(ALL_SURV_EXPRESSION$OS) #not normal
median(ALL_SURV_EXPRESSION$OS)
min(ALL_SURV_EXPRESSION$OS)
max(ALL_SURV_EXPRESSION$OS)
ALL_SURV_EXPRESSION %>%
  filter(tumor != "Benign") %>%
  pull(OS) %>%
  median(na.rm = TRUE)
# t test
benign_os <- ALL_SURV_EXPRESSION %>%
  filter(tumor == "Benign") %>%
  pull(OS)
malignant_os <- ALL_SURV_EXPRESSION %>%
  filter(tumor == "OC") %>%
  pull(OS)
wilcox.test(benign_os, malignant_os)
#SURVIVAL OBJECT #####################################################
# Create a survival object
surv_obj <- Surv(time = ALL_SURV_EXPRESSION$OS, event = ALL_SURV_EXPRESSION$STATUS)

# Kaplan-Meier survival fit for the whole dataset
km_fit <- survfit(surv_obj ~ 1)

#SURVIVAL OF GENE EXPRESSIONS, ALL####################################################
#list of genes
##COX REGRESION, UNIVARIATE, using continuous variables######################
univ_results <- lapply(genes, function(gene) {
  formula <- as.formula(paste("Surv(OS, STATUS) ~", gene))
  cox_model <- coxph(formula, data = ALL_SURV_EXPRESSION)
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

#show significant 
sig_genes <- univ_df[univ_df$pvalue < 0.1, ]
print(sig_genes) #only HES1

##LONG RANK, KM PLOTS, with univariable COX HR and CI######################

plots <- list()

for (gene in genes_f) {
  # Clean up gene name (remove "_f")
  gene_clean <- sub("_f$", "", gene)
  
  # Fit KM survival curve
  fit <- survfit(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                 data = ALL_SURV_EXPRESSION)
  
  # Extract Cox HR + CI
  hr_row <- univ_df[univ_df$Gene == gene, ]
  hr_text <- sprintf("HR = %.2f (95%% CI: %.2f–%.2f)", 
                     hr_row$HR, hr_row$lower95, hr_row$upper95)
  
  # KM log-rank p-value
  survdiff_res <- survdiff(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                           data = ALL_SURV_EXPRESSION)
  pval_km <- 1 - pchisq(survdiff_res$chisq, length(survdiff_res$n) - 1)
  pval_text <- sprintf("Log-rank p = %.3f", pval_km)
  
  # Subtitle (HR + CI + p-value)
  subtitle_text <- paste0(hr_text, "\n", pval_text)
  
  # Make legend labels with italic gene name
  legend_labels <- c(
    bquote(italic(.(gene_clean)) ~ " Low"),
    bquote(italic(.(gene_clean)) ~ " High")
  )
  
  # Plot
  p <- ggsurvplot(
    fit,
    data = ALL_SURV_EXPRESSION,
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
      subtitle = subtitle_text
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 10, hjust = 0.5, lineheight = 1.1)
    )
  
  plots[[gene]] <- p
}

# Combine plots
combined_plot <- wrap_plots(plots, ncol = 4)
combined_plot
ggsave(
  filename = "KM_combined_plot_w_HR.png",  # output file name
  plot = combined_plot,               # the patchwork plot object
  width = 16,                         # width in inches
  height = 20,                        # height in inches
  dpi = 300                           # resolution
)

##separate the 10 gene expression plot form the notch and wnt ##########
# Subset plots for the two gene sets
plots_genes10 <- plots[names(plots) %in% genes10_f]
plots_genes_notch <- plots[names(plots) %in% genes_notch_f]

# Combine separately
combined_plot_genes10 <- wrap_plots(plots_genes10, ncol = 5)
combined_plot_genes_notch <- wrap_plots(plots_genes_notch, ncol = 5)

# Display
combined_plot_genes10
combined_plot_genes_notch
#save 10 gene
ggsave(
  filename = "KM_combined_plot_w_HR_10_gene_20250915.png",  # output file name
  plot = combined_plot_genes10,               # the patchwork plot object
  width = 20,                         # width in inches
  height = 10,                        # height in inches
  dpi = 300                           # resolution
)
#save notch
ggsave(
  filename = "KM_combined_plot_w_HR_notch_20250915.png",  # output file name
  plot = combined_plot_genes_notch,               # the patchwork plot object
  width = 20,                         # width in inches
  height = 10,                        # height in inches
  dpi = 300                           # resolution
)


##COX REGRESION, UNIVARIATE, only OC#################################
univ_results2 <- lapply(genes, function(gene) {
  formula <- as.formula(paste("Surv(OS, STATUS) ~", gene))
  cox_model <- coxph(formula, data = OC_SURV_EXPRESSION)
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
univ_df2 <- do.call(rbind, univ_results2)
univ_df2$Gene <- genes_f
univ_df2 <- univ_df2[, c("Gene", "HR", "lower95", "upper95", "pvalue")]
print(univ_df2)

##LONG RANK, KM PLOTS, only OC#############
# Create an empty list to store plots
plots2 <- list()

for (gene in genes_f) {
  # Clean up gene name (remove "_f")
  gene_clean <- sub("_f$", "", gene)
  
  # Fit KM survival curve
  fit <- survfit(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                 data = OC_SURV_EXPRESSION)
  
  # Extract Cox HR + CI
  hr_row <- univ_df2[univ_df2$Gene == gene, ]
  hr_text <- sprintf("HR = %.2f (95%% CI: %.2f–%.2f)", 
                     hr_row$HR, hr_row$lower95, hr_row$upper95)
  
  # KM log-rank p-value
  survdiff_res <- survdiff(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                           data = OC_SURV_EXPRESSION)
  pval_km <- 1 - pchisq(survdiff_res$chisq, length(survdiff_res$n) - 1)
  
  # Format p-value: bold if significant
  if (pval_km < 0.05) {
    pval_expr <- bquote(bold("Log-rank p = " ~ .(sprintf("%.3f", pval_km))))
  } else {
    pval_expr <- bquote("Log-rank p = " ~ .(sprintf("%.3f", pval_km)))
  }
  
  # Combine HR line and p-value line in subtitle
  subtitle_expr <- bquote(.(hr_text) ~ "\n" ~ .(pval_expr))
  
  # Make legend labels with italic gene name
  legend_labels <- c(
    bquote(italic(.(gene_clean)) ~ " Low"),
    bquote(italic(.(gene_clean)) ~ " High")
  )
  
  # Plot
  p <- ggsurvplot(
    fit,
    data = OC_SURV_EXPRESSION,
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
  
  plots2[[gene]] <- p
}

# Combine plots
combined_plot2 <- wrap_plots(plots2, ncol = 4)
combined_plot2

# Save to file
ggsave(
  filename = "KM_combined_plot_OC_only_w_HR.png",
  plot = combined_plot2,
  width = 16,
  height = 20,
  dpi = 300
)

##COX REGRESION, UNIVARIATE, HGSOC only ###################
univ_results3 <- lapply(genes_f, function(gene) {
  formula <- as.formula(paste("Surv(OS, STATUS) ~", gene))
  cox_model <- coxph(formula, data = HGSOC_SURV_EXPRESSION)
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
univ_df3 <- do.call(rbind, univ_results3)
univ_df3$Gene <- genes_f
univ_df3 <- univ_df3[, c("Gene", "HR", "lower95", "upper95", "pvalue")]
print(univ_df3)

##LONG RANK, KM PLOTS, only HGSOC#############
# Create an empty list to store plots
plots3 <- list()

for (gene in genes_f) {
  # Clean up gene name (remove "_f")
  gene_clean <- sub("_f$", "", gene)
  
  # Fit KM survival curve
  fit <- survfit(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                 data = HGSOC_SURV_EXPRESSION)
  
  # Extract Cox HR + CI
  hr_row <- univ_df3[univ_df3$Gene == gene, ]
  hr_text <- sprintf("HR = %.2f (95%% CI: %.2f–%.2f)", 
                     hr_row$HR, hr_row$lower95, hr_row$upper95)
  
  # KM log-rank p-value
  survdiff_res <- survdiff(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                           data = HGSOC_SURV_EXPRESSION)
  pval_km <- 1 - pchisq(survdiff_res$chisq, length(survdiff_res$n) - 1)
  
  # Format p-value: bold if significant
  if (pval_km < 0.05) {
    pval_expr <- bquote(bold("Log-rank p = " ~ .(sprintf("%.3f", pval_km))))
  } else {
    pval_expr <- bquote("Log-rank p = " ~ .(sprintf("%.3f", pval_km)))
  }
  
  # Combine HR line and p-value line in subtitle
  subtitle_expr <- bquote(.(hr_text) ~ "\n" ~ .(pval_expr))
  
  # Make legend labels with italic gene name
  legend_labels <- c(
    bquote(italic(.(gene_clean)) ~ " Low"),
    bquote(italic(.(gene_clean)) ~ " High")
  )
  
  # Plot
  p <- ggsurvplot(
    fit,
    data = HGSOC_SURV_EXPRESSION,
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
  
  plots3[[gene]] <- p
}

# Combine plots
combined_plot3 <- wrap_plots(plots3, ncol = 4)

# save
ggsave(
  filename = "KM_combined_plot_HGSOC_only_w_HR.png",
  plot = combined_plot3,
  width = 16,
  height = 20,
  dpi = 300
)

#MUlTIVARIABLE COX REGRESSION ########################################
table(ALL_SURV_EXPRESSION$Stage4, ALL_SURV_EXPRESSION$STATUS, useNA = "a") #3 na in stage have status 0
table(ALL_SURV_EXPRESSION$Grade2, ALL_SURV_EXPRESSION$STATUS, useNA = "a") #grade 11, has status
table(ALL_SURV_EXPRESSION$CA125_f, ALL_SURV_EXPRESSION$STATUS, useNA = "a") #8 na
hist(HGSOC_SURV_EXPRESSION$Age) 

##Multivariable cox, HGSOC only CA125 and AGE#############################
multi_results <- lapply(genes, function(gene) {
  
  # Select relevant columns: OS, STATUS, gene, Amžius, CA125
  df_sub <- HGSOC_SURV_EXPRESSION[, c("OS", "STATUS", gene, "Age", "CA125")]
  
  # Remove rows with missing values for this gene or covariates
  df_sub <- na.omit(df_sub)
  
  # If too few patients, return NA
  if(nrow(df_sub) < 2) return(data.frame(HR=NA, lower95=NA, upper95=NA, pvalue=NA, N=0))
  
  # Fit Cox
  formula <- as.formula(paste("Surv(OS, STATUS) ~", gene, "+ Age + CA125"))
  cox_model <- coxph(formula, data = df_sub)
  sum_cox <- summary(cox_model)
  
  # Extract gene effect
  gene_row <- which(rownames(sum_cox$coefficients) == gene)
  data.frame(
    HR = sum_cox$conf.int[gene_row, "exp(coef)"],
    lower95 = sum_cox$conf.int[gene_row, "lower .95"],
    upper95 = sum_cox$conf.int[gene_row, "upper .95"],
    pvalue = sum_cox$coefficients[gene_row, "Pr(>|z|)"],
    N = nrow(df_sub)
  )
})

multi_df <- do.call(rbind, multi_results)
multi_df$Gene <- genes
multi_df <- multi_df[, c("Gene", "HR", "lower95", "upper95", "pvalue", "N")]
multi_df$Gene <- paste0(multi_df$Gene, "_f") #fix for plot, does not mean I use factor variables
print(multi_df)

#plot with multi
plots3x <- list()

for (gene in genes_f) {
  gene_clean <- sub("_f$", "", gene)
  
  # Fit KM survival curve
  fit <- survfit(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                 data = HGSOC_SURV_EXPRESSION)
  
  ## --- Univariate HR ---
  hr_row <- univ_df3[univ_df3$Gene == gene, ]
  if (nrow(hr_row) == 0 || any(is.na(hr_row$HR))) {
    hr_text_expr <- bquote("Univ HR = NA")
  } else {
    hr_text_expr <- bquote("Univ HR = " * .(sprintf("%.2f (95%% CI: %.2f–%.2f)", 
                                                    hr_row$HR, hr_row$lower95, hr_row$upper95)))
  }
  
  ## --- KM log-rank p-value ---
  survdiff_res <- survdiff(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                           data = HGSOC_SURV_EXPRESSION)
  pval_km <- 1 - pchisq(survdiff_res$chisq, length(survdiff_res$n) - 1)
  if (pval_km < 0.05) {
    pval_expr <- bquote(bold("Log-rank p = " ~ .(sprintf("%.3f", pval_km))))
  } else {
    pval_expr <- bquote("Log-rank p = " ~ .(sprintf("%.3f", pval_km)))
  }
  
  ## --- Multivariable HR (Age + CA125) ---
  multi_row <- multi_df[multi_df$Gene == gene, ]
  if (nrow(multi_row) == 0 || any(is.na(multi_row$HR))) {
    hr_multi_expr <- bquote("Multiv HR = NA")
  } else {
    # color red if p < 0.05
    if (!is.na(multi_row$pvalue) && multi_row$pvalue < 0.05) {
      hr_multi_expr <- bquote(bold("Multiv HR = " * .(sprintf("%.2f (95%% CI: %.2f–%.2f, N=%d)", 
                                                              multi_row$HR, multi_row$lower95, 
                                                              multi_row$upper95, multi_row$N))))
    } else {
      hr_multi_expr <- bquote("Multiv HR = " * .(sprintf("%.2f (95%% CI: %.2f–%.2f, N=%d)", 
                                                         multi_row$HR, multi_row$lower95, 
                                                         multi_row$upper95, multi_row$N)))
    }
  }
  
  ## --- Combine 3 lines in subtitle: Univ HR | Log-rank p | Multiv HR ---
  subtitle_expr <- bquote(.(hr_text_expr) ~ "\n" ~ .(pval_expr) ~ "\n" ~ .(hr_multi_expr))
  
  ## --- Legend labels ---
  legend_labels <- c(
    bquote(italic(.(gene_clean)) ~ " Low"),
    bquote(italic(.(gene_clean)) ~ " High")
  )
  
  ## --- KM plot ---
  p <- ggsurvplot(
    fit,
    data = HGSOC_SURV_EXPRESSION,
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
  
  plots3x[[gene]] <- p
}

# Combine plots and add overall title
combined_plot3x <- wrap_plots(plots3x, ncol = 4) +
  plot_annotation(
    title = "Survival analysis of HGSOC",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )

# Save
ggsave(
  filename = "KM_combined_plot_HGSOC_only_w_HR_MULTI.png",
  plot = combined_plot3x,
  width = 30,
  height = 25,
  dpi = 300
)
###separate the 10 gene expression plot form the notch and wnt ##########
# Subset plots for the two gene sets
plots_genes10_3x <- plots3x[names(plots3x) %in% genes10_f]
plots_genes_notch_3x <- plots3x[names(plots3x) %in% genes_notch_f]

# Combine separately
combined_plot_genes10_3x <- wrap_plots(plots_genes10_3x, ncol = 5)+
  plot_annotation(
    title = "Survival analysis of HGSOC",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )

combined_plot_genes_notch_3x <- wrap_plots(plots_genes_notch_3x, ncol = 5)+
  plot_annotation(
    title = "Survival analysis of HGSOC",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )

#save 10 gene
ggsave(
  filename = "KM_combined_plot_w_HR_HGSOC_10_gene_20250915.png",  # output file name
  plot = combined_plot_genes10_3x,               # the patchwork plot object
  width = 35,                         # width in inches
  height = 10,                        # height in inches
  dpi = 300                           # resolution
)
#save notch
ggsave(
  filename = "KM_combined_plot_w_HR_HGSOC_notch_20250915.png",  # output file name
  plot = combined_plot_genes_notch_3x,               # the patchwork plot object
  width = 35,                         # width in inches
  height = 10,                        # height in inches
  dpi = 300                           # resolution
)

##Multivariable cox, OC only CA125 and AGE#############################
multi_results2 <- lapply(genes, function(gene) {
  
  # Select relevant columns: OS, STATUS, gene, Amžius, CA125
  df_sub <- OC_SURV_EXPRESSION[, c("OS", "STATUS", gene, "Age", "CA125")]
  
  # Remove rows with missing values for this gene or covariates
  df_sub <- na.omit(df_sub)
  
  # If too few patients, return NA
  if(nrow(df_sub) < 2) return(data.frame(HR=NA, lower95=NA, upper95=NA, pvalue=NA, N=0))
  
  # Fit Cox
  formula <- as.formula(paste("Surv(OS, STATUS) ~", gene, "+ Age + CA125"))
  cox_model <- coxph(formula, data = df_sub)
  sum_cox <- summary(cox_model)
  
  # Extract gene effect
  gene_row <- which(rownames(sum_cox$coefficients) == gene)
  data.frame(
    HR = sum_cox$conf.int[gene_row, "exp(coef)"],
    lower95 = sum_cox$conf.int[gene_row, "lower .95"],
    upper95 = sum_cox$conf.int[gene_row, "upper .95"],
    pvalue = sum_cox$coefficients[gene_row, "Pr(>|z|)"],
    N = nrow(df_sub)
  )
})

multi_df2 <- do.call(rbind, multi_results2)
multi_df2$Gene <- genes
multi_df2 <- multi_df2[, c("Gene", "HR", "lower95", "upper95", "pvalue", "N")]
multi_df2$Gene <- paste0(multi_df2$Gene, "_f") #fix for plot, does not mean I use factor variables
print(multi_df2)

#plot with multi
plots2x <- list()

for (gene in genes_f) {
  gene_clean <- sub("_f$", "", gene)
  
  # Fit KM survival curve
  fit <- survfit(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                 data = OC_SURV_EXPRESSION)
  
  ## --- Univariate HR ---
  hr_row <- univ_df2[univ_df2$Gene == gene, ]
  if (nrow(hr_row) == 0 || any(is.na(hr_row$HR))) {
    hr_text_expr <- bquote("Univ HR = NA")
  } else {
    hr_text_expr <- bquote("Univ HR = " * .(sprintf("%.2f (95%% CI: %.2f–%.2f)", 
                                                    hr_row$HR, hr_row$lower95, hr_row$upper95)))
  }
  
  ## --- KM log-rank p-value ---
  survdiff_res <- survdiff(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                           data = OC_SURV_EXPRESSION)
  pval_km <- 1 - pchisq(survdiff_res$chisq, length(survdiff_res$n) - 1)
  if (pval_km < 0.05) {
    pval_expr <- bquote(bold("Log-rank p = " ~ .(sprintf("%.3f", pval_km))))
  } else {
    pval_expr <- bquote("Log-rank p = " ~ .(sprintf("%.3f", pval_km)))
  }
  
  ## --- Multivariable HR (Age + CA125) ---
  multi_row <- multi_df2[multi_df2$Gene == gene, ]
  if (nrow(multi_row) == 0 || any(is.na(multi_row$HR))) {
    hr_multi_expr <- bquote("Multiv HR = NA")
  } else {
    # Bold if p < 0.05
    if (!is.na(multi_row$pvalue) && multi_row$pvalue < 0.05) {
      hr_multi_expr <- bquote(bold("Multiv HR = " * .(sprintf("%.2f (95%% CI: %.2f–%.2f, N=%d)", 
                                                              multi_row$HR, multi_row$lower95, 
                                                              multi_row$upper95, multi_row$N))))
    } else {
      hr_multi_expr <- bquote("Multiv HR = " * .(sprintf("%.2f (95%% CI: %.2f–%.2f, N=%d)", 
                                                         multi_row$HR, multi_row$lower95, 
                                                         multi_row$upper95, multi_row$N)))
    }
  }
  
  ## --- Combine 3 lines in subtitle ---
  subtitle_expr <- bquote(.(hr_text_expr) ~ "\n" ~ .(pval_expr) ~ "\n" ~ .(hr_multi_expr))
  
  ## --- Legend labels ---
  legend_labels <- c(
    bquote(italic(.(gene_clean)) ~ " Low"),
    bquote(italic(.(gene_clean)) ~ " High")
  )
  
  ## --- KM plot ---
  p <- ggsurvplot(
    fit,
    data = OC_SURV_EXPRESSION,
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
  
  plots2x[[gene]] <- p
}

# Combine plots and add overall title
combined_plotx2 <- wrap_plots(plots2x, ncol = 4) +
  plot_annotation(
    title = "Survival analysis of OC",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )

# Save to file
ggsave(
  filename = "KM_combined_plot_OC_only_w_HR_MULTI.png",
  plot = combined_plotx2,
  width = 30,
  height = 25,
  dpi = 300
)

###separate the 10 gene expression plot form the notch and wnt ##########
# Subset plots for the two gene sets
plots_genes10_2x <- plots2x[names(plots2x) %in% genes10_f]
plots_genes_notch_2x <- plots2x[names(plots2x) %in% genes_notch_f]

# Combine separately
combined_plot_genes10_2x <- wrap_plots(plots_genes10_2x, ncol = 5)+
  plot_annotation(
    title = "Survival analysis of OC",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )
combined_plot_genes_notch_2x <- wrap_plots(plots_genes_notch_2x, ncol = 5)+
  plot_annotation(
    title = "Survival analysis of OC",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )

#save 10 gene
ggsave(
  filename = "KM_combined_plot_w_HR_OC_10_gene_20250925.png",  # output file name
  plot = combined_plot_genes10_2x,               # the patchwork plot object
  width = 35,                         # width in inches
  height = 10,                        # height in inches
  dpi = 400                           # resolution
)
#save notch
ggsave(
  filename = "KM_combined_plot_w_HR_OC_notch_20250925.png",  # output file name
  plot = combined_plot_genes_notch_2x,               # the patchwork plot object
  width = 35,                         # width in inches
  height = 10,                        # height in inches
  dpi = 400                           # resolution
)


##Multivariable cox, ALL CASES,  CA125 and AGE#############################
multi_results3 <- lapply(genes, function(gene) {
  
  # Select relevant columns: OS, STATUS, gene, Amžius, CA125
  df_sub <- ALL_SURV_EXPRESSION[, c("OS", "STATUS", gene, "Age", "CA125")]
  
  # Remove rows with missing values for this gene or covariates
  df_sub <- na.omit(df_sub)
  
  # If too few patients, return NA
  if(nrow(df_sub) < 2) return(data.frame(HR=NA, lower95=NA, upper95=NA, pvalue=NA, N=0))
  
  # Fit Cox
  formula <- as.formula(paste("Surv(OS, STATUS) ~", gene, "+ Age + CA125"))
  cox_model <- coxph(formula, data = df_sub)
  sum_cox <- summary(cox_model)
  
  # Extract gene effect
  gene_row <- which(rownames(sum_cox$coefficients) == gene)
  data.frame(
    HR = sum_cox$conf.int[gene_row, "exp(coef)"],
    lower95 = sum_cox$conf.int[gene_row, "lower .95"],
    upper95 = sum_cox$conf.int[gene_row, "upper .95"],
    pvalue = sum_cox$coefficients[gene_row, "Pr(>|z|)"],
    N = nrow(df_sub)
  )
})

multi_df3 <- do.call(rbind, multi_results2)
multi_df3$Gene <- genes
multi_df3 <- multi_df3[, c("Gene", "HR", "lower95", "upper95", "pvalue", "N")]
multi_df3$Gene <- paste0(multi_df3$Gene, "_f") #fix for plot, does not mean I use factor variables
print(multi_df3)

#plot 
plots <- list()

for (gene in genes_f) {
  gene_clean <- sub("_f$", "", gene)
  
  # Fit KM survival curve
  fit <- survfit(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                 data = ALL_SURV_EXPRESSION)
  
  ## --- Univariate HR ---
  hr_row <- univ_df[univ_df$Gene == gene, ]
  if (nrow(hr_row) == 0 || any(is.na(hr_row$HR))) {
    hr_text_expr <- bquote("Univ HR = NA")
  } else {
    hr_text_expr <- bquote("Univ HR = " * .(sprintf("%.2f (95%% CI: %.2f–%.2f)", 
                                                    hr_row$HR, hr_row$lower95, hr_row$upper95)))
  }
  
  ## --- KM log-rank p-value ---
  survdiff_res <- survdiff(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                           data = ALL_SURV_EXPRESSION)
  pval_km <- 1 - pchisq(survdiff_res$chisq, length(survdiff_res$n) - 1)
  if (pval_km < 0.05) {
    pval_expr <- bquote(bold("Log-rank p = " ~ .(sprintf("%.3f", pval_km))))
  } else {
    pval_expr <- bquote("Log-rank p = " ~ .(sprintf("%.3f", pval_km)))
  }
  
  ## --- Multivariable HR (Age + CA125) from multi_df3 ---
  multi_row <- multi_df3[multi_df3$Gene == gene, ]
  if (nrow(multi_row) == 0 || any(is.na(multi_row$HR))) {
    hr_multi_expr <- bquote("Multiv HR = NA")
  } else {
    if (!is.na(multi_row$pvalue) && multi_row$pvalue < 0.05) {
      hr_multi_expr <- bquote(bold("Multiv HR = " * .(sprintf("%.2f (95%% CI: %.2f–%.2f, N=%d)", 
                                                              multi_row$HR, multi_row$lower95, 
                                                              multi_row$upper95, multi_row$N))))
    } else {
      hr_multi_expr <- bquote("Multiv HR = " * .(sprintf("%.2f (95%% CI: %.2f–%.2f, N=%d)", 
                                                         multi_row$HR, multi_row$lower95, 
                                                         multi_row$upper95, multi_row$N)))
    }
  }
  
  ## --- Combine 3 lines in subtitle ---
  subtitle_expr <- bquote(.(hr_text_expr) ~ "\n" ~ .(pval_expr) ~ "\n" ~ .(hr_multi_expr))
  
  ## --- Legend labels ---
  legend_labels <- c(
    bquote(italic(.(gene_clean)) ~ " Low"),
    bquote(italic(.(gene_clean)) ~ " High")
  )
  
  ## --- KM plot ---
  p <- ggsurvplot(
    fit,
    data = ALL_SURV_EXPRESSION,
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
  
  plots[[gene]] <- p
}

# Combine plots and add overall title
combined_plot <- wrap_plots(plots, ncol = 4) +
  plot_annotation(
    title = "Survival analysis of all ovarian samples",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )

# Save to file
ggsave(
  filename = "KM_combined_plot_all_OC_w_HR_MULTI.png",
  plot = combined_plot,
  width = 30,
  height = 25,
  dpi = 300
)
###separate the 10 gene expression plot form the notch and wnt ##########
# Subset plots for the two gene sets
plots_genes10_1x <- plots[names(plots) %in% genes10_f]
plots_genes_notch_1x <- plots[names(plots) %in% genes_notch_f]

# Combine separately
combined_plot_genes10_1x <- wrap_plots(plots_genes10_1x, ncol = 5)+
  plot_annotation(
    title = "Survival analysis of all ovarian samples",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )
combined_plot_genes_notch_1x <- wrap_plots(plots_genes_notch_1x, ncol = 5)+
  plot_annotation(
    title = "Survival analysis of all ovarian samples",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )

#save 10 gene
ggsave(
  filename = "KM_combined_plot_w_HR_ALL_10_gene_20250915.png",  # output file name
  plot = combined_plot_genes10_1x,               # the patchwork plot object
  width = 35,                         # width in inches
  height = 10,                        # height in inches
  dpi = 300                           # resolution
)
#save notch
ggsave(
  filename = "KM_combined_plot_w_HR_ALL_notch_20250915.png",  # output file name
  plot = combined_plot_genes_notch_1x,               # the patchwork plot object
  width = 35,                         # width in inches
  height = 10,                        # height in inches
  dpi = 300                           # resolution
)

#METHYLATION, all cases#################################
str(ALL_SURV_EXPRESSION[, colnames(ALL_SURV_EXPRESSION) %in% methylation])
#change methylation back to factor variable
ALL_SURV_EXPRESSION[, methylation] <- lapply(ALL_SURV_EXPRESSION[, methylation], function(x) {
  factor(x, levels = c(0, 1), labels = c("not methylated", "methylated"))
})

## METHYLATION COX REGRESION, UNIVARIATE######################
univ_results_met <- lapply(methylation, function(gene) {
  formula <- as.formula(paste("Surv(OS, STATUS) ~", gene))
  cox_model <- coxph(formula, data = ALL_SURV_EXPRESSION)
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
met_univ_df <- do.call(rbind, univ_results_met)
met_univ_df$Gene <- methylation
met_univ_df <- met_univ_df[, c("Gene", "HR", "lower95", "upper95", "pvalue")]
print(met_univ_df)

##METHYLATION Multivariable cox, all cases CA125 and AGE#############################
# Multivariable Cox regression including a factor variable
cox_model1 <- coxph(Surv(OS, STATUS) ~ ARID1A_met  + Age + CA125, data = ALL_SURV_EXPRESSION)
cox_model2 <- coxph(Surv(OS, STATUS) ~ CDX2  + Age + CA125, data = ALL_SURV_EXPRESSION)
cox_model3 <- coxph(Surv(OS, STATUS) ~ ALX4  + Age + CA125, data = ALL_SURV_EXPRESSION)
cox_model4 <- coxph(Surv(OS, STATUS) ~ HOPX  + Age + CA125, data = ALL_SURV_EXPRESSION)

summary(cox_model1)
summary(cox_model2)
summary(cox_model3)
summary(cox_model4)
# Run Cox models for all methylation genes
cox_results <- lapply(methylation, function(gene) {
  formula <- as.formula(paste0("Surv(OS, STATUS) ~ ", gene, " + Age + CA125"))
  model <- coxph(formula, data = ALL_SURV_EXPRESSION)
  s <- summary(model)
  
  # Find all rows that correspond to this gene
  idx <- grep(paste0("^", gene), rownames(s$coefficients))
  
  # Skip if gene not found (safeguard)
  if(length(idx) == 0) return(NULL)
  
  # Extract info for each row (use lapply in case factor has multiple levels)
  do.call(rbind, lapply(idx, function(i) {
    data.frame(
      Gene = rownames(s$coefficients)[i],
      HR = s$coefficients[i, "exp(coef)"],
      Lower95 = s$conf.int[i, "lower .95"],
      Upper95 = s$conf.int[i, "upper .95"],
      p.value = s$coefficients[i, "Pr(>|z|)"],
      n = model$n
    )
  }))
})

cox_results_df <- do.call(rbind, cox_results)
cox_results_df

#fix names for later
cox_results_df$Gene <- methylation
#plot with multivariable, methylation analysis
plots_met <- list()

for (gene in methylation) {
  # Clean gene name
  gene_clean <- gene
  
  # Fit KM survival curve
  fit <- survfit(as.formula(paste("Surv(OS, STATUS) ~", gene)), data = ALL_SURV_EXPRESSION)
  
  ## --- Extract univariable Cox HR + CI from met_univ_df ---
  hr_row_uni <- met_univ_df[met_univ_df$Gene == gene, ]
  if (nrow(hr_row_uni) == 0 || any(is.na(hr_row_uni$HR))) {
    hr_text_uni <- bquote("Uni HR = NA")
  } else {
    hr_text_uni <- bquote("Uni HR = " * .(sprintf("%.2f (95%% CI: %.2f–%.2f)", 
                                                  hr_row_uni$HR, hr_row_uni$lower95, hr_row_uni$upper95)))
  }
  
  ## --- Extract multivariable Cox HR + CI from cox_results_df ---
  hr_row_multi <- cox_results_df[cox_results_df$Gene == gene, ]
  if (nrow(hr_row_multi) == 0 || any(is.na(hr_row_multi$HR))) {
    hr_text_multi <- bquote("Multi HR = NA")
  } else {
    hr_text_multi <- bquote("Multi HR = " * .(sprintf("%.2f (95%% CI: %.2f–%.2f)", 
                                                      hr_row_multi$HR, hr_row_multi$Lower95, hr_row_multi$Upper95)))
  }
  
  ## --- KM log-rank p-value ---
  survdiff_res <- survdiff(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                           data = ALL_SURV_EXPRESSION)
  pval_km <- 1 - pchisq(survdiff_res$chisq, length(survdiff_res$n) - 1)
  if (pval_km < 0.05) {
    pval_expr <- bquote(bold("Log-rank p = " ~ .(sprintf("%.3f", pval_km))))
  } else {
    pval_expr <- bquote("Log-rank p = " ~ .(sprintf("%.3f", pval_km)))
  }
  
  ## --- Combine Uni HR, Multi HR, and log-rank p-value in subtitle ---
  subtitle_expr <- bquote(.(hr_text_uni) ~ "; " ~ .(hr_text_multi) ~ "\n" ~ .(pval_expr))
  
  ## --- Legend labels ---
  legend_labels <- c("not methylated", "methylated")
  
  ## --- KM plot ---
  p <- ggsurvplot(
    fit,
    data = ALL_SURV_EXPRESSION,
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
  
  plots_met[[gene]] <- p
}

# Combine plots with overall title
combined_met_plot <- wrap_plots(plots_met, ncol = 4) +
  plot_annotation(
    title = "Survival analysis of promoter methylatiom, all ovarian tissue samples",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

# Save to file
ggsave(
  filename = "KM_combined_methylation_w_HR_multi.png",
  plot = combined_met_plot,
  width = 25,
  height = 5,
  dpi = 300
)
#METHYLATION, HGSOC cases#################################
str(HGSOC_SURV_EXPRESSION[, colnames(HGSOC_SURV_EXPRESSION) %in% methylation])
#change methylation back to factor variable
HGSOC_SURV_EXPRESSION[, methylation] <- lapply(HGSOC_SURV_EXPRESSION[, methylation], function(x) {
  factor(x, levels = c(0, 1), labels = c("not methylated", "methylated"))
})

##HGSOC METHYLATION COX REGRESION, UNIVARIATE######################
univ_results_met2 <- lapply(methylation, function(gene) {
  formula <- as.formula(paste("Surv(OS, STATUS) ~", gene))
  cox_model <- coxph(formula, data = HGSOC_SURV_EXPRESSION)
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
met_univ_df2 <- do.call(rbind, univ_results_met2)
met_univ_df2$Gene <- methylation
met_univ_df2 <- met_univ_df2[, c("Gene", "HR", "lower95", "upper95", "pvalue")]
print(met_univ_df2)

##HGSOC METHYLATION MULTIVARIABLE REGRESSION############## 
cox_results2 <- lapply(methylation, function(gene) {
  formula <- as.formula(paste0("Surv(OS, STATUS) ~ ", gene, " + Age + CA125"))
  model <- coxph(formula, data = HGSOC_SURV_EXPRESSION)
  s <- summary(model)
  
  # Find all rows that correspond to this gene
  idx <- grep(paste0("^", gene), rownames(s$coefficients))
  
  # Skip if gene not found (safeguard)
  if(length(idx) == 0) return(NULL)
  
  # Extract info for each row (use lapply in case factor has multiple levels)
  do.call(rbind, lapply(idx, function(i) {
    data.frame(
      Gene = rownames(s$coefficients)[i],
      HR = s$coefficients[i, "exp(coef)"],
      Lower95 = s$conf.int[i, "lower .95"],
      Upper95 = s$conf.int[i, "upper .95"],
      p.value = s$coefficients[i, "Pr(>|z|)"],
      n = model$n
    )
  }))
})

cox_results_df2 <- do.call(rbind, cox_results2)
cox_results_df2

#fix names for later
cox_results_df2$Gene <- methylation


#plot HGSOC METHYLATION< with multivariable and univariable
plots_met2 <- list()

for (gene in methylation) {
  # Clean gene name
  gene_clean <- gene
  
  # Fit KM survival curve
  fit <- survfit(as.formula(paste("Surv(OS, STATUS) ~", gene)), data = HGSOC_SURV_EXPRESSION)
  
  ## --- Extract univariable Cox HR + CI from met_univ_df ---
  hr_row_uni <- met_univ_df2[met_univ_df2$Gene == gene, ]
  if (nrow(hr_row_uni) == 0 || any(is.na(hr_row_uni$HR))) {
    hr_text_uni <- bquote("Uni HR = NA")
  } else {
    hr_text_uni <- bquote("Uni HR = " * .(sprintf("%.2f (95%% CI: %.2f–%.2f)", 
                                                  hr_row_uni$HR, hr_row_uni$lower95, hr_row_uni$upper95)))
  }
  
  ## --- Extract multivariable Cox HR + CI from cox_results_df ---
  hr_row_multi <- cox_results_df2[cox_results_df2$Gene == gene, ]
  if (nrow(hr_row_multi) == 0 || any(is.na(hr_row_multi$HR))) {
    hr_text_multi <- bquote("Multi HR = NA")
  } else {
    hr_text_multi <- bquote("Multi HR = " * .(sprintf("%.2f (95%% CI: %.2f–%.2f)", 
                                                      hr_row_multi$HR, hr_row_multi$Lower95, hr_row_multi$Upper95)))
  }
  
  ## --- KM log-rank p-value ---
  survdiff_res <- survdiff(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                           data = HGSOC_SURV_EXPRESSION)
  pval_km <- 1 - pchisq(survdiff_res$chisq, length(survdiff_res$n) - 1)
  if (pval_km < 0.05) {
    pval_expr <- bquote(bold("Log-rank p = " ~ .(sprintf("%.3f", pval_km))))
  } else {
    pval_expr <- bquote("Log-rank p = " ~ .(sprintf("%.3f", pval_km)))
  }
  
  ## --- Combine Uni HR, Multi HR, and log-rank p-value in subtitle ---
  subtitle_expr <- bquote(.(hr_text_uni) ~ "; " ~ .(hr_text_multi) ~ "\n" ~ .(pval_expr))
  
  ## --- Legend labels ---
  legend_labels <- c("not methylated", "methylated")
  
  ## --- KM plot ---
  p <- ggsurvplot(
    fit,
    data = HGSOC_SURV_EXPRESSION,
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
  
  plots_met2[[gene]] <- p
}

# Combine plots with overall title
combined_met_plot2 <- wrap_plots(plots_met2, ncol = 4) +
  plot_annotation(
    title = "Survival analysis of promoter methylation, HGSOC tissue samples",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

# Save to file
ggsave(
  filename = "KM_combined_methylation_HGSOC_w_HR_multi.png",
  plot = combined_met_plot2,
  width = 25,
  height = 5,
  dpi = 300
)

#METHYLATION, OC cases#################################
str(OC_SURV_EXPRESSION[, colnames(OC_SURV_EXPRESSION) %in% methylation])
#change methylation back to factor variable
OC_SURV_EXPRESSION[, methylation] <- lapply(OC_SURV_EXPRESSION[, methylation], function(x) {
  factor(x, levels = c(0, 1), labels = c("not methylated", "methylated"))
})

##OC METHYLATION COX REGRESION, UNIVARIATE######################
univ_results_met3 <- lapply(methylation, function(gene) {
  formula <- as.formula(paste("Surv(OS, STATUS) ~", gene))
  cox_model <- coxph(formula, data = OC_SURV_EXPRESSION)
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
met_univ_df3 <- do.call(rbind, univ_results_met3)
met_univ_df3$Gene <- methylation
met_univ_df3 <- met_univ_df3[, c("Gene", "HR", "lower95", "upper95", "pvalue")]
print(met_univ_df3)

##OC METHYLATION MULTIVARIABLE REGRESSION############## 
cox_results3 <- lapply(methylation, function(gene) {
  formula <- as.formula(paste0("Surv(OS, STATUS) ~ ", gene, " + Age + CA125"))
  model <- coxph(formula, data = OC_SURV_EXPRESSION)
  s <- summary(model)
  
  # Find all rows that correspond to this gene
  idx <- grep(paste0("^", gene), rownames(s$coefficients))
  
  # Skip if gene not found (safeguard)
  if(length(idx) == 0) return(NULL)
  
  # Extract info for each row (use lapply in case factor has multiple levels)
  do.call(rbind, lapply(idx, function(i) {
    data.frame(
      Gene = rownames(s$coefficients)[i],
      HR = s$coefficients[i, "exp(coef)"],
      Lower95 = s$conf.int[i, "lower .95"],
      Upper95 = s$conf.int[i, "upper .95"],
      p.value = s$coefficients[i, "Pr(>|z|)"],
      n = model$n
    )
  }))
})

cox_results_df3 <- do.call(rbind, cox_results3)
cox_results_df3

#fix names for later
cox_results_df3$Gene <- methylation

#plot OC METHYLATION, with multivariable and univariable
plots_met3 <- list()

for (gene in methylation) {
  # Clean gene name
  gene_clean <- gene
  
  # Fit KM survival curve
  fit <- survfit(as.formula(paste("Surv(OS, STATUS) ~", gene)), data = OC_SURV_EXPRESSION)
  
  ## --- Extract univariable Cox HR + CI from met_univ_df ---
  hr_row_uni <- met_univ_df3[met_univ_df3$Gene == gene, ]
  if (nrow(hr_row_uni) == 0 || any(is.na(hr_row_uni$HR))) {
    hr_text_uni <- bquote("Uni HR = NA")
  } else {
    hr_text_uni <- bquote("Uni HR = " * .(sprintf("%.2f (95%% CI: %.2f–%.2f)", 
                                                  hr_row_uni$HR, hr_row_uni$lower95, hr_row_uni$upper95)))
  }
  
  ## --- Extract multivariable Cox HR + CI from cox_results_df ---
  hr_row_multi <- cox_results_df3[cox_results_df3$Gene == gene, ]
  if (nrow(hr_row_multi) == 0 || any(is.na(hr_row_multi$HR))) {
    hr_text_multi <- bquote("Multi HR = NA")
  } else {
    hr_text_multi <- bquote("Multi HR = " * .(sprintf("%.2f (95%% CI: %.2f–%.2f)", 
                                                      hr_row_multi$HR, hr_row_multi$Lower95, hr_row_multi$Upper95)))
  }
  
  ## --- KM log-rank p-value ---
  survdiff_res <- survdiff(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                           data = OC_SURV_EXPRESSION)
  pval_km <- 1 - pchisq(survdiff_res$chisq, length(survdiff_res$n) - 1)
  if (pval_km < 0.05) {
    pval_expr <- bquote(bold("Log-rank p = " ~ .(sprintf("%.3f", pval_km))))
  } else {
    pval_expr <- bquote("Log-rank p = " ~ .(sprintf("%.3f", pval_km)))
  }
  
  ## --- Combine Uni HR, Multi HR, and log-rank p-value in subtitle ---
  subtitle_expr <- bquote(.(hr_text_uni) ~ "; " ~ .(hr_text_multi) ~ "\n" ~ .(pval_expr))
  
  ## --- Legend labels ---
  legend_labels <- c("not methylated", "methylated")
  
  ## --- KM plot ---
  p <- ggsurvplot(
    fit,
    data = OC_SURV_EXPRESSION,
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
  
  plots_met3[[gene]] <- p
}

# Combine plots with overall title
combined_met_plot3 <- wrap_plots(plots_met3, ncol = 4) +
  plot_annotation(
    title = "Survival analysis of promoter methylation, OC tissue samples",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

# Save to file
ggsave(
  filename = "KM_combined_methylation_OC_w_HR_multi.png",
  plot = combined_met_plot3,
  width = 25,
  height = 5,
  dpi = 300
)


#LITHUANIAN PLOT, GENES#################################
#plot with multi
plots2x <- list()

for (gene in genes_f) {
  gene_clean <- sub("_f$", "", gene)
  
  # Fit KM survival curve
  fit <- survfit(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                 data = OC_SURV_EXPRESSION)
  
  ## --- Univariate HR ---
  hr_row <- univ_df2[univ_df2$Gene == gene, ]
  if (nrow(hr_row) == 0 || any(is.na(hr_row$HR))) {
    hr_text_expr <- bquote("Uni PR = NA")
  } else {
    hr_text_expr <- bquote("Uni PR = " * .(sprintf("%.2f (95 %% PI: %.2f–%.2f)", 
                                                   hr_row$HR, hr_row$lower95, hr_row$upper95)))
  }
  
  ## --- KM log-rank p-value ---
  survdiff_res <- survdiff(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                           data = OC_SURV_EXPRESSION)
  pval_km <- 1 - pchisq(survdiff_res$chisq, length(survdiff_res$n) - 1)
  if (pval_km < 0.05) {
    pval_expr <- bquote(bold("Log-rank p = " ~ .(sprintf("%.3f", pval_km))))
  } else {
    pval_expr <- bquote("Log-rank p = " ~ .(sprintf("%.3f", pval_km)))
  }
  
  ## --- Multivariable HR (Age + CA125) ---
  multi_row <- multi_df2[multi_df2$Gene == gene, ]
  if (nrow(multi_row) == 0 || any(is.na(multi_row$HR))) {
    hr_multi_expr <- bquote("Multiv PR = NA")
  } else {
    # Bold if p < 0.05
    if (!is.na(multi_row$pvalue) && multi_row$pvalue < 0.05) {
      hr_multi_expr <- bquote(bold("Multi PR = " * .(sprintf("%.2f (95 %% PI: %.2f–%.2f, N = %d)", 
                                                             multi_row$HR, multi_row$lower95, 
                                                             multi_row$upper95, multi_row$N))))
    } else {
      hr_multi_expr <- bquote("Multi PR = " * .(sprintf("%.2f (95 %% PI: %.2f–%.2f, N = %d)", 
                                                        multi_row$HR, multi_row$lower95, 
                                                        multi_row$upper95, multi_row$N)))
    }
  }
  
  ## --- Combine 3 lines in subtitle ---
  subtitle_text <- paste0(
    "Uni PR = ", ifelse(is.na(hr_row$HR), "NA",
                        sprintf("%.2f (95%% PI: %.2f–%.2f)", hr_row$HR, hr_row$lower95, hr_row$upper95)
    ),
    "\nLog-rank p = ", sprintf("%.3f", pval_km),
    "\nMulti PR = ", ifelse(is.na(multi_row$HR), "NA",
                            sprintf("%.2f (95%% PI: %.2f–%.2f, N = %d)", 
                                    multi_row$HR, multi_row$lower95, multi_row$upper95, multi_row$N)
    )
  )
  
  ## --- Legend labels ---
  legend_labels <- c(
    bquote(italic(.(gene_clean)) ~ " Maža raiška"),
    bquote(italic(.(gene_clean)) ~ " Didelė raiška")
  )
  
  ## --- KM plot ---
  p <- ggsurvplot(
    fit,
    data = OC_SURV_EXPRESSION,
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
      subtitle = subtitle_text,
      y = "Išgyvenamumas",
      x = "Laikas, mėnesiais"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, lineheight = 1.3)
    )
  
  plots2x[[gene]] <- p
}

# Combine plots and add overall title
combined_plotx2 <- wrap_plots(plots2x, ncol = 4) +
  plot_annotation(
    title = "Išgyvenamumo analizė KV grupėje",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )
###separate the 10 gene expression plot form the notch and wnt ##########
# Subset plots for the two gene sets
plots_genes10_2x <- plots2x[names(plots2x) %in% genes10_f]
plots_genes_notch_2x <- plots2x[names(plots2x) %in% genes_notch_f]

# Combine separately
combined_plot_genes10_2x <- wrap_plots(plots_genes10_2x, ncol = 5)+
  plot_annotation(
    title = "Išgyvenamumo analizė KV grupėje",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )

combined_plot_genes_notch_2x <- wrap_plots(plots_genes_notch_2x, ncol = 5)+
  plot_annotation(
    title = "Išgyvenamumo analizė KV grupėje",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )

#save 10 gene
ggsave(
  filename = "KM_combined_plot_w_HR_OC_10_gene_20250925LT2.png",  # output file name
  plot = combined_plot_genes10_2x,               # the patchwork plot object
  width = 25,                         # width in inches
  height = 12,                        # height in inches
  dpi = 400                           # resolution
)
#save notch
ggsave(
  filename = "KM_combined_plot_w_HR_OC_notch_20250925LT2.png",  # output file name
  plot = combined_plot_genes_notch_2x,               # the patchwork plot object
  width = 25,                         # width in inches
  height = 12,                        # height in inches
  dpi = 400                           # resolution
)

#LITHUANIAN plot OC METHYLATION ##########################
plots_met3 <- list()

for (gene in methylation) {
  gene_clean <- gene
  
  # Fit KM
  fit <- survfit(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                 data = OC_SURV_EXPRESSION)
  
  ## --- Univariate ---
  hr_row_uni <- met_univ_df3[met_univ_df3$Gene == gene, ]
  if (nrow(hr_row_uni) == 0 || any(is.na(hr_row_uni$HR))) {
    uni_text <- "Uni PR = NA"
  } else {
    uni_text <- sprintf("Uni PR = %.2f (95%% PI: %.2f–%.2f)",
                        hr_row_uni$HR, hr_row_uni$lower95, hr_row_uni$upper95)
  }
  
  ## --- Multivariable ---
  hr_row_multi <- cox_results_df3[cox_results_df3$Gene == gene, ]
  if (nrow(hr_row_multi) == 0 || any(is.na(hr_row_multi$HR))) {
    multi_text <- "Multi PR = NA"
  } else {
    multi_text <- sprintf("Multi PR = %.2f (95%% CI: %.2f–%.2f), n = 48",
                          hr_row_multi$HR, hr_row_multi$Lower95, hr_row_multi$Upper95)
  }
  
  ## --- KM log-rank p-value ---
  survdiff_res <- survdiff(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                           data = OC_SURV_EXPRESSION)
  pval_km <- 1 - pchisq(survdiff_res$chisq, length(survdiff_res$n) - 1)
  pval_text <- sprintf("Log-rank p = %.3f", pval_km)
  
  ## --- Combine lines ---
  subtitle_text <- paste0(uni_text, "\n", multi_text, "\n", pval_text)
  
  ## --- KM plot ---
  p <- ggsurvplot(
    fit,
    data = OC_SURV_EXPRESSION,
    pval = FALSE,
    risk.table = TRUE,
    legend.title = "",
    palette = c("blue", "red")
  )$plot +
    scale_color_manual(
      values = c("blue", "red"),
      labels = c("nemetilintas", "metilintas")
    ) +
    labs(
      title = bquote(italic(.(gene_clean))),
      subtitle = subtitle_text,
      y = "Išgyvenamumas",
      x = "Laikas, mėnesiais"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, lineheight = 1.3) # bigger & spaced
    )
  
  plots_met3[[gene]] <- p
}

# Combine plots with overall title
combined_met_plot3 <- wrap_plots(plots_met3, ncol = 4) +
  plot_annotation(
    title = "Išgyvenamumo analizė kiaušidžių vėžio grupėje",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

# Save to file
ggsave(
  filename = "KM_combined_methylation_OC_w_HR_multi0925.png",
  plot = combined_met_plot3,
  width = 25,
  height = 6,
  dpi = 600
  
)

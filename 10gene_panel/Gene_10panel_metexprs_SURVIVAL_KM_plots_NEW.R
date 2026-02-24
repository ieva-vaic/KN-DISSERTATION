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
library(ggtext)
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
  survdiff_res <- survdiff(
    as.formula(paste("Surv(OS, STATUS) ~", gene)),
    data = OC_SURV_EXPRESSION
  )
  
  pval_km <- 1 - pchisq(survdiff_res$chisq,
                        length(survdiff_res$n) - 1)
  
  logrank_text <- if (pval_km < 0.05) {
    paste0("**Log-rank p = ", sprintf("%.3f", pval_km), "**")
  } else {
    paste0("Log-rank p = ", sprintf("%.3f", pval_km))
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
  subtitle_text <- paste(
    paste0(
      "Uni PR = ",
      ifelse(is.na(hr_row$HR), "NA",
             sprintf("%.2f (95%% PI: %.2f–%.2f)",
                     hr_row$HR, hr_row$lower95, hr_row$upper95))
    ),
    logrank_text,
    paste0(
      "Multi PR = ",
      ifelse(is.na(multi_row$HR), "NA",
             sprintf("%.2f (95%% PI: %.2f–%.2f, N = %d)",
                     multi_row$HR, multi_row$lower95,
                     multi_row$upper95, multi_row$N))
    ),
    sep = "<br>"
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
    palette = c("blue", "red"),
    font.main = 17,        # title
    font.x = 17,
    font.y = 17,
    font.tickslab = 12,
    font.legend = 17
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
      plot.title = element_text(
        hjust = 0.5, face = "bold", size = 20
      ),
      plot.subtitle = element_markdown(
        size = 17, hjust = 0.5, lineheight = 1.3
      ),
      axis.title = element_text(size = 17),
      axis.text = element_text(size = 17),
      legend.text = element_text(size = 17)
    )
  
  plots2x[[gene]] <- p
}

# Combine plots and add overall title
combined_plotx2 <- wrap_plots(plots2x, ncol = 4)
###separate the 10 gene expression plot form the notch and wnt ##########
# Subset plots for the two gene sets
plots_genes10_2x <- plots2x[names(plots2x) %in% genes10_f]
plots_genes_notch_2x <- plots2x[names(plots2x) %in% genes_notch_f]

# Combine separately
combined_plot_genes10_2x <- wrap_plots(plots_genes10_2x, ncol = 2)

combined_plot_genes_notch_2x <- wrap_plots(plots_genes_notch_2x, ncol = 2)

#save 10 gene############################################################
# 5 rows, 2 columns
nrow <- 5
ncol <- 2

panel_width  <- 8  # in per panel
panel_height <- 4  # in per panel

ggsave(
  filename = "KM_combined_plot_w_HR_OC_10_gene_20260130LT.png",  # output file name
  plot = combined_plot_genes10_2x,
  width  = ncol * panel_width,
  height = nrow * panel_height,
  dpi    = 400  # higher DPI for publication quality
)
#save notch
ggsave(
  filename = "KM_combined_plot_w_HR_OC_notch_20260130LT.png",  # output file name
  plot = combined_plot_genes_notch_2x,               # the patchwork plot object
  width  = ncol * panel_width,
  height = nrow * panel_height,
  dpi    = 400                           # resolution
)

#ENGLISH PLOT, GENES#################################
#plot with multi
plots2xEN <- list()

for (gene in genes_f) {
  gene_clean <- sub("_f$", "", gene)
  
  # Fit KM survival curve
  fit <- survfit(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                 data = OC_SURV_EXPRESSION)
  
  ## --- Univariate HR ---
  hr_row <- univ_df2[univ_df2$Gene == gene, ]
  if (nrow(hr_row) == 0 || any(is.na(hr_row$HR))) {
    hr_text_expr <- bquote("Uni HR = NA")
  } else {
    hr_text_expr <- bquote("Uni HR = " * .(sprintf("%.2f (95 %% CI: %.2f–%.2f)", 
                                                   hr_row$HR, hr_row$lower95, hr_row$upper95)))
  }
  
  ## --- KM log-rank p-value ---
  survdiff_res <- survdiff(
    as.formula(paste("Surv(OS, STATUS) ~", gene)),
    data = OC_SURV_EXPRESSION
  )
  
  pval_km <- 1 - pchisq(survdiff_res$chisq,
                        length(survdiff_res$n) - 1)
  
  logrank_text <- if (pval_km < 0.05) {
    paste0("**Log-rank p = ", sprintf("%.3f", pval_km), "**")
  } else {
    paste0("Log-rank p = ", sprintf("%.3f", pval_km))
  }
  ## --- Multivariable HR (Age + CA125) ---
  multi_row <- multi_df2[multi_df2$Gene == gene, ]
  if (nrow(multi_row) == 0 || any(is.na(multi_row$HR))) {
    hr_multi_expr <- bquote("Multiv HR = NA")
  } else {
    # Bold if p < 0.05
    if (!is.na(multi_row$pvalue) && multi_row$pvalue < 0.05) {
      hr_multi_expr <- bquote(bold("Multi HR = " * .(sprintf("%.2f (95 %% CI: %.2f–%.2f, N = %d)", 
                                                             multi_row$HR, multi_row$lower95, 
                                                             multi_row$upper95, multi_row$N))))
    } else {
      hr_multi_expr <- bquote("Multi HR = " * .(sprintf("%.2f (95 %% CI: %.2f–%.2f, N = %d)", 
                                                        multi_row$HR, multi_row$lower95, 
                                                        multi_row$upper95, multi_row$N)))
    }
  }
  
  ## --- Combine 3 lines in subtitle ---
  subtitle_text <- paste(
    paste0(
      "Uni HR = ",
      ifelse(is.na(hr_row$HR), "NA",
             sprintf("%.2f (95%% CI: %.2f–%.2f)",
                     hr_row$HR, hr_row$lower95, hr_row$upper95))
    ),
    logrank_text,
    paste0(
      "Multi HR = ",
      ifelse(is.na(multi_row$HR), "NA",
             sprintf("%.2f (95%% CI: %.2f–%.2f, N = %d)",
                     multi_row$HR, multi_row$lower95,
                     multi_row$upper95, multi_row$N))
    ),
    sep = "<br>"
  )
  ## --- Legend labels ---
  legend_labels <- c(
    bquote(italic(.(gene_clean)) ~ " low expression"),
    bquote(italic(.(gene_clean)) ~ " high expression")
  )
  
  ## --- KM plot ---
  p <- ggsurvplot(
    fit,
    data = OC_SURV_EXPRESSION,
    pval = FALSE,
    risk.table = TRUE,
    legend.title = "",
    palette = c("blue", "red"),
    font.main = 17,        # title
    font.x = 17,
    font.y = 17,
    font.tickslab = 12,
    font.legend = 17
  )$plot +
    scale_color_manual(
      values = c("blue", "red"),
      labels = legend_labels
    ) +
    labs(
      title = bquote(italic(.(gene_clean))),
      subtitle = subtitle_text,
      #y = "Išgyvenamumas",
      x = "Time, months"
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5, face = "bold", size = 20
      ),
      plot.subtitle = element_markdown(
        size = 18, hjust = 0.5, lineheight = 1.3
      ),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 18),
      legend.text = element_text(size = 18)
    )
  
  plots2xEN[[gene]] <- p
}

# Combine plots and add overall title
combined_plotx2EN <- wrap_plots(plots2xEN, ncol = 4)
###separate the 10 gene expression plot form the notch and wnt ##########
# Subset plots for the two gene sets
plots_genes10_2xEN <- plots2xEN[names(plots2xEN) %in% genes10_f]
plots_genes_notch_2xEN <- plots2xEN[names(plots2xEN) %in% genes_notch_f]

# Combine separately
combined_plot_genes10_2xEN <- wrap_plots(plots_genes10_2xEN, ncol = 2)

combined_plot_genes_notch_2xEN <- wrap_plots(plots_genes_notch_2xEN, ncol = 2)

##save 10 gene#######################################################
ggsave(
  filename = "KM_combined_plot_w_HR_OC_10_gene_20260130EN.png",  # output file name
  plot = combined_plot_genes10_2xEN,               # the patchwork plot object
  width  = ncol * panel_width,
  height = nrow * panel_height,
  dpi    = 400                              # resolution
)
#save notch
ggsave(
  filename = "KM_combined_plot_w_HR_OC_notch_20260130EN.png",  # output file name
  plot = combined_plot_genes_notch_2xEN,               # the patchwork plot object
  width  = ncol * panel_width,
  height = nrow * panel_height,
  dpi    = 400                            # resolution
)
## Combine separately EN HOROZONTAL##########################################
combined_plot_genes10_2xEN <- wrap_plots(plots_genes10_2xEN, ncol = 5)

combined_plot_genes_notch_2xEN <- wrap_plots(plots_genes_notch_2xEN, ncol = 5)

# 5 rows, 2 columns
nrow <- 2
ncol <- 5

panel_width  <- 7  # in per panel
panel_height <- 6 # in per panel
#save 10 gene
ggsave(
  filename = "KM_combined_plot_w_HR_OC_10_gene_20260218EN.png",  # output file name
  plot = combined_plot_genes10_2xEN,               # the patchwork plot object
  width  = ncol * panel_width,
  height = nrow * panel_height,
  dpi    = 300                              # resolution
)
#save notch
ggsave(
  filename = "KM_combined_plot_w_HR_OC_notch_20260218EN.png",  # output file name
  plot = combined_plot_genes_notch_2xEN,               # the patchwork plot object
  width  = ncol * panel_width,
  height = nrow * panel_height,
  dpi    = 300                            # resolution
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
    multi_text <- sprintf("Multi PR = %.2f (95%% PI: %.2f–%.2f), n = 48",
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
      plot.title = element_text(hjust = 0.5, face = "bold", size = 25),
      plot.subtitle = element_text(size = 23, hjust = 0.5, lineheight = 1.3), # bigger & spaced
      legend.text = element_text(size = 20),
      axis.title.x = element_text(size = 22),      # bigger x-axis title
      axis.title.y = element_text(size = 22)      # bigger y-axis title
      )
  
  plots_met3[[gene]] <- p
}

# Combine plots with overall title
combined_met_plot3 <- wrap_plots(plots_met3, ncol = 2)

# Save to file methylation plots ###########################################
# 2 rows, 2 columns
nrow2 <- 2
ncol2 <- 2

panel_width2  <- 7  # in per panel
panel_height2 <- 6  # in per panel

ggsave(
  filename = "KM_combined_methylation_OC_w_HR_multi20260130.png",
  plot = combined_met_plot3,
  width  = ncol2 * panel_width2,
  height = nrow2 * panel_height2,
  dpi    = 400  
)
#ENGLISH plot OC METHYLATION ##########################
plots_met3EN <- list()

for (gene in methylation) {
  gene_clean <- gene
  
  # Fit KM
  fit <- survfit(as.formula(paste("Surv(OS, STATUS) ~", gene)), 
                 data = OC_SURV_EXPRESSION)
  
  ## --- Univariate ---
  hr_row_uni <- met_univ_df3[met_univ_df3$Gene == gene, ]
  if (nrow(hr_row_uni) == 0 || any(is.na(hr_row_uni$HR))) {
    uni_text <- "Uni HR = NA"
  } else {
    uni_text <- sprintf("Uni HR = %.2f (95%% CI: %.2f–%.2f)",
                        hr_row_uni$HR, hr_row_uni$lower95, hr_row_uni$upper95)
  }
  
  ## --- Multivariable ---
  hr_row_multi <- cox_results_df3[cox_results_df3$Gene == gene, ]
  if (nrow(hr_row_multi) == 0 || any(is.na(hr_row_multi$HR))) {
    multi_text <- "Multi HR = NA"
  } else {
    multi_text <- sprintf("Multi HR = %.2f (95%% CI: %.2f–%.2f), n = 48",
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
      labels = c("not methylated", "methylated")
    ) +
    labs(
      title = bquote(italic(.(gene_clean))),
      subtitle = subtitle_text,
      #y = "Išgyvenamumas",
      x = "Time, months"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 25),
      plot.subtitle = element_text(size = 23, hjust = 0.5, lineheight = 1.3), # bigger & spaced
      legend.text = element_text(size = 20),
      axis.title.x = element_text(size = 22),      # bigger x-axis title
      axis.title.y = element_text(size = 22)      # bigger y-axis title
    )
  
  plots_met3EN[[gene]] <- p
}

# Combine plots with overall title
combined_met_plot3EN <- wrap_plots(plots_met3EN, ncol = 2)

# Save to file
ggsave(
  filename = "KM_combined_methylation_OC_w_HR_multiEN20260130.png",
  plot = combined_met_plot3EN,
  width  = ncol2 * panel_width2,
  height = nrow2 * panel_height2,
  dpi    = 400  
)
#HGSOC for supplements LITHUANIAN######################################
##COX REGRESION, UNIVARIATE, HGSOC only ###################
univ_results3 <- lapply(genes, function(gene) {
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
    hr_text_expr <- bquote("Univ PR = NA")
  } else {
    hr_text_expr <- bquote("Univ PR = " * .(sprintf("%.2f (95%% PI: %.2f–%.2f)", 
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
      hr_multi_expr <- bquote(bold("Multiv PR = " * .(sprintf("%.2f (95%% PI: %.2f–%.2f, N=%d)", 
                                                              multi_row$HR, multi_row$lower95, 
                                                              multi_row$upper95, multi_row$N))))
    } else {
      hr_multi_expr <- bquote("Multiv PR = " * .(sprintf("%.2f (95%% PI: %.2f–%.2f, N=%d)", 
                                                         multi_row$HR, multi_row$lower95, 
                                                         multi_row$upper95, multi_row$N)))
    }
  }
  
  ## --- Combine 3 lines in subtitle: Univ HR | Log-rank p | Multiv HR ---
  subtitle_text <- paste0(
    ifelse(is.na(hr_row$HR), "Univ PR = NA", 
           sprintf("Univ PR = %.2f (95%% PI: %.2f–%.2f)",
                   hr_row$HR, hr_row$lower95, hr_row$upper95)),
    "<br>",
    paste0("Log-rank p = ", sprintf("%.3f", pval_km)),
    "<br>",
    ifelse(is.na(multi_row$HR), "Multiv PR = NA", 
           sprintf("Multiv PR = %.2f (95%% PI: %.2f–%.2f, N=%d)",
                   multi_row$HR, multi_row$lower95,
                   multi_row$upper95, multi_row$N))
  )
  
  ## --- Legend labels ---
  legend_labels <- c(
    bquote(italic(.(gene_clean)) ~ " maža raiška"), #for lithuanian
    bquote(italic(.(gene_clean)) ~ " didelė raiška")
  )
  
  ## --- KM plot ---
  p <- ggsurvplot(
    fit,
    data = HGSOC_SURV_EXPRESSION,
    pval = FALSE,
    risk.table = TRUE,
    legend.title = "",
    palette = c("blue", "red"),
    font.main = 17,        # title
    font.x = 17,
    font.y = 17,
    font.tickslab = 12,
    font.legend = 17
  )$plot +
    scale_color_manual(
      values = c("blue", "red"),
      labels = legend_labels
    ) +
    labs(
      title = gene_clean,
      subtitle = subtitle_text,
      y = "Išgyvenamumas",
      x = "Laikas, mėnesiais"
    ) +
    theme(
      plot.title = element_text(
        hjust = 0.5, face = "italic", size = 20
      ),
      plot.subtitle = element_markdown(
        size = 18, hjust = 0.5, lineheight = 1.3
      ),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 18),
      legend.text = element_text(size = 18)
    )
  
  plots3x[[gene]] <- p
}

# Combine plots and add overall title
 combined_plot3x <- wrap_plots(plots3x, ncol = 4) 
###separate the 10 gene expression plot form the notch and wnt ##########
# Subset plots for the two gene sets
plots_genes10_3x <- plots3x[names(plots3x) %in% genes10_f]
plots_genes_notch_3x <- plots3x[names(plots3x) %in% genes_notch_f]

# Combine separately
combined_plot_genes10_3x <- wrap_plots(plots_genes10_3x, ncol = 2)+
  plot_annotation(
    title = "Išgyvenamumas HGSOC grupėje",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )

combined_plot_genes_notch_3x <- wrap_plots(plots_genes_notch_3x, ncol = 2)+
  plot_annotation(
    title = "Išgyvenamumas HGSOC grupėje",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  )

#save 10 gene SUPPLEMENTS########################################
ggsave(
  filename = "KM_combined_plot_w_HR_HGSOC_10_gene_20260130.png",  # output file name
  plot = combined_plot_genes10_3x,               # the patchwork plot object
  width  = ncol * panel_width,
  height = nrow * panel_height,
  dpi    = 400                           # resolution
)
#save notch
ggsave(
  filename = "KM_combined_plot_w_HR_HGSOC_notch_20260130.png",  # output file name
  plot = combined_plot_genes_notch_3x,               # the patchwork plot object
  width  = ncol * panel_width,
  height = nrow * panel_height,
  dpi    = 400                          # resolution
)

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
    hr_text_uni <- bquote("Uni PR = NA")
  } else {
    hr_text_uni <- bquote("Uni PR = " * .(sprintf("%.2f (95%% PI: %.2f–%.2f)", 
                                                  hr_row_uni$HR, hr_row_uni$lower95, hr_row_uni$upper95)))
  }
  
  ## --- Extract multivariable Cox HR + CI from cox_results_df ---
  hr_row_multi <- cox_results_df2[cox_results_df2$Gene == gene, ]
  if (nrow(hr_row_multi) == 0 || any(is.na(hr_row_multi$HR))) {
    hr_text_multi <- bquote("Multi PR = NA")
  } else {
    hr_text_multi <- bquote("Multi PR = " * .(sprintf("%.2f (95%% PI: %.2f–%.2f)", 
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
  legend_labels <- c("nemetilintas", "metilintas")
  
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
      subtitle = subtitle_expr,
      y = "Išgyvenamumas",
      x = "Laikas, mėnesiais"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, lineheight = 1.1)
    )
  
  plots_met2[[gene]] <- p
}

# Combine plots with overall title
combined_met_plot2 <- wrap_plots(plots_met2, ncol = 2) +
  plot_annotation(
    title = "Išgyvenamumumas HGSOC grupėje",
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

# Save to file

# 2 rows, 2 columns
nrow2 <- 2
ncol2 <- 2

panel_width2  <- 8  # in per panel
panel_height2 <- 4  # in per panel
ggsave(
  filename = "KM_combined_methylation_HGSOC_w_HR_multi20260130.png",
  plot = combined_met_plot2,
  width  = ncol2 * panel_width2,
  height = nrow2 * panel_height2,
  dpi    = 400 
)

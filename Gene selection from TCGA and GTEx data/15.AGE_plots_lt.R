#KN-DISSERTATION project. Genes selection via statistical analysis of TCGA and GTEx data 
#(only scripts producing lithuanian figures, included in the dissertation)
#Plot correlation with age
# Load packages ##########################################
Sys.setenv(LANG = "en")
library(tidyverse)
library(FSA)
library(grid)
library(ggpubr)
library(magick)
library(gt)
library(multcomp)
library(FSA)
library(cowplot)
library(patchwork)
#set directory of the data
setwd("../TCGA-OV-RISK-PROJECT/Public data RDSs/")
#set genes of interest
expression <- c("EXO1" ,  "RAD50" , "PPT2" ,  "LUC7L2" ,"PKP3"  , "CDCA5" ,
                "ZFPL1" , "VPS33B", "GRB7" ,  "TCEAL4")
#TEST####################################################################
# Load test data ###################################
full_test <- read.csv("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/full_test_clinicals20250619.csv")
full_test$AGE2 <- factor(full_test$AGE2, ordered = T)
#pivot long
age_long <- full_test %>%
  dplyr::select(all_of(expression), AGE2) %>%
  pivot_longer(cols = all_of(expression), names_to = "gene", values_to = "expression") %>%
  filter(!is.na(AGE2))
#BOXPLOT: Age TEST ########################################
# Step 1: Compute Kruskal-Wallis test for each gene (number version)
# kw_df <- lapply(expression, function(gene) {
#   test <- kruskal.test(full_test[[gene]] ~ full_test$AGE2)
#   data.frame(
#     gene = gene,
#     statistic = test$statistic,
#     p_value = test$p.value
#   )
# }) %>%
#   do.call(rbind, .) %>%
#   mutate(
#     label = ifelse(
#       p_value < 0.001,
#       "italic(p) < 0.001",
#       paste0("italic(p) == ", signif(p_value, 3))
#     ),
#     gene_label = paste0("italic(", gene, ")")
#   )

# Step 1: Compute Kruskal-Wallis test for each gene and convert to stars
kw_df <- lapply(expression, function(gene) {
  test <- kruskal.test(full_test[[gene]] ~ full_test$AGE2)
  data.frame(
    gene = gene,
    statistic = test$statistic,
    p_value = test$p.value
  )
}) %>%
  do.call(rbind, .) %>%
  mutate(
    # Convert p-values to significance stars
    label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    # Keep gene label for facet
    gene_label = paste0("italic(", gene, ")")
  )


# Step 2: Add labels to long-format data
age_long <- age_long %>%
  mutate(gene_label = paste0("italic(", gene, ")")) %>%
  left_join(kw_df[, c("gene", "label")], by = "gene")
#plot
age_test_plot <- 
  ggplot(age_long, aes(x = AGE2, y = expression)) +
  geom_boxplot(aes(fill = AGE2), outlier.shape = NA, alpha = 0.6) +
  geom_jitter(aes(color = AGE2), width = 0.2, alpha = 0.8, size = 1) +
  geom_text(
    data = kw_df,
    aes(x = 1, y = Inf, label = label),
    inherit.aes = FALSE,
    parse = TRUE,
    hjust = -0.1, vjust = 1.2
  ) +
  facet_wrap(~ gene_label, scales = "free_y", labeller = label_parsed, nrow = 2, ncol = 5) +
  theme_minimal() +
  labs(
    title = "Genų raiškos sąsaja su amžiumi testavimo imtyje",
    x = "Amžiaus intervalai", y = "Genų raiška",
    fill = "Amžius", color = "Amžius"
  ) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))+
  scale_fill_brewer(palette = "Pastel2") +
  scale_color_brewer(palette = "Pastel2") 
#save plot
png("C:/Users/Ieva/rprojects/outputs_all/DISS/age_groups_vs_expressiontest20260127.png",
    width = 15, height = 10, res = 300, units = "cm") # width and height in pixels, resolution in dpi
age_test_plot #
grid.text("B", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
          gp = gpar(fontsize = 18, fontface = "bold"))
dev.off() # Close the PNG device

#Count age intervals #########
full_test$id <-  substr(trimws(full_test$barcode), 1, 4)
table(full_test$AGE2,full_test$id, useNA = "a")
table(full_test$id, useNA = "a")

#TRAIN ##################################################################
# Load train data ###################################
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot_2025.RDS")
#filter for lasso genes
gtex_genes <- readRDS("gtcga_elastic_2025.RDS")
gtex_filtered_counts_train <- gtex_counts_train[colnames(gtex_counts_train) %in% gtex_genes] 
#filter for TCGA cases
gtex_filtered_counts_train2 <- gtex_filtered_counts_train %>%
  dplyr::filter(grepl("^TCGA", rownames(gtex_filtered_counts_train)))
dim(gtex_filtered_counts_train2)

#upload age
full_clinicals <- read.csv("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/pheno_best202520619.csv")
rownames(full_clinicals) <- full_clinicals$X
train_data <- gtex_filtered_counts_train[, colnames(gtex_filtered_counts_train) %in% c(expression)]
rownames(full_clinicals) == rownames(train_data)
full_age_df <- merge(full_clinicals, train_data, by = "row.names", all = TRUE)
rownames(full_age_df) <- full_age_df$Row.names
full_age_df$Row.names <- NULL

#pivot long
age_long1 <- full_age_df %>%
  dplyr::select(all_of(expression), AGE2) %>%
  pivot_longer(cols = all_of(expression), names_to = "gene", values_to = "expression") %>%
  filter(!is.na(AGE2))
#BOXPLOT: Age TRAIN ########################################
# Step 1: Compute Kruskal-Wallis test for each gene
kw_df1 <- lapply(expression, function(gene) {
  test <- kruskal.test(full_age_df[[gene]] ~ full_age_df$AGE2)
  data.frame(
    gene = gene,
    statistic = test$statistic,
    p_value = test$p.value
  )
}) %>%
  do.call(rbind, .) %>%
  mutate(
    # round to 3 decimals
    label = ifelse(
      p_value < 0.001,
      "italic(p) < 0.001",
      paste0("italic(p) == ", format(round(p_value, 3), nsmall = 3))
    ),
    gene_label = paste0("italic(", gene, ")")
  )

# Step 2: Add labels to long-format data
age_long1 <- age_long1 %>%
  mutate(gene_label = paste0("italic(", gene, ")")) %>%
  left_join(kw_df1[, c("gene", "label")], by = "gene")
#plot
# Step 3: Plot with Kruskal-Wallis p-values
ggplot(age_long1, aes(x = AGE2, y = expression)) +
  geom_boxplot() +
  geom_text(
    data = kw_df1,
    aes(x = 1, y = Inf, label = label),
    inherit.aes = FALSE,
    parse = TRUE,
    hjust = -0.1, vjust = 1.2
  ) +
  facet_wrap(~ gene_label, scales = "free_y", labeller = label_parsed) +
  theme_minimal() +
  labs(title = "Genų raiškos sąsaja su amžiumi",
       x = "Amžiaus intervalai", y = "Genų raiška")
#make significant comparison df
full_age_df$AGE2 <- factor(full_age_df$AGE2, ordered = FALSE)
dunn_results <- lapply(expression, function(gene) {
  formula <- as.formula(paste0("`", gene, "` ~ AGE2"))
  test <- dunnTest(formula, data = full_age_df, method = "bh")
  result <- test$res
  result$gene <- gene
  result
})
dunn_df <- do.call(rbind, dunn_results)
dunn_sig <- dunn_df %>%
  dplyr::filter(P.adj < 0.05)
dunn_sig
#make plot with age and significant results, test
# Step 1: Split the 'Comparison' column into two groups
sig_pairs <- dunn_sig %>%
  mutate(
    group1 = sub(" - .*", "", Comparison),
    group2 = sub(".* - ", "", Comparison)
  )

# Step 2: Add y.position for each gene to avoid label overlap
# Use max expression per gene from age_long to space labels properly
y_positions <- age_long1 %>%
  group_by(gene) %>%
  summarize(max_expr = max(expression, na.rm = TRUE))

sig_pairs <- sig_pairs %>%
  left_join(y_positions, by = "gene") %>%
  group_by(gene) %>%
  mutate(
    y.position = max_expr + 0.3 * row_number()  # spacing out labels vertically
  ) %>%
  ungroup()

# Step 3: Format p-values as strings for labels
sig_pairs <- sig_pairs %>%
  mutate(
    label = ifelse(P.adj < 0.001, "p < 0.001", paste0("p = ", signif(P.adj, 2)))
  )

# Step 4: Add gene_label for facet matching (italic gene names)
sig_pairs <- sig_pairs %>%
  mutate(gene_label = paste0("italic(", gene, ")"))

# Step 5: change label
sig_pairs <- sig_pairs %>%
  mutate(
    label = ifelse(
      P.adj < 0.001,
      "p < 0.001",
      paste0("p = ", format(round(P.adj, 3), nsmall = 3))
    )
  )
# Step 6: change label to *
sig_pairs <- sig_pairs %>%
  mutate(
    label = case_when(
      P.adj < 0.001 ~ "***",
      P.adj < 0.01  ~ "**",
      P.adj < 0.05  ~ "*",
      TRUE          ~ ""  # no star for non-significant
    )
  )
# FINAL plot
age_train_plot <- ggplot(age_long1, aes(x = AGE2, y = expression)) +
  geom_boxplot(outlier.shape = NA, aes(fill = AGE2), alpha = 0.6) +   # colored boxes
  geom_jitter(aes(color = AGE2), width = 0.2, alpha = 0.8, size = 1) +  # colored points
  facet_wrap(~ gene_label, scales = "free_y", labeller = label_parsed, , nrow = 2, ncol = 5) +
  theme_minimal() +
  labs(title = "Genų raiškos sąsaja su amžiumi mokymosi imtyje",
       x = "Amžiaus intervalai", y = "Genų raiška",
       fill = "Amžius",   
       color = "Amžius")  +
  stat_pvalue_manual(
    data = sig_pairs,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    parse = FALSE,
    tip.length = 0.01,
    size = 5
  ) +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))+  # put legend below for clarity
  scale_fill_brewer(palette = "Pastel2") +
  scale_color_brewer(palette = "Pastel2") 
#save
png("C:/Users/Ieva/rprojects/outputs_all/DISS/age_groups_vs_expressiontrain20260127.png",
    width = 15, height = 15, res = 300, units = "cm") # width and height in pixels, resolution in dpi
age_train_plot #
grid.text("A", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
          gp = gpar(fontsize = 18, fontface = "bold"))
dev.off() # Close the PNG device

#PLOT TOGETHER##########################
# Read images
img1 <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/age_groups_vs_expressiontrain20260127.png")
img2 <- image_read("C:/Users/Ieva/rprojects/outputs_all/DISS/age_groups_vs_expressiontest20260127.png")

# Stack vertically
stacked <- image_append(c(img1, img2), stack = TRUE)  # stack=TRUE → vertical

# Save result
image_write(stacked, "C:/Users/Ieva/rprojects/outputs_all/DISS/combined_age_plots20260127.png")

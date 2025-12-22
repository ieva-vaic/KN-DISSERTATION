#KN-DISSERTATION project. 10 gene panel - genes selected via statistical analysis of TCGA and GTEx data
#KN 10 gene analysis data (cv 2 %, one threshold, 65 cases!, restored KN-59) 2025-02-17
#t.tests between sample groups, OC vs benign
Sys.setenv(LANG = "en")
#libraries
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
#set wd for plots
setwd("C:/Users/Ieva/rprojects/outputs_all/DISS/")
#data
OC_full <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/OC_10_genes_clean_2025_02_14.RDS")
#expression genes list
expression <- c("EXO1", "RAD50","PPT2", "LUC7L2","PKP3", "CDCA5","ZFPL1","VPS33B", "GRB7","TCEAL4")

#SHAPIRO TESTS###############################################
table(OC_full$tumor) 
# OVCa vs benign
shapiro_results <- OC_full[, c(3:13 )] %>%
  pivot_longer(cols = -tumor, names_to = "gene", values_to = "value") %>%
  group_by(tumor, gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop") %>%
  filter(p_value < 0.05)
shapiro_results #not notmal for OC PPT2, LUC7L2, PKP3, TCEAL4

#var.test####################################################
#var test OVCa vs benign
var_results <- OC_full[, c(3:13 )] %>%
  pivot_longer(cols = -tumor, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = var.test(value[tumor == unique(tumor)[1]], 
                       value[tumor == unique(tumor)[2]])$p.value,
    .groups = "drop"
  ) %>%
  filter(p_value < 0.05)
var_results #luc7l2 not equal variances
#2GROUPS#####################
#test 2 groups according to variance and normalcy
#MANN/ WILCOX for "PPT2", "LUC7L2", "PKP3", "TCEAL4"
#STJUDENTS T test  "EXO1","RAD50","CDCA5","ZFPL1","VPS33B","GRB7" 

#melt table for expression
OC_table <- melt(OC_full[, c(3:13 )], id.vars="tumor",  measure.vars=expression)
#Mann-whiney / wilcox test: not normal 
wilcox.test_OC_2 <- OC_table %>%
  filter(variable %in% c("PPT2", "LUC7L2", "PKP3", "TCEAL4")) %>%
  group_by(variable) %>%
  summarise(
    p_value = wilcox.test(value ~ tumor)$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%  # Adjust p-values
  filter(p_value < 0.1)  # Keep significant results
wilcox.test_OC_2

#Stujents t: normal
stat.test_OC_2 <- OC_table %>%
  filter(variable %in% c( "EXO1","RAD50", "CDCA5","ZFPL1","VPS33B","GRB7")) %>%
  group_by(variable) %>%
  summarise(
    p_value = t.test(value ~ tumor, var.equal = TRUE)$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%  # Adjust p-values
  filter(p_value < 0.1)  # Keep significant results
stat.test_OC_2

##tribble all tests together 2 groups################
each.vs.ref_sig_tumor <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "Gerybiniai",   "KV", 0.0246       , -5, "EXO1", #stjudents
  "Gerybiniai",   "KV", 0.00197       ,  -2, "RAD50",#stjudents
  "Gerybiniai",   "KV", 0.000301 ,  -4, "PPT2",#WILCOX
  "Gerybiniai",   "KV",  0.000925 , 1, "LUC7L2", #wilcox
  "Gerybiniai",   "KV",  0.000301 , -4, "PKP3", #wilcox
  "Gerybiniai",   "KV", 0.00213      , -2, "CDCA5",#stjudents
  "Gerybiniai",   "KV", 0.00282      , -3, "ZFPL1",#stjudents
  "Gerybiniai",   "KV", 0.00253      , -4, "VPS33B", #stjudents
  "Gerybiniai",   "KV", 0.00000000242, -1, "GRB7",# stjudents
  "Gerybiniai",   "KV", 0.0000337, 3, "TCEAL4",# wilcox
)
# 3 digits after the ,
each.vs.ref_sig_tumor$p.adj_custom <- ifelse(each.vs.ref_sig_tumor$p.adj < 0.001, 
                                             "p < 0.001", 
                                             paste0("p = ", sprintf("%.3f",
                                                                    each.vs.ref_sig_tumor$p.adj)))
#change to lithuanin
OC_table <- OC_table %>% mutate(tumor = case_when(
  tumor  == "OC" ~ "KV",
  tumor  == "Benign" ~ "Gerybiniai"
))

##boxplot 2 groups###############
custom_colors <- c("KV" = "pink2", "Gerybiniai" = "blue")
tumor_plot <- ggplot(OC_table, aes(x=tumor, y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = tumor)) +
  geom_jitter(aes(color = tumor), size=2, alpha=0.5) +
  ylab(label = expression("Santykinė genų raiška, normalizuota pagal  " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(each.vs.ref_sig_tumor, label = "p.adj_custom") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  #labs(tag = "A")+ 
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(labels = function(x) gsub("-", "\u2212", x))+
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x))) #add long "-" signs

tumor_plot #

##save plot 2 groups ###########################
png("10genes_boxplot2groups_2025-06-20.png"
    , width = 3000, height = 2300, res = 300) # width and height in pixels, resolution in dpi
tumor_plot #
dev.off() # Close the PNG device

##FC 2 groups ##################
expression_df_2gr <- OC_full[, c("tumor", expression, "KN")]
rownames(expression_df_2gr) <- expression_df_2gr$KN
expression_df_2gr <- expression_df_2gr[, -12]
#REMOVE NA
expression_df2_2gr <- expression_df_2gr %>%
  filter(complete.cases(.))

exp_df2_2gr <- expression_df2_2gr %>%
  melt(id.vars="tumor",  measure.vars=expression) %>%
  rename(gene = variable) %>%
  rename(expression = value)%>%
  rename(condition = tumor) %>%
  group_by(condition, gene)%>%
  summarise(mean_expression = mean(expression, drop_na= T)) %>%
  spread(condition, mean_expression)

mean_expression_tumor2 <- exp_df2_2gr %>%
  mutate(fold_change_HB = log2(2^(`OC` - `Benign`)))
mean_expression_tumor2

## EN boxplot 2 groups###############
OC_tableEN <- melt(OC_full[, c(3:13 )], id.vars="tumor",  measure.vars=expression)
##tribble all tests together 2 groups################
each.vs.ref_sig_tumorEN <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "Benign",   "OC", 0.0246       , -5, "EXO1", #stjudents
  "Benign",   "OC", 0.00197       ,  -2, "RAD50",#stjudents
  "Benign",   "OC", 0.000301 ,  -4, "PPT2",#WILCOX
  "Benign",   "OC",  0.000925 , 1, "LUC7L2", #wilcox
  "Benign",   "OC",  0.000301 , -4, "PKP3", #wilcox
  "Benign",   "OC", 0.00213      , -2, "CDCA5",#stjudents
  "Benign",   "OC", 0.00282      , -3, "ZFPL1",#stjudents
  "Benign",   "OC", 0.00253      , -4, "VPS33B", #stjudents
  "Benign",   "OC", 0.00000000242, -1, "GRB7",# stjudents
  "Benign",   "OC", 0.0000337, 3, "TCEAL4",# wilcox
)
# 3 digits after the ,
each.vs.ref_sig_tumorEN$p.adj_custom <- ifelse(each.vs.ref_sig_tumorEN$p.adj < 0.001, 
                                               "p < 0.001", 
                                               paste0("p = ", sprintf("%.3f",
                                                                      each.vs.ref_sig_tumorEN$p.adj)))


custom_colorsEN <- c("OC" = "pink2", "Benign" = "blue")
tumor_plotEN <- ggplot(OC_tableEN, aes(x=tumor, y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = tumor)) +
  geom_jitter(aes(color = tumor), size=2, alpha=0.5) +
  ylab(label = expression("Relative expression, normalized to  " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(each.vs.ref_sig_tumorEN, label = "p.adj_custom") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  #labs(tag = "A")+ 
  scale_fill_manual(values = custom_colorsEN) +
  scale_color_manual(values = custom_colorsEN) +
  scale_y_continuous(labels = function(x) gsub("-", "\u2212", x))+
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x))) #add long "-" signs

tumor_plotEN #

##save plot 2 groups ###########################
png("10genes_boxplot2groupsEN_20251217.png"
    , width = 3000, height = 2300, res = 350) # width and height in pixels, resolution in dpi
tumor_plotEN #
dev.off() # Close the PNG device

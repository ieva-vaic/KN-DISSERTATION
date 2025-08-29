#KN-DISSERTATION project. 10 gene panel - genes selected via statistical analysis of TCGA and GTEx data
#KN 10 gene analysis data (cv 2 %, one threshold, 65 cases!, restored KN-59) 2025-02-17
#t.tests between sample groups, all variations
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

#3 GROUPS:HGSOC, OTHER, BENIGN
shapiro_results2 <- OC_full[, c(2:12 )] %>%
  pivot_longer(cols = -type, names_to = "gene", values_to = "value") %>%
  group_by(type, gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop") %>%
  filter(p_value < 0.05)
shapiro_results2 #not notmal for OC LUC7L2, PKP3, TCEAL4
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

#var.test to 3 groups#############################
#hgsoc vs benign
var_results_benign_hgsoc <- OC_full[, c(2:12 )] %>%
  filter(type != "Other")%>%
  pivot_longer(cols = -type, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = var.test(value[type == unique(type)[1]], 
                       value[type == unique(type)[2]])$p.value,
    .groups = "drop"
  ) %>%
  filter(p_value < 0.05)
var_results_benign_hgsoc #luc7l2 not equal variances
#hgsoc vs other
var_results_other_hgsoc <- OC_full[, c(2:12 )] %>%
  filter(type != "Benign")%>%
  pivot_longer(cols = -type, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = var.test(value[type == unique(type)[1]], 
                       value[type == unique(type)[2]])$p.value,
    .groups = "drop"
  ) %>%
  filter(p_value < 0.05)
var_results_other_hgsoc # pkp3 unequal
#other vs benign
var_results_other_benign <- OC_full[, c(2:12 )] %>%
  filter(type != "HGSOC")%>%
  pivot_longer(cols = -type, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = var.test(value[type == unique(type)[1]], 
                       value[type == unique(type)[2]])$p.value,
    .groups = "drop"
  ) %>%
  filter(p_value < 0.05)
var_results_other_benign # luc7l2 

#TEST EXAMPLES ####################################
#NORMAL DISTRIBUTION AND VARIANCE: STUJENT T.TEST
#NOT NORMAL DISTRIBUTION, NORMAL VARIANCE: Mann-Whitney U test (Wilcoxon rank-sum test)
#NOT NORMAL DISTRIBUTION AND VARIANCE:  Welch’s t-test (if near normal) OR Mann-Whitney U test
#NORMAL DISTRIBUTION BUT NOT VARIANCE: Welch’s t-test 

##STJUDENT TEST ###################################
#melt table for expression
OC_table <- melt(OC_full[, c(3:13 )], id.vars="tumor",  measure.vars=expression)
#Stujents t test for two groups
stat.test_OC <- OC_table %>%
  group_by(variable) %>%
  summarise(
    p_value = t.test(value ~ tumor, var.equal = TRUE)$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%  # Adjust p-values
  filter(p_value < 0.1)  # Keep significant results
stat.test_OC

##WELCH TEST ###################################
#Welch test (normal, unequal variances), for two groups
welch.test_OC <- OC_table %>%
  group_by(variable) %>%
  summarise(
    p_value = t.test(value ~ tumor, var.equal = FALSE)$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%  # Adjust p-values
  filter(p_value < 0.1)  # Keep significant results
welch.test_OC

##MANN-WHITNEY (Wilcoxon) TEST ###################################
#Wilcoxon test (not normal), for two groups
wilcox.test_OC <- OC_table %>%
  group_by(variable) %>%
  summarise(
    p_value = wilcox.test(value ~ tumor)$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%  # Adjust p-values
  filter(p_value < 0.1)  # Keep significant results
wilcox.test_OC

#I have multiple groups
##PAIRWISE STJUDENTS T TEST ###################################
#stjundents test (normal, equal variances)
t.test_3groups <- Group3_table %>%
  group_by(variable) %>%
  t_test(value ~ type,
         p.adjust.method = "BH", 
         var.equal = TRUE, #stjudents
         paired = FALSE, 
         #detailed=TRUE 
  )%>%
  filter(p.adj < 0.05) %>%  # Keep significant results
  mutate(across(c(p.adj), ~ format(., scientific = FALSE)))  # Format p-values to remove scientific notation
t.test_3groups

##PAIRWISE WELCH TEST ###################################
#Welch test (normal, unequal variances)
welch.test_3groups <- Group3_table %>%
  group_by(variable) %>%
  t_test(value ~ type,
         p.adjust.method = "BH",
         var.equal = FALSE, #welch
         paired = FALSE) %>% 
  filter(p.adj < 0.1)
welch.test_3groups

##PAIRWISE MANN-WHITNEY (Wilcoxon) TEST ###################################
#melt table for expression
Group3_table <- melt(OC_full[, c(2:12 )], id.vars="type",  measure.vars=expression)
#Wilcoxon test (not normal)
wilcox.test_3groups <- Group3_table %>%
  group_by(variable) %>%
  pairwise_wilcox_test(value ~ type,
                       #ref.group = "Benign" , #only if one group is needed
                       p.adjust.method = "BH") %>%
  filter(p.adj < 0.1)
wilcox.test_3groups

#3GROUPS#################################################
#tests 3 groups according to variance and normalcy:
#MANN/ WILCOX for "LUC7L2", "PKP3", "TCEAL4"
#WELCH for "PPT2" with others group
#STJUDENTS T test (without others) "PPT2","EXO1","RAD50","CDCA5","ZFPL1","VPS33B"

#MANN WITHNY / Wilcox test (not normal)
Group3_table <- melt(OC_full[, colnames(OC_full) %in%
                               c("type", "LUC7L2", "PKP3", "TCEAL4" )],
                     id.vars="type",
                     measure.vars= c( "LUC7L2", "PKP3", "TCEAL4"))
#test wilcox
wilcox.test_3groups <- Group3_table %>%
  group_by(variable) %>%
  pairwise_wilcox_test(value ~ type,
                       p.adjust.method = "BH") %>%
  filter(p.adj < 0.1)
wilcox.test_3groups

#Welch's t test (normal dist not variance)
#for ppt2 with others
welch.test_others<- OC_full[, colnames(OC_full) %in% c("type", "PPT2")] %>%
  t_test(PPT2 ~ type,
                  p.adjust.method = "BH",
                  var.equal = FALSE, #welch
                  paired = FALSE) 
welch.test_others #only for PPT2 comparisons to others, ALL NS

#t tests stjudents (normal dist. + equal var.)
Group3_table_normal <- melt(OC_full[, colnames(OC_full) %in%
                               c("type", "PPT2" , "EXO1","RAD50", "CDCA5", "ZFPL1","GRB7",  "VPS33B" )],
                     id.vars="type",
                     measure.vars= c("PPT2" , "EXO1", "RAD50", "CDCA5", "ZFPL1","GRB7", "VPS33B"))
#stjudents tests
t.test_3groups <- Group3_table_normal %>%
  group_by(variable) %>%
  t_test(value ~ type,
         p.adjust.method = "BH", 
         var.equal = TRUE, #stjudents
         paired = FALSE, 
         #detailed=TRUE 
  )%>%
  filter(p.adj < 0.05) %>%  # Keep significant results
  mutate(across(c(p.adj), ~ format(., scientific = FALSE)))  # Format p-values to remove scientific notation
t.test_3groups

##tribble all tests together 3 groups#####################
each.vs.ref_sig <- 
  each.vs.ref_sig2 <- tibble::tribble(
    ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
    "Gerybiniai",   "HGSOC", 0.000231, -4, "PPT2", #stjudent's
    "Gerybiniai",   "HGSOC", 0.03, -5, "EXO1", #stjudent's
    "Gerybiniai",   "HGSOC", 0.003, -2, "RAD50", #stjudent's
    "Gerybiniai",   "Kiti KV", 0.006, -1, "RAD50", #stjudent's
    "Gerybiniai",   "HGSOC", 0.000354 , -3, "CDCA5",  #stjudent's
    "HGSOC",   "Kiti KV", 0.047, -2, "CDCA5",  #stjudent's
    "Gerybiniai",   "HGSOC", 0.008, -4, "ZFPL1", #stjudent's
    "Gerybiniai",   "Kiti KV", 0.008, -3, "ZFPL1",#stjudent's
    "Gerybiniai",   "HGSOC", 0.00000000972, -2, "GRB7", #stjudent's
    "Gerybiniai",   "Kiti KV", 0.0000126, -1, "GRB7", #stjudent's
    "Gerybiniai",   "HGSOC", 0.001, -5, "VPS33B",#stjudent's
    
    
    "Gerybiniai",   "HGSOC", 0.000606, 0, "LUC7L2", #mann
    "Gerybiniai",   "HGSOC", 0.000666, -5, "PKP3", #mann
    "Gerybiniai",   "Kiti KV", 0.009, -4, "PKP3", #mann
    "Gerybiniai",   "HGSOC", 0.0000000591, 4, "TCEAL4",#mann
    "Gerybiniai",   "Kiti KV", 0.000569, 5, "TCEAL4", #mann
    "HGSOC",   "Kiti KV", 0.000569, 3, "TCEAL4", #mann
  )
# 3 digits after the ,
each.vs.ref_sig$p.adj_custom <- ifelse(each.vs.ref_sig$p.adj < 0.001, 
                                          "p < 0.001", 
                                          paste0("p = ", sprintf("%.3f",
                                          each.vs.ref_sig$p.adj)))
##boxplot 3 groups ########################################
#melt one table for expression
Group3_table_full <- melt(OC_full[, colnames(OC_full) %in%
                               c("type", expression )],
                     id.vars="type",
                     measure.vars= expression)
#rename to lt 
Group3_table_full$type
Group3_table_full <- Group3_table_full %>%
  mutate(type = case_when(
    type == "HGSOC" ~ "HGSOC",
    type == "Benign" ~ "Gerybiniai",
    type == "Other" ~ "Kiti KV"
  ))

custom_colors <- c("HGSOC" = "deeppink","Kiti KV" = "lightpink", "Gerybiniai" = "blue") 
OC_plot <- ggplot(Group3_table_full, aes(x=type , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = type )) +
  geom_jitter(aes(color = type ), size=1, alpha=0.5) +
  ylab(label = expression("Santykinė genų raiška, normalizuota pagal  " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(each.vs.ref_sig, label = "p.adj_custom") + #pvalue
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
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x))) #add long "-" signs

OC_plot

##save png 3 groups####################################
png("10genes_boxplot3groups_2025-06-20.png",
    width = 3200, height = 2300, res = 300) # width and height in pixels, resolution in dpi
OC_plot #
dev.off() # Close the PNG device

##FC 3 groups###########################################
expression_df <- OC_full[, c("type", expression, "KN")]
rownames(expression_df) <- expression_df$KN
expression_df <- expression_df[, -12]
#REMOVE NA
expression_df2 <- expression_df %>%
  filter(complete.cases(.))

exp_df2 <- expression_df2 %>%
  melt(id.vars="type",  measure.vars=expression) %>%
  rename(gene = variable) %>%
  rename(expression = value)%>%
  rename(condition = type) %>%
  group_by(condition, gene)%>%
  summarise(mean_expression = mean(expression, drop_na= T)) %>%
  spread(condition, mean_expression)

mean_expression_OC2 <- exp_df2 %>%
  mutate(fold_change_HB = log2(2^(`HGSOC` - `Benign`))) %>% #benign vs hgsoc 
  mutate(fold_change_HO = log2(2^(`HGSOC` - `Other`))) %>% #others vs hgsoc 
  mutate(fold_change_OB = log2(2^(`Other` / `Benign`))) #benign vs others
mean_expression_OC2 

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


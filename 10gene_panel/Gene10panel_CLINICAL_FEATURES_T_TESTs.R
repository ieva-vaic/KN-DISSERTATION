#KN-DISSERTATION project. 10 gene panel - genes selected via statistical analysis of TCGA and GTEx data
#KN 10 gene analysis data (cv 2 %, one threshold, 65 cases) 2025-02-12
#t.tests, with clinical 2025-02-12
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
library(ggpubr)
library(patchwork)
#save wd for plots
setwd("C:/Users/Ieva/rprojects/outputs_all/DISS/")
#data###################################################################
OC_full <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/OC_10_genes_clean_2025_02_14.RDS")
#expression genes list
expression <- c("EXO1", "RAD50","PPT2", "LUC7L2","PKP3", "CDCA5","ZFPL1","VPS33B", "GRB7","TCEAL4")

#Choose clinical features of intrest#########################################
#Age
hist(OC_full$Age)
#STAGE
#chek witch groups to test for stage
fisher.test(table(OC_full$Stage4, OC_full$type)) #not viable because significant difference
fisher.test(table(OC_full$Stage2, OC_full$type)) #not viable because significant difference

#HGSOC only
HGSOC_only <- OC_full %>%
  filter(type == "HGSOC") 
table(HGSOC_only$Stage4, useNA = "a")

table(HGSOC_only$type, HGSOC_only$Stage2) 
table(HGSOC_only$type, HGSOC_only$Stage4) 
#OC only
OC_only <- OC_full %>%
  filter(tumor == "OC") 
table(OC_only$Stage4, useNA = "a")

fisher.test(table(OC_only$type, OC_only$Stage2)) #significant
fisher.test(table(OC_only$type, OC_only$Stage4)) #extra significant

#GRADE
fisher.test(table(OC_full$Grade2, OC_full$type)) #extra significant
#basicaly the same as stage vs grade, but only 3 other is highgrade

#CA125
fisher.test(table(OC_full$CA125_f, OC_full$type)) #extra significant
#only 5 NOT increased: 4 in norm AND 1 IN OTHER
hist(OC_full$CA125, breaks  = 15) #some very high values
sum(!is.na(OC_full$CA125_f)) #55 values
sum(!is.na(OC_full$CA125_post_treatment)) #35 values
hist(OC_full$CA125_post_treatment) #35 values
#in OC only
OC_full <- OC_full %>%
  mutate(CA125_post_op_f = ifelse(CA125_post_treatment > 35, "CA125 increase", "Norm"))
table(OC_full$CA125_post_op_f, useNA = "a") #11 vs 24

#shapiro tests###############################################
colnames(OC_full)
#just the gene expression
shapiro_results <- lapply(OC_full[,3:12], shapiro.test)
shapiro_results #not normal overall TCEAL4, GRB7, PKP3, LUC7L2, 
#age
shapiro.test(OC_full$Age) # p =0.3321
#grade
shapiro_results_grade <- OC_full[, c(18, 3:12 )] %>%
  pivot_longer(cols = -Grade2, names_to = "gene", values_to = "value") %>%
  group_by(Grade2, gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop") %>%
  filter(p_value < 0.05)
shapiro_results_grade #not normal LUC7L2, PKP3, TCEAL4
#CA125 pre op
shapiro_results_CA125_f <- OC_full[, c(17, 3:12 )] %>%
  pivot_longer(cols = -CA125_f, names_to = "gene", values_to = "value") %>%
  group_by(CA125_f, gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop") %>%
  filter(p_value < 0.05)
shapiro_results_CA125_f #not normal GRB7, LUC7L2, PKP3, PPT2, TCEAL4
#Ca125 post op
shapiro_results_CA125_f2 <- OC_full[, c(27, 3:12 )] %>%
  pivot_longer(cols = -CA125_post_op_f, names_to = "gene", values_to = "value") %>%
  group_by(CA125_post_op_f, gene) %>%
  summarise(p_value = shapiro.test(value)$p.value, .groups = "drop") %>%
  filter(p_value < 0.05)
shapiro_results_CA125_f2 #not normal LUC7L2, TCEAL4
#var.tests#############################################################

#GRADE, all OC from this point on
table(OC_full$Grade2)
var_results_grade <- OC_full[, c(18, 3:12 )] %>%
  filter(Grade2 != is.na(Grade2))%>%
  pivot_longer(cols = -Grade2, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = var.test(value[Grade2 == unique(Grade2)[1]], 
                       value[Grade2 == unique(Grade2)[2]])$p.value,
    .groups = "drop"
  ) %>%
  filter(p_value < 0.05)
var_results_grade #luc7l2 

#CA125 pre op
table(OC_full$CA125_f)
var_results_ca125 <- OC_full[, c(17, 3:12 )] %>%
  filter(CA125_f != is.na(CA125_f))%>%
  pivot_longer(cols = -CA125_f, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = var.test(value[CA125_f == unique(CA125_f)[1]], 
                       value[CA125_f == unique(CA125_f)[2]])$p.value,
    .groups = "drop"
  ) %>%
  filter(p_value < 0.05)
var_results_ca125 #all same variances

#CA125 post op
table(OC_full$CA125_post_op_f)
var_results_ca125_post <- OC_full[, c(27, 3:12 )] %>%
  filter(CA125_post_op_f != is.na(CA125_post_op_f))%>%
  pivot_longer(cols = -CA125_post_op_f, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    p_value = var.test(value[CA125_post_op_f == unique(CA125_post_op_f)[1]], 
                       value[CA125_post_op_f == unique(CA125_post_op_f)[2]])$p.value,
    .groups = "drop"
  ) %>%
  filter(p_value < 0.05)
var_results_ca125_post #all same variances


#t tests GRADE#########################################
#not normal LUC7L2, PKP3, TCEAL4
Grade_table <- melt(OC_only[, c(18, 3:12)], id.vars="Grade2", 
                    measure.vars=expression)
wilcox.test_grade <- Grade_table %>%
group_by(variable) %>%
  filter(variable %in% c("LUC7L2", "PKP3", "TCEAL4")) %>%
  summarise(
    p_value = wilcox.test(value ~ Grade2)$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%  # Adjust p-values
  filter(p_adj < 0.1)  # Keep significant results
wilcox.test_grade #TCEAL4    
#normall all else
grade_normal <- setdiff(expression, c("LUC7L2", "PKP3", "TCEAL4"))
stjudents.test_grade <- Grade_table %>%
  filter(variable %in% grade_normal) %>%
  group_by(variable) %>%
  summarise(
    p_value = t.test(value ~ Grade2, var.equal = TRUE)$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%  # Adjust p-values
  filter(p_adj < 0.1)  # Keep significant results
stjudents.test_grade #none

##tribble GRADE ################
each.vs.ref_sig_grade <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "G1",   "G3", 0.0385, 1, "TCEAL4" #wilcox
)

##boxplot GRADE #############################
Grade_table_OC <- melt(OC_only[, colnames(OC_only) %in%
                                       c("Grade2", expression )],
                          id.vars="Grade2",
                          measure.vars= expression) %>%
                          drop_na(Grade2)


custom_colors <- c("G1" = "blue","G3" = "deeppink") 
Grade_plot <- ggplot(Grade_table_OC, aes(x=Grade2 , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = Grade2 )) +
  geom_jitter(aes(color = Grade2 ), size=1, alpha=0.5) +
  ylab(label = expression("Santykinė genų raiška, normalizuota pagal  " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(each.vs.ref_sig_grade, label = "p.adj") + #pvalue
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
  scale_y_continuous(labels = function(x) gsub("-", "\u2212", x))


Grade_plot

##save plot####################################################
png("10gene_grade_20250623.png"
    , width = 3000, height = 2300, res = 300) # width and height in pixels, resolution in dpi
Grade_plot 
dev.off() # Close the PNG device


##FC GRADE ##################
expression_df_GRADE <- OC_only[, c("Grade2", expression, "KN")]
rownames(expression_df_GRADE) <- expression_df_GRADE$KN
expression_df_GRADE <- expression_df_GRADE[, -12]
#REMOVE NA
#expression_df_GRADE <- expression_df_GRADE %>%
#  filter(complete.cases(.))

exp_df_GRADE <- expression_df_GRADE %>%
  melt(id.vars="Grade2",  measure.vars=expression) %>%
  rename(gene = variable) %>%
  rename(expression = value)%>%
  rename(condition = Grade2) %>%
  group_by(condition, gene)%>%
  summarise(mean_expression = mean(expression, drop_na= T)) %>%
  spread(condition, mean_expression)

mean_expression_GRADE <- exp_df_GRADE %>%
  mutate(fold_change_3_1 = log2( 2^( `G3` -`G1`) ) )
mean_expression_GRADE 
#coreliation with CA125 at diagnosis#################
shapiro_results_CA125 <- OC_full[, c(20, 3:12 )] %>%
  map(~ shapiro.test(.x)$p.value) 
shapiro_results_CA125 # all not normal
CA125_test <- OC_full[, colnames(OC_full) %in% c(expression, "CA125")]

results_nn1 <- lapply(OC_full[, colnames(OC_full) %in% expression], 
                      function(x) cor.test(x, OC_full$CA125, method = "spearman"))
results_nn1


# Spearman results (nn)
df_spearman1 <- data.frame(
  variable = names(results_nn1),
  R = sapply(results_nn1, function(x) if (is.list(x)) x$estimate else NA),
  p_value = sapply(results_nn1, function(x) if (is.list(x)) x$p.value else NA),
  method = "Spearman"
)


##plot ca125 corr #######################
# Create the plot list
plot_list1 <- lapply(df_spearman1$variable, function(gene) {
  method <- df_spearman1$method[df_spearman1$variable == gene]
  r_val <- round(df_spearman1$R[df_spearman1$variable == gene], 2)
  p_val <- signif(df_spearman1$p_value[df_spearman1$variable == gene], 2)
  
  ggplot(CA125_test, aes_string(x = "CA125", y = gene)) +
    geom_point(color = "steelblue") +
    geom_smooth(method = "lm", se = FALSE, color = "darkred") +
    labs(
      #title = bquote("CA125 koncentracija vs"~italic(.(gene))),
      subtitle = paste0( 
        "  R = ", r_val, 
        " | p = ", p_val),
      x = "CA125 U/ml",
      y = bquote("Raiška"~italic(.(gene)))
    ) +
    theme_minimal() +
    theme(plot.subtitle = element_text(size = 10, face = "italic", color = "gray30"))
})
# Show plots one by one (optional)
for (p in plot_list1) print(p)

# Combine into one figure
combined_plot1 <- wrap_plots(plot_list1, ncol = 2) +
  plot_annotation(title = "Koreliacija tarp CA125 koncentracijos ir genų raiškos")
combined_plot1



png("10gene_CA125_20251222.png"
    , width = 3000, height = 5000, res = 300) # width and height in pixels, resolution in dpi
combined_plot1
dev.off() # 


#t tests CA125 pre op ####################
#not normal GRB7, LUC7L2, PKP3, PPT2, TCEAL4
ca125_table <- melt(OC_full[, c(17, 3:12)], id.vars="CA125_f", 
                    measure.vars=expression)
wilcox.test_ca125 <- ca125_table %>%
  group_by(variable) %>%
  filter(variable %in% c("GRB7", "LUC7L2", "PKP3", "PPT2", "TCEAL4")) %>%
  summarise(
    p_value = wilcox.test(value ~ CA125_f)$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%  # Adjust p-values
  filter(p_adj < 0.05)  # Keep significant results
wilcox.test_ca125 #TCEAL4 LUC7L2      
#normall all else
ca125_normal <- setdiff(expression, c("GRB7", "LUC7L2", "PKP3","PPT2", "TCEAL4"))
stjudents.test_ca125 <- ca125_table %>%
  filter(variable %in% ca125_normal) %>%
  group_by(variable) %>%
  summarise(
    p_value = t.test(value ~ CA125_f, var.equal = TRUE)$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%  # Adjust p-values
  filter(p_adj < 0.05) #none
stjudents.test_ca125
##tribble CA125 pre op ################
each.vs.ref_sig_ca125 <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "CA125 padidėjimas",   "Norma" , 0.0499, 2, "TCEAL4", #wilcox
  "CA125 padidėjimas",   "Norma" , 0.0259, 0, "LUC7L2" #wilcox
)

CA125_table_OC <- melt(OC_full[, colnames(OC_full) %in%
                                 c("CA125_f", expression )],
                       id.vars="CA125_f",
                       measure.vars= expression) %>%
  drop_na(CA125_f) %>%
  mutate(CA125_f = case_when(
    CA125_f == "CA125 increase" ~ "CA125 padidėjimas",
    CA125_f == "Norm" ~ "Norma"
  ))

##boxplot CA125 pre op ####################
custom_colors <- c("CA125 padidėjimas" = "blue","Norma" = "lightpink") 
Ca125_plot <- ggplot(CA125_table_OC, aes(x=CA125_f , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = CA125_f )) +
  geom_jitter(aes(color = CA125_f ), size=1, alpha=0.5) +
  ylab(label = expression("Santykinė genų raiška, normalizuota " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(each.vs.ref_sig_ca125, label = "p.adj") + #pvalue
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
  scale_y_continuous(labels = function(x) gsub("-", "\u2212", x))

Ca125_plot

##save plot ca125 pre op######################################
png("10_gene_boxplot_ca125_2025_06_23.png"
    , width = 3500, height = 2300, res = 300) # width and height in pixels, resolution in dpi
Ca125_plot #
dev.off() # Close the PNG device

##FC CA125 pre op ########################
expression_df_ca125 <- OC_full[, c("CA125_f", expression, "KN")]
rownames(expression_df_ca125) <- expression_df_ca125$KN
expression_df_ca125 <- expression_df_ca125[, -12]
#REMOVE NA
expression_df_ca125 <- expression_df_ca125 %>%
  filter(complete.cases(.))

exp_df_ca125 <- expression_df_ca125 %>%
  melt(id.vars="CA125_f",  measure.vars=expression) %>%
  rename(gene = variable) %>%
  rename(expression = value)%>%
  rename(condition = CA125_f) %>%
  group_by(condition, gene)%>%
  summarise(mean_expression = mean(expression, drop_na= T)) %>%
  spread(condition, mean_expression)

mean_expression_ca125 <- exp_df_ca125 %>%
  mutate(fold_change_ca125 = log2(2^(`CA125 increase`- `Norm`) ))
mean_expression_ca125 

#t tests CA125 post op ###############################
#not normal LUC7L2, TCEAL4
ca125_table2 <- melt(OC_full[, c(27, 3:12)], id.vars="CA125_post_op_f", 
                    measure.vars=expression)
wilcox.test_ca125_2 <- ca125_table2 %>%
  group_by(variable) %>%
  filter(variable %in% c( "LUC7L2", "TCEAL4")) %>%
  summarise(
    p_value = wilcox.test(value ~ CA125_post_op_f)$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%  # Adjust p-values
  filter(p_adj < 0.1)  # Keep significant results
wilcox.test_ca125_2 # LUC7L2      
#normall all else
ca125_normal2 <- setdiff(expression, c( "LUC7L2", "TCEAL4"))
stjudents.test_ca125_2 <- ca125_table2 %>%
  filter(variable %in% ca125_normal2) %>%
  group_by(variable) %>%
  summarise(
    p_value = t.test(value ~ CA125_post_op_f, var.equal = TRUE)$p.value
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%  # Adjust p-values
  filter(p_adj < 0.1) #none
stjudents.test_ca125_2 

##tribble CA125 post op ################
each.vs.ref_sig_ca125_2 <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "CA125 increase",   "Norm" , 0.0565, 0, "LUC7L2" #wilcox
)

CA125_table_OC2 <- melt(OC_full[, colnames(OC_full) %in%
                                 c("CA125_post_op_f", expression )],
                       id.vars="CA125_post_op_f",
                       measure.vars= expression) %>%
  drop_na(CA125_post_op_f)

##boxplot CA125 pre op ####################
custom_colors <- c("CA125 increase" = "darkviolet","Norm" = "darkgreen") 
Ca125_plot2 <- ggplot(CA125_table_OC2, aes(x=CA125_post_op_f , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = CA125_post_op_f )) +
  geom_jitter(aes(color = CA125_post_op_f ), size=1, alpha=0.5) +
  ylab(label = expression("Relative expression, normalized to  " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(each.vs.ref_sig_ca125_2, label = "p.adj") + #pvalue
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
  ggtitle("Gene expression correlation with CA125 post op. levels")+
  scale_y_continuous(labels = function(x) gsub("-", "\u2212", x))

Ca125_plot2

##save plot#############################
png("10gene_boxplot_ca125_post_op_2025_02_17.png"
    , width = 3000, height = 2300, res = 300) # width and height in pixels, resolution in dpi
Ca125_plot2 #
dev.off() # Close the PNG device

##FC CA125 post op ########################
expression_df_ca125_2 <- OC_full[, c("CA125_post_op_f", expression, "KN")]
rownames(expression_df_ca125_2) <- expression_df_ca125_2$KN
expression_df_ca125_2 <- expression_df_ca125_2[, -12]
#REMOVE NA
expression_df_ca125_2 <- expression_df_ca125_2 %>%
  filter(complete.cases(.))

exp_df_ca125_2 <- expression_df_ca125_2 %>%
  melt(id.vars="CA125_post_op_f",  measure.vars=expression) %>%
  rename(gene = variable) %>%
  rename(expression = value)%>%
  rename(condition = CA125_post_op_f) %>%
  group_by(condition, gene)%>%
  summarise(mean_expression = mean(expression, drop_na= T)) %>%
  spread(condition, mean_expression)

mean_expression_ca125_2 <- exp_df_ca125_2 %>%
  mutate(fold_change_ca125_post_op = `CA125 increase` / `Norm`) 
mean_expression_ca125_2

#coreliation with age#################
table(OC_full$Age, useNA ="a")
age_test <- OC_full[, colnames(OC_full) %in% c(expression, "Age")]
#not normal TCEAL4, GRB7, PKP3, LUC7L2
nn <- c("TCEAL4", "GRB7", "PKP3", "LUC7L2")
on <- expression[!expression %in% nn]
results_nn <- lapply(OC_full[, colnames(OC_full) %in% nn], 
                     function(x) cor.test(x, OC_full$Age, method = "spearman"))
results_nn

results_n <- lapply(OC_full[, colnames(OC_full) %in% on], 
                     function(x) cor.test(x, OC_full$Age, method = "pearson"))
results_n

#make dfs
df_pearson <- data.frame(
  variable = names(results_n),
  R = sapply(results_n, function(x) if (is.list(x)) x$estimate else NA),
  p_value = sapply(results_n, function(x) if (is.list(x)) x$p.value else NA),
  method = "Pearson"
)
# Spearman results (nn)
df_spearman <- data.frame(
  variable = names(results_nn),
  R = sapply(results_nn, function(x) if (is.list(x)) x$estimate else NA),
  p_value = sapply(results_nn, function(x) if (is.list(x)) x$p.value else NA),
  method = "Spearman"
)

# Combine
df_all <- rbind(df_pearson, df_spearman)

##plot age #######################
# Create the plot list
plot_list <- lapply(df_all$variable, function(gene) {
  method <- df_all$method[df_all$variable == gene]
  r_val <- round(df_all$R[df_all$variable == gene], 2)
  p_val <- signif(df_all$p_value[df_all$variable == gene], 2)
  
  ggplot(age_test, aes_string(x = "Age", y = gene)) +
    geom_point(color = "steelblue") +
    geom_smooth(method = "lm", se = FALSE, color = "darkred") +
    labs(
      title = bquote("Amžius vs"~italic(.(gene))),
      subtitle = paste0( 
                        "  R = ", r_val, 
                        " | p = ", p_val),
      x = "Amžius metais",
      y = bquote("Raiška"~italic(.(gene)))
    ) +
    theme_minimal() +
    theme(plot.subtitle = element_text(size = 10, face = "italic", color = "gray30"))
})
# Show plots one by one (optional)
for (p in plot_list) print(p)

# Combine into one figure
combined_plot <- wrap_plots(plot_list, ncol = 2) +
  plot_annotation(title = "Ryšys tarp amžiaus ir genų raiškos")
combined_plot
##save age plot to PNG####################
ggsave(
  filename = "10_gene_age10genes20250623.png",
  plot = combined_plot,
  width = 8,       # adjust width as needed
  height = 11,       # adjust height as needed
  dpi = 300         # high resolution
)
#ENGLISH PLOTS #####################################
#EN CA125 plot################################################
##tribble CA125 pre op\
each.vs.ref_sig_ca125EN <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "CA125 increase",   "Norm" , 0.0499, 2, "TCEAL4", #wilcox
  "CA125 increase",   "Norm" , 0.0259, 0, "LUC7L2" #wilcox
)

CA125_table_OCEN <- melt(OC_full[, colnames(OC_full) %in%
                                   c("CA125_f", expression )],
                         id.vars="CA125_f",
                         measure.vars= expression) %>%
  drop_na(CA125_f) 
##boxplot CA125 pre op ####################
custom_colors <- c("CA125 increase" = "blue","Norm" = "lightpink") 
Ca125_plot <- ggplot(CA125_table_OCEN, aes(x=CA125_f , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = CA125_f )) +
  geom_jitter(aes(color = CA125_f ), size=1, alpha=0.5) +
  ylab(label = expression("Relative gene expression, normalised to " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(each.vs.ref_sig_ca125EN, label = "p.adj") + #pvalue
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
  scale_y_continuous(labels = function(x) gsub("-", "\u2212", x))

Ca125_plot

##save plot ca125 pre op######################################
png("10_gene_boxplot_ca125_2025_0927EN.png"
    , width = 3500, height = 2300, res = 300) # width and height in pixels, resolution in dpi
Ca125_plot #
dev.off() # Close the PNG device

table(OC_full$CA125_f, useNA = "a")

#EN coreliation with age#################
table(OC_full$Age, useNA ="a")
age_test <- OC_full[, colnames(OC_full) %in% c(expression, "Age")]
#not normal TCEAL4, GRB7, PKP3, LUC7L2
nn <- c("TCEAL4", "GRB7", "PKP3", "LUC7L2")
on <- expression[!expression %in% nn]
results_nn <- lapply(OC_full[, colnames(OC_full) %in% nn], 
                     function(x) cor.test(x, OC_full$Age, method = "spearman"))
results_nn

results_n <- lapply(OC_full[, colnames(OC_full) %in% on], 
                    function(x) cor.test(x, OC_full$Age, method = "pearson"))
results_n

#make dfs
df_pearson <- data.frame(
  variable = names(results_n),
  R = sapply(results_n, function(x) if (is.list(x)) x$estimate else NA),
  p_value = sapply(results_n, function(x) if (is.list(x)) x$p.value else NA),
  method = "Pearson"
)
# Spearman results (nn)
df_spearman <- data.frame(
  variable = names(results_nn),
  R = sapply(results_nn, function(x) if (is.list(x)) x$estimate else NA),
  p_value = sapply(results_nn, function(x) if (is.list(x)) x$p.value else NA),
  method = "Spearman"
)

# Combine
df_all <- rbind(df_pearson, df_spearman)

##plot age #######################
# Create the plot list
plot_list <- lapply(df_all$variable, function(gene) {
  method <- df_all$method[df_all$variable == gene]
  r_val <- round(df_all$R[df_all$variable == gene], 2)
  p_val <- signif(df_all$p_value[df_all$variable == gene], 2)
  
  ggplot(age_test, aes_string(x = "Age", y = gene)) +
    geom_point(color = "steelblue") +
    geom_smooth(method = "lm", se = FALSE, color = "darkred") +
    labs(
      #title = bquote("Age vs"~italic(.(gene))),
      subtitle = paste0( 
        "  R = ", r_val, 
        " | p = ", p_val),
      x = "Age, years ",
      y = bquote("Gene Expression"~italic(.(gene)))
    ) +
    theme_minimal() +
    theme(plot.subtitle = element_text(size = 10, face = "italic", color = "gray30"))
})
# Show plots one by one (optional)
for (p in plot_list) print(p)

# Combine into one figure
combined_plot <- wrap_plots(plot_list, ncol = 2) 
combined_plot
#save in english
ggsave(
  filename = "10_gene_age10genes202501023EN.png",
  plot = combined_plot,
  width = 8,       # adjust width as needed
  height = 11,       # adjust height as needed
  dpi = 500         # high resolution
)
#EN grade plot #############################

Grade_plotEN <- ggplot(Grade_table_OC, aes(x=Grade2 , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = Grade2 )) +
  geom_jitter(aes(color = Grade2 ), size=1, alpha=0.5) +
  ylab(label = expression("Relative expression, normalized to  " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(each.vs.ref_sig_grade, label = "p.adj") + #pvalue
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
  scale_y_continuous(labels = function(x) gsub("-", "\u2212", x))


Grade_plotEN

png("10gene_gradeEN_20251217.png"
    , width = 3000, height = 2300, res = 300) # width and height in pixels, resolution in dpi
Grade_plotEN
dev.off() # Cl

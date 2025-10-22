#KN-DISSERTATION project. MET_EXPRS - methylation, expression, mutation data
#EXPRESSION / BOXPLOT
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
#DATA - FROM THE PUBLICATION MET_EXPRS
KN_data <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/KN_data1114_essential.rds")
colnames(KN_data)
#biomarker groups
#expression
raiska <- colnames(KN_data[18:27])
#methylation
metilinimas <- colnames(KN_data[28:31])
biomarkers <- c(raiska, metilinimas)
#change tumor
KN_data$tumor <- recode(KN_data$tumor, OvCa = "OC", Benign ="Benign")

#EXPR BOXPLOT OVca##############################################
#melt table for expression
tumor_table <- melt(KN_data, id.vars="tumor",  measure.vars=raiska)
tumor_table <- tumor_table %>%
  mutate(tumor = recode(tumor, "OC" = "KV", "Benign" = 
                          "Gerybiniai")) 
#STATISTICS: wilcox
stat.test_tumor <- tumor_table %>%
  group_by(variable) %>%
  pairwise_wilcox_test(value ~ tumor,
                       p.adjust.method = "BH",
                       pool.sd = TRUE) #equal variance
print(stat.test_tumor) #results

each.vs.ref_sig_tumor <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "Gerybiniai",   "KV", 0.013, -5, "NOTCH1", 
  "Gerybiniai",   "KV", 0.000238,  -2, "NOTCH2",
  "Gerybiniai",   "KV", 0.026,  -2, "NOTCH3",
  "Gerybiniai",   "KV", 0.000426, -3, "NOTCH4",
  "Gerybiniai",   "KV", 0.059, -2, "ARID1A",
  "Gerybiniai",   "KV", 0.0000201, 0, "CTNNB1",
  "Gerybiniai",   "KV", 0.000176, -4, "FBXW7",
  #"Gerybiniai",   "KV", 0.842, -4, "JAG2",
  "Gerybiniai",   "KV", 0.00044, -5, "DLL1",
  "Gerybiniai",   "KV", 0.002, -2, "HES1",
)
# iki 3 skaiciu po kablelio numazinimas
each.vs.ref_sig_tumor$p.adj_custom <- ifelse(each.vs.ref_sig_tumor$p.adj < 0.001, 
                                             "p < 0.001", 
                                             paste0("p = ", sprintf("%.3f",
                                                                    each.vs.ref_sig_tumor$p.adj)))
#plot
custom_colors <- c("KV" = "#cf5784", "Gerybiniai" = "#929cef")
tumor_plot <- ggplot(tumor_table, aes(x=tumor, y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = tumor)) +
  geom_jitter(aes(color = tumor), size=2, alpha=0.5) +
  ylab(label = expression("Santykinė raiška, normalizuota pagal  " * italic("GAPDH"))) + 
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

#SAVE PNG
#set directory for saving
setwd("C:/Users/Ieva/rprojects/outputs_all/DISS/")
#save
png("met_exprs_boxplot_OVca_output1020.png", width = 3000, height = 2700, res = 400) # width and height in pixels, resolution in dpi
tumor_plot #
dev.off() # Close the PNG device

#fc all cases
raiska_df <- KN_data[, c("tumor", raiska, "patient_id_aud")]
rownames(raiska_df) <- raiska_df$patient_id_aud
raiska_df <- raiska_df[, -12]
#exp df
exp_df <- raiska_df %>%
  melt(id.vars="tumor",  measure.vars=raiska) %>%
  rename(gene = variable) %>%
  rename(expression = value)%>%
  rename(condition = tumor) %>%
  group_by(condition, gene)%>%
  summarise(mean_expression = mean(expression, drop_na= T)) %>%
  spread(condition, mean_expression)

mean_expression_tumor2 <- exp_df %>%
  mutate(fold_change_HB = log2(2^`OC` / 2^`Benign`))
mean_expression_tumor2#biger values

#fc 64 cases, because DLL1 has a missing data point
raiska_df <- KN_data[, c("tumor", raiska, "patient_id_aud")]
#remove one dll1 value
raiska_df <-  raiska_df %>% filter(!is.na(DLL1))
dim(raiska_df) #64 cases left
#continue
rownames(raiska_df) <- raiska_df$patient_id_aud
raiska_df <- raiska_df[, -12]
#exp df
exp_df <- raiska_df %>%
  melt(id.vars="tumor",  measure.vars=raiska) %>%
  rename(gene = variable) %>%
  rename(expression = value)%>%
  rename(condition = tumor) %>%
  group_by(condition, gene)%>%
  summarise(mean_expression = mean(expression, drop_na= T)) %>%
  spread(condition, mean_expression)

mean_expression_tumor2 <- exp_df %>%
  mutate(fold_change_HB = log2(2^`OC` / 2^`Benign`))
mean_expression_tumor2#biger values

#EXPR BOXPLOT 3 groups ########################
#melt table for expression
HOB_table <- melt(KN_data, id.vars="Grupė_Ieva",  measure.vars=raiska)

HOB_table <- HOB_table %>%
  mutate(tumor = recode(Grupė_Ieva, "Other" = "Kiti KV", "Benign" = 
                          "Gerybiniai")) 
#wilcox
stat.test_HOB <- HOB_table %>%
  group_by(variable) %>%
  pairwise_wilcox_test(value ~ tumor,
                       #ref.group = "Benign" , #only if one group is needed
                       p.adjust.method = "BH")%>%
  filter(group1 == "Gerybiniai")%>%
  filter(group2 == "Kiti KV")

stat.test_HOB
#fill in the values
each.vs.ref_sig_HOB <- 
  each.vs.ref_sig2 <- tibble::tribble(
    ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
    "Gerybiniai",   "HGSOC",    0.036, -5, "NOTCH1",
    "Gerybiniai",   "HGSOC",    0.000333,  -2, "NOTCH2",
    #"Gerybiniai",   "HGSOC",    0.068, -2, "NOTCH3", 
    "Gerybiniai",   "HGSOC", 0.0000852, -3, "NOTCH4", 
    #  "Gerybiniai",   "HGSOC", 0.184, -2, "ARID1A", 
    "Gerybiniai",   "HGSOC", 0.000000522, -0.5, "CTNNB1",
    "Gerybiniai",   "HGSOC",0.00000357, -4, "FBXW7",
    #  "Gerybiniai",   "HGSOC", 1, -4, "JAG2",
    "Gerybiniai",   "HGSOC", 0.0000852, -4, "DLL1",
    "Gerybiniai",   "HGSOC",0.00000357, -3, "HES1",
    
    #"Kiti KV",   "HGSOC",    0.786, -3, "NOTCH1",
    #"Kiti KV",   "HGSOC",    1,  -3, "NOTCH2",
    #"Kiti KV",   "HGSOC",    0.297, -2.5, "NOTCH3",
    #  "Kiti KV",   "HGSOC", 0.178, -2.5, "NOTCH4",
    # "Kiti KV",   "HGSOC", 0.715, -2.5, "ARID1A",  
    #"Kiti KV",   "HGSOC", 0.052, -0.5, "CTNNB1",
    #"Kiti KV",   "HGSOC", 0.062, -3.5, "FBXW7", #lost, BUT NONAJUDTED SIGNIFICANT
    #  "Kiti KV",   "HGSOC", 1, -3.5, "JAG2",
    "Kiti KV",   "HGSOC", 0.02, -5, "DLL1",
    "Kiti KV",   "HGSOC",0.000214, -2, "HES1",
    
    #"Kiti KV",   "Gerybiniai", 0.068 , -4, "NOTCH1", #non-adjusted significnt
    "Kiti KV",   "Gerybiniai", 0.002,  -1, "NOTCH2",
    #"Kiti KV",   "Gerybiniai", 0.057, -1, "NOTCH3", 
    #"Kiti KV",   "Gerybiniai", 0.068 , -2, "NOTCH4", #non-adjusted significnt
    #"Kiti KV",   "Gerybiniai", 0.212, -1, "ARID1A",  
    "Kiti KV",   "Gerybiniai", 0.002, 0, "CTNNB1",
    #"Kiti KV",   "Gerybiniai",0.083, -3, "FBXW7",
    #"Kiti KV",   "Gerybiniai",  1, -3, "JAG2",
    #"Kiti KV",   "Gerybiniai", 0.06  , -3, "DLL1",
    #"Kiti KV",   "Gerybiniai", 0.877 , -2, "HES1"
  )

# iki 3 skaiciu po kablelio numazinimas
each.vs.ref_sig_HOB$p.adj_custom <- ifelse(each.vs.ref_sig_HOB$p.adj < 0.001, 
                                           "p < 0.001", 
                                           paste0("p = ", sprintf("%.3f", each.vs.ref_sig_HOB$p.adj)))
custom_colors <- c("HGSOC" = "#cf5784","Kiti KV" = "pink", "Gerybiniai" = "#929cef") 
HOB_x <- ggplot(HOB_table, aes(x=tumor , y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = tumor )) +
  geom_jitter(aes(color = tumor ), size=2, alpha=0.5) +
  ylab(label = expression("Santykinė raiška, normalizuota pagal " * italic("  GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(each.vs.ref_sig_HOB, label = "p.adj_custom") + #pvalue
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 11, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x))) #add long "-" signs

HOB_x

#SAVE PNG
png("met_exprs_boxplot_3groups_output1020.png",
    width = 4200, height = 3000, res = 400) # width and height in pixels, resolution in dpi
HOB_x #
dev.off() # Close the PNG device

#fc all cases 3 groups
raiska_df2 <- KN_data[, c("Grupė_Ieva", raiska, "patient_id_aud")]
rownames(raiska_df2) <- raiska_df2$patient_id_aud
raiska_df2 <- raiska_df2[, -12]
#exp df
exp_df2 <- raiska_df2 %>%
  melt(id.vars="Grupė_Ieva",  measure.vars=raiska) %>%
  rename(gene = variable) %>%
  rename(expression = value)%>%
  rename(condition = Grupė_Ieva) %>%
  group_by(condition, gene)%>%
  summarise(mean_expression = mean(expression, drop_na= T)) %>%
  spread(condition, mean_expression)

mean_expression_HOB <- exp_df2 %>%
  mutate(fold_change_HB = log2(2^(`HGSOC` - `Benign`)) ) %>% #benign vs hgsoc
  mutate(fold_change_OB = log2(2^(`Other` - `Benign`) ))%>% #benign vs other
  mutate(fold_change_HO = log2(2^(`HGSOC` - `Other`) ) )
mean_expression_HOB

#fc 64 cases 3 groups - for DLL1
raiska_df2 <- KN_data[, c("Grupė_Ieva", raiska, "patient_id_aud")]
rownames(raiska_df2) <- raiska_df2$patient_id_aud
raiska_df2 <- raiska_df2[, -12]
#remove one dll1 value
raiska_df2 <-  raiska_df2 %>% filter(!is.na(DLL1))
dim(raiska_df2) #64 cases left
#exp df
exp_df2 <- raiska_df2 %>%
  melt(id.vars="Grupė_Ieva",  measure.vars=raiska) %>%
  rename(gene = variable) %>%
  rename(expression = value)%>%
  rename(condition = Grupė_Ieva) %>%
  group_by(condition, gene)%>%
  summarise(mean_expression = mean(expression, drop_na= T)) %>%
  spread(condition, mean_expression)

mean_expression_HOB <- exp_df2 %>%
  mutate(fold_change_HB = log2(2^(`HGSOC` - `Benign`)) ) %>% #benign vs hgsoc
  mutate(fold_change_OB = log2(2^(`Other` - `Benign`) ))%>% #benign vs other
  mutate(fold_change_HO = log2(2^(`HGSOC` - `Other`) ) )
mean_expression_HOB

#EXPR BOXPLOT GRADE##############################################
table(KN_data$Grade2, useNA = "a")
table(KN_data$Grade2, KN_data$Grupė_Ieva, useNA = "a")
Grade_data <- KN_data[!is.na(KN_data$Grade2),]

#melt table for expression
grade_table <- melt(Grade_data, id.vars="Grade2",  measure.vars=raiska)
#t test
stat.test_grade2 <- grade_table %>%
  group_by(variable) %>%
  pairwise_wilcox_test(value ~ Grade2, p.adjust.method = "BH")
stat.test_grade2

each.vs.ref_sig_grade <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  #"G1",   "G3", 0.573, -5, "NOTCH1",
  #"G1",   "G3", 0.988,  -2, "NOTCH2",
  #"G1",   "G3", 0.637, -2, "NOTCH3",
  #"G1",   "G3", 0.401, -3, "NOTCH4",
  #"G1",   "G3", 0.573, -3, "ARID1A",
  "G1",   "G3", 0.05, -1.5, "CTNNB1",
  "G1",   "G3", 0.054, -4, "FBXW7", #not significant anymore
  #"G1",   "G3", 0.335, -3, "JAG2",
  "G1",   "G3", 0.054, -6, "DLL1",
  "G1",   "G3", 0.000364, -2, "HES1",
  
)
# iki 3 skaiciu po kablelio numazinimas
each.vs.ref_sig_grade$p.adj_custom <- ifelse(each.vs.ref_sig_grade$p.adj < 0.001, 
                                             "p < 0.001", 
                                             paste0("p = ", sprintf("%.3f",
                                                                    each.vs.ref_sig_grade$p.adj)))
#plot
custom_colors <- c("G3" = "#CF5784", "G1" = "#929cef", "NA" = "grey")
grade_plot <- ggplot(grade_table, aes(x=Grade2, y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = Grade2)) +
  geom_jitter(aes(color = Grade2), size=2, alpha=0.5) +
  ylab(label = expression("Santykinė raiška, normalizuota pagal  " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(each.vs.ref_sig_grade, label = "p.adj_custom") + #pvalue
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
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x))) #add long "-" signs

grade_plot

#SAVE PNG AT THE SAME SIZE

png("metexprs_BOXPLOT_grade_output1020.png",
    width = 3000, height = 2600, res = 400) # width and height in pixels, resolution in dpi
grade_plot #
dev.off() # Close the PNG device

#fc all cases grade
raiska_df3 <- KN_data[, c("Grade2", raiska, "patient_id_aud")]
rownames(raiska_df3) <- raiska_df3$patient_id_aud
raiska_df3 <- raiska_df3[, -12]
#exp df
exp_df3 <- raiska_df3 %>%
  melt(id.vars="Grade2",  measure.vars=raiska) %>%
  rename(gene = variable) %>%
  rename(expression = value)%>%
  rename(condition = Grade2) %>%
  group_by(condition, gene)%>%
  summarise(mean_expression = mean(expression, drop_na= T)) %>%
  spread(condition, mean_expression)

mean_expression_GRADE <- exp_df3 %>%
  mutate(fold_change_G3_1 = log2(2^(`G3` - `G1`) ))  #G3 vs G1
mean_expression_GRADE

#EXPR BOXPLOT STAGE IN HGSOC ONLY ###############
KN_HGSOC <- KN_data[KN_data$Grupė_Ieva == "HGSOC", ] #42 people
KN_HGSOC$Grupė_Ieva <- droplevels(KN_HGSOC$Grupė_Ieva)
table(KN_HGSOC$Grupė_Ieva) #42 left
table(KN_HGSOC$Stage4, useNA = "a")

stage4_table <- melt(KN_HGSOC, id.vars = "Stage4", measure.vars= raiska)
#wilcox test
stat.test_stage4_w <- stage4_table %>%
  group_by(variable) %>%
  pairwise_wilcox_test(value ~ Stage4,
                       p.adjust.method = "BH")%>%
  filter(p.adj <0.1)

stat.test_stage4_w

each.vs.ref_sig_stage4 <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "II",   "IV", 0.026   , -3, "NOTCH2",
  "II",   "IV", 0.046     , 0, "CTNNB1",
  "III",   "IV", 0.042    , -1, "CTNNB1",
  "II",   "IV", 0.046     , -3, "FBXW7",
  "III",   "IV", 0.007     , -4, "FBXW7",
  "II",   "III", 0.075     , -1, "ARID1A",
  "III",   "IV", 0.053     , -2, "ARID1A",
  "II",   "IV", 0.053     , -4, "NOTCH1",
  "II",   "III", 0.083     , -3.5, "NOTCH2",
  "III",   "IV", 0.086     , -2, "NOTCH2",
)
# iki 3 skaiciu po kablelio numazinimas
each.vs.ref_sig_stage4$p.adj_custom <- ifelse(each.vs.ref_sig_stage4$p.adj < 0.001, 
                                              "p < 0.001", 
                                              paste0("p = ", sprintf("%.3f",
                                                                     each.vs.ref_sig_stage4$p.adj)))
#plot
custom_colors <- c("IV" = "#cf5784","III" = "pink", "II" = "#929cef")
stage4_plot <- ggplot(stage4_table, aes(x=Stage4, y=value, fill = variable)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = Stage4)) +
  geom_jitter(aes(color = Stage4), size=2, alpha=0.5) +
  ylab(label = expression("Santykinė raiška, normalizuota pagal " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(each.vs.ref_sig_stage4, label = "p.adj_custom") + #pvalue
  theme_minimal()+
  xlab("") +
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  #labs(x=NULL, title ="Gene expression in relation to FIGO stage in HGSOC")+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x))) #add long "-" signs

stage4_plot #

#SAVE PNG
png("met_exprs_BOXPLOT_HGSOC_stage_output1020.png", 
    width = 3000, height = 2700, res = 400) # width and height in pixels, resolution in dpi
stage4_plot #
dev.off() # Close the PNG device

#fc all cases, stage
raiska_df4 <- KN_HGSOC[, c("Stage4", raiska, "patient_id_aud")]
rownames(raiska_df4) <- raiska_df4$patient_id_aud
raiska_df4 <- raiska_df4[, -12]
#exp df
exp_df4 <- raiska_df4 %>%
  melt(id.vars="Stage4",  measure.vars=raiska) %>%
  rename(gene = variable) %>%
  rename(expression = value)%>%
  rename(condition = Stage4) %>%
  group_by(condition, gene)%>%
  summarise(mean_expression = mean(expression, drop_na= T)) %>%
  spread(condition, mean_expression)

mean_expression_hgsoc_stage<- exp_df4 %>%
  mutate(fold_change_III_II = log2(2^(`III` - `II`) ))%>%
  mutate(fold_change_IV_II = log2(2^(`IV` - `II`) )) %>%
  mutate(fold_change_IV_III =log2(2^( `IV` - `III`) ))
mean_expression_hgsoc_stage

#EXPR BOXPLOT CA125 ###########################
table(KN_data$CA125_f, useNA = "a")#factor of ca125

ca_table <- melt(KN_data, id.vars="CA125_f",  measure.vars=raiska)
ca_table$CA125_f <- factor(ca_table$CA125_f, levels = c("Norm", "CA125 increase") )
ca_table <- ca_table %>%
  drop_na() %>%
  mutate(CA125_f = recode(CA125_f, "Norm" = "Norma", "CA125 increase" = 
                          "CA125 padidėjimas")) 

stat.test_ca125_w<- ca_table %>%
  group_by(variable) %>%
  pairwise_wilcox_test(value ~ CA125_f, p.adjust.method = "BH")
stat.test_ca125_w #not normal 

each.vs.ref_sig_ca <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  # "Norma",   "CA125 padidėjimas",0.193, -5, "NOTCH1",
  "Norma",   "CA125 padidėjimas",0.042,  -1, "NOTCH2",
  "Norma",   "CA125 padidėjimas",0.0494, -2, "NOTCH3",
  "Norma",   "CA125 padidėjimas",0.016, -4, "NOTCH4", 
  # "Norma",   "CA125 padidėjimas",0.372, -2, "ARID1A",
  "Norma",   "CA125 padidėjimas",0.029, -0.5, "CTNNB1", #not normal
  "Norma",   "CA125 padidėjimas",0.023, -4, "FBXW7",
  # "Norma",   "CA125 padidėjimas",0.588, -4, "JAG2",
  # "Norma",   "CA125 padidėjimas",0.123 , -6, "DLL1",
  "Norma",   "CA125 padidėjimas",0.031, -1.5, "HES1"
)

# 3 digtits over the dot
each.vs.ref_sig_ca$p.adj_custom <- ifelse(each.vs.ref_sig_ca$p.adj < 0.001, 
                                          "p < 0.001", 
                                          paste0("p = ", sprintf("%.3f", each.vs.ref_sig_ca$p.adj)))
custom_colors4 <- c("Norma" = "#929cef", "CA125 padidėjimas" = "#cf5784")
ca_expr_plot <- ggplot(ca_table, aes(x=CA125_f, y=value)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = CA125_f)) +
  geom_jitter(aes(color = CA125_f), size=2, alpha=0.5) +
  ylab(label = expression("Santykinė raiška, normalizuota pagal " * italic("GAPDH"))) + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  add_pvalue(each.vs.ref_sig_ca, label = "p.adj_custom") + #pvalue
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), #turn
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  scale_colour_manual(values = custom_colors4)+
  scale_fill_manual(values = custom_colors4)+
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x))) #add long "-" signs

ca_expr_plot

#SAVE PNG 
png("met_exprs_BOXPLOT__ca125_output0530.png", width = 3000, height = 2300, res = 300) # width and height in pixels, resolution in dpi
ca_expr_plot #
dev.off() # Close the PNG device

#fc
raiska_df5 <- KN_data[, c("CA125_f", raiska, "patient_id_aud")]
rownames(raiska_df5) <- raiska_df5$patient_id_aud
raiska_df5 <- raiska_df5[, -12]
#exp df
exp_df5 <- raiska_df5 %>%
  melt(id.vars="CA125_f",  measure.vars=raiska) %>%
  rename(gene = variable) %>%
  rename(expression = value)%>%
  rename(condition = CA125_f) %>%
  group_by(condition, gene)%>%
  summarise(mean_expression = mean(expression, drop_na= T)) %>%
  spread(condition, mean_expression)

mean_expression_ca125<- exp_df5 %>%
  mutate(fold_change_ca125 = log2(2^(`CA125 increase` - `Norm`) ))
mean_expression_ca125

#fc 64 cases for DLL1
raiska_df5 <- KN_data[, c("CA125_f", raiska, "patient_id_aud")]
rownames(raiska_df5) <- raiska_df5$patient_id_aud
raiska_df5 <- raiska_df5[, -12]
#remove one dll1 value
raiska_df5 <-  raiska_df5 %>% filter(!is.na(DLL1))
dim(raiska_df5) #64 cases left
#exp df
exp_df5 <- raiska_df5 %>%
  melt(id.vars="CA125_f",  measure.vars=raiska) %>%
  rename(gene = variable) %>%
  rename(expression = value)%>%
  rename(condition = CA125_f) %>%
  group_by(condition, gene)%>%
  summarise(mean_expression = mean(expression, drop_na= T)) %>%
  spread(condition, mean_expression)

mean_expression_ca125<- exp_df5 %>%
  mutate(fold_change_ca125 = `CA125 increase` / `Norm`) 
mean_expression_ca125

#EXPR age test########################
shapiro_results <- lapply(KN_data[,colnames(KN_data)%in% c(raiska, "Amžius")], shapiro.test)
shapiro_results #NOTCH4 NOTCH1

age_test <- KN_data[, colnames(KN_data) %in% c(raiska, "Amžius")]
#not normal TCEAL4, GRB7, PKP3, LUC7L2
nn <- c("NOTCH4", "NOTCH1")
on <- raiska[!raiska %in% nn]
#SPEARMAN CORR
results_nn <- lapply(KN_data[, colnames(KN_data) %in% nn], 
                     function(x) cor.test(x, KN_data$Amžius, method = "spearman"))
results_nn
#PREARSON CORR
#get a df because more genes are tested
results_n <- sapply(KN_data[, colnames(KN_data) %in% on], 
                  function(x) {
                    test <- cor.test(x, KN_data$Amžius, method = "pearson")
                    c(cor = unname(test$estimate), p = test$p.value)
                  })

results_df <- as.data.frame(t(results_n))
results_df$variable <- rownames(results_df)
rownames(results_df) <- NULL
results_df

#plot all correlations
# reshape to long format
plot_df <- KN_data[, c("Amžius", on)]
plot_df <- melt(plot_df, id.vars = "Amžius")

# compute correlations + label text
stats <- plot_df %>%
  group_by(variable) %>%
  summarise(
    cor = cor(Amžius, value, use = "complete.obs"),
    p   = cor.test(Amžius, value)$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = paste0("r = ", round(cor, 2),
                        ", p = ", signif(p, 3)))

# plot
corr_plot <- ggplot(plot_df, aes(x = Amžius, y = value)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  facet_wrap(~ variable, scales = "free_y") +
  # use annotate per facet: match variable names
  geom_text(
    data = stats,
    mapping = aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1, vjust = 1.2,
    inherit.aes = FALSE
  ) +
  theme_minimal() +
  labs(y = "Expression", x = "Amžius")

#SAVE PNG 
png("met_exprs_CORELATION_age_output0530.png", width = 3000, height = 2300, res = 300) # width and height in pixels, resolution in dpi
corr_plot #
dev.off() # Close the PNG device

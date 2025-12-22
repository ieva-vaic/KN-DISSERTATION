##KN-DISSERTATION project. Mutation data (ARID1A and CTNNB1 genes)
#Plots with mutations
#libraries
library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(ggpubr)
library(scales)  
library(reshape2)
#set wd for plots###########################################
setwd("C:/Users/Ieva/rprojects/outputs_all/DISS/")
#upload data##############################################################
KN_data <- read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/audiniu_mut_exprs_met20250709.xlsx")
#mutation genes
expression3 <- c("ARID1A","CTNNB1", "FBXW7")
#expression genes
raiska <- c("ARID1A", "CTNNB1")
#get arid1a
colnames(KN_data)
Arid1a_columns <- c("Histology", "Grupė_Ieva", "ARID1A","CTNNB1", "FBXW7",
                    "ARID1A_met", "ARID1A_tumor_mut","ARID1A_tumor_VUS",  "ARID1A_tumor_type",
                    "CTNNB1_tumor_mut",  "TP53_tumor_mut", "KN" )
#make arid1a df
ARID1A_df <- KN_data[, colnames(KN_data) %in% Arid1a_columns]
#fix arid1a df
str(ARID1A_df)
ARID1A_df <- data.frame(ARID1A_df)
ARID1A_df$Grupė_Ieva <- factor(ARID1A_df$Grupė_Ieva)
levels(ARID1A_df$Grupė_Ieva) <- c("Gerybinis", "HGSOC", "Kiti KV")
ARID1A_df$ARID1A_met <- factor(ARID1A_df$ARID1A_met)
ARID1A_df$ARID1A_met <- factor(ARID1A_df$ARID1A_met, levels = c(0, 1), 
                               labels = c("Nemetilintas", "Metilintas"))


#ARID1A_df$ARID1A_tumor_VUS <- factor(ARID1A_df$ARID1A_tumor_VUS)
#ARID1A_df$ARID1A_tumor_type <- factor(ARID1A_df$ARID1A_tumor_type)
#fix mutations
cols <- c("ARID1A_tumor_mut", "CTNNB1_tumor_mut",  "TP53_tumor_mut")
ARID1A_df[cols] <- lapply(ARID1A_df[cols], function(x) ifelse(is.na(x), "Be mutacijų", "Mutacija"))
##only up to kn-95 the samples had mutation data, others should be NA
rownames(ARID1A_df) <- ARID1A_df$KN
not_mut <- c("KN-96" , "KN-97",  "KN-99",  "KN-100" ,"KN-101", "KN-103", "KN-104" ,
             "KN-105", "KN-106", "KN-107", "KN-108", "KN-109", "KN-110", "KN-111", "KN-112")
ARID1A_df[row.names(ARID1A_df) %in% not_mut, cols] <- NA
#histology to lt version
ARID1A_df$Histology <- ifelse(ARID1A_df$Histology %in% c("Endometrial", "Endometriod"), "Endometrioidinis",
                              ARID1A_df$Histology )
ARID1A_df$Histology <- ifelse(ARID1A_df$Histology == "Endometriois", "Endometriozė",
                              ARID1A_df$Histology )
ARID1A_df$Histology <- ifelse(ARID1A_df$Histology == "Cystis", "Cista",
                              ARID1A_df$Histology )
ARID1A_df$Histology <- ifelse(ARID1A_df$Histology == "Mucinous", "Mucininis",
                              ARID1A_df$Histology )
ARID1A_df$Histology <- ifelse(ARID1A_df$Histology == "Serous", "Serozinis",
                              ARID1A_df$Histology )
ARID1A_df$Histology <- ifelse(ARID1A_df$Histology == "Clear cell", "Šviesių lastelių",
                              ARID1A_df$Histology )
ARID1A_df$Histology <- ifelse(ARID1A_df$Histology == "RSS", "Riziką mažinanti operacija",
                              ARID1A_df$Histology )
ARID1A_df$Histology <- ifelse(ARID1A_df$Histology == "Granuloza", "Granulosa",
                              ARID1A_df$Histology )
ARID1A_df$Histology <- factor(ARID1A_df$Histology)
#rename arid1a VUS
ARID1A_df$ARID1A_tumor_VUS[is.na(ARID1A_df$ARID1A_tumor_VUS)] <- "Be mutaciju"
ARID1A_df$ARID1A_tumor_VUS <- ifelse(ARID1A_df$ARID1A_tumor_VUS == "Pathogenic", "Patogeninė, VUS",
                                     ARID1A_df$ARID1A_tumor_VUS )
ARID1A_df$ARID1A_tumor_VUS <- ifelse(ARID1A_df$ARID1A_tumor_VUS == "VUS, benign", "Gerybinė, VUS",
                                         ARID1A_df$ARID1A_tumor_VUS )
ARID1A_df$ARID1A_tumor_VUS

#Methylation vs expression ARID1A ##########################
shapiro.test(ARID1A_df$ARID1A)
var.test(ARID1A_df$ARID1A ~ ARID1A_df$ARID1A_met)
arid1a_p <- t.test(ARID1A_df$ARID1A ~ ARID1A_df$ARID1A_met, alternative = "two.sided",  var.equal = T ) # 0.1649
arid1a_p
p_text <- paste0("p = ", signif(arid1a_p$p.value, 3))
formatted_p <- sprintf("p = %.3f", arid1a_p$p.value)
p_df <- data.frame(
  group1 = "Nemetilintas",
  group2 = "Metilintas",
  y.position = max(ARID1A_df$ARID1A, na.rm = TRUE) + 0.5,
  p.value = formatted_p,
  inherit.aes = FALSE
)
custom_gray_colors = c("Metilintas" = "lightpink", "Nemetilintas" = "lightblue4")
ARID1A_met <- ggplot(ARID1A_df, aes(x = ARID1A_met, y = ARID1A, fill = ARID1A_met))+
geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = ARID1A_met )) +
  geom_jitter(aes(color = ARID1A_met ), size=1, alpha=0.5) +
  stat_pvalue_manual(p_df, label = "p.value",tip.length = 0.01,
                     inherit.aes = FALSE) +
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1.5)+
  scale_fill_manual(values = custom_gray_colors) +
  scale_color_manual(values = custom_gray_colors)+
  ggtitle(expression(italic("ARID1A") *" promotoriaus metilinimas " 
                     * italic("vs") * " raiška"))+
  ylab(label =
         expression("Santykinė "* italic("ARID1A") *" raiška, normalizuota pagal  " * italic("GAPDH")))+
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x))) #add long "-" signs
ARID1A_met
#save plot
png("ARID1A_MET_RAISKA_BOXPLOT_20251009.png", width = 400, height = 400,
    res =  100, units = "px", pointsize = 14) # width and height in pixels, resolution in dpi
ARID1A_met
dev.off()

##EN Methylation vs expression ARID1A ##########################
p_df <- data.frame(
  group1 = "Not methylated",
  group2 = "Methylated",
  y.position = max(ARID1A_df$ARID1A, na.rm = TRUE) + 0.5,
  p.value = formatted_p,
  inherit.aes = FALSE
)
ARID1A_dfEN <- ARID1A_df %>%
  mutate(ARID1A_met = recode(ARID1A_met,
                             "Nemetilintas" = "Not methylated",
                             "Metilintas"   = "Methylated"))

custom_gray_colorsEN = c("Methylated" = "lightpink", "Not methylated" = "lightblue4")
ARID1A_metEN <- ggplot(ARID1A_dfEN, aes(x = ARID1A_met, y = ARID1A, fill = ARID1A_met))+
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = ARID1A_met )) +
  geom_jitter(aes(color = ARID1A_met ), size=1, alpha=0.5) +
  stat_pvalue_manual(p_df, label = "p.value",tip.length = 0.01,
                     inherit.aes = FALSE) +
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  annotate("text", x = -Inf, y = Inf, label = "A", hjust = -0.5, vjust = 1.5)+
  scale_fill_manual(values = custom_gray_colorsEN) +
  scale_color_manual(values = custom_gray_colorsEN)+
  ggtitle(expression(italic("ARID1A") *" promoter methylation " 
                     * italic("vs") * " expression"))+
  ylab(label =
         expression("Relative "* italic("ARID1A") *" expression, normalized to  " * italic("GAPDH")))+
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x))) #add long "-" signs
ARID1A_metEN

#save plot
png("ARID1A_MET_RAISKA_BOXPLOT_20251212.png", width = 400, height = 400,
    res =  100, units = "px", pointsize = 14) # width and height in pixels, resolution in dpi
ARID1A_metEN
dev.off()

#FC
raiska_df <- ARID1A_df[, c("ARID1A_met", raiska, "KN")]
rownames(raiska_df) <- raiska_df$KN
raiska_df <- raiska_df[, -5]
#exp df
exp_df <- raiska_df %>%
  melt(id.vars="ARID1A_met",  measure.vars=raiska) %>%
  rename(gene = variable) %>%
  rename(expression = value)%>%
  rename(condition = ARID1A_met) %>%
  group_by(condition, gene)%>%
  summarise(mean_expression = mean(expression, drop_na= T)) %>%
  spread(condition, mean_expression)
mean_expression_tumor2 <- exp_df %>%
  mutate(fold_change_HB = log2(2^`Metilintas` / 2^`Nemetilintas`))
mean_expression_tumor2

#t test for mutations################################################
##mut  vs raiska arid1a##################################
ARID1A_df_mut <- ARID1A_df %>% filter(!is.na(ARID1A_tumor_mut))
var.test(ARID1A_df_mut$ARID1A ~ ARID1A_df_mut$ARID1A_tumor_mut)
arid1a_p2 <- t.test(ARID1A_df_mut$ARID1A ~ ARID1A_df_mut$ARID1A_tumor_mut, alternative = "two.sided",  var.equal = T ) # 0.1649
arid1a_p2
p_text2 <- paste0("T-test p = ", signif(arid1a_p2$p.value, 3))

formatted_p2 <- sprintf("p = %.3f", arid1a_p2$p.value)
p_df2 <- data.frame(
  group1 = "Be mutacijų",
  group2 = "Mutacija",
  y.position = max(ARID1A_df$ARID1A, na.rm = TRUE) + 0.5,
  p.value = formatted_p2
)
custom_gray_colors2 = c("Mutacija" = "lightpink", "Be mutacijų" = "lightblue4")
ARID1A_mut <- ggplot(ARID1A_df_mut, aes(x = ARID1A_tumor_mut, y = ARID1A, fill = ARID1A_tumor_mut))+
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = ARID1A_tumor_mut )) +
  geom_jitter(aes(color = ARID1A_tumor_mut ), size=1, alpha=0.5) +
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1.5)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values = custom_gray_colors2) +
  scale_color_manual(values = custom_gray_colors2)+
  stat_pvalue_manual(p_df2, label = "p.value",tip.length = 0.01,
                     inherit.aes = FALSE) +
  ggtitle(expression(italic("ARID1A") *" sekos pokyčiai " 
                     * italic("vs") * " raiška"))+
  ylab(label =
         expression("Santykinė "* italic("ARID1A") *" raiška, normalizuota pagal  " * italic("GAPDH")))+
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x))) #add long "-" signs
ARID1A_mut
#save png
png("ARID1A_MUT_RAISKA_BOXPLOT_20251009.png", width = 400, height = 400,
    res =  100, units = "px", pointsize = 14) # width and height in pixels, resolution in dpi
ARID1A_mut# Render the heatmap
dev.off() # Close the PNG device

### EN mut  vs raiska arid1a##################################
ARID1A_df_mutEN <- ARID1A_df_mut %>%
  mutate(ARID1A_tumor_mut = recode(ARID1A_tumor_mut,
                                   "Be mutacijų" = "No mutations",
                                   "Mutacija"   = "Mutation"))

p_df2EN <- data.frame(
    group1 = "No mutations",
    group2 = "Mutation",
    y.position = max(ARID1A_df_mutEN$ARID1A, na.rm = TRUE) + 0.5,
    p.value = formatted_p2
  )
custom_gray_colors2EN = c("Mutation" = "lightpink", "No mutations" = "lightblue4")
ARID1A_mutEN <- ggplot(ARID1A_df_mutEN, aes(x = ARID1A_tumor_mut, y = ARID1A, fill = ARID1A_tumor_mut))+
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = ARID1A_tumor_mut )) +
  geom_jitter(aes(color = ARID1A_tumor_mut ), size=1, alpha=0.5) +
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  annotate("text", x = -Inf, y = Inf, label = "B", hjust = -0.5, vjust = 1.5)+
  stat_boxplot(geom ='errorbar')+
  scale_fill_manual(values = custom_gray_colors2EN) +
  scale_color_manual(values = custom_gray_colors2EN)+
  stat_pvalue_manual(p_df2EN, label = "p.value",tip.length = 0.01,
                     inherit.aes = FALSE) +
  ggtitle(expression(italic("ARID1A") *" mutations " 
                     * italic("vs") * " expression"))+
  ylab(label =
         expression("Relative "* italic("ARID1A") *" expression, normalized to  " * italic("GAPDH")))+
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x))) #add long "-" signs
ARID1A_mutEN

#save png
png("ARID1A_MUT__RAISKA_BOXPLOT_EN20251212.png", width = 400, height = 400,
    res =  100, units = "px", pointsize = 14) 
ARID1A_mutEN# 
dev.off() 


#FC
raiska_df <- ARID1A_df_mut[, c("ARID1A_tumor_mut", raiska, "KN")]
rownames(raiska_df) <- raiska_df$KN
raiska_df <- raiska_df[, -5]
#exp df
exp_df <- raiska_df %>%
  melt(id.vars="ARID1A_tumor_mut",  measure.vars=raiska) %>%
  rename(gene = variable) %>%
  rename(expression = value)%>%
  rename(condition = ARID1A_tumor_mut) %>%
  group_by(condition, gene)%>%
  summarise(mean_expression = mean(expression, drop_na= T)) %>%
  spread(condition, mean_expression)
mean_expression_tumor2 <- exp_df %>%
  mutate(fold_change_HB = log2(2^`Mutacija` / 2^`Be mutacijų`))
mean_expression_tumor2#biger values


#simplify mut + mutacija #########################################
table(ARID1A_df_mut$met_mut_arid1a)
ARID1A_df_mut <- ARID1A_df_mut %>%
  mutate(arid1a_met_mut2 = case_when(
    met_mut_arid1a %in% c("Metilintas_Mutacija", "Metilintas_Be mutacijų", "Nemetilintas_Mutacija") ~ "Metilinimas arba mutacija",
    TRUE ~ "Nemetilintas, Be mutacijų"  # keep other values unchanged
  )) 
table(ARID1A_df_mut$arid1a_met_mut2)
#test distribution and variance
by(ARID1A_df_mut$ARID1A, ARID1A_df_mut$arid1a_met_mut2, shapiro.test)
var.test(ARID1A_df_mut$ARID1A ~ ARID1A_df_mut$arid1a_met_mut2)
#t test
arid1a_p3 <- t.test(ARID1A_df_mut$ARID1A ~ ARID1A_df_mut$arid1a_met_mut2, alternative = "two.sided",  var.equal = T ) # 0.1649
arid1a_p3
#format p value
formatted_p3 <- sprintf("p = %.3f", arid1a_p3$p.value)
p_df3 <- data.frame(
  group1 = "Metilinimas arba mutacija",
  group2 = "Nemetilintas, Be mutacijų",
  y.position = max(ARID1A_df_mut$ARID1A, na.rm = TRUE) + 0.5,
  p.value = formatted_p3
)
#plot
custom_gray_colors3 = c("Nemetilintas, Be mutacijų" = "lightblue4", "Metilinimas arba mutacija" = "lightpink")
ARID1A_met4 <- ggplot(ARID1A_df_mut, aes(x = arid1a_met_mut2, y = ARID1A, fill = arid1a_met_mut2))+
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = arid1a_met_mut2 )) +
  geom_jitter(aes(color = arid1a_met_mut2 ), size=1, alpha=0.5) +
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.5, vjust = 1.5)+
  scale_fill_manual(values = custom_gray_colors3) +
  scale_color_manual(values = custom_gray_colors3)+
  ggtitle(expression(italic("ARID1A") *" raiškos sąsaja su mutacijų ir promotoriaus metilinimo kombinacija" ))+
  stat_pvalue_manual(p_df3, label = "p.value",tip.length = 0.01,
                     inherit.aes = FALSE) +
  ylab(label =
         expression("Santykinė "* italic("ARID1A") *" raiška, normalizuota pagal  " * italic("GAPDH")))+
  scale_y_continuous(labels = function(x) 
    gsub("-", "\u2212", as.character(x))) #add long "-" signs
ARID1A_met4
#save png
png("ARID1A_MUT_MET_RAISKA_BOXPLOT_20251009.png", width = 650, height = 400,
    res =  100, units = "px", pointsize = 14) 
ARID1A_met4# 
dev.off() 

#FC
raiska_df <- ARID1A_df_mut[, c("arid1a_met_mut2", raiska, "KN")]
rownames(raiska_df) <- raiska_df$KN
raiska_df <- raiska_df[, -5]
#exp df
exp_df <- raiska_df %>%
  melt(id.vars="arid1a_met_mut2",  measure.vars=raiska) %>%
  rename(gene = variable) %>%
  rename(expression = value)%>%
  rename(condition = arid1a_met_mut2) %>%
  group_by(condition, gene)%>%
  summarise(mean_expression = mean(expression, drop_na= T)) %>%
  spread(condition, mean_expression)
mean_expression_tumor2 <- exp_df %>%
  mutate(fold_change_HB = log2(2^`Metilinimas arba mutacija` / 2^`Nemetilintas, Be mutacijų`))
mean_expression_tumor2

#BOXPLOT mutation vs methyaltion################
#median arid1a expression
ARID1A_df_mut <- ARID1A_df_mut %>%
  mutate(ARID1A_expression_group = ifelse(ARID1A > median(ARID1A, na.rm = TRUE), "big", "small"))
#fisher test expression group with mutation status
fisher.test(table(ARID1A_df_mut$ARID1A_expression_group, ARID1A_df_mut$ARID1A_tumor_mut))
fisher.test(table(ARID1A_df_mut$ARID1A_expression_group, ARID1A_df_mut$arid1a_met_mut2))
#PLOT barplots with methylation status and mutation status
# Prepare data
tbl <- table(ARID1A_df_mut$ARID1A_met, ARID1A_df_mut$ARID1A_tumor_mut)
fisher_test <- fisher.test(tbl)
df <- as.data.frame(tbl)
# Rename columns for clarity
colnames(df) <- c("Methylation", "Mutation", "Count")
# Calculate percentages by methylation status 
df <- df %>%
  group_by(Methylation) %>%
  mutate(Percent = Count / sum(Count))
#plot
barplot_arid <- ggplot(df, aes(x = Methylation, y = Percent, fill = Mutation)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = scales::percent(Percent, accuracy = 1)),
            position = position_stack(vjust = 0.5), size = 5,  color = "white") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("pink", "hotpink")) +
  labs(y = expression(italic("ARID1A") * " sekos pokyčiai, %"),
       x = NULL,
       fill = expression(italic("ARID1A") * " mutacija"),
       title = expression(italic("ARID1A") * " sekos pokyčiai " *
                            italic("vs") * " promotoriaus metilinimas"),
       subtitle = paste("Fisherio testo p =", signif(fisher_test$p.value, 3))) +
  theme_minimal(base_size = 14)
#show
barplot_arid
#save
png("ARID1A_MUT_met_fishers_20250930.png", width = 700, height = 500,
    res =  100, units = "px", pointsize = 12) # width and height in pixels, resolution in dpi
barplot_arid# Render the heatmap
dev.off() # Close the PNG device


#EN BOXPLOT MUTATIONS VS METHTLATION#################
ARID1A_df_mutEN <- ARID1A_df_mutEN %>%
  mutate(ARID1A_met = recode(ARID1A_met,
                             "Nemetilintas" = "Not methylated",
                             "Metilintas"   = "Methylated"))
#median arid1a expression
ARID1A_df_mutEN <- ARID1A_df_mutEN %>%
  mutate(ARID1A_expression_group = ifelse(ARID1A > median(ARID1A, na.rm = TRUE), "big", "small"))
# Prepare data
tblEN <- table(ARID1A_df_mutEN$ARID1A_met, ARID1A_df_mutEN$ARID1A_tumor_mut)
fisher_test <- fisher.test(tblEN)
dfEN <- as.data.frame(tblEN)
# Rename columns for clarity
colnames(dfEN) <- c("Methylation", "Mutation", "Count")
# Calculate percentages by methylation status 
dfEN <- dfEN %>%
  group_by(Methylation) %>%
  mutate(Percent = Count / sum(Count))
#plot
barplot_aridEN <- ggplot(dfEN, aes(x = Methylation, y = Percent, fill = Mutation)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = scales::percent(Percent, accuracy = 1)),
            position = position_stack(vjust = 0.5), size = 5,  color = "white") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("pink", "hotpink")) +
  labs(y = expression(italic("ARID1A") * " mutations, %"),
       x = NULL,
       fill = expression(italic("ARID1A") * " mutation"),
       title = expression(italic("ARID1A") * " mutation " *
                            italic("vs") * " promoter methylation"),
       subtitle = paste("Fisher test p =", signif(fisher_test$p.value, 3))) +
  theme_minimal(base_size = 14)
#show
barplot_aridEN

#save
png("ARID1A_MUT_met_fishers_EN20251212.png", width = 700, height = 500,
    res =  100, units = "px", pointsize = 12) # width and height in pixels, resolution in dpi
barplot_aridEN# Render the heatmap
dev.off() # Close the PNG device
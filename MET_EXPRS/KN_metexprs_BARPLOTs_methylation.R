#KN-DISSERTATION project. MET_EXPRS - methylation, expression, mutation data
#METHYLATIO / BARPLOTS
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
#methylation data must be 0 1
KN_data$HOPX <- as.numeric(KN_data$HOPX) -1 
KN_data$ALX4 <- as.numeric(KN_data$ALX4) -1 
KN_data$CDX2 <- as.numeric(KN_data$CDX2) -1 
KN_data$ARID1A_met <- as.numeric(KN_data$ARID1A_met) -1
#biomarker groups
#expression
raiska <- colnames(KN_data[18:27])
#methylation
metilinimas <- colnames(KN_data[28:31])
biomarkers <- c(raiska, metilinimas)
#change tumor
KN_data$tumor <- recode(KN_data$tumor, OvCa = "KV", Benign ="Gerybinis")
#separate groups of samples
# HGSOC vs others = OC only
KN_OTHER_HGSOC <- KN_data[KN_data$Grupė_Ieva != "Benign", ] 
KN_OTHER_HGSOC$Grupė_Ieva <- droplevels(KN_OTHER_HGSOC$Grupė_Ieva)
table(KN_OTHER_HGSOC$Grupė_Ieva) #57 left 42 vs 15
# HGSOC vs benign
KN_BENIGN_HGSOC <- KN_data[KN_data$Grupė_Ieva != "Other", ] 
KN_BENIGN_HGSOC$Grupė_Ieva <- droplevels(KN_BENIGN_HGSOC$Grupė_Ieva)
table(KN_BENIGN_HGSOC$Grupė_Ieva) #51 left 42 vs 9
# Other vs benign
KN_BENIGN_OTHER <- KN_data[KN_data$Grupė_Ieva != "HGSOC", ] 
KN_BENIGN_OTHER$Grupė_Ieva <- droplevels(KN_BENIGN_OTHER$Grupė_Ieva)
table(KN_BENIGN_OTHER$Grupė_Ieva) #24 left 9 vs 15
# HGSOC only
KN_HGSOC <- KN_data[KN_data$Grupė_Ieva == "HGSOC", ] #42 people
KN_HGSOC$Grupė_Ieva <- droplevels(KN_HGSOC$Grupė_Ieva)
table(KN_HGSOC$Grupė_Ieva) #42 left

#Methylation barplots OC vs benign###################################
#number of samples with methylation in each group
table(KN_data$ARID1A_met)
table(KN_data$CDX2)
table(KN_data$ALX4)
table(KN_data$HOPX)
#fisher test
fisher.test(table(KN_data$tumor, KN_data$ARID1A_met)) #1
fisher.test(table(KN_data$tumor, KN_data$HOPX)) #0.0229
fisher.test(table(KN_data$tumor, KN_data$ALX4)) #0.03554
fisher.test(table(KN_data$tumor, KN_data$CDX2)) #0.02684
#barplot
metilinimas_tumor_table <- KN_data[, colnames(KN_data) %in% c(metilinimas, "tumor")]

met_long_tumor <- pivot_longer(metilinimas_tumor_table, 
                               cols = -tumor, 
                               names_to = "biomarker", 
                               values_to = "value")

# Rename biomarker column
met_long_tumor <- met_long_tumor %>%
  # Rename
  mutate(biomarker = recode(biomarker, 
                            "HOPX_met" = "HOPX",      
                            "ALX4_met" = "ALX4",   
                            "CDX2_met" = "CDX2",
                            "ARID1A_met" = "ARID1A")) %>%
  mutate(
    value = case_when(
      value == "1" ~ "metilintas",
      value == "0" ~ "nemetilintas",
      TRUE ~ as.character(value)  # Keep other values unchanged
    ),
    value = factor(value, levels = c( "nemetilintas", "metilintas"))
  )

# Calculate counts
df_perc_tumor <- met_long_tumor %>%
  group_by(tumor, biomarker, value) %>%
  summarise(count = n()) %>%
  group_by(tumor, biomarker) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup() %>%
  mutate(value = factor(value, levels = c( "nemetilintas", "metilintas")))  
# Ensuring factor levels for consistent ordering

# Plot stacked barplots
p1 <- ggplot(df_perc_tumor, aes(x = tumor, y = percentage, fill = value)) +
  geom_bar(stat = "identity")  +
  facet_wrap(~ biomarker, scales = "free_y", nrow = 1, ncol = 4) +
  labs(x = "Tumor type", y = "Promotorių metilinimas, %", fill = "Value") +
  theme_minimal()+
  labs(x=NULL,  
       #title = "Promoter methylation in gynecologic tumors"
  )+
  theme(legend.position = "bottom",
    strip.text.x = element_text(
      size = 15, face = "bold.italic"
    ),
    legend.text = element_text(face = "italic"),
    plot.title = element_text(hjust = 0.5))+
  labs(fill = "Promotorių metilinimo statusas") 
p1

#add p values
anno <- data.frame(x1 = c(1, 1, 1, 1), x2 = c(2, 2, 2, 2), 
                   y1 = c(101, 101, 101, 101), y2 = c(102, 102, 102, 102), 
                   xstar = c(1.5, 1.5, 1.5, 1.5), ystar = c(104, 104, 104, 104),
                   lab = c("p = 0,023", "p = 0,027" , "p = 0,036", "p = 1"),
                   biomarker = c("HOPX", "CDX2", "ALX4", "ARID1A"), 
                   value = c("nemetilintas", "nemetilintas", "nemetilintas", "nemetilintas"))
anno

new_value_colors <- c("metilintas" = "#cf5784", "nemetilintas" = "#dcbeff") #spalvos BUVO #CBC3E3 ; #AFE1AF

pp1 <- p1+
  geom_text(data = anno, aes(x = xstar,  y = ystar, label = lab))+ #pvalues
  geom_text(aes(label = paste0(round(percentage), "%")), #balti skaiciai
            position = position_stack(vjust = 0.5), 
            size = 5, 
            color = "white", 
            fontface = "bold") +
  geom_segment(data = anno, aes(x = x1, xend = x1,  #visu 3 reikia nubraizyti brakets
                                y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno, aes(x = x2, xend = x2, 
                                y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno, aes(x = x1, xend = x2, 
                                y = y2, yend = y2),
               colour = "black") +
  scale_fill_manual(values = new_value_colors, labels = c( 'nemetilintas', 'metilintas'))
pp1


#SAVE PNG 
#set directory for saving plots
setwd("C:/Users/Ieva/rprojects/outputs_all/DISS/")
#save
png("met_BARPLOT_OVca_output20251020.png",
    width = 3000, height = 2200, res = 360) # width and height in pixels, resolution in dpi
pp1# Render the heatmap
dev.off() # Close the PNG device


#Methylation barplots 3 groups#################
#fisher tests
#HGSOC vs benign
fisher.test(table(KN_BENIGN_HGSOC$Grupė_Ieva, KN_BENIGN_HGSOC$ARID1A_met)) #1
fisher.test(table(KN_BENIGN_HGSOC$Grupė_Ieva, KN_BENIGN_HGSOC$HOPX)) #0.043
fisher.test(table(KN_BENIGN_HGSOC$Grupė_Ieva, KN_BENIGN_HGSOC$ALX4)) #0.06
fisher.test(table(KN_BENIGN_HGSOC$Grupė_Ieva, KN_BENIGN_HGSOC$CDX2)) #0.06
#HGSOC vs other
fisher.test(table(KN_OTHER_HGSOC$Grupė_Ieva, KN_OTHER_HGSOC$ARID1A_met)) #0.55
fisher.test(table(KN_OTHER_HGSOC$Grupė_Ieva, KN_OTHER_HGSOC$HOPX)) #0.76
fisher.test(table(KN_OTHER_HGSOC$Grupė_Ieva, KN_OTHER_HGSOC$ALX4)) #1
fisher.test(table(KN_OTHER_HGSOC$Grupė_Ieva, KN_OTHER_HGSOC$CDX2)) #0.22
#benign vs other
fisher.test(table(KN_BENIGN_OTHER$Grupė_Ieva, KN_BENIGN_OTHER$ARID1A_met)) #1
fisher.test(table(KN_BENIGN_OTHER$Grupė_Ieva, KN_BENIGN_OTHER$HOPX)) #0.04
fisher.test(table(KN_BENIGN_OTHER$Grupė_Ieva, KN_BENIGN_OTHER$ALX4)) #0.09
fisher.test(table(KN_BENIGN_OTHER$Grupė_Ieva, KN_BENIGN_OTHER$CDX2)) #0.009

#barplot 3  groups
metilinimas_group_table <- KN_data[, colnames(KN_data) %in% c(metilinimas, "Grupė_Ieva")]

met_long_group <- pivot_longer(metilinimas_group_table, 
                               cols = -Grupė_Ieva, 
                               names_to = "biomarker", 
                               values_to = "value")

# Rename biomarker column
met_long_group <- met_long_group %>%
  # Rename
  mutate(biomarker = recode(biomarker, 
                            "HOPX_met" = "HOPX",      
                            "ALX4_met" = "ALX4",   
                            "CDX2_met" = "CDX2",
                            "ARID1A_met" = "ARID1A")) %>%
  mutate(
    value = case_when(
      value == "0" ~ "nemetilintas",
      value == "1" ~ "metilintas",
      TRUE ~ as.character(value)  # Keep other values unchanged
    ),
    value = factor(value, levels = c("metilintas", "nemetilintas"))
  )%>%
  # Rename
  mutate(Grupė_Ieva = recode(Grupė_Ieva, 
                             "Benign" = "Gerybiniai",      
                             "HGSOC" = "HGSOC",   
                             "Other" = "Kiti KV"))

# Calculate counts
df_perc_group <- met_long_group %>%
  group_by(Grupė_Ieva, biomarker, value) %>%
  summarise(count = n()) %>%
  group_by(Grupė_Ieva, biomarker) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup() %>%
  mutate(value = factor(value, levels = c("nemetilintas", "metilintas")))  
# Ensuring factor levels for consistent ordering

# Plot stacked barplots
p3 <- ggplot(df_perc_group, aes(x = Grupė_Ieva, y = percentage, fill = value)) +
  geom_bar(stat = "identity")  +
  facet_wrap(~ biomarker, scales = "free_y", nrow = 1, ncol = 4) +
  labs(x = "Tumor type", y = "Promotorių metilinimas, procentais", fill = "Value") +
  theme_minimal()+
  labs(x=NULL,  
       #title = "Promoter methylation in gynecologic tumors"
  )+
  theme(legend.position = "bottom",  
        strip.text.x = element_text(
          size = 15, face = "bold.italic"
        ),
        legend.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5))+
  labs(fill = "Promotorių metilinimo statusas") 
p3

#add p values HGSOC VS BENIGN
anno <- data.frame(x1 = c(1, 1, 1, 1), x2 = c(2, 2, 2, 2), 
                   y1 = c(101, 101, 101, 101), y2 = c(102, 102, 102, 102), 
                   xstar = c(1.5, 1.5, 1.5, 1.5), ystar = c(104, 104, 104, 104),
                   lab = c("p = 0,043", "p = 0,06" , "p = 0,06", "p = 1"),
                   biomarker = c("HOPX", "CDX2", "ALX4", "ARID1A"), 
                   value = c("nemetilintas", "nemetilintas", "nemetilintas", "nemetilintas"))
anno

pp3 <- p3+
  geom_text(data = anno, aes(x = xstar,  y = ystar, label = lab))+ #pvalues
  geom_text(aes(label = paste0(round(percentage), "%")), #balti skaiciai
            position = position_stack(vjust = 0.5), 
            size = 5, 
            color = "white", 
            fontface = "bold") +
  geom_segment(data = anno, aes(x = x1, xend = x1,  #visu 3 reikia nubraizyti brakets
                                y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno, aes(x = x2, xend = x2, 
                                y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno, aes(x = x1, xend = x2, 
                                y = y2, yend = y2),
               colour = "black") +
  scale_fill_manual(values = new_value_colors, labels = c('nemetilintas', 'metilintas'))
pp3

#benign vs others
anno2 <- data.frame(x1 = c(1, 1,1,1), x2 = c(3, 3, 3, 3), 
                    y1 = c(105, 105, 105, 105), y2 = c(107, 107, 107, 107), 
                    xstar = c(2, 2, 2,2), ystar = c(109, 109, 109, 109),
                    lab = c("p = 0,048", "p = 0,01" , "p = 0,08", "p = 1"),
                    biomarker = c("HOPX", "CDX2", "ALX4", "ARID1A"), 
                    value = c("nemetilintas", "nemetilintas", "nemetilintas", "nemetilintas"))
anno2

ppp3 <- pp3 +
  geom_text(data = anno2, aes(x = xstar,  y = ystar, label = lab))+ #pvalues
  geom_segment(data = anno2, aes(x = x1, xend = x1,  #visu 3 reikia nubraizyti brakets
                                 y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno2, aes(x = x2, xend = x2, 
                                 y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno2, aes(x = x1, xend = x2, 
                                 y = y2, yend = y2),
               colour = "black") +
  scale_fill_manual(values = new_value_colors, labels = c('nemetilintas', 'metilintas'))

ppp3

# #HGSOC vs others
# anno3 <- data.frame(x1 = c(2, 2,2,2), x2 = c(3, 3, 3, 3), 
#                     y1 = c(103, 103, 103, 103), y2 = c(104, 104, 104, 104), 
#                     xstar = c(2.5,2.5,2.5, 2.5), ystar = c(106, 106, 106, 106), #skaiciai neveikia 
#                     lab = c("p = 0,762", "p = 0,220" , "p = 1", "p = 0,551"),
#                     biomarker = c("HOPX", "CDX2", "ALX4", "ARID1A"), 
#                     value = c("nemetilintas", "nemetilintas", "nemetilintas", "nemetilintas"))
# anno3
# 
# ppp3_final <- ppp3 +
#   geom_text(data = anno3, aes(x = xstar,  y = ystar, label = lab))+ #pvalues
#   geom_segment(data = anno3, aes(x = x1, xend = x1,  #visu 3 reikia nubraizyti brakets
#                                  y = y1, yend = y2),
#                colour = "black") +
#   geom_segment(data = anno3, aes(x = x2, xend = x2, 
#                                  y = y1, yend = y2),
#                colour = "black") +
#   geom_segment(data = anno3, aes(x = x1, xend = x2, 
#                                  y = y2, yend = y2),
#                colour = "black") 
# 
# ppp3_final

#SAVE PNG 
png("met_BARPLOT_3GR_output20251020.png", 
    width = 3000, height = 2600, res = 370) # width and height in pixels, resolution in dpi
ppp3# Render the heatmap
dev.off() # Close the PNG device

#methylation, count % for text####################################################
#arid1a
ARD_TAB <- table(KN_data$ARID1A_met)
prop.table(ARD_TAB) * 100
#cdx
cdx_TAB <- table(KN_data$CDX2)
prop.table(cdx_TAB) * 100
#alx
alx_TAB <- table(KN_data$ALX4)
prop.table(alx_TAB) * 100
#hopx
hoxp_TAB <- table(KN_data$HOPX)
prop.table(hoxp_TAB) * 100
#all or nothing, methylation
KN_data <- KN_data %>%
  mutate(METHYLATION = if_else(rowSums(across(c(ARID1A_met, CDX2, ALX4, HOPX))) > 0, 1, 0))
all_met_TAB <- table(KN_data$METHYLATION)
prop.table(all_met_TAB) * 100

#Chek clinical features against methylation###################################
#stage, all 4 stages
fisher.test(table(KN_data$Stage4, KN_data$ARID1A_met)) #0.1765
fisher.test(table(KN_data$Stage4, KN_data$HOPX)) #0.8597
fisher.test(table(KN_data$Stage4, KN_data$ALX4)) #0.7301
fisher.test(table(KN_data$Stage4, KN_data$CDX2)) #0.4025

#stage, grouped into 2 stages
fisher.test(table(KN_data$Stage2, KN_data$ARID1A_met)) #0.5254
fisher.test(table(KN_data$Stage2, KN_data$HOPX)) #1
fisher.test(table(KN_data$Stage2, KN_data$ALX4)) #0.7458
fisher.test(table(KN_data$Stage2, KN_data$CDX2)) #0.1913

#grade, grouped into 2 stages
fisher.test(table(KN_data$Grade2, KN_data$ARID1A_met)) #1
fisher.test(table(KN_data$Grade2, KN_data$HOPX)) #0.669
fisher.test(table(KN_data$Grade2, KN_data$ALX4)) #1
fisher.test(table(KN_data$Grade2, KN_data$CDX2)) #0.2115

#CA125
fisher.test(table(KN_data$CA125_f, KN_data$ARID1A_met)) #1
fisher.test(table(KN_data$CA125_f, KN_data$HOPX)) #0.16
fisher.test(table(KN_data$CA125_f, KN_data$ALX4)) #1
fisher.test(table(KN_data$CA125_f, KN_data$CDX2)) #1

#AGE
sum(!is.na(KN_data$Amžius)) #I have 65 cases
shapiro.test(KN_data$Amžius) #normal, not gonna chek further in groups

# Perform the Mann-Whitney U test for each factor variable 
results_age <- lapply(metilinimas, function(var) {
  test <- wilcox.test(KN_data$Amžius ~ KN_data[[var]], data = KN_data)
  list(variable = var, p_value = test$p.value)
})
results_age <- do.call(rbind, lapply(results_age, as.data.frame))
print(results_age) #nothing is significant

#boxplot age vs methylation
age_table <- melt(KN_data, id.vars="Amžius",  measure.vars=metilinimas)
age_table$value <- as.factor(age_table$value)
age_table$value <- recode(age_table$value, "0" = "nemetilintas", "1" = "metilintas" )
age_table$variable <- recode(age_table$variable, "ARID1A_met" = "ARID1A")

each.vs.ref_age <- tibble::tribble(
  ~group1, ~group2, ~p.adj,   ~y.position, ~variable,
  "metilintas",   "nemetilintas", 0.96, 100, "HOPX",
  "metilintas",   "nemetilintas", 0.57,  100, "ALX4",
  "metilintas",   "nemetilintas", 0.35, 100, "CDX2",
  "metilintas",   "nemetilintas", 0.54, 100, "ARID1A"
)

age_plot <- ggplot(age_table, aes(x=value, y=Amžius)) +
  geom_boxplot( outlier.shape = NA , alpha=0.3, aes(fill = value)) +
  geom_jitter(aes(color = value), size=1, alpha=0.5) +
  ylab(label = "Amžius, metais") + 
  facet_wrap(.~ variable, nrow = 2, scales = "free") +
  theme_minimal()+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold.italic"
    ),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5))+
  labs(x=NULL)+
  stat_boxplot(geom ='errorbar')+
  add_pvalue(each.vs.ref_age, label = "p.adj")  #pvalue
age_plot #

png("met_BOXPLOT_age_20250602.png", width = 2000, height = 2000, res = 300) # width and height in pixels, resolution in dpi
age_plot# Render the heatmap
dev.off() # Close the PNG device


#HGSOC only: Chek clinical features against methylation###########
KN_HGSOC #HGSOC DF

#HGSOC stage, all 4 stages
fisher.test(table(KN_HGSOC$Stage4, KN_HGSOC$ARID1A_met)) #0.04763 *
set.seed(123) 
fisher.test(table(KN_HGSOC$Stage4, KN_HGSOC$ARID1A_met), simulate.p.value = TRUE, B = 10000) #0.0433 *
fisher.test(table(KN_HGSOC$Stage4, KN_HGSOC$HOPX)) #0.8793
fisher.test(table(KN_HGSOC$Stage4, KN_HGSOC$ALX4)) #1
fisher.test(table(KN_HGSOC$Stage4, KN_HGSOC$CDX2)) #0.2885

#HGSOC stage, grouped into 2 stages
fisher.test(table(KN_HGSOC$Stage2, KN_HGSOC$ARID1A_met)) #0.2596
fisher.test(table(KN_HGSOC$Stage2, KN_HGSOC$HOPX)) #1
fisher.test(table(KN_HGSOC$Stage2, KN_HGSOC$ALX4)) #1
fisher.test(table(KN_HGSOC$Stage2, KN_HGSOC$CDX2)) #1

#grade, grouped into 2 stages
fisher.test(table(KN_OTHER_HGSOC$Grade2, KN_OTHER_HGSOC$ARID1A_met)) #1
fisher.test(table(KN_OTHER_HGSOC$Grade2, KN_OTHER_HGSOC$HOPX)) #0.669
fisher.test(table(KN_OTHER_HGSOC$Grade2, KN_OTHER_HGSOC$ALX4)) #1
fisher.test(table(KN_OTHER_HGSOC$Grade2, KN_OTHER_HGSOC$CDX2)) #0.2115

#Barplot HGSOC 4 FIGO stages################################
table(KN_HGSOC$Stage4, KN_HGSOC$ARID1A_met, useNA = "a")
#barplot 3  groups
metilinimas_stage_table <- KN_HGSOC[, colnames(KN_HGSOC) %in%
                                      c(metilinimas, "Stage4")]

met_long_stage <- pivot_longer(metilinimas_stage_table, 
                               cols = -Stage4, 
                               names_to = "biomarker", 
                               values_to = "value")

# Rename biomarker column
met_long_stage <- met_long_stage %>%
  # Rename
  mutate(biomarker = recode(biomarker, 
                            "HOPX_met" = "HOPX",      
                            "ALX4_met" = "ALX4",   
                            "CDX2_met" = "CDX2",
                            "ARID1A_met" = "ARID1A")) %>%
  mutate(
    value = case_when(
      value == "0" ~ "nemetilintas",
      value == "1" ~ "metilintas",
      TRUE ~ as.character(value)  # Keep other values unchanged
    ),
    value = factor(value, levels = c("metilintas", "nemetilintas"))
  )

# Calculate counts
df_perc_stage <- met_long_stage %>%
  group_by(Stage4, biomarker, value) %>%
  summarise(count = n()) %>%
  group_by(Stage4, biomarker) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup() %>%
  mutate(value = factor(value, levels = c("nemetilintas", "metilintas")))  
df_perc_stage
# Ensuring factor levels for consistent ordering

# Plot stacked barplots
plot3 <- ggplot(df_perc_stage, aes(x = Stage4, y = percentage, fill = value)) +
  geom_bar(stat = "identity")  +
  facet_wrap(~ biomarker, scales = "free_y", nrow = 1, ncol = 4) +
  labs(x = "Tumor FIGO stage", y = "Promotorių metilinimas, procentais", fill = "Value") +
  theme_minimal()+
  labs(x=NULL,  
       #title = "Promoter methylation in HGSOC according to FIGO stage"
  )+
  theme(legend.position = "bottom",
    strip.text.x = element_text(
      size = 15, face = "bold.italic"
    ),
    legend.text = element_text(face = "italic"),
    plot.title = element_text(hjust = 0.5))+
  labs(fill = "Promotorių metilinimo statusas") 
plot3

#add p values HGSOC VS BENIGN
anno_plot <- data.frame(x1 = c(1, 1, 1, 1), x2 = c(3, 3, 3, 3), 
                   y1 = c(101, 101, 101, 101), y2 = c(102, 102, 102, 102), 
                   xstar = c(2, 2, 2, 2), ystar = c(104, 104, 104, 104),
                   lab = c("p = 0,9", "p = 0,3" , "p = 1", "p = 0,048"),
                   biomarker = c("HOPX", "CDX2", "ALX4", "ARID1A"), 
                   value = c("nemetilintas", "nemetilintas", "nemetilintas", "nemetilintas"))
anno_plot

plot33 <- plot3+
  geom_text(data = anno_plot, aes(x = xstar,  y = ystar, label = lab))+ #pvalues
  geom_text(aes(label = paste0(round(percentage), "%")), #balti skaiciai
            position = position_stack(vjust = 0.5), 
            size = 5, 
            color = "white", 
            fontface = "bold") +
  geom_segment(data = anno_plot, aes(x = x1, xend = x1,  #visu 3 reikia nubraizyti brakets
                                y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno_plot, aes(x = x2, xend = x2, 
                                y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno_plot, aes(x = x1, xend = x2, 
                                y = y2, yend = y2),
               colour = "black") +
  scale_fill_manual(values = new_value_colors, labels = c('nemetilintas', 'metilintas'))
plot33

#SAVE PNG 

png("met_BARPLOT_HGSOC_FIGO_output20251020.png",
    width = 4000, height = 2500, res = 400) # width and height in pixels, resolution in dpi
plot33# Render the heatmap
dev.off() # Close the PNG device

#ENGLISH PLOTS ####################################
##EN OC barplot ########
metilinimas_tumor_tableEN <- metilinimas_tumor_table  %>%
  mutate(tumor = recode(tumor, "KV" = "OC", "Gerybinis" = 
                          "Benign")) 

met_long_tumorEN <- pivot_longer(metilinimas_tumor_tableEN, 
                                 cols = -tumor, 
                                 names_to = "biomarker", 
                                 values_to = "value")

# Rename biomarker column
met_long_tumorEN <- met_long_tumorEN %>%
  # Rename
  mutate(biomarker = recode(biomarker, 
                            "HOPX_met" = "HOPX",      
                            "ALX4_met" = "ALX4",   
                            "CDX2_met" = "CDX2",
                            "ARID1A_met" = "ARID1A")) %>%
  mutate(
    value = case_when(
      value == "1" ~ "methylated",
      value == "0" ~ "not methylated",
      TRUE ~ as.character(value)  # Keep other values unchanged
    ),
    value = factor(value, levels = c( "not methylated", "methylated"))
  )

# Calculate counts
df_perc_tumorEN <- met_long_tumorEN %>%
  group_by(tumor, biomarker, value) %>%
  summarise(count = n()) %>%
  group_by(tumor, biomarker) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup() %>%
  mutate(value = factor(value, levels = c( "not methylated", "methylated")))  
# Ensuring factor levels for consistent ordering

# Plot stacked barplots
p1EN <- ggplot(df_perc_tumorEN, aes(x = tumor, y = percentage, fill = value)) +
  geom_bar(stat = "identity")  +
  facet_wrap(~ biomarker, scales = "free_y", nrow = 1, ncol = 4) +
  labs(x = "Tumor type", y = "Promoter methylation, %", fill = "Value") +
  theme_minimal()+
  labs(x=NULL,  
       #title = "Promoter methylation in gynecologic tumors"
  )+
  theme(legend.position = "bottom",
        strip.text.x = element_text(
          size = 15, face = "bold.italic"
        ),
        legend.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5))+
  labs(fill = "Promoter methylation") 
p1EN

#add p values
anno <- data.frame(x1 = c(1, 1, 1, 1), x2 = c(2, 2, 2, 2), 
                   y1 = c(101, 101, 101, 101), y2 = c(102, 102, 102, 102), 
                   xstar = c(1.5, 1.5, 1.5, 1.5), ystar = c(104, 104, 104, 104),
                   lab = c("p = 0,023", "p = 0,027" , "p = 0,036", "p = 1"),
                   biomarker = c("HOPX", "CDX2", "ALX4", "ARID1A"), 
                   value = c("not methylated", "not methylated", "not methylated", "not methylated"))
anno

new_value_colors <- c("methylated" = "#cf5784", "not methylated" = "#dcbeff") #spalvos BUVO #CBC3E3 ; #AFE1AF

pp1EN <- p1EN+
  geom_text(data = anno, aes(x = xstar,  y = ystar, label = lab))+ #pvalues
  geom_text(aes(label = paste0(round(percentage), "%")), #balti skaiciai
            position = position_stack(vjust = 0.5), 
            size = 5, 
            color = "white", 
            fontface = "bold") +
  geom_segment(data = anno, aes(x = x1, xend = x1,  #visu 3 reikia nubraizyti brakets
                                y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno, aes(x = x2, xend = x2, 
                                y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno, aes(x = x1, xend = x2, 
                                y = y2, yend = y2),
               colour = "black") +
  scale_fill_manual(values = new_value_colors, labels = c( 'not methylated', 'methylated'))
pp1EN

#save
#set directory for saving
setwd("C:/Users/Ieva/rprojects/outputs_all/DISS/")
#save
png("met_exprs_boxplot_OVca_ENoutput20251212.png", width = 3000, height = 2700, res = 400) # width and height in pixels, resolution in dpi
pp1EN #
dev.off() # Close the PNG device

## 3 group EN plot ########################
metilinimas_group_tableEN <- KN_data[, colnames(KN_data) %in% c(metilinimas, "Grupė_Ieva")]

met_long_groupEN <- pivot_longer(metilinimas_group_tableEN, 
                                 cols = -Grupė_Ieva, 
                                 names_to = "biomarker", 
                                 values_to = "value")

# Rename biomarker column
met_long_groupEN <- met_long_groupEN %>%
  # Rename
  mutate(biomarker = recode(biomarker, 
                            "HOPX_met" = "HOPX",      
                            "ALX4_met" = "ALX4",   
                            "CDX2_met" = "CDX2",
                            "ARID1A_met" = "ARID1A")) %>%
  mutate(
    value = case_when(
      value == "0" ~ "not methylated",
      value == "1" ~ "methylated",
      TRUE ~ as.character(value)  # Keep other values unchanged
    ),
    value = factor(value, levels = c("methylated", "not methylated"))
  )

# Calculate counts
df_perc_groupEN <- met_long_groupEN %>%
  group_by(Grupė_Ieva, biomarker, value) %>%
  summarise(count = n()) %>%
  group_by(Grupė_Ieva, biomarker) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup() %>%
  mutate(value = factor(value, levels = c("not methylated", "methylated")))  
# Ensuring factor levels for consistent ordering

# Plot stacked barplots
p3EN <- ggplot(df_perc_groupEN, aes(x = Grupė_Ieva, y = percentage, fill = value)) +
  geom_bar(stat = "identity")  +
  facet_wrap(~ biomarker, scales = "free_y", nrow = 1, ncol = 4) +
  labs(x = "Tumor type", y = "Promoter methylation, %", fill = "Value") +
  theme_minimal()+
  labs(x=NULL,  
       #title = "Promoter methylation in gynecologic tumors"
  )+
  theme(legend.position = "bottom",  
        strip.text.x = element_text(
          size = 15, face = "bold.italic"
        ),
        legend.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5))+
  labs(fill = "Promoter methylation") 
p3EN

#add p values HGSOC VS BENIGN
annoEN <- data.frame(x1 = c(1, 1, 1, 1), x2 = c(2, 2, 2, 2), 
                     y1 = c(101, 101, 101, 101), y2 = c(102, 102, 102, 102), 
                     xstar = c(1.5, 1.5, 1.5, 1.5), ystar = c(104, 104, 104, 104),
                     lab = c("p = 0,043", "p = 0,06" , "p = 0,06", "p = 1"),
                     biomarker = c("HOPX", "CDX2", "ALX4", "ARID1A"), 
                     value = c("not methylated", "not methylated", "not methylated", "not methylated"))
annoEN

pp3EN<- p3EN+
  geom_text(data = annoEN, aes(x = xstar,  y = ystar, label = lab))+ #pvalues
  geom_text(aes(label = paste0(round(percentage), "%")), #balti skaiciai
            position = position_stack(vjust = 0.5), 
            size = 5, 
            color = "white", 
            fontface = "bold") +
  geom_segment(data = annoEN, aes(x = x1, xend = x1,  #visu 3 reikia nubraizyti brakets
                                y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = annoEN, aes(x = x2, xend = x2, 
                                y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = annoEN, aes(x = x1, xend = x2, 
                                y = y2, yend = y2),
               colour = "black") +
  scale_fill_manual(values = new_value_colors, labels = c('not methylated', 'methylated'))
pp3EN

#benign vs others
anno2EN <- data.frame(x1 = c(1, 1,1,1), x2 = c(3, 3, 3, 3), 
                      y1 = c(105, 105, 105, 105), y2 = c(107, 107, 107, 107), 
                      xstar = c(2, 2, 2,2), ystar = c(109, 109, 109, 109),
                      lab = c("p = 0,048", "p = 0,01" , "p = 0,08", "p = 1"),
                      biomarker = c("HOPX", "CDX2", "ALX4", "ARID1A"), 
                      value = c("not methylated", "not methylated", "not methylated", "not methylated"))
anno2EN

ppp3EN <- pp3EN +
  geom_text(data = anno2, aes(x = xstar,  y = ystar, label = lab))+ #pvalues
  geom_segment(data = anno2, aes(x = x1, xend = x1,  #visu 3 reikia nubraizyti brakets
                                 y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno2, aes(x = x2, xend = x2, 
                                 y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno2, aes(x = x1, xend = x2, 
                                 y = y2, yend = y2),
               colour = "black") +
  scale_fill_manual(values = new_value_colors, labels = c('not methylated', 'methylated'))

ppp3EN

#save
png("met_exprs_boxplot_HGSOC_ENoutput20251212.png", width = 4000, height = 2700, res = 400) # width and height in pixels, resolution in dpi
ppp3EN #
dev.off() # Close the PNG device

##EN Barplot HGSOC 4 FIGO stages################################
table(KN_HGSOC$Stage4, KN_HGSOC$ARID1A_met, useNA = "a")
#barplot 3  groups
metilinimas_stage_table <- KN_HGSOC[, colnames(KN_HGSOC) %in%
                                      c(metilinimas, "Stage4")]

met_long_stageEN <- pivot_longer(metilinimas_stage_table, 
                                 cols = -Stage4, 
                                 names_to = "biomarker", 
                                 values_to = "value")

# Rename biomarker column
met_long_stageEN <- met_long_stage %>%
  # Rename
  mutate(biomarker = recode(biomarker, 
                            "HOPX_met" = "HOPX",      
                            "ALX4_met" = "ALX4",   
                            "CDX2_met" = "CDX2",
                            "ARID1A_met" = "ARID1A")) %>%
  mutate(
    value = case_when(
      value == "0" ~ "not methylated",
      value == "1" ~ "methylated",
      TRUE ~ as.character(value)  # Keep other values unchanged
    ),
    value = factor(value, levels = c("methylated", "not methylated"))
  )

# Calculate counts
df_perc_stageEN <- met_long_stageEN %>%
  group_by(Stage4, biomarker, value) %>%
  summarise(count = n()) %>%
  group_by(Stage4, biomarker) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup() %>%
  mutate(value = factor(value, levels = c("not methylated", "methylated")))  
df_perc_stageEN
# Ensuring factor levels for consistent ordering

# Plot stacked barplots
plot3EN <- ggplot(df_perc_stageEN, aes(x = Stage4, y = percentage, fill = value)) +
  geom_bar(stat = "identity")  +
  facet_wrap(~ biomarker, scales = "free_y", nrow = 1, ncol = 4) +
  labs(x = "Tumor FIGO stage", y = "Promoter methylation, %", fill = "Value") +
  theme_minimal()+
  labs(x=NULL,  
       #title = "Promoter methylation in HGSOC according to FIGO stage"
  )+
  theme(legend.position = "bottom",
        strip.text.x = element_text(
          size = 15, face = "bold.italic"
        ),
        legend.text = element_text(face = "italic"),
        plot.title = element_text(hjust = 0.5))+
  labs(fill = "Promoter methylation") 
plot3EN

#add p values HGSOC VS BENIGN
anno_plotEN <- data.frame(x1 = c(1, 1, 1, 1), x2 = c(3, 3, 3, 3), 
                          y1 = c(101, 101, 101, 101), y2 = c(102, 102, 102, 102), 
                          xstar = c(2, 2, 2, 2), ystar = c(104, 104, 104, 104),
                          lab = c("p = 0,9", "p = 0,3" , "p = 1", "p = 0,048"),
                          biomarker = c("HOPX", "CDX2", "ALX4", "ARID1A"), 
                          value = c("not methylated", "not methylated", "not methylated", "not methylated"))
anno_plotEN

plot33EN <- plot3EN+
  geom_text(data = anno_plotEN, aes(x = xstar,  y = ystar, label = lab))+ #pvalues
  geom_text(aes(label = paste0(round(percentage), "%")), #balti skaiciai
            position = position_stack(vjust = 0.5), 
            size = 5, 
            color = "white", 
            fontface = "bold") +
  geom_segment(data = anno_plotEN, aes(x = x1, xend = x1,  #visu 3 reikia nubraizyti brakets
                                       y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno_plotEN, aes(x = x2, xend = x2, 
                                       y = y1, yend = y2),
               colour = "black") +
  geom_segment(data = anno_plotEN, aes(x = x1, xend = x2, 
                                       y = y2, yend = y2),
               colour = "black") +
  scale_fill_manual(values = new_value_colors, labels = c('not methylated', 'methylated'))
plot33EN

#save
png("met_exprs_boxplot_HGSOCSTAGE_ENoutput20251212.png", width = 3000, height = 2000, res = 400) # width and height in pixels, resolution in dpi
plot33EN #
dev.off() # Close the PNG device


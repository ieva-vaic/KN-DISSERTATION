#KN-DISSERTATION project. MET_EXPRS - methylation, expression, mutation data
#HEATMAPs
Sys.setenv(LANG = "en")
#libraries
library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(circlize)
#data upload##################################################################
merged_df <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/merged_all_data_diss.RDS")
expression <- c("EXO1", "RAD50","PPT2", "LUC7L2","PKP3", "CDCA5","ZFPL1","VPS33B", "GRB7","TCEAL4", 
                "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "ARID1A", "CTNNB1", "FBXW7", "JAG2", "DLL1", "HES1")

# All expression heatmap #######################################################
Heat_data <- merged_df[, c("KN",
                           expression 
)]
Heat_data <- as.data.frame(Heat_data)
rownames(Heat_data) <- Heat_data[, 1]
Heat_data <- Heat_data[, -1]
heatmap_raiska <- Heatmap(as.matrix(Heat_data))
#Clinical data###################################################################
clinical2 <- merged_df[, c("KN", "Histology", "Tumor", "CA125_f" ,"STAGE4", "Grade.x", "Age", "CA125 po gydymo", "patient_id_aud")]
clinical2 <- as.data.frame(clinical2)
rownames(clinical2) <- clinical2$KN
clinical2 <- clinical2[, -1]
head(clinical2)
#fix na so it shows on the legend
clinical2$CA125_f <- factor(clinical2$CA125_f, levels = c("Padidėjimas", "Norma", "NA"), exclude = NULL)
clinical2$STAGE4 <- replace(clinical2$STAGE4, is.na(clinical2$STAGE4), "NA" )
clinical2$Grade.x <-replace(clinical2$Grade.x, is.na(clinical2$Grade.x), "NA" )

#FIX CA125 post op
clinical2 <- clinical2 %>%
  mutate(CA125_status_post_op = case_when(
    `CA125 po gydymo` > 10 ~ "Padidėjimas",
    `CA125 po gydymo` <= 10 ~ "Norma",
    TRUE ~ NA_character_   # keep NA values
  ))
clinical2$CA125_status_post_op
clinical2$`CA125 po gydymo`
clinical2$CA125_status_post_op <-replace(clinical2$CA125_status_post_op,
                                        is.na(clinical2$CA125_status_post_op), "NA" )

#get survival data as of 2025-09-11
SURVIVAL_KN <- openxlsx::read.xlsx("../../OTHER DATA/KN-DISSERTATION FILES/KN_MIRTIES_FAILAS_20250911.xlsx")
#make only surv df
SURV <- SURVIVAL_KN[, c(2, 3,20, 21)]
head(SURV)

#join with main data
clinical <- left_join(clinical2, SURV, by = "patient_id_aud")
#fix status
clinical <- clinical %>%
  mutate(STATUS  = case_when(
    STATUS  == 0 ~ "Gyvos",
    STATUS  == 1 ~ "Mirę",
    TRUE ~ NA_character_   # keep NA values
  ))
clinical$STATUS

clinical2$STATUS <-replace(clinical$STATUS,is.na(clinical$STATUS), "NA" )

# clinical data annotation
col_age <- colorRamp2(c(40, 90), c( "#9cd4c4", "#3c402f")) #age colors
row_ha = rowAnnotation(Histologija = clinical$Histology,
                       Navikas = clinical$Tumor, Stadija = clinical$STAGE4,
                       `Diferenciacijos laipsnis` = clinical$Grade.x,
                       Amžius = clinical$Age, CA125 = clinical$CA125_f, 
                       `CA125 po gydymo` = clinical$CA125_status_post_op, 
                       `Mirties statusas` = clinical$STATUS,
                       col = list(Histologija = c("Šviesių lastelių" = "lightblue",
                                                  "Cista" = "lightgreen", 
                                                  "Endometrioidinis" = "green",
                                                  "Endometriozė" = "darkgreen", 
                                                  "Granulosa" = "turquoise",
                                                  "HGSOC" = "deeppink",
                                                  "Mioma" = "red", 
                                                  "Mucininis" = "yellow",
                                                  "Riziką mažinanti operacija" = "orange", 
                                                  "Serozinis" = "lightpink"), 
                                  Navikas = c("HGSOC" = "#a89cd4",
                                              "Gerybinis" = "#d49cac", 
                                              "Kiti KV" = "darkblue"), 
                                  Stadija = c("I" = "#9cd4c4",  
                                              "II" = "#c8d49c", 
                                              "III" = "#d49cac", 
                                              "IV" = "#a89cd4", 
                                              "NA" = "grey"), 
                                  `Diferenciacijos laipsnis` = c("G1"="#9cd4c4", 
                                                                 "G3" = "#a89cd4", 
                                                                 "NA" = "grey"),
                                  CA125 = c("Norma"="#9cd4c4", 
                                            "Padidėjimas" = "#3c402f",
                                            "NA" = "grey"),
                                  `CA125 po gydymo` = c("Norma"="#9cd4c4", 
                                            "Padidėjimas" = "#3c402f",
                                            "NA" = "grey"),
                                  `Mirties statusas` = c("Gyvos"="#9cd4c4", 
                                                        "Mirę" = "#3c402f",
                                                        "NA" = "grey"),
                                  Amžius = col_age
                       ))

#All expression heatmap with clinical data#####################################
#expression colors
col_fun = colorRamp2(c(2, -5, -10, -15), c("#8564fb",  "#64b3fb","#e088bd", "#af2745"))
#col_fun = colorRamp2(c(2, 0, -2, -4, -6, -8, -10, -12, -14, -16),
#                     c("#e7e0fe", "#cec1fd", "#8564fb", "#64b3fb","#93cafc","#325a7e", "#e088bd", "#af2745", "#9e233e", "#4f121f"))
labels <- c("\u221215", "\u221210", "\u22125", "2" ) #THIS IS NEEDED FOR MDPI AT LEAST - LONG MINUS SIGNS
#order if histology
order_hist <- c( "Riziką mažinanti operacija", "Mioma",  "Cista", "Endometriozė", "Endometrioidinė",
                 "Granuloza", "Mucininis", "Šviesių lastelių", "Serozinis", "HGSOC" )
merged_df$Histology <- factor(merged_df$Histology, levels = order_hist)
merged_df <- merged_df[order(merged_df$Histology, merged_df$Tumor), ]
order_list <- merged_df$KN
#final heatmap
heatmap_raiska <- Heatmap(as.matrix(Heat_data),cluster_columns = TRUE,
                          name = "Santykinė genų raiška",  
                          right_annotation = row_ha, col = col_fun, 
                          row_order = order_list,
                          row_split = clinical$Tumor, 
                          column_names_gp = gpar(fontface = "italic"),
                          column_title = "Santykinė genų raiška", 
                          row_names_gp = gpar(fontsize = 8), 
                          heatmap_legend_param = list( #THIS IS FOR THE LONG MINUS SIGNS
                            at = c(2, -5, -10, -15),   # Legend positions
                            labels = labels
                          )     # Adjusted labels
)
heatmap_raiska

# Methylation data heatmap ##################################################
met_data <- merged_df[, c("patient_id_aud", "ALX4" ,"CDX2","ARID1A_met", "HOPX")]
met_data <- as.matrix(met_data)
rownames(met_data) <- met_data[, 1]
met_data <- met_data[, -1]
colnames(met_data) <- c("ALX4 ", "CDX2", "ARID1A met", "HOPX")
met_data <- ifelse(met_data == 1, "metilintas", "nemetilintas") # recode for matrix
heatmap_metilinimas <- Heatmap(met_data, name = "Promoter methylation status")
#metilinimo spalvos
col_fun_met <- colorRampPalette(c("#e088bd", "#70a1fb"))(2)
heatmap_metilinimas <- Heatmap(met_data, name = "Promotorių metilinimo statusas",
                               col = col_fun_met, 
                               column_names_gp = gpar(fontface = "italic"), 
                               column_title = "Promotorių metilinimo statusas", 
                               #row_order = order_list
)
heatmap_metilinimas   

# list of heatmap: add the two together
htlist <- heatmap_metilinimas + heatmap_raiska

htlist 

draw(htlist, row_km = 1, row_split = clinical$Tumor)

#SAVE PNG
png("C:/Users/Ieva/rprojects/outputs_all/DISS/HEATMAP_ALL_data_20250825.png", width = 8000, height = 5500,
    res = 510, units = "px", pointsize = 12) # width and height in pixels, resolution in dpi
draw(htlist, row_km = 1, row_split = clinical$Tumor)# Render the heatmap
dev.off() # Close the PNG device

#Methylation only ####################################################
met_data1 <- met_data
colnames(met_data1)[3] <- "ARID1A"
heatmap_metilinimas_only <- Heatmap(met_data1, name = "Promotorių metilinimo statusas",
                                    col = col_fun_met, 
                                    right_annotation = row_ha,
                                    row_split = clinical$Tumor, 
                                    column_names_gp = gpar(fontface = "italic")
)
heatmap_metilinimas_only

#SAVE PNG
png("C:/Users/Ieva/rprojects/outputs_all/DISS/heatmap_output_met20251001.png", width =4000, height = 5500,
    res = 510, units = "px", pointsize = 12) # width and height in pixels, resolution in dpi
draw(heatmap_metilinimas_only)# Render the heatmap
dev.off() # Close the PNG device

#Notch, wnt, arid1a gene expression only ###################################################
colnames(Heat_data)
heatmap_raiska1 <- Heatmap(as.matrix(Heat_data[, 11:20]),cluster_columns = TRUE,
                           name = "Santykinė genų raiška",  
                           right_annotation = row_ha, col = col_fun, 
                           row_order = order_list,
                           row_split = clinical$Tumor, 
                           column_names_gp = gpar(fontface = "italic"),
                           column_title = "Santykinė genų raiška", 
                           row_names_gp = gpar(fontsize = 8), 
                           heatmap_legend_param = list( #THIS IS FOR THE LONG MINUS SIGNS
                             at = c(2, -5, -10, -15),   # Legend positions
                             labels = labels
                           )     # Adjusted labels
)
heatmap_raiska1

#SAVE PNG
png("C:/Users/Ieva/rprojects/outputs_all/DISS/heatmap_output_expr1_20251001.png", width =7000, height = 5500,
    res = 510, units = "px", pointsize = 12) # width and height in pixels, resolution in dpi
draw(heatmap_raiska1)# Render the heatmap
dev.off() # Close the PNG device

#Other gene expression only##########################################################
colnames(Heat_data)
heatmap_raiska2 <- Heatmap(as.matrix(Heat_data[, 1:10]),cluster_columns = TRUE,
                           name = "Santykinė genų raiška",  
                           right_annotation = row_ha, col = col_fun, 
                           row_order = order_list,
                           row_split = clinical$Tumor, 
                           column_names_gp = gpar(fontface = "italic"),
                           column_title = "Santykinė genų raiška", 
                           row_names_gp = gpar(fontsize = 8), 
                           heatmap_legend_param = list( #THIS IS FOR THE LONG MINUS SIGNS
                             at = c(2, -5, -10, -15),   # Legend positions
                             labels = labels
                           )     # Adjusted labels
)
heatmap_raiska2

#SAVE PNG
png("C:/Users/Ieva/rprojects/outputs_all/DISS/heatmap_output_expr2_20251001.png", width =7000, height = 5500,
    res = 510, units = "px", pointsize = 12) # width and height in pixels, resolution in dpi
draw(heatmap_raiska2)# Render the heatmap
dev.off() # Close the PNG device

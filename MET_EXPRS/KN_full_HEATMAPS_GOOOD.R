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
#fix clinicals
#fix na so it shows on the legend
merged_df$CA125_f <- factor(merged_df$CA125_f, levels = c("Padidėjimas", "Norma", "NA"), exclude = NULL)
merged_df$STAGE4 <- replace(merged_df$STAGE4, is.na(merged_df$STAGE4), "NA" )
merged_df$Grade.x <-replace(merged_df$Grade.x, is.na(merged_df$Grade.x), "NA" )
merged_df <- merged_df %>%
  mutate(CA125_status_post_op = case_when(
    `CA125 po gydymo` > 10 ~ "Padidėjimas",
    `CA125 po gydymo` <= 10 ~ "Norma",
    TRUE ~ NA_character_   # keep NA values
  ))
merged_df$CA125_status_post_op
merged_df$`CA125 po gydymo`
merged_df$CA125_status_post_op <-replace(merged_df$CA125_status_post_op,
                                         is.na(merged_df$CA125_status_post_op), "NA" )
#get survival data as of 2025-09-11
SURVIVAL_KN <- openxlsx::read.xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/KN_MIRTIES_FAILAS_20250911.xlsx")
#make only surv df
SURV <- SURVIVAL_KN[, c(2, 3,20, 21)]
head(SURV)

#join with main data
merged_df <- left_join(merged_df, SURV, by = "patient_id_aud")
#fix status
merged_df <- merged_df %>%
  mutate(STATUS  = case_when(
    STATUS  == 0 ~ "Gyvos",
    STATUS  == 1 ~ "Mirę",
    TRUE ~ NA_character_   # keep NA values
  ))
merged_df$STATUS

merged_df$STATUS <-replace(merged_df$STATUS,is.na(merged_df$STATUS), "NA" )

# clinical data annotation
col_age <- colorRamp2(c(40, 90), c( "#9cd4c4", "#3c402f")) #age colors
row_ha = rowAnnotation(Histologija = merged_df$Histology,
                       Navikas = merged_df$Tumor, Stadija = merged_df$STAGE4,
                       `Diferenciacijos laipsnis` = merged_df$Grade.x,
                       Amžius = merged_df$Age, CA125 = merged_df$CA125_f, 
                       `CA125 po gydymo` = merged_df$CA125_status_post_op, 
                       `Mirties statusas` = merged_df$STATUS,
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
# Methylation data heatmap ##################################################
met_data <- merged_df[, c("patient_id_aud", "ALX4" ,"CDX2","ARID1A_met", "HOPX")]
met_data <- as.matrix(met_data)
rownames(met_data) <- met_data[, 1]
met_data <- met_data[, -1]
colnames(met_data) <- c("ALX4 ", "CDX2", "ARID1A met", "HOPX")
met_data <- ifelse(met_data == 1, "metilintas", "nemetilintas") # recode for matrix
#metilinimo spalvos
col_fun_met <- colorRampPalette(c("#e088bd", "#70a1fb"))(2)
heatmap_metilinimas <- Heatmap(met_data, name = "Promotorių metilinimas",
                               right_annotation = row_ha, 
                               col = col_fun_met, 
                               column_names_gp = gpar(fontface = "italic"), 
                               row_split = merged_df$Tumor)
heatmap_metilinimas


# All expression heatmap #######################################################
Heat_data <- merged_df[, c("patient_id_aud",
                           expression 
)]
Heat_data <- as.data.frame(Heat_data)
rownames(Heat_data) <- Heat_data[, 1]
Heat_data <- Heat_data[, -1]
heatmap_raiska <- Heatmap(as.matrix(Heat_data))

#expression colors
col_fun = colorRamp2(c(2, -5, -10, -15), c("#8564fb",  "#64b3fb","#e088bd", "#af2745"))
#col_fun = colorRamp2(c(2, 0, -2, -4, -6, -8, -10, -12, -14, -16),
#                     c("#e7e0fe", "#cec1fd", "#8564fb", "#64b3fb","#93cafc","#325a7e", "#e088bd", "#af2745", "#9e233e", "#4f121f"))
labels <- c("\u221215", "\u221210", "\u22125", "2" ) #THIS IS NEEDED FOR MDPI AT LEAST - LONG MINUS SIGNS
#order if histology
order_hist <- c( "Riziką mažinanti operacija", "Mioma",  "Cista", "Endometriozė", "Endometrioidinė",
                 "Granuloza", "Mucininis", "Šviesių lastelių", "Serozinis", "HGSOC" )
#final heatmap
heatmap_raiska <- Heatmap(as.matrix(Heat_data),cluster_columns = TRUE,
                          name = "Santykinė genų raiška",  
                          #right_annotation = row_ha,
                          col = col_fun, 
                          row_split = merged_df$Tumor, 
                          column_names_gp = gpar(fontface = "italic"),
                          column_title = "Santykinė genų raiška", 
                          row_names_gp = gpar(fontsize = 8), 
                          heatmap_legend_param = list( #THIS IS FOR THE LONG MINUS SIGNS
                            at = c(2, -5, -10, -15),   # Legend positions
                            labels = labels
                          )     # Adjusted labels
)
heatmap_raiska


# list of heatmap: add the two together
htlist <- heatmap_raiska +heatmap_metilinimas 

htlist 

draw(htlist, row_km = 1, row_split = merged_df$Tumor)

#SAVE PNGs#########################################
png("C:/Users/Ieva/rprojects/outputs_all/DISS/heatmap_output_ALL20251204.png", width =4000, height = 5500,
    res = 510, units = "px", pointsize = 12) # width and height in pixels, resolution in dpi
draw(htlist, row_km = 1, row_split = merged_df$Tumor)
dev.off() # Close the PNG device

png("C:/Users/Ieva/rprojects/outputs_all/DISS/heatmap_output_METHILATION20251204.png", width =3000, height = 5500,
    res = 510, units = "px", pointsize = 12) # width and height in pixels, resolution in dpi
draw(heatmap_metilinimas)
dev.off() # Close the PNG device

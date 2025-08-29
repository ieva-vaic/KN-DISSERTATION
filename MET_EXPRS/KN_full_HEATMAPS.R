#KN-DISSERTATION project. MET_EXPRS - methylation, expression, mutation data
#HEATMAPs
Sys.setenv(LANG = "en")
#libraries
library(tidyverse)
library(readxl)
library(ComplexHeatm)
library(circlize)
#data upload##################################################################
merged_df <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/merged_all_data_diss.RDS")
# All expression heatmap #######################################################
Heat_data <- merged_df[, c("KN",
                           expression        
                           
)]
Heat_data <- as.data.frame(Heat_data)
rownames(Heat_data) <- Heat_data[, 1]
Heat_data <- Heat_data[, -1]
heatmap_raiska <- Heatmap(as.matrix(Heat_data))
#Clinical data###################################################################
clinical <- merged_df[, c("KN", "Histology", "Tumor", "CA125_f" ,"STAGE4", "Grade.x", "Age")]
clinical <- as.data.frame(clinical)
rownames(clinical) <- clinical$KN
clinical <- clinical[, -1]
head(clinical)
#fix na so it shows on the legend
clinical$CA125_f <- factor(clinical$CA125_f, levels = c("Padidėjimas", "Norma", "NA"), exclude = NULL)
clinical$STAGE4 <- replace(clinical$STAGE4, is.na(clinical$STAGE4), "NA" )
clinical$Grade.x <-replace(clinical$Grade.x, is.na(clinical$Grade.x), "NA" )

# clinical data annotation
col_age <- colorRamp2(c(40, 90), c( "#9cd4c4", "#3c402f")) #age colors
row_ha = rowAnnotation(Histologija = clinical$Histology,
                       Navikas = clinical$Tumor, Stadija = clinical$STAGE4,
                       `Diferenciacijos laipsnis` = clinical$Grade.x,
                       Amžius = clinical$Age, CA125 = clinical$CA125_f,
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
met_data <- KN_data[, c("patient_id_aud", "ALX4" ,"CDX2","ARID1A_met", "HOPX")]
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
png("C:/Users/Ieva/rprojects/outputs_all/DISS/heatmap_output_met20250605.png", width =4000, height = 5500,
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
png("C:/Users/Ieva/rprojects/outputs_all/DISS/heatmap_output_expr1_20250604.png", width =7000, height = 5500,
    res = 510, units = "px", pointsize = 12) # width and height in pixels, resolution in dpi
draw(heatmap_raiska1)# Render the heatmap
dev.off() # Close the PNG device

#Other gene expression only##########################################################
colnames(Heat_data)
col_fun2 = colorRamp2(c(2, -5, -10, -15), c("#64fbb3","#88e0bd","#7EC8E3", "#2745af"))
heatmap_raiska2 <- Heatmap(as.matrix(Heat_data[, 1:10]),cluster_columns = TRUE,
                           name = "Santykinė genų raiška",  
                           right_annotation = row_ha, col = col_fun2, 
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
png("C:/Users/Ieva/rprojects/outputs_all/DISS/heatmap_output_expr2_20250604.png", width =7000, height = 5500,
    res = 510, units = "px", pointsize = 12) # width and height in pixels, resolution in dpi
draw(heatmap_raiska2)# Render the heatmap
dev.off() # Close the PNG device
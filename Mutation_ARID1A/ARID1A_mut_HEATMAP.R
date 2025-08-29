##KN-DISSERTATION project. Mutation data (ARID1A and CTNNB1 genes)
#HEATMAP with mutation data
#libraries
library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(ggpubr)
#set wd for ploting
setwd("C:/Users/Ieva/rprojects/outputs_all/DISS/")
#upload data #############################################################
KN_data <- read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/audiniu_mut_exprs_met20250709.xlsx")
#get arid1a data and other mutation data
colnames(KN_data)
Arid1a_columns <- c("Histology", "Grupė_Ieva", "ARID1A","CTNNB1",
                    "ARID1A_met", "ARID1A_tumor_mut","ARID1A_tumor_VUS", "ARID1A_tumor_type",
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
# make "NA" 
ARID1A_df$ARID1A_tumor_mut[is.na(ARID1A_df$ARID1A_tumor_mut)] <- "NA"
ARID1A_df$CTNNB1_tumor_mut[is.na(ARID1A_df$CTNNB1_tumor_mut)] <- "NA"
#fix mutations - types
ARID1A_df$ARID1A_tumor_type[is.na(ARID1A_df$ARID1A_tumor_type)] <- "NA"
ARID1A_df$CTNNB1_tumor_mut2 <- ifelse(ARID1A_df$CTNNB1_tumor_mut == "Mutacija", "Missense",
                              ARID1A_df$CTNNB1_tumor_mut )
ARID1A_df$ARID1A_tumor_type2 <- ifelse(ARID1A_df$ARID1A_tumor_type == "frameshiftDeletion", "Rėmelio poslinkio iškrita",
                                      ARID1A_df$ARID1A_tumor_type )
ARID1A_df$ARID1A_tumor_type2 <- ifelse(ARID1A_df$ARID1A_tumor_type2 == "NA", "Be mutacijų",
                                       ARID1A_df$ARID1A_tumor_type2 )
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

ARID1A_df[row.names(ARID1A_df) %in% not_mut, "ARID1A_tumor_type2"] <- NA
#heatmap of expression###########################################
expression3 <- c("ARID1A","CTNNB1")
Heat_data <- KN_data[, c("KN",
                         expression3        
                         
)]
Heat_data <- as.data.frame(Heat_data)
rownames(Heat_data) <- Heat_data[, 1]
Heat_data <- Heat_data[, -1]
heatmap_raiska <- Heatmap(as.matrix(Heat_data))
#show expression
heatmap_raiska

#heatmap clinical data#########################################
colnames(ARID1A_df)
# clinical data annotation
row_ha = rowAnnotation(`ARID1A mutacija` = ARID1A_df$ARID1A_tumor_type2,
                       `CTNNB1 mutacija` = ARID1A_df$CTNNB1_tumor_mut2,
                       `ARID1A metilinimas` = ARID1A_df$ARID1A_met,
                       Histologija = ARID1A_df$Histology,
                       Navikas = ARID1A_df$Grupė_Ieva, 
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
                                  `ARID1A metilinimas` = c("Nemetilintas" = "#9cd4c4",  
                                                           "Metilintas" = "#c8d49c"), 
                                  `ARID1A mutacija` = c("Be mutacijų"="#9cd4c4", 
                                                        "missense" = "#a89cd4", 
                                                        "nonsense" = "darkblue",
                                                        "Rėmelio poslinkio iškrita" = "maroon",
                                                        "NA" = "grey"),
                                  `CTNNB1 mutacija` = c("Be mutacijų"="#9cd4c4", 
                                                        "Missense" = "#a89cd4",
                                                        "NA" = "grey")
                       ),
                       annotation_legend_param = list(
                         Histologija = list(title_gp = gpar(fontface = "italic")),
                         Navikas = list(title_gp = gpar(fontface = "italic")),
                         `ARID1A metilinimas` = list(title_gp = gpar(fontface = "italic")),
                         `ARID1A mutacija` = list(title_gp = gpar(fontface = "italic")),
                         `CTNNB1 mutacija` = list(title_gp = gpar(fontface = "italic"))
                       ),
                       annotation_name_gp = gpar(fontface = "italic")
)
#expression colors
col_fun = colorRamp2(c(2, -5, -10, -15), c("#8564fb",  "#64b3fb","#e088bd", "#af2745"))
#col_fun = colorRamp2(c(2, 0, -2, -4, -6, -8, -10, -12, -14, -16),
#                     c("#e7e0fe", "#cec1fd", "#8564fb", "#64b3fb","#93cafc","#325a7e", "#e088bd", "#af2745", "#9e233e", "#4f121f"))
labels <- c("\u221215", "\u221210", "\u22125", "2" ) #THIS IS NEEDED FOR MDPI AT LEAST - LONG MINUS SIGNS

#final heatmap#####################################
heatmap_raiska <- Heatmap(as.matrix(Heat_data),cluster_columns = TRUE,
                          name = "Santykinė genų raiška",  
                          right_annotation = row_ha, col = col_fun, 
                          row_split = ARID1A_df$ARID1A_tumor_mut, 
                          column_names_gp = gpar(fontface = "italic"),
                          column_title = "Santykinė genų raiška", 
                          row_names_gp = gpar(fontsize = 8), 
                          heatmap_legend_param = list( #THIS IS FOR THE LONG MINUS SIGNS
                            at = c(2, -5, -10, -15),   # Legend positions
                            labels = labels, title_gp = gpar(fontface = "italic")
                          )     # Adjusted labels
)
heatmap_raiska
#save png
png("heatmap_mut20250711_longwise.png", width = 2000, height = 2500,
    res = 300, units = "px", pointsize = 12) # width and height in pixels, resolution in dpi
draw(heatmap_raiska)# Render the heatmap
dev.off() # Close the PNG device

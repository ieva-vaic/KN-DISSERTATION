#KN-DISSERTATION project. 10 gene panel - genes selected via statistical analysis of TCGA and GTEx data
#KN 10 gene analysis data (cv 2 %, one threshold, 65 cases!, restored KN-59) 2025-02-17
#HETAMAP
Sys.setenv(LANG = "en")
#libraries
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
#upload data#########################################################################
OC_full <- readRDS("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/OC_10_genes_clean_2025_02_14.RDS")
dim(OC_full) #65 cases
expression <- c("EXO1", "RAD50","PPT2", "LUC7L2","PKP3", "CDCA5","ZFPL1","VPS33B", "GRB7","TCEAL4")
#heatmap form expression data####################################
Heat_data <- OC_full[, c("KN",
                        expression        
                        
)]
Heat_data <- as.data.frame(Heat_data)
rownames(Heat_data) <- Heat_data[, 1]
Heat_data <- Heat_data[, -1]
heatmap_raiska <- Heatmap(as.matrix(Heat_data))
# clinical data as df##########################################################
clinical <- OC_full[, c("KN", "type", "tumor", "CA125" ,"Stage4", "Grade2", "Age", "CA125_f")]
clinical <- as.data.frame(clinical)
rownames(clinical) <- clinical$KN
head(clinical)
#fix na so it shows on the legend
clinical$CA125_f <- factor(clinical$CA125_f, levels = c("CA125 increase", "Norm", "NA"), exclude = NULL)
clinical$CA125_f <- recode(clinical$CA125_f, `CA125 increase`= "CA125 padidėjimas", Norm = "Norma")
clinical$Stage4 <- replace(clinical$Stage4, is.na(clinical$Stage4), "NA" )
clinical$Grade2 <-replace(clinical$Grade2, is.na(clinical$Grade2), "NA" )
clinical$type <- factor(clinical$type)
clinical$type <- recode(clinical$type, `Benign`= "Gerybiniai", Other = "Kiti KV")

##add more clinical and experimental data form another df (mostly for histology info)
OC_clinical_full <- readxl::read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/KN_2024-11-04-data.xlsx")
OC_clinical_full <- OC_clinical_full[ !(OC_clinical_full$patient_id_aud %in% "KN-034"), ] #remove uneeded case

#put together the dfs
rownames(OC_clinical_full) <-OC_clinical_full$patient_id_aud
OC_clinical_full$KN <- OC_clinical_full$patient_id_aud
OC_clinical_full$KN <- gsub("^(KN-)0+", "\\1", OC_clinical_full$KN) 
OC_clinical_full$KN %in% clinical$KN
OC_clinical_full <- OC_clinical_full[, colnames(OC_clinical_full) %in% c("KN", "Histology")]
merged_df <- left_join(clinical, OC_clinical_full, by = "KN")
head(merged_df)

#remane histology info to lithuanian
merged_df$Histology <- ifelse(merged_df$Histology %in% c("Endometrial", "Endometriod"), "Endometrioidinis",
                              merged_df$Histology )
merged_df$Histology <- ifelse(merged_df$Histology == "Endometriois", "Endometriozė",
                              merged_df$Histology )
merged_df$Histology <- ifelse(merged_df$Histology == "Cystis", "Cista",
                              merged_df$Histology )
merged_df$Histology <- ifelse(merged_df$Histology == "Mucinous", "Mucininis",
                              merged_df$Histology )
merged_df$Histology <- ifelse(merged_df$Histology == "Serous", "Serozinis",
                              merged_df$Histology )
merged_df$Histology <- ifelse(merged_df$Histology == "Clear cell", "Šviesių lastelių",
                              merged_df$Histology )
merged_df$Histology <- ifelse(merged_df$Histology == "RSS", "Riziką mažinanti operacija",
                              merged_df$Histology )
merged_df$Histology <- ifelse(merged_df$Histology == "Granuloza", "Granulosa",
                              merged_df$Histology )

# clinical data annotation###########################
col_age <- colorRamp2(c(40, 90), c( "#9cd4c4", "#3c402f")) #age colors
row_ha = rowAnnotation(Histologija = merged_df$Histology,
                       Navikas = merged_df$type, Stadija = merged_df$Stage4,
                       `Diferenciacijos laipsnis` = merged_df$Grade2,
                        Amžius = merged_df$Age, CA125 = merged_df$CA125_f,
                        
                       col = list(Navikas = c("HGSOC" = "#a89cd4",
                                                "Gerybiniai" = "#d49cac", 
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
                                            "CA125 padidėjimas" = "#3c402f",
                                            "NA" = "grey"),
                                  Amžius = col_age,
                                  Histologija = c("Šviesių lastelių" = "lightblue",
                                                  "Cista" = "lightgreen", 
                                                  "Endometrioidinis" = "green",
                                                  "Endometriozė" = "darkgreen", 
                                                  "Granulosa" = "turquoise",
                                                  "HGSOC" = "deeppink",
                                                  "Mioma" = "red", 
                                                  "Mucininis" = "yellow",
                                                  "Riziką mažinanti operacija" = "orange", 
                                                  "Serozinis" = "lightpink")
                       ))


# disease type, row order
clinical$type <- factor(clinical$type, levels = c("Gerybiniai",  "Kiti KV", "HGSOC"))

#expression colors
col_fun = colorRamp2(c(2, -5, -10, -15), c("#8564fb",  "#64b3fb","#e088bd", "#af2745"))
#col_fun = colorRamp2(c(2, 0, -2, -4, -6, -8, -10, -12, -14, -16),
#                     c("#e7e0fe", "#cec1fd", "#8564fb", "#64b3fb","#93cafc","#325a7e", "#e088bd", "#af2745", "#9e233e", "#4f121f"))
labels <- c("\u221215", "\u221210", "\u22125", "2" ) #THIS IS NEEDED FOR MDPI AT LEAST - LONG MINUS SIGNS

#final heatmap####################################################
heatmap_raiska <- Heatmap(as.matrix(Heat_data), 
                          cluster_rows = F,
                          name = "Santykinė genų raiška",  
                          right_annotation = row_ha,
                          col = col_fun, 
                          row_split = clinical$type, 
                          column_names_gp = gpar(fontface = "italic"),
                          column_title = "Santykinė genų raiška", 
                          row_names_gp = gpar(fontsize = 8), 
                          #row_order = rows_order,
                          heatmap_legend_param = list( #THIS IS FOR THE LONG MINUS SIGNS
                            at = c(2, -5, -10, -15),   # Legend positions
                            labels = labels)     # Adjusted labels
)
heatmap_raiska

#save png######################################3
png("C:/Users/Ieva/rprojects/outputs_all/DISS/heatmap_10genes_65_cases_2025-06-23.png", width = 3000, height = 2500, res = 300) # width and height in pixels, resolution in dpi
heatmap_raiska
dev.off() # Close the PNG device

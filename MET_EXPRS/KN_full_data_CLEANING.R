#KN-DISSERTATION project. MET_EXPRS - methylation, expression, mutation data
#CLEAN data
Sys.setenv(LANG = "en")
#libraries
library(tidyverse)
library(readxl)
#data upload##################################################################
KN_data <- read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/KN_2024-11-04-data.xlsx")
#remove the single endometrial cancer case
KN_data <- KN_data[ !(KN_data$patient_id_aud %in% "KN-034"), ]
#remove the t4 form kn-76 (RSS)case
idx <- which(KN_data$patient_id_aud == "KN-076" & KN_data$`FIGO STAGE` == "IV")
KN_data[idx, c("FIGO STAGE", "Node", "Metastasis")] <- NA

#expression genes list
expression <- c("EXO1", "RAD50","PPT2", "LUC7L2","PKP3", "CDCA5","ZFPL1","VPS33B", "GRB7","TCEAL4", 
                "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "ARID1A", "CTNNB1", "FBXW7", "JAG2", "DLL1", "HES1")

#second df with clinical data###############################################
OC_clinical_full <- readxl::read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/OC_10_genes_ALTERATIVE_w_brca_2025_03-13.xlsx")
#remove the t4 form kn-76 (RSS)case
idx <- which(OC_clinical_full$KN == "KN-76" & OC_clinical_full$STAGE == "IV")
OC_clinical_full[idx, c("STAGE", "L/M", "MTS")] <- NA

#put together the dfs#######################################################
OC_clinical_full$KN
KN_data$KN <- gsub("^(KN-)0+", "\\1", KN_data$patient_id_aud) 
KN_data$KN %in% OC_clinical_full$KN

merged_df <- left_join(KN_data, OC_clinical_full, by = "KN")
colnames(merged_df)

# change names to lithuanian#################################################
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
table(merged_df$Histology, useNA = "a")
merged_df$Tumor <- factor(merged_df$Tumor)
levels(merged_df$Tumor) <- c("Gerybinis", "HGSOC", "Kiti KV")

#reorder levels#####################################################################
merged_df$Tumor <- factor(merged_df$Tumor, levels = c("Gerybinis", "Kiti KV", "HGSOC" ))
merged_df$CA125_f <- factor(merged_df$CA125_f)
levels(merged_df$CA125_f) <- c("Padidėjimas", "Norma")
merged_df$STAGE4 <- ifelse(merged_df$STAGE == "IA&IA", "IA", merged_df$STAGE)
merged_df$STAGE4 <- gsub("^([IVX]+)[A-Z]$", "\\1", merged_df$STAGE4)
merged_df$Grade.x <- ifelse(merged_df$Grade.x == "G2&G1", "G1", merged_df$Grade.x)
merged_df$Grade.x <- ifelse(merged_df$Grade.x == "G1&G1", "G1", merged_df$Grade.x)
merged_df$Grade.x <- ifelse(merged_df$Grade.x %in% c("GB", "GL"), NA, merged_df$Grade.x)
table(merged_df$Grade.x, useNA = "a")

#remove uneeded info #############################################################
colnames(merged_df)
remove_rows <- c("Vardas ir pavarde", "Operacijos data", "Pacientas", "Laboratorinis kodas",                
                 "Komentarai", "Kraujas" , "Šlapimas" , "Nuoplovos", "Audinys" , 
                 "Mėginio gavimo data", "CA125.y", "Ca 125 po gydymo", "Grade.y")
merged_df <- merged_df[, !c(colnames(merged_df) %in% remove_rows)]

#save###########################################################################
saveRDS(merged_df, "C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/merged_all_data_diss.RDS")

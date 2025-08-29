#KN-DISSERTATION project. 10 gene panel - genes selected via statistical analysis of TCGA and GTEx data
#KN 10 gene analysis data (cv 2 %, one threshold) 2025-03-13 (full brca data)
#FIX CLINICAL data, alterantive features
#libraries
library(tidyverse)
library(readxl)
#upload data##################################################################
#data -> only gene expression of 10 genes
OC_data <- read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/OC_data_IEVA_one_thr_cvs_fixed59_20250214.xlsx")
#expression genes list
expression <- c("EXO1", "RAD50", "PPT2", "LUC7L2", "PKP3", "CDCA5", "ZFPL1", "VPS33B", "GRB7", "TCEAL4")
colnames(OC_data)[4:13] <- expression #rename columns
#clinical data alternative
OC_clinical_full <- readxl::read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/KN_clinical_2025-03-13.xlsx")
dim(OC_clinical_full)
colnames(OC_clinical_full)

#merge expression data with clinical ###################
#chek the names
OC_clinical_full$KN #clinical
OC_data$KN #no 0
diff <- setdiff( OC_data$KN, OC_clinical_full$KN ) #61 witch clinical
print(diff)
intersection <- intersect( OC_data$KN, OC_clinical_full$KN ) #61 witch clinical
print(intersection)
OC_clinical_slimmer <- OC_clinical_full %>%
  filter(OC_clinical_full$KN %in% intersection)
dim(OC_clinical_slimmer) #66
#perform join
OC_alt_full <- merge(OC_data, OC_clinical_slimmer,by = "KN" )
OC_alt_full <- OC_alt_full[, -2] #remove nr.
dim(OC_alt_full)
rownames(OC_alt_full) <-OC_alt_full$KN 
#need to remove the endometrial carcinoma case (KN- 34)
OC_alt_full <- OC_alt_full[!rownames(OC_alt_full) %in% "KN-34", ] #65 left
dim(OC_alt_full)
#fix some clinicals######################################
OC_alt_full$Tumor <- as.factor(OC_alt_full$Tumor)
OC_alt_full$OC <- fct_recode(OC_alt_full$Tumor, "OC" = "HGSOC", "OC" = "Others")
#brca fix remove 3
OC_alt_full$BRCA_mut_klinikinė_kodas <- factor(OC_alt_full$BRCA_mut_klinikinė_kodas)
table(OC_alt_full$BRCA_mut_klinikinė_kodas)
OC_alt_full$BRCA_mut_klinikinė_kodas[OC_alt_full$BRCA_mut_klinikinė_kodas == 3] <- NA
OC_alt_full$BRCA_mut_klinikinė_kodas <- droplevels(OC_alt_full$BRCA_mut_klinikinė_kodas)
table(OC_alt_full$BRCA_mut_klinikinė_kodas, useNA = "a")
#fix ca125
OC_alt_full$CA125[OC_alt_full$CA125 == "Neatlikta"] <- NA
OC_alt_full$CA125_f <- factor(ifelse(as.numeric(OC_alt_full$CA125) > 35, "Increase", "Norm"))
table(OC_alt_full$CA125_f )
OC_alt_full <- OC_alt_full %>% rename(Age = `Amžius diagnozės metu`)
OC_alt_full$Age

#save df############################################################
saveRDS(OC_alt_full, "C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/OC_10_genes_ALTERATIVE_w_brca_2025_03-13.RDS")
openxlsx::write.xlsx(OC_alt_full, "C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/OC_10_genes_ALTERATIVE_w_brca_2025_03-13.xlsx")


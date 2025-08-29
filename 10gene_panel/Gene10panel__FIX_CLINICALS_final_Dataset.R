#KN-DISSERTATION project. 10 gene panel - genes selected via statistical analysis of TCGA and GTEx data
#KN 10 gene analysis data (cv 2 %, one threshold) from 2025-02-14
#FIX CLINICAL data
#libraries
library(tidyverse)
library(readxl)
#Add files###########################################
#data -> only gene expression of 10 genes
OC_data <- read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/OC_data_IEVA_one_thr_cvs_fixed59_20250214.xlsx")
#expression genes list
expression <- c("EXO1", "RAD50", "PPT2", "LUC7L2", "PKP3", "CDCA5", "ZFPL1", "VPS33B", "GRB7", "TCEAL4")
colnames(OC_data)[4:13] <- expression #rename columns
#clinical data
OC_clinical <- readxl::read_xlsx("C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/OCclinical20250207.xlsx")
dim(OC_clinical)
#Clean clinical
colnames(OC_clinical)
#fix patients#################################################################
#KN-76 is benign but it has stage - remove rss clin data
OC_clinical <- OC_clinical %>% mutate(
  `FIGO STAGE` = if_else(Grupe_full == "RSS", NA, `FIGO STAGE`),
  Node = if_else(Grupe_full == "RSS", NA, Node),
  Metastasis = if_else(Grupe_full == "RSS", NA, Metastasis)
)
#fix clinical#################################################################
#Stage factors 
table(OC_clinical$`FIGO STAGE`, useNA = "a") #9 na
#stage 4 groups
OC_clinical <- OC_clinical %>%
  mutate(Stage4 = recode(`FIGO STAGE`, "IA&IA" = "I", "IA" = "I", "IB" = "I",
                         "IIA" = "II", "IIB" = "II", "IIIA" = "III",
                         "IIIB" = "III", "IIIC" = "III", "IVB" = "IV"))
table(OC_clinical$Stage4, useNA = "a")
#stage 2 groups
OC_clinical <- OC_clinical %>%
  mutate(Stage2 = recode(Stage4, "I" = "I&II", "II" = "I&II", "III" = "III&IV",
                         "IV" = "III&IV"))
table(OC_clinical$Stage2, OC_clinical$Grupė_Ieva, useNA = "a")
OC_clinical$Stage2 <- as.factor(OC_clinical$Stage2)
OC_clinical$Stage2 <- relevel(OC_clinical$Stage2, ref = "I&II")
table(OC_clinical$Stage2, useNA = "a")
#CA125 factor
OC_clinical <- OC_clinical %>%
  mutate(CA125_f = ifelse(CA125 > 35, "CA125 increase", "Norm"))
table(OC_clinical$CA125_f, useNA = "a") #6 norm
#Grade 2 groups#
OC_clinical <- OC_clinical %>%
  mutate(Grade2 = recode(Grade, "G1&G1" = "G1", "G2&G1" = "G1", "GB" = NA_character_, "GL" = NA_character_))
table(OC_clinical$Grade2, useNA = "a") #G1 7 vs G3 42
#create OC vs benign factor
OC_clinical$tumor <- ifelse(OC_clinical$Grupė_Ieva %in% c("HGSOC", "Other"),
                            "OC", OC_clinical$Grupė_Ieva )
OC_clinical$tumor <- factor(OC_clinical$tumor)
levels(OC_clinical$tumor) <- c("Benign", "OC")
table(OC_clinical$tumor, useNA = "a")
#remove random columns#
OC_clinical <- OC_clinical[, -c(1:6)]
#rename original columns #######################################
OC_clinical <- OC_clinical %>%
  rename(Age = Amžius, CA125_post_treatment = `CA125 po gydymo`,
         Stage_original = `FIGO STAGE`, Grade_original = Grade,
         Group3 = Grupė_Ieva) %>% #reorder
  select("KN","tumor", "Group3", "Age", "Stage4", "Stage2","CA125_f", "Grade2",
         "Grupe_full", "CA125","CA125_post_treatment", "Stage_original",
         "Grade_original", "Node", "Metastasis", "Histology")
#MERGE clinical with gene expression #####################################
OC_full <- merge(OC_data, OC_clinical,by = "KN" )
OC_full <- OC_full[, -2] #remove nr.
OC_full$type == OC_full$Group3 #chek type for good merge
OC_full <- OC_full[, !names(OC_full)  %in% "Group3"] #remove the equal column

#remove Endometrial cancer case#################################################
OC_full <- OC_full %>%
  filter(Grupe_full!= "Endometrial carcinoma") #65 cases left

#save df#####################################################################
saveRDS(OC_full, "C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/OC_10_genes_clean_2025_02_14.RDS")
openxlsx::write.xlsx(OC_full, "C:/Users/Ieva/rprojects/OTHER DATA/KN-DISSERTATION FILES/OC_10_genes_clean_2025_02_14.xlsx")
                     
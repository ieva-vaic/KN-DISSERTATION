#KN-DISSERTATION project. Genes selection via statistical analysis of TCGA and GTEx data 
#(only scripts producing lithuanian figures, included in the dissertation)
#THIS IS TCGA-OV-RISK-GENES project script No. 11
#HEATMAP 
# Load packages ##########################################
Sys.setenv(LANG = "en")
library(ComplexHeatmap)
library(tidyverse)
library(RColorBrewer) 
library(circlize)
#set directory of the data
setwd("C:/Users/Ieva/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/")
# Load train data ###################################
gtex_counts_train <- readRDS("train_gtcga_normcounts_prot_2025.RDS")
gtex_genes <- readRDS("gtcga_elastic_2025.RDS") #only the genes left after lasso
gtex_filtered_counts_train <- gtex_counts_train[colnames(gtex_counts_train) %in% gtex_genes] 
dim(gtex_filtered_counts_train) #489 samples
#get the rank 
res_coef_gtex_R <- readRDS("C:/Users/Ieva/rprojects/TCGA-OV-RISK-PROJECT/Public data RDSs/elastic_net_model_gtex_2025.RDS" )
RANK <- names(res_coef_gtex_R)
#get TCGA df
gtex_filtered_counts_train2 <- as.data.frame(t(gtex_filtered_counts_train))
gtex_filtered_counts_train2 <- gtex_filtered_counts_train2 %>%  dplyr::select(starts_with("TCGA")) 
dim(gtex_filtered_counts_train2)
#336 samples #214 genes - leaves only the genes left at the lasso step
#get GTEX df
gtex_filtered_counts_train3 <- as.data.frame(t(gtex_filtered_counts_train))
gtex_filtered_counts_train3 <- gtex_filtered_counts_train3 %>%  dplyr::select(starts_with("GTEX")) 
dim(gtex_filtered_counts_train3)
#336 samples #214 genes - leaves only the genes left at the lasso step
#Load full clinical data#######################################
pheno_full <- readRDS("joinedTCGA_XENA_clinical.RDS")
#reduce to train data
tcga_cases <- colnames(gtex_filtered_counts_train2)
sum(pheno_full$barcode %in% tcga_cases) #336
pheno_full <- pheno_full[pheno_full$barcode %in% tcga_cases, ]
dim(pheno_full)

#Leave only features of interest in clinical data##############################
pheno_best <- pheno_full[, (names(pheno_full)
                            %in% c("barcode", "STAGE", "race", "vital_status", 
                                   "ageatinitialpathologicdiagnosis", "neoplasmhistologicgrade", 
                                   "lymphaticinvasion", "treatment_type"  ))] # left
#fix up grade by removing GB, GX, and singular Grade 4 case
pheno_best$neoplasmhistologicgrade <- recode(pheno_best$neoplasmhistologicgrade,
                                             "GB" = NA_character_,
                                             "GX" = NA_character_,
                                             "G4" = NA_character_)
table(pheno_best$neoplasmhistologicgrade)
#rename clinical features to more conventional names
pheno_best <- pheno_best %>%
  rename(grade = neoplasmhistologicgrade,
         stage = STAGE,
         `vital status` = vital_status,
         `treatment type` = treatment_type,
         age = ageatinitialpathologicdiagnosis,
         `lymphatic invasion` = lymphaticinvasion)
colnames(pheno_best)
#fix up lymphatic invasion
pheno_best$`lymphatic invasion` <- factor(pheno_best$`lymphatic invasion`, levels = c(0 , 1 , NA ),
                                          labels = c("No invasion",
                                                     "Lymphatic invasion"
                                          ))

pheno_best$`lymphatic invasion` <- as.character(pheno_best$`lymphatic invasion`)
#add gtex to clinical data as NAs
pheno_best <- bind_rows(pheno_best, tibble(barcode = colnames(gtex_filtered_counts_train3)))
rownames(pheno_best) <- pheno_best$barcode
pheno_best <- pheno_best %>%
  arrange(match(barcode, rownames(gtex_filtered_counts_train)))
pheno_best <- pheno_best %>%
  mutate(across(everything(), ~replace(., is.na(.), "NA")))
#fix the numeric variable
pheno_best$age <- as.numeric(pheno_best$age)  

#CREATE WITH HEATMAP CLINICAL DATA FEATURES ######################
col_age <- colorRamp2(c(30, 100), c( "#9cd4c4", "#3c402f"))
row_ha = rowAnnotation(Race  = pheno_best$race, 
                       Age = pheno_best$age, 
                       `Vital status` = pheno_best$`vital status`,
                       `Treatment type` = pheno_best$`treatment type`,
                       Grade = pheno_best$grade,
                       Stage =pheno_best$stage,
                       `Lymphatic invasion` = pheno_best$`lymphatic invasion`,
                       #choose colors
                       col = list(Race = c("american indian or alaska native" = "pink",
                                           "asian" = "#9cd4c4",
                                           "black or african american" = "darkgreen",
                                           "native hawaiian or other pacific islander" = "turquoise",
                                           "white" = "lightgrey",
                                           "not reported" = "darkgrey", 
                                           "NA" = "grey"), 
                                  Age = col_age,
                                  `Vital status` = c("Alive" = "#9cd4c4",  
                                                     "Dead" = "#a89cd4", 
                                                     "NA" = "grey"),
                                  `Treatment type` = c("Pharmaceutical Therapy, NOS" = "#9cd4c4",  
                                                       "Radiation Therapy, NOS" = "#a89cd4", 
                                                       "NA" = "grey"),
                                  Grade = c("G2" = "#9cd4c4",  
                                            "G3" = "#a89cd4", 
                                            "NA" = "grey"),
                                  Stage = c("Stage I" = "#9cd4c4",  
                                            "Stage II" = "#a89cd4",
                                            "Stage III" = "pink",  
                                            "Stage IV" = "turquoise", 
                                            "NA" = "grey"),
                                  `Lymphatic invasion` = c("No invasion" = "#9cd4c4",  
                                                           "Lymphatic invasion" = "#a89cd4", 
                                                           "NA" = "grey")
                       ))

#HIGHLITE THE CHOSEN GENES####################################
#gene names:
expression <- c("EXO1", "RAD50","PPT2", "LUC7L2","PKP3", "CDCA5",
                "ZFPL1","VPS33B", "GRB7","TCEAL4")
expression %in% colnames(gtex_counts_train) #check if they are in the data
#chosen genes highlated in red
col_colors <- sapply(colnames(gtex_counts_train), function(x) {
  if (x %in% expression) {
    "red"  # Color selected columns red
  } else {
    "black"  # Default color
  }
})
#create groupings according to clinical data 
snames = rownames(gtex_counts_train)
group = substr(snames, 1, 4)
group = as.factor(group)
levels(group) <- c("GTEx", "TCGA-OV")

#GENERATE TCGA HEATMAP##############################################
# Generate heatmap with colored column names
gtex_counts_train <- data.matrix(gtex_filtered_counts_train)
heatmap_tcga <- Heatmap(as.matrix(gtex_filtered_counts_train), 
                        row_split = group,   
                        show_row_names = F,
                        column_names_gp = gpar(fontsize = 6, col = col_colors, fontface = "italic"), 
                        row_names_gp = gpar(fontsize = 2), # 
                        heatmap_legend_param = list(title = "Gene Expression"),
                        right_annotation = row_ha,
                        cluster_rows = F,
                        cluster_columns = F,
                        column_order = RANK)

#save tcga heatmap
png("C:/Users/Ieva/rprojects/outputs_all/DISS/gtextcga_heatmap2025-04-010.png",
    width = 6500, height = 4000, res = 400) # width and height in pixels, resolution in dpi
heatmap_tcga #
dev.off() # Close the PNG device

#CREATE SMALLER HEATMAP OF CHOSEN GENES#######################################
small_train_df <- gtex_filtered_counts_train[colnames(gtex_filtered_counts_train) %in% expression]
small_train_matrix <- data.matrix(small_train_df)
heatmap_tcga2 <- Heatmap(as.matrix(small_train_matrix), 
                         row_split = group,   
                         show_row_names = F,
                         column_names_gp = gpar(fontsize = 6, col = col_colors, fontface = "italic"), 
                         row_names_gp = gpar(fontsize = 2), # 
                         heatmap_legend_param = list(title = "Gene Expression"),
                         right_annotation = row_ha,
                         cluster_rows = F)
#save smaller heatmap
png("C:/Users/Ieva/rprojects/outputs_all/DISS/gtextcga_heatmap_small2025-04-010.png",
    width = 2000, height = 2500, res = 300) 
heatmap_tcga2 
dev.off()

#ADD GTEx AGE###########################################################
# Load gtex sample data ###################################
#dowloaded https://gtexportal.org/home/downloads/adult-gtex/metadata on 2025-04-18
metadata <- read.table("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# view column names: id, sex (2= Female),
#DTHHRDY = hardy sclase of death (1 = trauma (sudden), 2 = suden death (<1h) e.g. miocardal infactrion, 
#3 = patients that vere ill but death was unexpected (>24h),
#4 = long illness such as cancer, 0 = ventilator cases)
colnames(metadata)

#Filter for test data#################
#ids
ids <- colnames(gtex_filtered_counts_train3)
metadata$SUBJID
#trim ids in gtex_df_train
# Keep only the first two segments
donor_ids <- sub("^([^-]+-[^-]+).*", "\\1", ids)
# Result
donor_ids %in% metadata$SUBJID #all represented in metadata
#filter metadata
filtered_metadata <- metadata[metadata$SUBJID %in% donor_ids, ]
#look at age of gtex data
filtered_metadata$AGE
#make the same factors for tcga data:
pheno_best$AGE <- cut(pheno_best$age,
                      breaks = seq(30, 80, by = 10),   # From 30 to 80 in 10-year steps
                      right = FALSE,                  # e.g., 60 is in "60-69"
                      labels = c("30-39", "40-49", "50-59", "60-69", "70-79"))

pheno_best$AGE
#fix pheno_best
pheno_best$SUBJID <- sub("^([^-]+-[^-]+).*", "\\1", pheno_best$barcode)
#merge
pheno_best$AGE2 <- filtered_metadata$AGE[match(pheno_best$SUBJID, filtered_metadata$SUBJID)]
pheno_best$AGE2 <- factor(pheno_best$AGE2)
pheno_best$AGE2 <- ifelse(!is.na(pheno_best$AGE),
                          as.character(pheno_best$AGE),
                          as.character(pheno_best$AGE2))

# Convert back to factor with consistent levels (optional)
pheno_best$AGE2 <- factor(pheno_best$AGE2, levels = levels(pheno_best$AGE))
pheno_best$AGE2
#SAVE
#write.csv(pheno_best, "C:/Users/Ieva/rprojects/outputs_all/DISS/pheno_best202520619.csv")

#CREATE HEATMAP CLINICAL DATA FEATURES with GTEX AGE######################
row_ha2 = rowAnnotation(Race  = pheno_best$race, 
                        Age = pheno_best$AGE2, 
                        `Vital status` = pheno_best$`vital status`,
                        `Treatment type` = pheno_best$`treatment type`,
                        Grade = pheno_best$grade,
                        Stage =pheno_best$stage,
                        `Lymphatic invasion` = pheno_best$`lymphatic invasion`,
                        #choose colors
                        col = list(Race = c("american indian or alaska native" = "pink",
                                            "asian" = "#9cd4c4",
                                            "black or african american" = "darkgreen",
                                            "native hawaiian or other pacific islander" = "turquoise",
                                            "white" = "lightgrey",
                                            "not reported" = "darkgrey", 
                                            "NA" = "grey"), 
                                   Age = c("30-39" = "#E6D6F5",  
                                           "40-49" = "#C49DDE", 
                                           "50-59" = "#9C6CD3",  
                                           "60-69" = "#713AB6", 
                                           "70-79" = "#4B1C74", 
                                           "NA" = "grey"),
                                   `Vital status` = c("Alive" = "#9cd4c4",  
                                                      "Dead" = "#a89cd4", 
                                                      "NA" = "grey"),
                                   `Treatment type` = c("Pharmaceutical Therapy, NOS" = "#9cd4c4",  
                                                        "Radiation Therapy, NOS" = "#a89cd4", 
                                                        "NA" = "grey"),
                                   Grade = c("G2" = "#9cd4c4",  
                                             "G3" = "#a89cd4", 
                                             "NA" = "grey"),
                                   Stage = c("Stage I" = "#9cd4c4",  
                                             "Stage II" = "#a89cd4",
                                             "Stage III" = "pink",  
                                             "Stage IV" = "turquoise", 
                                             "NA" = "grey"),
                                   `Lymphatic invasion` = c("No invasion" = "#9cd4c4",  
                                                            "Lymphatic invasion" = "#a89cd4", 
                                                            "NA" = "grey")
                        ))

#MAKE BIG HEATMAP WITH GETX AGE################
heatmap_tcga4 <- Heatmap(as.matrix(gtex_filtered_counts_train), 
                         row_split = group,   
                         show_row_names = F,
                         column_names_gp = gpar(fontsize = 6, col = col_colors, fontface = "italic"), 
                         row_names_gp = gpar(fontsize = 2), # 
                         heatmap_legend_param = list(title = "Gene Expression"),
                         right_annotation = row_ha2,
                         cluster_rows = F)
#save big heatmap2 (with gtex age)
png("C:/Users/Ieva/rprojects/outputs_all/DISS/gtextcga_heatmap_GTEXAGE_train_BIG2025-04-23.png",
    width = 6500, height = 4000, res = 300) 
heatmap_tcga4 
dev.off() 

#lithuanize heatmap################################
#rename groups to lithuanian
pheno_best2 <- pheno_best
#race - rasė
pheno_best2$race <- as.factor(pheno_best2$race )
print(levels((pheno_best2$race)))
pheno_best2$race <- recode(pheno_best2$race,
                           "american indian or alaska native" = "Indėnai arba Aliaskos vietiniai",
                           "asian" = "Azijiečiai",
                           "black or african american" = "Juodaodžiai",
                           "native hawaiian or other pacific islander" = "Vietiniai Havajiečiai arba kiti Ramiojo vandenyno salų gyventojai",
                           "not reported" =   "neatsakė",  "white" = "Baltieji") 
#survival - išgyvenamumas
pheno_best2$`vital status` <- as.factor(pheno_best2$`vital status` )
print(levels((pheno_best2$`vital status`)))
pheno_best2$`vital status` <- recode(pheno_best2$`vital status`,
                                     "Alive" = "Gyvas",
                                     "Dead" = "Miręs")
#treatment type - gydymo tipas
pheno_best2$`treatment type`<- as.factor(pheno_best2$`treatment type`)
print(levels((pheno_best2$`treatment type`)))
pheno_best2$`treatment type` <- recode(pheno_best2$`treatment type`,
                                     "Pharmaceutical Therapy, NOS" = "Terapija vaistais",
                                     "Radiation Therapy, NOS"    = "Radioterapija")
#stage - stadija
pheno_best2$stage <- as.factor(pheno_best2$stage )
print(levels((pheno_best2$stage)))
pheno_best2$stage <- recode(pheno_best2$stage,
                            "Stage I" = "I Stadija",
                            "Stage II" = "II Stadija",
                            "Stage III" = "III Stadija",
                            "Stage IV" = "IV Stadija")
#lymphatic invasion - limfmazgiai
pheno_best2$`lymphatic invasion` <- as.factor(pheno_best2$`lymphatic invasion` )
print(levels((pheno_best2$`lymphatic invasion`)))
pheno_best2$`lymphatic invasion` <- recode(pheno_best2$`lymphatic invasion`,
                            "Lymphatic invasion" = "Invazija į limfmazgius",
                            "No invasion" = "Nėra invazijos į limfmazgius")
#FIX ANOTATION OF CLINICAL DATA - pirmiausia klinikinių anotaciją susitvarkyt reikia
row_ha2 = rowAnnotation(Rasė  = pheno_best2$race, 
                        Amžius = pheno_best2$AGE2, 
                        `Išgyvenamumo statusas` = pheno_best2$`vital status`,
                        `Gydymo tipas` = pheno_best2$`treatment type`,
                        `Diferenciacijos laipsnis` = pheno_best2$grade,
                        Stadija =pheno_best2$stage,
                        `Invazija į limfmazgius` = pheno_best2$`lymphatic invasion`,
                        #choose colors
                        col = list(Rasė = c("Indėnai arba Aliaskos vietiniai" = "pink",
                                            "Azijiečiai" = "#9cd4c4",
                                            "Juodaodžiai" = "darkgreen",
                                            "Vietiniai Havajiečiai arba kiti Ramiojo vandenyno salų gyventojai" = "turquoise",
                                            "Baltieji" = "lightgrey",
                                            "neatsakė" = "darkgrey", 
                                            "NA" = "grey"), 
                                   Amžius = c("30-39" = "#E6D6F5",  
                                           "40-49" = "#C49DDE", 
                                           "50-59" = "#9C6CD3",  
                                           "60-69" = "#713AB6", 
                                           "70-79" = "#4B1C74", 
                                           "NA" = "grey"),
                                   `Išgyvenamumo statusas` = c("Gyvas" = "#9cd4c4",  
                                                      "Miręs" = "#a89cd4", 
                                                      "NA" = "grey"),
                                   `Gydymo tipas` = c("Terapija vaistais" = "#9cd4c4",  
                                                        "Radioterapija" = "#a89cd4", 
                                                        "NA" = "grey"),
                                   `Diferenciacijos laipsnis` = c("G2" = "#9cd4c4",  
                                             "G3" = "#a89cd4", 
                                             "NA" = "grey"),
                                   Stadija = c("I Stadija" = "#9cd4c4",  
                                             "II Stadija" = "#a89cd4",
                                             "III Stadija" = "pink",  
                                             "IV Stadija" = "turquoise", 
                                             "NA" = "grey"),
                                   `Invazija į limfmazgius` = c("Nėra invazijos į limfmazgius" = "#9cd4c4",  
                                                            "Invazija į limfmazgius" = "#a89cd4", 
                                                            "NA" = "grey")
                        ))

#LT: MAKE BIG HEATMAP WITH GETX AGE################
#gene names:
expression <- c("EXO1", "RAD50","PPT2", "LUC7L2","PKP3", "CDCA5",
                "ZFPL1","VPS33B", "GRB7","TCEAL4")
expression %in% colnames(gtex_filtered_counts_train) #check if they are in the data
#chosen genes highlated in red
col_colors <- sapply(colnames(gtex_filtered_counts_train), function(x) {
  if (x %in% expression) {
    "red"  # Color selected columns red
  } else {
    "black"  # Default color
  }
})
#GENERATE LT HEATMAP ###############################
heatmap_tcga4 <- Heatmap(as.matrix(gtex_filtered_counts_train), 
                         row_split = group,   
                         show_row_names = F,
                         column_names_gp = gpar(fontsize = 6, col = col_colors, fontface = "italic"), 
                         row_names_gp = gpar(fontsize = 2), # 
                         heatmap_legend_param = list(title = "Genų raiška"),
                         right_annotation = row_ha2,
                         cluster_rows = F,
                         cluster_columns = F,
                         column_order = RANK)
#save lithuanian heatmap
png("C:/Users/Ieva/rprojects/outputs_all/DISS/big_elastic_net_heatmap20251208.png",
    width = 7000, height = 3000, res = 350) 
heatmap_tcga4 
dev.off()

#BONUS: count age groups###########################################
pheno_best$id <-  substr(trimws(pheno_best$barcode), 1, 4)

table(pheno_best$AGE2, pheno_best$id, useNA = "a")
table(pheno_best$id, useNA = "a")

#BONUS: chek if my other genes are in the 214 gene list
full214 <- rownames(gtex_filtered_counts_train2)
expression <- c("EXO1", "RAD50","PPT2", "LUC7L2","PKP3", "CDCA5","ZFPL1","VPS33B", "GRB7","TCEAL4", 
                "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "ARID1A", "CTNNB1", "FBXW7", "JAG2", "DLL1", "HES1")
expression[expression %in% full214] #not in the list

library(readxl)
# 定义路径
folder_path <- "F:/3_Gene/3_GWAS/lung_GWAS/our_data/"


file_paths <- list.files(folder_path, pattern = "\\.raw$", full.names = TRUE)


resnps_list <- lapply(file_paths, read.csv, sep="")

list2env(setNames(resnps_list, gsub(".raw", "", basename(file_paths))), envir = .GlobalEnv)

UKB_results_cn_fev1 <- read_excel("F:/3_Gene/3_GWAS/lung_GWAS/UKB_results_cn_fev1.xlsx")

merged_SNP <- Reduce(function(x, y) merge(x, y, by = "IID", all = TRUE), list(resnps3, resnps4, resnps6, resnps17, resnps22))

merged_SNP <- merged_SNP[, !grepl("\\.x", colnames(merged_SNP))]
merged_SNP <- merged_SNP[, !grepl("\\.y", colnames(merged_SNP))]

merged_SNP <- merged_SNP[, !colnames(merged_SNP) %in% c("FID","PAT","MAT","SEX","PHENOTYPE")]

SNP_deletGT = SNP_deletGT[, !colnames(SNP_deletGT) %in% c("IID")]


colnames(SNP_deletGT) <- sub("^(rs\\d+).*", "\\1", colnames(SNP_deletGT))


SNP_deletGT_colnames <- colnames(SNP_deletGT)


fev1_results_reorder <- UKB_results_cn_fev1[match(SNP_deletGT_colnames, UKB_results_cn_fev1$snp), ]

betas <- fev1_results_reorder$beta

prs <- rowSums(sweep(SNP_deletGT, 2, betas, FUN = "*"), na.rm = TRUE)

SNP_FEV1 = merged_SNP

SNP_FEV1$PRS <- prs
library(rcompanion)

PRS_表型 <- read_excel("F:/3_Gene/3_GWAS/lung_GWAS/our_data//PRS_表型.xlsx")
FamilyID <- read_excel("F:/3_Gene/3_GWAS/lung_GWAS/our_data//FamilyID.xlsx")

PRS_表型$height_BLOM = blom(PRS_表型$height,method = "blom")
PRS_表型$FEV1_BLOM = blom(PRS_表型$FEV1,method = "blom")
PRS_表型$FVC_BLOM = blom(PRS_表型$FVC1,method = "blom")

colnames(PRS_表型)[colnames(PRS_表型) == "no"] <- "IID"

merged_FEV1 <- merge(SNP_FEV1, PRS_表型, by = "IID", all.x = TRUE)
merged_FEV1 = merge(merged_FEV1, FamilyID, by = "IID", all.x = TRUE)

library(geepack)

fev1_result = geeglm(FEV1_BLOM ~ age  + gender + smoking + height_BLOM + PRS, id = Family_ID, data = merged_FEV1, corstr = "exchangeable" )
summary(fev1_result)

#######################fvc
UKB_results_cn_fvc <- read_excel("F:/3_Gene/3_GWAS/lung_GWAS/UKB_results_cn_fvc.xlsx")

snp_list <- UKB_results_cn_fvc$snp

SNP_FVC <- SNP_deletGT[, colnames(SNP_deletGT) %in% snp_list]

SNP_FVC

SNP_FVC_colnames <- colnames(SNP_FVC)

FVC_results_reorder <- UKB_results_cn_fvc[match(SNP_FVC_colnames, UKB_results_cn_fvc$snp), ]

betas_FVC <- FVC_results_reorder$beta

prs_FVC <- rowSums(sweep(SNP_FVC, 2, betas_FVC, FUN = "*"), na.rm = TRUE)

SNP_FVC$PRS <- prs_FVC

merged_FVC <- merge(SNP_FVC, PRS_表型, by = "IID", all.x = TRUE)
merged_FVC = merge(merged_FVC, FamilyID, by = "IID", all.x = TRUE)

library(geepack)

fVC_result = geeglm(FVC_BLOM ~ age  + gender + smoking + height_BLOM + PRS, id = Family_ID, data = merged_FVC, corstr = "exchangeable" )
summary(fVC_result)

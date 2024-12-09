##########combine all chr######
library(dplyr)
library(readr)

chr23 <- read.delim("chr23_012.raw",sep = " ")


list_file <- list.files(pattern = "*.raw") %>% 
  lapply(read.delim, stringsAsFactors=F,sep = " ") %>% 
  bind_cols

save(list_file,file = "chr1_22.RData") 

library(dplyr)

coln <- colnames(list_file)
coln

tests <- subset(list_file,select = -c(FID...17,IID...18,PAT...19,MAT...20,SEX...21,
                                      PHENOTYPE...22,FID...48,IID...49,PAT...50,MAT...51,
                                      SEX...52,PHENOTYPE...53,FID...62,IID...63,PAT...64,
                                      MAT...65,SEX...66,PHENOTYPE...67,FID...94,IID...95,
                                      PAT...96,MAT...97,SEX...98,PHENOTYPE...99,FID...112,
                                      IID...113,PAT...114,MAT...115,SEX...116,
                                      PHENOTYPE...117,FID...123,IID...124,PAT...125,
                                      MAT...126,SEX...127,PHENOTYPE...128,FID...131,
                                      IID...132,PAT...133,MAT...134,SEX...135,PHENOTYPE...136,
                                      FID...140,IID...141,PAT...142,MAT...143,SEX...144,
                                      PHENOTYPE...145,FID...169,IID...170,PAT...171,
                                      MAT...172,SEX...173,PHENOTYPE...174,FID...176,
                                      IID...177,PAT...178,MAT...179,SEX...180,PHENOTYPE...181,
                                      FID...188,IID...189,PAT...190,MAT...191,SEX...192,
                                      PHENOTYPE...193,FID...228,IID...229,PAT...230,MAT...231,
                                      SEX...232,PHENOTYPE...233,FID...236,IID...237,PAT...238,
                                      MAT...239,SEX...240,PHENOTYPE...241,FID...245,IID...246,
                                      PAT...247,MAT...248,SEX...249,PHENOTYPE...250,FID...333,
                                      IID...334,PAT...335,MAT...336,SEX...337,PHENOTYPE...338,
                                      FID...346,IID...347,PAT...348,MAT...349,SEX...350,
                                      PHENOTYPE...351,FID...369,IID...370,PAT...371,MAT...372,
                                      SEX...373,PHENOTYPE...374,FID...413,IID...414,PAT...415,
                                      MAT...416,SEX...417,PHENOTYPE...418,FID...439,IID...440,
                                      PAT...441,MAT...442,SEX...443,PHENOTYPE...444,FID...446,
                                      IID...447,PAT...448,MAT...449,SEX...450,PHENOTYPE...451))

colnames(tests)[1] <- "FID"

chr1_23=dplyr::left_join(tests,chr23,by="FID")

snp <- subset(chr1_23,select = -c(IID...2,PAT...3,MAT...4,SEX...5,PHENOTYPE...6,IID,PAT,MAT,SEX,PHENOTYPE))

save(snp,file = "SNP.RData")


####ukb basic data####
bd <- read.table("F:/3_Gene/3_GWAS/lung_GWAS/ukb_data/base data/ukb672468.tab", header=TRUE, sep="\t")
colnames(bd)
#############FVC##############

meancol <- c("f.3062.0.0","f.3062.0.1","f.3062.0.2")
bd$fvc_1 <- rowMeans ( bd[meancol],na.rm = T )

meancol <- c("f.3062.1.0","f.3062.1.1","f.3062.1.2")
bd$fvc_2 <- rowMeans ( bd[meancol],na.rm = T )

meancol <- c("f.3062.2.0","f.3062.2.1","f.3062.2.2")
bd$fvc_3 <- rowMeans ( bd[meancol],na.rm = T )

meancol <- c("f.3062.3.0","f.3062.3.1","f.3062.3.2")
bd$fvc_4 <- rowMeans ( bd[meancol],na.rm = T )


nrow(bd[bd$fvc_1 == 'NaN' & (bd$fvc_2 != 'NaN'| bd$fvc_3 != 'NaN' |bd$fvc_4 != 'NaN'), ])

rows_ok <- subset(bd,!(bd$fvc_1 == 'NaN' & (bd$fvc_2 != 'NaN'| bd$fvc_3 != 'NaN' |bd$fvc_4 != 'NaN')))
rows_na <- subset(bd,(bd$fvc_1 == 'NaN' & (bd$fvc_2 != 'NaN'| bd$fvc_3 != 'NaN' |bd$fvc_4 != 'NaN')))


rows_na$fvc <- ifelse(rows_na$fvc_1 != "NaN",rows_na$fvc_1,
                      ifelse(rows_na$fvc_2 != "NaN",rows_na$fvc_2,
                             ifelse(rows_na$fvc_3 != "NaN",rows_na$fvc_3,
                                    ifelse(rows_na$fvc_4 != "NaN",rows_na$fvc_4,NA))))
rows_ok$fvc <- ifelse(!is.na(rows_ok$fvc_1),rows_ok$fvc_1,NA)

rows_fvc <- rbind(rows_ok,rows_na)

################FEV1############

meancol <- c("f.3063.0.0","f.3063.0.1","f.3063.0.2")
rows_fvc$FEV1_1 <- rowMeans ( rows_fvc[meancol],na.rm = T )

meancol <- c("f.3063.1.0","f.3063.1.1","f.3063.1.2")
rows_fvc$FEV1_2 <- rowMeans ( rows_fvc[meancol],na.rm = T )

meancol <- c("f.3063.2.0","f.3063.2.1","f.3063.2.2")
rows_fvc$FEV1_3 <- rowMeans ( rows_fvc[meancol],na.rm = T )

meancol <- c("f.3063.3.0","f.3063.3.1","f.3063.3.2")
rows_fvc$FEV1_4 <- rowMeans ( rows_fvc[meancol],na.rm = T )

nrow(rows_fvc[rows_fvc$FEV1_1 == 'NaN' & (rows_fvc$FEV1_2 != 'NaN'| rows_fvc$FEV1_3 != 'NaN' |rows_fvc$FEV1_4 != 'NaN'), ])

rows_ok <- subset(rows_fvc,!(rows_fvc$FEV1_1 == 'NaN' & (rows_fvc$FEV1_2 != 'NaN'| rows_fvc$FEV1_3 != 'NaN' |rows_fvc$FEV1_4 != 'NaN')))
rows_na <- subset(rows_fvc,(rows_fvc$FEV1_1 == 'NaN' & (rows_fvc$FEV1_2 != 'NaN'| rows_fvc$FEV1_3 != 'NaN' |rows_fvc$FEV1_4 != 'NaN')))


rows_na$FEV1 <- ifelse(rows_na$FEV1_1 != "NaN",rows_na$FEV1_1,
                       ifelse(rows_na$FEV1_2 != "NaN",rows_na$FEV1_2,
                              ifelse(rows_na$FEV1_3 != "NaN",rows_na$FEV1_3,
                                     ifelse(rows_na$FEV1_4 != "NaN",rows_na$FEV1_4,NA))))
rows_ok$FEV1 <- ifelse(!is.na(rows_ok$FEV1_1),rows_ok$FEV1_1,NA)

rows_FEV1 <- rbind(rows_ok,rows_na)

rm(rows_fvc,rows_na,rows_ok)

#####FEV1/FVC######

nrow(rows_FEV1[!is.na(rows_FEV1$fvc) & !is.na(rows_FEV1$FEV1), ])

nrow(rows_FEV1[is.na(rows_FEV1$fvc) & is.na(rows_FEV1$FEV1), ]) ###都缺失的有44970

rows_FEV1$fev1fvc <- rows_FEV1$FEV1/rows_FEV1$fvc

save(rows_FEV1,file = "basic0828.RData")


names <- c("J12","J13","J14","J15","J16","J17","J18","J20","J21","J22",
           "J40","J41","J42","J43","J44","J45","J46","J47","J60","J61",
           "J62","J63","J64","J65","J66","J67","J68","J69","J70","J80",
           "J81","J82","J83","J84","J85","J86","J90","J91","J92","J93",
           "J94","J96","J98") 
library(tidyverse)

LUNG_disease <- rows_FEV1 %>%
  filter(if_any(everything(), ~ str_detect(.x, names))) 

library(dplyr)

no_lungdiease <- rows_FEV1 %>% anti_join(LUNG_disease) 人

save(no_lungdiease,file = "NO_LUNG_diease.RData")

colnames(rows_FEV1)


covariates_nodiease <- subset(no_lungdiease, select = c(f.eid,f.31.0.0,f.34.0.0,f.52.0.0,f.21022.0.0,
                                                        f.21000.0.0,f.21000.1.0,f.21000.2.0,f.21000.3.0,
                                                        f.20116.0.0,f.20116.1.0,f.20116.2.0,f.20116.3.0))
LF_nodiease <- subset(no_lungdiease, select = c(f.eid,fvc_1,fvc_2,fvc_3,fvc_4,fvc,FEV1_1,FEV1_2,FEV1_3,FEV1_4,
                                                FEV1,fev1fvc))

height <- read.table("D:/3_Gene/3_GWAS/lung_GWAS/ukb_data/base data/height/ukb672468.tab", header=TRUE, sep="\t")

covariates_a_nodiease <- left_join(covariates_nodiease,height, by = "f.eid")

cova_nodiease <- subset(covariates_a_nodiease,select = c(f.eid,f.31.0.0,f.34.0.0,f.52.0.0,f.50.0.0,
                                                         f.21022.0.0,f.21000.0.0,f.20116.0.0,f.21001.0.0))

colnames(cova_nodiease) <- c("f.eid","sex","birthyear","birthmonth","height_base",
                             "age_recruit","race_base","smoke_base","BMI_base")

save(cova_nodiease,file = "covariates_nodiease_all.RData")
save(LF_nodiease,file = "LF_nodisease.RData")

######validated analysis#####
install.packages("devtools")
#It may be necessary to install required as not all package dependencies are installed by devtools:
install.packages(c("dplyr","tidyr","ggplot2","MASS","meta","ggrepel","DT"))
devtools::install_github("PheWAS/PheWAS")
library(PheWAS)

#####
geno <- as.matrix(snp)
colnames(LF_ALL)
PCA <- read.table("D:/2_Gene/3_GWAS/lung_GWAS/ukb_data/base data/pca/ukb672468.tab", header=TRUE, sep="\t")
PCA10 <- PCA[,c(1,2,3,4,5,6,7,8,9,10,11)]
save(PCA10,file = "PCA10.RData")

cova_chinese <- subset(cova_nodiease,race_base == 5)
cova_chinese_adj <-subset(cova_chinese,select = c(f.eid,sex,height_base,age_recruit,smoke_base))
cov_cn_all <- left_join(cova_chinese_adj,PCA10,by = "f.eid")
LF <- subset(LF_nodiease,select = c(f.eid,FEV1,fvc,fev1fvc))
LF_CN <- LF %>% 
  filter(f.eid %in% cov_cn_all$f.eid) 
SNP_CN <- snp%>% 
  filter(FID %in% cov_cn_all$f.eid)

colnames(SNP_CN)[1] <- "f.eid"

table(cov_cn_all$smoke_base)
cov_cn_all <- cov_cn_all %>% mutate_at(c('smoke_base'), ~na_if(., -3))

library(PheWAS)
results_cn=phewas(phenotypes=LF_CN,genotypes=SNP_CN,covariates=cov_cn_all,cores=4)

results_cn_sig <- subset(results_cn,p<=0.05)
results_cn_sig$snp_check <- results_cn_sig$snp
results_cn_sig <- separate(results_cn_sig,col=snp_check,into=c("snp_c","alle"), sep = "_")####split columns into multiple columns

results_cn_sig$snp <- sub("_[ACGT]$", "", results_cn_sig$snp)

results_cn_fev1 <- subset(results_cn_sig,phenotype == "FEV1")
results_cn_fvc <- subset(results_cn_sig,phenotype == "fvc")
results_cn_fev1fvc <- subset(results_cn_sig,phenotype == "fev1fvc")

library(readxl)
FEV1_GWAS_sig <- read_excel("F:/3_Gene/3_GWAS/lung_GWAS/FEV1_GWAS_sig.xlsx")
FVC_GWAS_sig <- read_excel("F:/3_Gene/3_GWAS/lung_GWAS/FVC_GWAS_sig.xlsx")
FEV1FVC_GWAS_sig <- read_excel("F:/3_Gene/3_GWAS/lung_GWAS/FEV1FVC_GWAS_sig.xlsx")

filtered_results_cn_fev1 <- semi_join(results_cn_fev1, FEV1_GWAS_sig, by = c("snp" = "rs"))
filtered_results_cn_fvc <- semi_join(results_cn_fvc, FVC_GWAS_sig, by = c("snp" = "rs"))
filtered_results_cn_fev1fvc <- semi_join(results_cn_fev1fvc, FEV1FVC_GWAS_sig, by = c("snp" = "rs"))


library(xlsx)
write.xlsx(filtered_results_cn_fev1, file = "UKB_results_cn_fev1.xlsx",
           sheetName = "all", append = FALSE)  #######write xlsx. need library xlsx
write.xlsx(filtered_results_cn_fvc, file = "UKB_results_cn_fvc.xlsx",
           sheetName = "all", append = FALSE)  #######write xlsx. need library xlsx

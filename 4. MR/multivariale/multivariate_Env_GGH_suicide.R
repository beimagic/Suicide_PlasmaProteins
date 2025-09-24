rm(list=ls())
library('dplyr')
library("data.table")
library("TwoSampleMR")
library("stringr")
library('ieugwasr')
library('RadialMR')
library("readr")
library("vroom")
library("tidyr")
library("tibble")
#expo1
expo1 <- "bmi"
expo1_top_dat <-  as.data.frame(fread(paste0(path_env_top, "Expo_", expo1))) #1e-6/clumping
#expo2
expo2 <- "GGH"
expo2_top_dat <- as.data.frame(fread(paste0(path_protein_top, "Expo_", expo2))) #1e-6/clumping

# clean data for rbind
expo1_top_dat  <- expo1_top_dat[, colnames(expo1_top_dat) %in% colnames(expo2_top_dat)]
expo1_top_dat$chr.exposure <- as.integer(expo1_top_dat$chr.exposure)

#row bind top-dat--tips
tophits_list <- list(expo1_top_dat, expo2_top_dat)
tophits <- bind_rows(tophits_list) %>% pull(SNP)

#row bind gwas
expo1_gwas_dat <-  as.data.frame(fread(paste0(path_env_gwas, "result_", expo1, ".chrall.txt")))
expo1_gwas_dat <- expo1_gwas_dat[,c(1:12)]
expo1_gwas_dat$phenotype <- expo1
expo1_gwas <- format_data(
  dat = expo1_gwas_dat,
  type='outcome',
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col ="A1",
  other_allele_col = "REF",
  pval_col = "P",
  chr_col = "CHR",
  phenotype_col = "phenotype")
expo1_gwas$outcome= expo1
expo2_gwas <-  as.data.frame(fread(paste0(path_protein_gwas, "GWAS_form_", expo2 ,".txt")))
full_gwas_list <- list(expo1_gwas, expo2_gwas)

# combine expo
exposure_dat <- get_mv_exposures(tophits_list, full_gwas_list)

# outcome
t2d  = as.data.table(fread(paste0(path_base, "EUR_noUKB_QC_nodup.noamb_MR_0.8.txt"), header=T,colClasses=NULL))
t2d$phenotype <- 'suicide'
outcome_dat <- format_data(
  dat=t2d,
  type = "outcome",
  snps = exposure_dat$SNP,
  header = TRUE,
  phenotype_col = "phenotype",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  chr_col = "CHR",
  pos_col = "BP" )
print(dim(outcome_dat))

#Once the data has been obtained, harmonise so that exposure_dat and outcome_dat are on the same reference allele
mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)

#Finally, perform the multivariable MR analysis
res_bmis <- mv_multiple(mvdat)

# Creating a tidy outcome
result_2smr <- res_bmis$result %>%
  split_outcome() %>%
  separate(outcome, "outcome", sep="[(]") %>% 
  mutate(outcome=stringr::str_trim(outcome))%>% 
  generate_odds_ratios() %>% 
  tidy_pvals()
write.table(result_2smr, paste0(path, expo1,"_", expo2, "_suicide.csv"),quote=F,row.names=T,col.names=T,sep="\t")
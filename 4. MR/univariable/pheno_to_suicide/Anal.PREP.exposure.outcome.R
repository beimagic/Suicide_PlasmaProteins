rm(list=ls())
library('dplyr')
library("data.table")
library("TwoSampleMR")
#load protein list
ls.IM <- read.csv("***.csv")
# outcome suicide GWAS
t2d  = as.data.table(fread(paste0(path_base, "EUR_noUKB.txt"), header=T,colClasses=NULL))
for (f in c(1:nrow(ls.IM))){
  print(f)
  traitA= ls.IM$Assay[f]
  traitB='suicide'
  pheno_dat <- as.data.frame(fread(paste0(path_snp, traitA, ".csv"))) # p < 5e-8
  ## read exposure data
  pheno_dat <- format_data(
    dat = pheno_dat,
    type='exposure',
    snp_col = "rsid",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col ="ALLELE1",
    other_allele_col = "ALLELE0",
    pval_col = "P",
    chr_col = "CHROM")
  pheno_dat$exposure = ls.IM$Assay[f]
  ## local clump 
  snp_clump <- ld_clump_local(
    dplyr::tibble(rsid=pheno_dat$SNP,pval=pheno_dat$pval.exposure),
    clump_kb = 1000, clump_r2 = 0.01,clump_p = 1,
    bfile='/***/g1000_eur',
    plink_bin="/***/plink")
  exp_dat <- pheno_dat[which(pheno_dat$SNP %in% snp_clump$rsid),]
  print(dim(exp_dat))
  write.table(exp_dat,file=paste0(path_exposure, "Expo_", traitA), col.names=TRUE,row.names =FALSE,sep="\t",quot=F)
  
  ## read outcome data
  t2d_out <- format_data(
    dat=t2d,
    type = "outcome",
    snps = exp_dat$SNP,
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
  print(dim(t2d_out))
  write.table(t2d_out,file=paste0(path_outcome, traitA, "_to_suicide"),col.names=TRUE,row.names =FALSE,sep="\t",quot=F)
  rm(exp_dat)
  rm(t2d_out)
  rm(outcome_data)
}

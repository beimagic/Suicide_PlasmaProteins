rm(list=ls())
library('dplyr')
library("data.table")
library("TwoSampleMR")
#load protein list
ls.IM <- read.csv("/***/.csv")
rsid <- read.delim("/***/rsid")
#exposure datasets
sui_exp <- as.data.table(fread(paste0(path_exposure, "suicide.expo_dat_1e-6"), header=T, colClasses=NULL))
#Prepare protein GWAS data path
protein_gwas_list <-  data.table(list.files(path_input, full.names = TRUE))
protein_gwas_list <- protein_gwas_list[!grep(".tar", protein_gwas_list$V1), ]
dim(protein_gwas_list)
for (f in c(1:nrow(ls.IM))){
  print(f)
  traitA='suicide'
  traitB= ls.IM$Assay[f]
  gwas.name=paste0(traitB, "_", ls.IM$UniProt[f])
  ind = grep(gwas.name, protein_gwas_list$V1)
  ## combine all chr for protein GWAS
  chr_list <- list.files(protein_gwas_list$V1[ind], full.names = T)
  chr_list <- chr_list[-23]
  #read all outcome datasets
  for (j in 1:22){
    print(paste0('outcome_chr',j))
    tmp <- fread(chr_list[j])
    tmp <- data.frame(tmp)
    if(j==1){
      outcome_data <- tmp
    }else{
      outcome_data <- rbind(outcome_data,tmp)
    }
  outcome_data <- merge(rsid,outcome_data,by='ID') 
  outcome_data$P <- 10^(-outcome_data$LOG10P)
  # 1e-6
  outcome_dat <- format_data(
    dat = outcome_data,
    type = "outcome",
    snps = sui_exp$SNP,
    header = TRUE,
    snp_col = "rsid",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col ="ALLELE1",
    other_allele_col = "ALLELE0",
    pval_col = "P",
    chr_col = "CHROM",
    phenotype_col = "phenotype")
  outcome_dat$outcome=traitB
  print(dim(outcome_dat))
  write.table(outcome_dat,file=paste0(path_outcome,"outcome_", traitB),col.names=T,row.names = F,sep="\t",quot=F)
  rm(outcome_dat)
}
}




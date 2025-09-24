# This script was used to perform protein-imaging association analysis: DWI
imaging_data <- as.data.frame(fread(paste0(path_wd,"DWI_imaging_all.csv")))
all_pheno_data <- as.data.frame(fread(paste0(path_wd, "**.csv")))
all_pheno_data <- merge(all_pheno_data, imaging_data, by="eid")
all_pheno_data$site <- factor(all_pheno_data$site)
all_pheno_data$protein_batch <- factor(all_pheno_data$protein_batch)
pheno_list <- as.data.frame(fread(paste0(path_wd, "sig_proteins.csv")))
lable <- as.data.frame(fread(paste0(path_wd,"DWI_Imaging.csv")))
TimeGap  <- as.data.frame(fread(paste0(path_cov,"TimeGap.csv")))
for (k in (1:nrow(lable))){
  ind <- grep(lable$FieldID[k], colnames(all_pheno_data))
  img_ID <- colnames(all_pheno_data)[ind]
  for (i in (1:nrow(pheno_list))){
    protein_ID <- pheno_list$Assay[i]
    data <- all_pheno_data[,c("eid", protein_ID, img_ID, cov_list_final)]
    colnames(data)[2:3] <-  c("protein", "img")
    time_gap <- TimeGap[,c("eid", pheno_list$Panel[i])]
    data <- merge(data, time_gap, by="eid")
    colnames(data)[ncol(data)] <-  "tg"
    data$age <- scale(data$age)
    data$TDI <- scale(data$TDI)
    data$BMI <- scale(data$BMI)
    data$protein <- scale(data$protein)
    data$img <- scale(data$img)
    data = data[complete.cases(data),]
    formula_str <-  paste("protein ~ img + ", paste(cov_list_final, collapse = " + "))
    formula <- as.formula(formula_str)
    pheno_lm <- lm(formula,data = data)
  }
}
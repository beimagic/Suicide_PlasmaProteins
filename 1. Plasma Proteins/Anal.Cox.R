# This script was used to perform cox analysis
all_pheno_data <- as.data.frame(fread(paste0(path_wd,"**.csv")))
TimeGap  <- as.data.frame(fread(paste0(path_cov,"TimeGap.csv")))
for (i in(1:nrow(sum_data))){
   IDs <- sum_data$Assay[i] 
   data <- all_pheno_data[,c("eid", "group", "diff",  IDs, "sex",  "age", "Ethnic", "TDI", "BMI", "smoking", "site",  "protein_batch", cov_list_final)]
   colnames(data)[4] <- "value"
   time_gap <- TimeGap[,c("eid", sum_data$Panel[i])]
   data <- merge(time_gap, data, by="eid")
   colnames(data)[2] <-  "tg"
   data = data[complete.cases(data),]
   formula_str <- paste("Surv(diff, group) ~ value + sex + age + Ethnic + TDI + BMI + tg + smoking + strata(site) + strata(protein_batch) +", 
                       paste(cov_list_final, collapse = " + "))
  formula <- as.formula(formula_str)
  pheno_cox <- coxph(formula,  id = eid, data = data)
  sum_data[i, 7:11] <- summary(pheno_cox)$coef[1,] 
}
write.table(sum_data, paste0(path_output, "proteins_cox.csv"), col.names = T, sep=',', row.names = F)

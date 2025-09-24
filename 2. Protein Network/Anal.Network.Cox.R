# This script was used to perform cox analysis
all_pheno_data <- as.data.frame(fread(paste0(path_wd,"*")))
net_data <- as.data.frame(fread(paste0(path_input, "net.value.txt")))
# merge 
all_pheno_data <- merge(all_pheno_data, net_data, by="eid")
for (i in(1:nrow(sum_data))){
   IDs <- sum_data$Modules[i] 
   data <- all_pheno_data[,c("eid", "group", "diff", IDs, cov_list_select)]
   colnames(data)[4] <- "value"
   data$age <- scale(data$age)
   data$TDI <- scale(data$TDI)
   data$BMI <- scale(data$BMI)
   data$value <- scale(data$value)
   data = data[complete.cases(data),]
  formula_str <- paste("Surv(diff, group) ~ value + sex + age + Ethnic + TDI + BMI + smoking + strata(site) +", 
                       paste(cov_list_select_2, collapse = " + "))
  formula <- as.formula(formula_str)
  pheno_cox <- coxph(formula, data = data)
  sum_data[i, 7:11] <- summary(pheno_cox)$coef[1,] 
}
write.table(sum_data, paste0(path_output, "proteins_modules_Cox.csv"), col.names = T, sep=',', row.names = F)

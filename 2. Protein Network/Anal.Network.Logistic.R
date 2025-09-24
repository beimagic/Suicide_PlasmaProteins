# This script was used to perform logistic-regression analysis
all_pheno_data <- as.data.frame(fread(paste0(path_wd,"**.csv")))
net_data <- as.data.frame(fread(paste0(path_input, "net.value.txt")))
# merge 
all_pheno_data <- merge(all_pheno_data, net_data, by="eid")
all_pheno_data$site <- factor(all_pheno_data$site)
all_pheno_data$protein_batch <- factor(all_pheno_data$protein_batch)
for (i in(1:nrow(sum_data))){
  IDs <- sum_data$Modules[i] 
  data <- all_pheno_data[,c("eid", "group", IDs, cov_list_select)]
  colnames(data)[3] <- "value"
  data$age <- scale(data$age)
  data$TDI <- scale(data$TDI)
  data$BMI <- scale(data$BMI)
  data$value <- scale(data$value)
  data = data[complete.cases(data),]
  formula_str <-  paste("group ~ value + ", paste(cov_list_select, collapse = " + "))
  formula <- as.formula(formula_str)
  pheno_logistic <- glm(formula, data = data, family = binomial(), control = glm.control(maxit = 200))
  sum_data[i, 7:10] <- summary(pheno_logistic)$coef[2,]
}
write.table(sum_data, paste0(path_output, "proteins_modules_logistic_regression.csv"), col.names = T, sep=',', row.names = F)
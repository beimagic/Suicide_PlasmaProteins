# This script was used to perform Logistic.regression for plasma proteins
# Prepare datasets 
all_pheno_data <- as.data.frame(fread(paste0(path_wd,"*")))
all_pheno_data$site <- factor(all_pheno_data$site)
all_pheno_data$protein_batch <- factor(all_pheno_data$protein_batch)
TimeGap  <- as.data.frame(fread(paste0(path_cov,"TimeGap.csv")))
# calculation
 for (i in(1:nrow(sum_data))){
  IDs <- sum_data$Assay[i] 
  data <- all_pheno_data[,c("eid", "group", IDs, cov_list_final)]
  colnames(data)[3] <- "value"
  time_gap <- TimeGap[,c("eid", sum_data$Panel[i])]
  data <- merge(time_gap, data, by="eid")
  colnames(data)[2] <-  "tg"
  data = data[complete.cases(data),]
  formula_str <-  paste("group ~ value + tg + ", paste(cov_list_final, collapse = " + "))
  formula <- as.formula(formula_str)
  pheno_logistic <- glm(formula, data = data, family = binomial(), control = glm.control(maxit = 200))

  #save statistics
  sum_data[i, 7:10] <- summary(pheno_logistic)$coef[2,]
  rm(pheno_logistic)
  rm(data)
}
colnames(sum_data)[7:10] <- c("beta", "se", "z","p")
write.table(sum_data, paste0(path_output, "proteins_logistic.csv"), col.names = T, sep=',', row.names = F)


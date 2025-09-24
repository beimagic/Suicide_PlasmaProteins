# This script was used to perform network clustering analysis
library(data.table)
library(netboost)
library(readxl)
pheno_data <- as.data.frame(fread(paste0(path_protein_data, "**.csv")))
#start network processing
results <- netboost(datan= pheno_data, 
                    stepno=20L,
                    filter_method = 'spearman',
                    progress = 1000L,
                    cores=10L,
                    soft_power=2, 
                    min_cluster_size=20L, 
                    n_pc=1, 
                    robust_PCs = TRUE,
                    scale=TRUE, 
                    ME_diss_thres=0.25, 
                    qc_plot=TRUE, 
                    method = 'spearman', 
                    verbose=3)
save(results,file = 'netboostReslult_all.RData')

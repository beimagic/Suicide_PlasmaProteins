get_mv_exposures <- function(tophits_list, full_gwas_list) {
  
  ### This is a modified version of `mv_extract_exposures` function in TwoSampleMR package.
  
  # Collapse list of exposures' tophits into a dataframe
  exposures <- bind_rows(tophits_list)
  
  # clump exposures: this will produce a list of instruments (shared and unique) of the given exposures
  temp <- exposures
  temp$id.exposure <- 1
 # temp <- clump_data(temp)
  snp_clump <- ld_clump_local(
    dplyr::tibble(rsid=temp$SNP,pval=temp$pval.exposure),
    clump_kb = 1000, clump_r2 = 0.01,clump_p = 1,
    bfile='/public/home2/nodecw_group/lizeyu/1210/ref_GRCh37/g1000_eur',
    plink_bin="/public/home/beizhang/Suicide/PRS/PRS_1st_R/plink1.9/plink")
  temp2 <- temp[which(temp$SNP %in% snp_clump$rsid),]
  exposures <- filter(exposures, SNP %in% temp2$SNP)
  
  # subset full gwas summary stats of each exposure to the list of SNPs (instruments) produced above
  for (i in 1:length(full_gwas_list)){
    full_gwas_list[[i]] <- full_gwas_list[[i]] %>% filter(SNP %in% exposures$SNP)
  }
  
  # Collapse lists of subset gwas into a dataframe
  d1 <- bind_rows(full_gwas_list) %>%
    distinct()
  
  ###  The logic of next steps is largely unchanged from the original function `mv_extract_exposures`
  
  # get auto-generated ids
  id_exposure <- unique(d1$id.outcome) 
  
  # convert first trait to exposure format  -- exp1 is exposure
  tmp_exposure <- d1 %>% filter(id.outcome == id_exposure[1]) %>% convert_outcome_to_exposure()
  # keep other traits (n>=2) as outcome -- exp2+ are outcomes
  tmp_outcome <- d1 %>% filter(id.outcome != id_exposure[1])
  
  # Harmonise against the first trait
  d <- harmonise_data(exposure_dat = tmp_exposure, 
                      outcome_dat = tmp_outcome, action=2)
  
  # Only keep SNPs that are present in all
  snps_not_in_all <- d %>% 
    count(SNP)  %>% 
    filter(n < length(tophits_list)-1) %>%
    pull(SNP)
  d <- filter(d, !SNP %in% snps_not_in_all)
  
  # Subset and concat data
  
  # for exp1 get exposure cols
  # dh1x <- d %>% filter(id.outcome == id.outcome[1]) %>% 
  #   select("SNP", contains("exposure"))
  dh1x <- d %>% filter(id.outcome == id.outcome[1])
  dh1x <-  dh1x[  ,c("SNP", grep("exposure", colnames(dh1x), value = TRUE))]
  
  # for exp2 get outcome cols
  # dh2x <-d %>%  select(SNP, contains("outcome"))
  dh2x <- d[  ,c("SNP", grep("outcome", colnames(d), value = TRUE))]
  
  # rename outcome to exposure in these
  names(dh2x) <- gsub("outcome", "exposure", names(dh2x) )
  if ("chr.exposure" %in% colnames(dh1x)){
  dh1x$chr.exposure <- as.integer(dh1x$chr.exposure) }
  if ("pos.exposure" %in% colnames(dh1x)){
    dh1x$pos.exposure <- as.integer(dh1x$pos.exposure)}
  # join together (drop not needed cols)
  # exposure_dat <- bind_rows(dh1x, dh2x) %>%  
  #   select(-c("samplesize.exposure" ,"mr_keep.exposure", "pval_origin.exposure")) %>% 
  #   distinct()
  exposure_dat <- bind_rows(dh1x, dh2x)
  # exposure_dat <- exposure_dat[ ,-c(15, 12, 13)]  %>% distinct()
  exposure_dat <- exposure_dat[ ,-which(colnames(exposure_dat) %in% c("samplesize.exposure", "mr_keep.exposure","pval_origin.exposure"))] %>% distinct()
  return(exposure_dat)
}  
  
  tidy_pvals<-function(df){
    # round up output values and keep p-vals in scientific notation
    df %>% 
      mutate(pval= as.character(pval)) %>% 
      mutate_if(is.numeric, round, digits=2) %>% 
      mutate(pval=as.numeric(pval),
             pval=scales::scientific(pval, digits = 2),
             pval=as.numeric(pval))
  }
  
  tidy_mvmr_output <- function(mvmr_res) {
    #  tidy up MVMR returned output
    mvmr_res %>%
      as.data.frame() %>% 
      rownames_to_column("exposure") %>% 
      rename(b=Estimate,
             se="Std. Error",
             pval="Pr(>|t|)") %>% 
      select(-c(`t value`)) %>% 
      TwoSampleMR::generate_odds_ratios()
  }

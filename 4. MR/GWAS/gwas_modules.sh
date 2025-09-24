for chr in {1..22};
do

/home1/beizhang/GWAS/plink/plink2 --bfile /***/ukb_imp_chr${chr} --linear firth-fallback hide-covar --variance-standardize --pheno /***/modules_value.txt --input-missing-phenotype -9 --covar /***/cov_gwas_module.txt --out /***/result_module.chr$chr

done


for j in {4}
do

for ((i=1;i<=22;i++)) 
do 
echo result_module.chr$i.Module$j.glm.linear;

done|xargs -i cat {} >>result_module.Module$j.chrall.txt

sed -i '/CHROM/d' result_module.Module$j.chrall.txt 

sed -i 'CHR	POS	SNP	REF	ALT	A1	A1_FREQ	OBS_CT	BETA	SE	T_STAT	P' result_module.Module$j.chrall.txt

done



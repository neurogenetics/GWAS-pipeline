# README

This is a script to quickly check for ancestry in genotype data based on hapmap-data

### also needs these.... but too big to upload to github
HAPMAP_hg19_new.log
HAPMAP_hg19_new.bed
HAPMAP_hg19_new.fam
HAPMAP_hg19_new.bim

## No ancestry outliers -> based on Hapmap3 PCA plot, should be near combined CEU/TSI

```
# note if you have NeuroX data then you need to convert snp-id's first since this only works with rs-ids
# neuroX_snps_for_hapmap.txt # extract these snps from neuroX
# neuroX_snps_for_hapmap_conversion.txt # conversion of IDs --update-map
# Keep in mind that this comparison with hapmap is based on the number of SNPs that overlap between your input dataset and hapmap

plink --bfile $FILENAME --bmerge HAPMAP_hg19_new --out hapmap3_bin_snplis --make-bed
plink --bfile $FILENAME --flip hapmap3_bin_snplis-merge.missnp --make-bed --out $FILENAME3
plink --bfile $FILENAME3 --bmerge HAPMAP_hg19_new --out hapmap3_bin_snplis --make-bed
plink --bfile $FILENAME3 --exclude hapmap3_bin_snplis-merge.missnp --out $FILENAME4 --make-bed
plink --bfile $FILENAME4 --bmerge HAPMAP_hg19_new --out hapmap3_bin_snplis --make-bed
plink --bfile hapmap3_bin_snplis --geno 0.01 --out pca --make-bed --pca 4

# then add some names here and there
grep "EUROPE" pca.eigenvec > eur.txt
grep "ASIA" pca.eigenvec > asia.txt
grep "AFRICA" pca.eigenvec > afri.txt
grep -v -f eur.txt pca.eigenvec | grep -v -f asia.txt | grep -v -f afri.txt > new_samples.txt
cut -d " " -f 3 after_gender.fam > new_samples_add.txt
paste new_samples_add.txt new_samples.txt > new_samples2.txt
paste eur_add.txt eur.txt > euro.txt
paste asia_add.txt asia.txt > asiao.txt
paste afri_add.txt afri.txt > afrio.txt
cat new_samples2.txt euro.txt asiao.txt afrio.txt > pca.eigenvec2

# R script for PCA plotting and filtering
R < PCA_in_R.R --no-save  

# then back to plink to remove outliers
plink --bfile $FILENAME --keep PCA_filtered_europeans.txt --make-bed --out after_gender_heterozyg_hapmap
cat PCA_filtered_asians.txt PCA_filtered_africans.txt PCA_filtered_mixed_race.txt > hapmap_outliers.txt
# this creates several plots and lists based on genetic ancestry
```

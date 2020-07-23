# GWAS-pipeline

Date: June 2019

#### NOTE!! see for a more updated version here => https://github.com/GP2code/GWAS

Authors: Cornelis Blauwendraat, Mike Nalls, Hirotaka Iwaki, Sara Bandres-Ciga, Mary Makarious, Ruth Chia, Frank Grenn, Hampton Leonard, Monica Diez-Fairen, Jeff Kim 

## General description and purpose
This is a general (somewhat comprehensive) description of the LNG GWAS pipeline which can be used to guide researchers on how-to run a GWAS.
Roughly the pipeline can be divided into five steps:

1. QC and data cleaning
2. Imputation
3. GWAS
4. Optional meta-analyses
5. Results visualization 


# 1. QC and data cleaning
The QC and data cleaning is very important prior imputation since the cleaner data you put in there the better and more accurate data you will get out. PLINK and GCTA are very easy program for doing this.

## Sample QC parameters


- Heterozygosity outliers (--het), F cut-off of -0.15 and <- 0.15 for inclusion 

High or low heterozygosity is an indication of failed experiments or contamination of DNA sample

```
plink --bfile $FILENAME --geno 0.01 --maf 0.05 --indep-pairwise 50 5 0.5 --out pruning
plink --bfile $FILENAME --extract pruning.prune.in --make-bed --out pruned_data
plink --bfile pruned_data --het --out prunedHet
awk '{if ($6 <= -0.15) print $0 }' prunedHet.het > outliers1.txt
awk '{if ($6 >= 0.15) print $0 }' prunedHet.het > outliers2.txt
cat outliers1.txt outliers2.txt > HETEROZYGOSITY_OUTLIERS.txt
cut -f 1,2 HETEROZYGOSITY_OUTLIERS.txt > all_outliers.txt
plink --bfile $FILENAME --remove all_outliers.txt --make-bed --out $FILENAME_after_heterozyg
```

- Call rate outliers (--mind), Call rate of >95% is preferred which is --mind 0.05

Low call rate is an indication of a failed experiment

```
plink --bfile $FILENAME --mind 0.05 --make-bed --out $FILENAME_after_call_rate
mv $FILENAME_after_call_rate.irem CALL_RATE_OUTLIERS.txt
```

- Genetic sex fails (--check-sex), F cut-off of 0.25 and 0.75 for inclusion --check-sex 0.25 0.75

Genetic sex fails are indication of a failed experiment or a sample-switch

Note use 2,781,479 and 155,701,383 for hg38

```
plink --bfile $FILENAME --check-sex 0.25 0.75 --maf 0.05 --out gender_check1
plink --bfile $FILENAME --chr 23 --from-bp 2699520 --to-bp 154931043 --maf 0.05 --geno 0.05 --hwe 1E-5 --check-sex  0.25 0.75 --out gender_check2 
grep "PROBLEM" gender_check1.sexcheck > problems1.txt
grep "PROBLEM" gender_check2.sexcheck > problems2.txt
cat problems1.txt problems2.txt > GENDER_FAILURES.txt
cut -f 1,2 GENDER_FAILURES.txt > samples_to_remove.txt
plink --bfile $FILENAME --remove samples_to_remove.txt --make-bed --out $FILENAME_after_gender
```

## optional steps depending general purpose of data cleaning
- No ancestry outliers -> based on Hapmap3 PCA data, should be near combined CEU/TSI, 6SD+/-

```
see script hapmap ancestry check folder this
```

- No relatedness (--grm-cutoff), Typically this would be set at 0.125 to remove cousins or more related individuals

Note you can also increase the grm-cut-off to e.g. 0.8 and remove duplicate samples. This would allow you keep more related individuals and use a linear mixed model.

```
gcta --bfile $FILENAME --make-grm --out GRM_matrix --autosome --maf 0.05 
gcta --grm-cutoff 0.125 --grm GRM_matrix --out GRM_matrix_0125 --make-grm
plink --bfile $FILENAME --keep GRM_matrix_0125.grm.id --make-bed --out $FILENAME_relatedness
```

## Variant QC parameters

- Variant missingness, basically if a variant has high missingness this suggests that the genotyping probe for this variant is not great 

```
plink --bfile $FILENAME --make-bed --out $FILENAME_geno --geno 0.05
```


- Missingness by case control (--test-missing), using P > 1E-4

```
plink --bfile $FILENAME --test-missing --out missing_snps 
awk '{if ($5 <= 0.0001) print $2 }' missing_snps.missing > missing_snps_1E4.txt
plink --bfile $FILENAME --exclude missing_snps_1E4.txt --make-bed --out $FILENAME_missing1
```

- Missing by haplotype (--test-mishap), using P > 1E-4
```
plink --bfile $FILENAME --test-mishap --out missing_hap 
awk '{if ($8 <= 0.0001) print $9 }' missing_hap.missing.hap > missing_haps_1E4.txt
sed 's/|/\n/g' missing_haps_1E4.txt > missing_haps_1E4_final.txt
plink --bfile $FILENAME --exclude missing_haps_1E4_final.txt --make-bed --out $FILENAME_missing2
```


- Hardy Weinberg SNP from controls only (--filter-controls --hwe 1E-4)
```
plink --bfile $FILENAME --filter-controls --hwe 1E-4 --write-snplist
plink --bfile $FILENAME --extract plink.snplist --make-bed --out $FILENAME_HWE
```

# optional steps depending general purpose of data cleaning 
- Minor allele frequency (--maf), In some instances this might not be used if specific rare variants are on your array that you do not want to exclude.

```
plink --bfile $FILENAME --maf 0.01 --make-bed --out $FILENAME_MAF
```

## Preparation of data prior to submission to Imputation server

Prepare data for HRC imputation using Will Rayner tool and continue with the plink file that is cleaned with the above described options.


```
# download file to check
wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim.v4.2.5.zip

# make your .frq file
plink --bfile $FILENAME --freq --out $FILENAME

perl HRC-1000G-check-bim.pl -b $FILENAME.bim -f $FILENAME.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h

# then run to fix your data
sh Run-plink.sh

# then make vcf files

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
  plink --bfile YOURFILE-updated-chr$chnum --recode vcf --chr $chnum --out YOURFILE$chnum 
done

## then sort and zip

module load vcftools

for chnum in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
  do
	vcf-sort YOURFILE$chnum.vcf | bgzip -c >  pre_impute_YOURFILE_$chnum.vcf.gz
done

# and then you are ready to submit to the imputation server
```


# 2. Imputation
Imputation will be done using the Michigan Imputation Server (BIG THANKS TO THEM FOR MAKING OUR LIVES SO MUCH EASIER!!)

https://imputationserver.sph.umich.edu/index.html

Using Reference panel: HRC 1.1 2016 and Eagle v2.3 Phasing

### More info on reference panel here :

https://imputationserver.readthedocs.io/en/latest/reference-panels/ and http://www.haplotype-reference-consortium.org/

### Output of Imputation server

Three files per input chromosomes e.g. for chromosome 21:

```
chr21.dose.vcf.gz.tbi
chr21.dose.vcf.gz
chr21.info.gz
```
The .dose.vcf.gz are your imputed genotypes with dosage information
The .dose.vcf.gz.tbi is the index file of your .vcf.gz file (which some programs require)
The .info.gz file is a file that provided some information on the included variants such as quality and frequency

Example of the .info.gz file (Note we copied and pasted some random variants here)
```
SNP	REF(0)	ALT(1)	ALT_Frq	MAF	AvgCall	Rsq	Genotyped	LooRsq	EmpR	EmpRsq	Dose0	Dose1
9:11739	G	A	0.00004	0.00004	0.99996	0.00023	Imputed	-	-	-	-	-
9:14665	G	A	0.20555	0.20555	0.82608	0.20321	Imputed	-	-	-	-	-
9:62179	G	A	0.22655	0.22655	0.99994	0.99976	Genotyped	0.998	0.999	0.99717	0.99920	0.00078
```
```
Most important columns are:
SNP = variant name and in this case Chromosome:Basepair
REF(0) = Reference allele
ALT(1) = Alternative allele
ALT_Frq = Alternative allele frequency
MAF = Minor Allele frequency
Rsq = This is the imputation quality, depending on who you ask, >0.3, >0.6 or >0.8 is good enough for GWAS.
Genotyped = This says either Imputed or Genotyped. Genotyped means that it was included in the input data.
```

### References:

Imputation server: https://www.ncbi.nlm.nih.gov/pubmed/27571263

HRC panel: https://www.ncbi.nlm.nih.gov/pubmed/27548312

Eagle Imputation: https://www.ncbi.nlm.nih.gov/pubmed/27694958


# 3. GWAS
For running GWAS many many tools and programs are available. We most commonly use either RVTESTS or PLINK. 

## Covariate generation
Covariates should always be included when running GWAS, we typically include at least AGE, SEX and PC1-5. AGE and SEX you (hopefully) already have. PC1-5 you can generate using several programs we typically use either PLINK or FlashPCA for this.

PLINK
more info here -> https://www.cog-genomics.org/plink2
```
# Filter data (you can speed up things by adding --memory 119500 --threads 19 in PLINK)
plink --bfile FILENAME --maf 0.01 --geno 0.01 --hwe 5e-6 --autosome --exclude exclusion_regions_hg19.txt 
--make-bed --out FILENAME_2  
# Prune snps 
plink --bfile FILENAME_2 --indep-pairwise 1000 10 0.02 --autosome --out pruned_data
# Extract pruned SNPs and only these variants will be used for PC calculation
plink --bfile FILENAME_2 --extract pruned_data.prune.in --make-bed --out FILENAME_3 
# Calculate/generate PCs based on pruned data set
plink --bfile FILENAME_3 --pca --out PCA
# then use the .eigenvec file
```

FlashPCA (you can speed up things by adding --numthreads 19 in FlashPCA)
more info here -> https://github.com/gabraham/flashpca
```
# Filter data
plink --bfile FILENAME --maf 0.01 --geno 0.01 --hwe 5e-6 --autosome --exclude exclusion_regions_hg19.txt 
--make-bed --out FILENAME_2  
# Prune snps 
plink --bfile FILENAME_2 --indep-pairwise 1000 10 0.02 --autosome --out pruned_data
# Extract pruned SNPs and only these variants will be used for PC calculation
plink --bfile FILENAME_2 --extract pruned_data.prune.in --make-bed --out FILENAME_3 
# Calculate/generate PCs based on pruned data set
flashpca --bfile FILENAME_3 --suffix _filter_pruned_forPCA.txt --numthreads 19
# then use the pcs_* file
```

Note that it is recommend (by FlashPCA authors and others) to exclude some regions prior creating PC's

hg19:
```
5 44000000 51500000 r1
6 25000000 33500000 r2
8 8000000 12000000 r3
11 45000000 57000000 r4
```
hg38:
```
1   47534328    51534328    r1
2   133742429   137242430   r2
2   182135273   189135274   r3
3   47458510    49962567    r4
3   83450849    86950850    r5
5   98664296    101164296   r6
5   129664307   132664308   r7
5   136164311   139164311   r8
6   24999772    35032223    r9
6   139678863   142178863   r10
8   7142478 13142491    r11
8   110987771   113987771   r12
11  87789108    90766832    r13
12  109062195   111562196   r14
20  33412194    35912078    r15
```

## Option1: RVTESTS
Why RVTESTS? It is very easy to use due to the similar file structure as PLINK takes, it has a good manual and it has very flexible options to use.

Files needed to run GWAS:
1. Phenotype file
2. Regions file
3. Imputed genotype data

### Phenotype file (and also covariate file)
Structure is very similar to PLINK. We generally add a couple of columns such as: SEX_cov,pheno_01,AGE,dataset,PC1-PC10

- SEX_cov is sex-1 this is because some programs prefer binairy coding of covariates over 1 and 2
- pheno_01 is pheno-1 this is again because some programs prefer binairy coding of phenotypes over 1 and 2
- AGE is for some phenotypes very important so would always recommend including this
- dataset is the data origin and can be an important covariate in some instances
- PC1-PC10 are the principal components created as above described.

```
FID	IID	patid	matid	sex	pheno	SEX_cov	pheno_01	AGE	dataset	PC1	PC2	PC3	PC4	PC5	PC6	PC7	PC8	PC9	PC10
sample1	sample1	0	0	1	2	0	1	56	dataset1	0.0128029	0.00134518	-0.0132907	0.00909855	-0.00238954	0.0327338	-0.00984243	-0.0295871	-0.00869277	0.0112763
sample2	sample2	0	0	2	1	1	0	48	dataset1	0.0149187	0.00623717	0.00486343	-0.0384283	-0.00951609	-0.0304118	0.0293294	0.00264849	-0.0339268	-0.00915441
sample3	sample3	0	0	2	1	1	0	78	dataset2	0.00596214	0.0154845	0.00780206	0.0832667	0.0516746	-0.0313716	0.00238327	-0.0361836	0.0244427	0.012373
```


### Regions file
This makes sure that you don't include low (imputation) quality and low frequency variants.
The code below creates two files, which are useful for when running GWAS and reformating the GWAS output

- maf001rsq03minimums_chr22.info, this is file looks like this:
```
SNP    ALT_Frq    Rsq
22:17029525    0.9407    0.5927
22:17030029    0.94073    0.59404
22:17030126    0.94039    0.59182
22:17030792    0.94112    0.60053
```

- maf001rsq03minimums_chr22.txt, this is file looks like this:
```
22:17029525-17029525
22:17030029-17030029
22:17030126-17030126
22:17030792-17030792
```

Code to run in R:

```
R
library(plyr)
# note1 you can speed this up with the data.table package and changing read.table and write.table to fread and fwrite
# note2 if you already unzipped your info.gz files update this to .info
# note3 you can change the MAF and Rsq filter to whatever you think is appropriate 

for(i in 1:22)
{
  input <- paste("chr",i,".info.gz", sep = "")
  data <- read.table(input, header = T)
  dat <- subset(data, MAF >= 0.001 & Rsq >= 0.30)
  dat$chr <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[1]]
  dat$bp <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[2]]
  da <- dat[,c("SNP","ALT_Frq","Rsq")]
  write.table(da, paste("maf001rsq03minimums_chr",i,".info",sep = ""), row.names = F, quote = F, sep = "\t")
}

for(i in 1:22)
{
  input <- paste("chr",i,".info.gz", sep = "")
  data <- read.table(input, header = T)
  dat <- subset(data, MAF >= 0.001 & Rsq >= 0.30)
  dat$chr <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[1]]
  dat$bp <- ldply(strsplit(as.character(dat$SNP), split = ":"))[[2]]
  dat$range <- paste(dat$chr, ":", dat$bp, "-", dat$bp, sep = "")
  da <- dat[,c("range")]
  write.table(da, paste("maf001rsq03minimums_chr",i,".txt",sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
}

```

### Imputed genotype data

This is the .vcf file that is the output from the imputation server.

They should be in this format chr$CHNUM.dose.vcf.gz

### Running the GWAS

Can start like this:
```
sbatch --cpus-per-task=10 --mem=48g --time=8:00:00 GWAS_per_chr.sh PHENO CHNUM KEEPFILE 
```
or simpler if not having a slurm submission system
```
sh GWAS_per_chr.sh PHENO CHNUM KEEPFILE 
```

```
#!/bin/bash
PHENO=$1
CHNUM=$2
KEEPFILE=$3

cd /the/location/of/your/files/

rvtest --noweb --hide-covar --rangeFile maf001rsq03minimums_chr$CHNUM.txt \
--out $PHENO.$KEEPFILE.chr$CHNUM --single wald \
--inVcf chr$CHNUM.dose.vcf.gz --dosage DS --pheno $PHENO.pheno.txt \
--pheno-name pheno --covar $PHENO.pheno.txt \
--covar-name SEX_COV,PC1,PC2,PC3,PC4,PC5 --peopleIncludeFile $KEEPFILE.txt
```

## option explanation
```
--noweb  == to ignore web-check
--hide-covar  == hide covariates 
--rangeFile maf001rsq03minimums_chr$CHNUM.txt  == the file made at STEP2
--out $PHENO.$KEEPFILE.chr$CHNUM  == the location and name of your output file(s)
--single wald  == type of test to use see for more tests-> http://zhanxw.github.io/rvtests/
--inVcf chr$CHNUM.dose.vcf.gz == VCF file as input
--dosage DS  == tells the program to use dosage format from your .vcf files
--pheno $PHENO.pheno.txt  == link to the phenotype file from STEP1
--pheno-name pheno   == name of the phenotype from the phenotype file from STEP1
--covar $PHENO.pheno.txt  == link to the phenotype file from STEP1
--covar-name SEX_COV,PC1,PC2,PC3,PC4,PC5  == name(s) of the covariates from the phenotype file from STEP1
--peopleIncludeFile $KEEPFILE.txt  == keep file from STEP1
```

## Option2: PLINK
Why PLINK? It is very easy to use, it has a good manual and it has very flexible options to use.
PLINK2 is often used, but typically that means version 1.9. The actual PLINK2 version is still under development and the first releases look very promising (https://www.cog-genomics.org/plink/2.0/) incredibly fast.


### Phenotype file
Structure is very similar to PLINK:
``` 
work in progress
```


### Imputed genotype data

This is the .vcf file that is the output from the imputation server
``` 
work in progress
```

### PLINK 1.9
Example 1
``` 
work in progress
```

### PLINK 2
Example 1
``` 
work in progress
```

## Option3: Linear mixed model
Why linear mixed model?


Example 1
``` 
work in progress
```

### References:

FlashPCA: https://www.ncbi.nlm.nih.gov/pubmed/28475694

RVTESTS: https://www.ncbi.nlm.nih.gov/pubmed/27153000

PLINK: https://www.ncbi.nlm.nih.gov/pubmed/25722852


# 4. Optional meta-analyses
Meta-analyses are used when you have multiple datasets.

## Prep of your GWAS results (based on RVTEST file)

merge all GWAS files per chromosome and merge all .info files
```
cat *.SingleWald.assoc | grep -v 'N_INFORMATIVE' > allChrs_FILE.assoc 
cat maf001rsq03minimums_chr*.info | grep -v 'Rsq' > allChrs_FILE.Info
```
Then in reformat in R
```
R
# note1 you can speed this up with the data.table package and changing read.table and write.table to fread and fwrite
# note2 we set here beta filtering at <5 and >-5 since typically these are unrealistic beta's coming from GWAS
infos <- read.table(paste("allChrs_FILE.Info"))
colnames(infos) <- c("SNP","ALT_Frq","Rsq")
assoc <- read.table(paste("allChrs_FILE.assoc"))
colnames(assoc) <- c("CHROM","POS","REF","ALT","N_INFORMATIVE","Test","Beta","SE","Pvalue")
data <- merge(infos, assoc, by.x = "SNP", by.y = "Test", all.y = T)
dat <- subset(data, Beta < 5 & Beta > -5 & !is.na(data$Pvalue))
dat$chr <- paste("chr",dat$CHROM, sep = "")
dat$markerID <- paste(dat$chr,dat$POS, sep = ":")
dat$minorAllele <- ifelse(dat$ALT_Frq <= 0.5, as.character(dat$ALT), as.character(dat$REF))
dat$majorAllele <- ifelse(dat$ALT_Frq <= 0.5, as.character(dat$REF), as.character(dat$ALT))
dat$beta <- ifelse(dat$ALT_Frq <= 0.5, dat$Beta, dat$Beta*-1)
dat$se <- dat$SE
dat$maf <- ifelse(dat$ALT_Frq <= 0.5, dat$ALT_Frq, 1 - dat$ALT_Frq)
dat$P <- dat$Pvalue
dat0 <- dat[,c("markerID","minorAllele","majorAllele","beta","se","maf","P")]
write.table(dat0, file=paste("toMeta.FILE.tab"), quote = F, sep = "\t", row.names = F)
```

## Create a metal file 
Create a metal file that looks like below and save it as my_metal_analysis.txt

```
#../generic-metal/metal metalAll.txt
#THIS SCRIPT EXECUTES AN ANALYSIS OF EIGHT STUDIES
#THE RESULTS FOR EACH STUDY ARE STORED IN FILES Inputfile1.txt THROUGH Inputfile8.txt
SCHEME  STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
LABEL TotalSampleSize as N # If input files have a column for the sample size labeled as 'N'
# LOAD THE FIRST SEVEN INPUT FILES

# UNCOMMENT THE NEXT LINE TO ENABLE GenomicControl CORRECTION
# GENOMICCONTROL ON

# === DESCRIBE AND PROCESS THE FIRST INPUT FILE ===
MARKER markerID
ALLELE minorAllele majorAllele
FREQ   maf
EFFECT beta
STDERR se
PVALUE P
WEIGHT N 
PROCESS toMeta.GWAS1.tab

# === DESCRIBE AND PROCESS THE SECOND INPUT FILE ===
MARKER markerID
ALLELE minorAllele majorAllele
FREQ   maf
EFFECT beta
STDERR se
PVALUE P
WEIGHT N
PROCESS toMeta.GWAS2.tab


# === DESCRIBE AND PROCESS THE SIXTH INPUT FILE ===
MARKER markerID
ALLELE minorAllele majorAllele
FREQ   maf
EFFECT beta
STDERR se
PVALUE P
WEIGHT N
PROCESS toMeta.GWAS3.tab


OUTFILE MY_THREE_GWAS_FILES_META .tbl
ANALYZE HETEROGENEITY

QUIT
``` 

## Run metal

Then run metal like this:

```
metal metal.txt
```

## Output will look something like this:
These are a couple of lines from one a recent GWAS, including 17 datasets. 

```
MarkerName      Allele1 Allele2 Freq1   FreqSE  MinFreq MaxFreq Effect  StdErr  P-value Direction       HetISq  HetChiSq        HetDf   HetPVal TotalSampleSize
chr4:90666041   t       c       0.6036  0.0229  0.5505  0.6295  0.6980  0.1169  2.348e-09       +-++++++-++++++++       40.4    26.827  16      0.04344 138030
chr3:8379719	a	g	0.2688	0.0101	0.2534	0.2884	0.0006	0.1582	0.997	-+-?+?-+??--++--+	0.0	8.536	12	0.742 140900

```
A couple of notes:

For the Direction column three options are possible +,- and ?. ? means the variant was not present in that dataset, + means positive beta and - means negative beta value.
The HetDf value is the number of included datasets for the variant -1 
Note the it is also recommend to filter for HetISq, anything higher than 80 is not reliable.

Also optional sorting and filtering:
```
# column 10 is the p-value column
sort -gk 10 META_ALL_DATA1.tbl > SORTED_META_ALL_DATA1.tbl
```

### References:

METAL: https://www.ncbi.nlm.nih.gov/pubmed/20616382


# 5. Results visualization 
Results visualization are crucial to explain and inspect the results. You can look at your results using Manhattan plot and QQ-plot. Also it is good to check the Lambda value from your GWAS.


## Manhattan plot

See here a great way of making a Manhattan plot:
https://github.com/ipdgc/Manhattan-Plotter

## QQ-plot

Making the QQ plot is a way to vizualize the potential inflation

```
data = read.table("YOUR_summaray_stats.TBL",header=T)
# optional subsetting e.g. number of datasets or MAF
# newdata <- subset(data, HetDf > 9)
# newdata <- subset(data, MAF > 0.05)
observed <- sort(p)
lobs <- -(log10(observed))
expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected)+1)))
# can update the name of output file if needed here, can also change to pdf("file.pdf")
png("qqplot.png")
# note that the range is here set to 10 on both X and Y axis you might want to change this if you have more significant hits
plot(c(0,10), c(0,10), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,10), ylim=c(0,10), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
dev.off()
```

## Lambda value
```
data = read.table("YOUR_summaray_stats.TBL",header=T)
# optional subsetting e.g. number of datasets or MAF
# newdata <- subset(data, HetDf > 9)
# newdata <- subset(data, MAF > 0.05)
p <- newdata$Pvalue
n <- length(newdata$Pvalue)
x2obs <- qchisq(as.numeric(as.character(p)), 1, lower.tail = FALSE)
x2exp <- qchisq((1:n)/n, 1, lower.tail = FALSE)
lambda <- median(x2obs)/median(x2exp)
lambda
# lambda1000 can be used to normalize unbalanced case-control numbers
# change the Ncases and Ncontrols
lambda1000 <- 1 + ( lambda -1 ) * (1/'Ncases' + 1/'Ncontrols')/(1/1000 + 1/1000)
lambda1000
```












# That's it.... Good luck!





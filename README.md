# GWAS-pipeline

Date: June 2019

Authors: Cornelis Blauwendraat, Mike Nalls

## General description and purpose
This is a general description of the LNG GWAS pipeline which can be used to guide researchers on how-to run a GWAS
Roughly the pipeline can be divided into four steps:

1. QC and data cleaning
2. Imputation
3. GWAS
4. Optional meta-analyses


# 1. QC and data cleaning
The QC and data cleaning is very important prior imputation since the cleaner data you put in there the better and more accurate data you will get out.

## Sample QC parameters



## Variant QC parameters


## preparation of data prior to submission to Imputation server

prepare data for HRC imputation Will Rayner tool


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

## Covariate files




## RVTESTS
Why RVTESTS? It is very easy to use due to the similar file structure as PLINK takes, it has a good manual and it has very flexible options to use.

Files needed to run GWAS:
1. Phenotype file
2. Regions file
3. Imputed genotype data

### Phenotype file
Structure is very similar to PLINK:
.....


### Regions file
This makes sure that you don't include low (imputation) quality and low frequency variants.
The code below creates two files:
- maf001rsq03minimums_chr22.info
- maf001rsq03minimums_chr22.txt


```
R
library(plyr)
# note1 you can speed this up with the data.table package and changing read.table and write.table to fread and fwrite
# note2 if you already unzipped your info.gz files update this to .info
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

This is the .vcf file that is the output from the imputation server


## PLINK
Why PLINK? It is very easy to use, it has a good manual and it has very flexible options to use.

### Phenotype file
Structure is very similar to PLINK:
.....


### Regions file
This makes sure that you don't include low (imputation) quality and low frequency variants.

```


some code



```

### Imputed genotype data

This is the .vcf file that is the output from the imputation server


### References:

RVTESTS:

PLINK:


# 4. Optional meta-analyses
Meta-analyses are used when you have multiple datasets 


### References:

METAL:

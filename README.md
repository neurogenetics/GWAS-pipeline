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


# QC and data cleaning
The QC and data cleaning is very important prior imputation since the cleaner data you put in there the better and more accurate data you will get out.

## Sample QC parameters



## Variant QC parameters


## preparation of data prior to submission to Imputation server

prepare data for HRC imputation Will Rayner tool


# Imputation
Imputation will be done using the Michigan Imputation Server (BIG THANKS TO THEM FOR MAKING OUR LIVES SO MUCH EASIER!!)

https://imputationserver.sph.umich.edu/index.html

Using Reference panel: HRC 1.1 2016 and Eagle v2.3 Phasing

### More info on reference panel here :

https://imputationserver.readthedocs.io/en/latest/reference-panels/ and http://www.haplotype-reference-consortium.org/

### Output of Imputation server

.vcf.gz files of all your input chromosomes e.g. for chromosome 21:
```
chr21.dose.vcf.gz.tbi
chr21.dose.vcf.gz
chr21.info.gz
```

### References:

Imputation server: https://www.ncbi.nlm.nih.gov/pubmed/27571263

HRC panel: https://www.ncbi.nlm.nih.gov/pubmed/27548312

Eagle Imputation: https://www.ncbi.nlm.nih.gov/pubmed/27694958


# GWAS
For running GWAS many many tools and programs are available. We most commonly use either RVTESTS or PLINK. 

### RVTESTS
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

```


some code



```

### Imputed genotype data

This is the .vcf file that is the output from the imputation server




### References:

RVTESTS:

PLINK:


# Optional meta-analyses



### References:

METAL:

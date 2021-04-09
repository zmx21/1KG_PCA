# 1KG PCA
## Description
Overlay samples onto 1KG reference in genetic principal component (PC) space. 

## Dependencies
 - GATK (version 1.93.2)
 - PLINK (v1.90)
 - PLINK2 (v2.00)
 - bcftools (v1.10.2)

## Data
### Download VCF Files:
`./src/download_1000genomes_vcf 1KG_data/`
### Concatenate VCF Files:
`bcftools concat ALL.chr{1..22}*vcf.gz -O z -o 1KG_data/joined.1000genomes.vcf.gz`
### Create VCF File Index:
`bcftools index -t 1KG_data/joined.1000genomes.vcf.gz`

## Usage
### Example 1: To overlay samples onto all 1KG superpopulations:
`./run_pca -vcf VCF_FILE -o OUT_DIR -p PREFIX`
 - `VCF_FILE` should be a bgzip zipped vcf file of the study cohort. hg19/GRCh37.
 - `OUT_DIR` specifies directory of output files/plots
 - `PREFIX` specifies file prefixes of output 
 
### Example 2: To overlay samples to a specific 1KG subpopulation:
`./run_pca -vcf VCF_FILE -o OUT_DIR -p PREFIX -pop POP`
 - `POP` should be EUR, EAS, SAS, AFR or AMR
### Additional options
`./run_pca -h` for full list of commands. 
 - `-d` path to directory with required 1KG data
 - `-s` directory to software executables 
 - `-maf` specifies minor allele frequency cutoff 
 - `-hwe` specified hardy-weinberg equilibrium p-val cutoff
 - `-t` number of CPU threads 

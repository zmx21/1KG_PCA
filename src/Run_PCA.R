if (!require('glue')) install.packages('glue'); library('glue')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
RunPCA1KG <- function(VCF_path,out_dir,software_dir = NULL,data_dir = './data/',prefix = '',maf_thresh = 0.01,hwe_thresh = 1e-6,n_cores = 20,sub_pop = NULL){

  if(!is.null(software_dir)){
    PLINK <- glue::glue("{software_dir}plink")
    bcftools <- glue::glue("{software_dir}bcftools")
    GCTA <-glue::glue("{software_dir}gcta")
  }else{
    PLINK <- glue::glue("plink")
    bcftools <- glue::glue("bcftools")
    GCTA <-glue::glue("gcta")
  }
  
  #Set Path and variables
  ref_file <- glue::glue("{data_dir}/joined.1000genomes.vcf.gz")
  excl_regions <- glue::glue("{data_dir}/exclusion_regions_hg19.txt") #https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium
  super_pop <- data.table::fread(glue::glue('{data_dir}/20131219.populations.tsv'))
  tbl_1KG <- data.table::fread(glue::glue('{data_dir}/20130606_g1k.ped')) %>% dplyr::left_join(super_pop,by = c('Population' = 'Population Code'))
  
  #Make dir to store tmp files
  tmp_dir <- glue::glue("{out_dir}/tmp/")
  system(glue::glue("mkdir -p {tmp_dir}"))
  
  ## TB_DAR ###
  #Remove Duplicate entries in VCF
  no_dup_vcf <- basename(gsub(x=VCF_path,pattern = '.vcf.gz',replacement = '.nodup.vcf.gz'))
  system(
    glue::glue("{bcftools} norm -d snps {VCF_path} | {bcftools} view -O z -o {tmp_dir}{no_dup_vcf}")
  )
  system(
    glue::glue("{bcftools} index -t {tmp_dir}{no_dup_vcf}")
  )
  
  bim_filt <- gsub(x=no_dup_vcf,pattern = '.vcf.gz',replacement = '.hwe')
  #Filter based on HWE and set SNP ID
  system(
    glue::glue(
      "{PLINK}2 --vcf {tmp_dir}{no_dup_vcf} --threads {n_cores} --hwe {hwe_thresh} --make-bed --keep-allele-order --set-all-var-ids @:#[b37]\\$r,\\$a  --out {tmp_dir}{bim_filt}"
    )
  )
  system(
    glue::glue(
      "{PLINK}2 --bfile {tmp_dir}{bim_filt} --threads {n_cores} --keep-allele-order --export vcf --const-fid --out {tmp_dir}{bim_filt}"
    )
  )
  
  # #Remove long range LD region
  system(
    glue::glue(
      "{PLINK}2 --bfile {tmp_dir}{bim_filt} --threads {n_cores} --exclude range {excl_regions} --make-bed --keep-allele-order --out {tmp_dir}{prefix}_tmp1 --autosome"
    )
  )

  ### 1KG ###
  #Solve SNP ID issue
  system(
    glue::glue(
      "{PLINK}2 --vcf {ref_file} --threads {n_cores} --make-bed --keep-allele-order --set-all-var-ids @:#[b37]\\$r,\\$a --out {tmp_dir}{prefix}_1KG --autosome"
    )
  )
  snps_1KG <- data.table::fread(glue::glue("{tmp_dir}{prefix}_1KG.bim"))
  snps_TB_DAR <- data.table::fread(glue::glue("{tmp_dir}{prefix}_tmp1.bim"))
  consensus_snps <- dplyr::inner_join(snps_TB_DAR,snps_1KG,by=c('V1'='V1','V4'='V4','V5'='V5','V6'='V6'))
  write(consensus_snps$V2.y,file = glue::glue("{tmp_dir}{prefix}_1KG_consensus_snps.txt"))
  system(
    glue::glue(
      "{PLINK}2 --bfile {tmp_dir}{prefix}_1KG --threads {n_cores} --keep-allele-order --make-bed --extract {tmp_dir}{prefix}_1KG_consensus_snps.txt --out {tmp_dir}{prefix}_1KG_consensus --autosome"
    )
  )

  #Generate Cohort with only consensus SNPs
  system(
    glue::glue(
      "{PLINK}2 --bfile {tmp_dir}{prefix}_tmp1 --threads {n_cores} --keep-allele-order --make-bed --extract {tmp_dir}{prefix}_1KG_consensus_snps.txt --out {tmp_dir}{prefix}_consensus --autosome"
    )
  )

  #Merge Cohort and 1KG
  system(
    glue::glue(
      "{PLINK} --bfile {tmp_dir}{prefix}_1KG_consensus --bmerge {tmp_dir}{prefix}_consensus --threads {n_cores} --keep-allele-order --make-bed --out {tmp_dir}{prefix}_1KG_Cohort_consensus"
    )
  )

  #Prune (r2 < 0.2)
  system(
    glue::glue(
      "{PLINK} --bfile {tmp_dir}{prefix}_1KG_Cohort_consensus --threads {n_cores} --keep-allele-order --indep-pairwise 200 100 0.2 --out {tmp_dir}{prefix}_1KG_Cohort_consensus --autosome"
    )
  )

  system(
    glue::glue(
      "{PLINK} --bfile {tmp_dir}{prefix}_1KG_Cohort_consensus --threads {n_cores} --keep-allele-order --extract {tmp_dir}{prefix}_1KG_Cohort_consensus.prune.in --make-bed --out {tmp_dir}{prefix}_1KG_Cohort_consensus_pruned --autosome"
    )
  )


  #Calculate PCA on the merged cohort
  n_pcs <- nrow(data.table::fread(glue::glue("{tmp_dir}{prefix}_1KG_Cohort_consensus_pruned.fam")))
  system(
    glue::glue(
      "{GCTA} --bfile {tmp_dir}{prefix}_1KG_Cohort_consensus_pruned --make-grm --out {tmp_dir}{prefix}_1KG_Cohort_consensus_pruned --thread-num {n_cores}"
    )
  )
  system(
    glue::glue(
      "{GCTA} --grm {tmp_dir}{prefix}_1KG_Cohort_consensus_pruned --pca 20 --out {out_dir}{prefix} --thread-num {n_cores}"
    )
  )

  #If subpopulation is specified, calculate PCA based on merge with subpop
  if(sub_pop!='NA'){
    #Calculate PCA on the merged AFR cohort
    Non_SubPop_Samples <- dplyr::filter(tbl_1KG,`Super Population` != sub_pop) %>% dplyr::select(`Individual ID`)
    consensus_fam_file <- data.table::fread(glue::glue("{tmp_dir}{prefix}_1KG_Cohort_consensus.fam")) %>% dplyr::filter(V2 %in% Non_SubPop_Samples$`Individual ID`) %>% dplyr::select(V1,V2)
    data.table::fwrite(consensus_fam_file,col.names = F,row.names = F,sep = ' ',file = glue::glue("{tmp_dir}{prefix}_1KG_Non_{sub_pop}.txt"))
    system(
      glue::glue(
        "{PLINK} --bfile {tmp_dir}{prefix}_1KG_Cohort_consensus_pruned --remove {tmp_dir}{prefix}_1KG_Non_{sub_pop}.txt --keep-allele-order --make-bed --out {tmp_dir}{prefix}_1KG_{sub_pop}_Cohort_consensus_pruned"
      )
    )
    
    system(
      glue::glue(
        "{GCTA} --bfile {tmp_dir}{prefix}_1KG_{sub_pop}_Cohort_consensus_pruned --make-grm --out {tmp_dir}{prefix}_1KG_{sub_pop}_Cohort_consensus_pruned --thread-num {n_cores}"
      )
    )
    system(
      glue::glue(
        "{GCTA} --grm {tmp_dir}{prefix}_1KG_{sub_pop}_Cohort_consensus_pruned --pca 20 --out {out_dir}{prefix}_{sub_pop} --thread-num {n_cores}"
      )
    )
  }
}
args = commandArgs(trailingOnly=TRUE)
VCF_path = args[[1]]
out_dir = args[[2]]
prefix = args[[3]]
sub_pop = args[[4]]
data_dir = args[[5]]
software_dir = args[[6]]
maf_thresh = as.numeric(args[[7]])
hwe_thresh = as.numeric(args[[8]])
n_cores = as.numeric(args[[9]])

RunPCA1KG(VCF_path = VCF_path,out_dir = out_dir,
          prefix = prefix,sub_pop=sub_pop,
          data_dir=data_dir,software_dir=software_dir,
          maf_thresh=maf_thresh,hwe_thresh=hwe_thresh,
          n_cores=n_cores)
super_pop <- data.table::fread(glue::glue('{data_dir}20131219.populations.tsv'))
tbl_1KG <- data.table::fread(glue::glue('{data_dir}20130606_g1k.ped')) %>% dplyr::left_join(super_pop,by = c('Population' = 'Population Code'))

#Load PCA results
pc_df <- data.table::fread(glue::glue("{out_dir}{prefix}.eigenvec"))
colnames(pc_df) <- c('FID','IID',paste0('PC',seq(1,ncol(pc_df) - 2)))
#Plot with all Superpopulations
merged_df <- pc_df %>%
  dplyr::left_join(tbl_1KG,c('IID'='Individual ID'))
merged_df$`Data Set` <- as.factor(ifelse(is.na(merged_df$Population),'Study Cohort','1000 Genomes'))
merged_df$`Super Population` <- as.factor(as.character(merged_df$`Super Population`))
merged_df$`Super Population` <-factor(merged_df$`Super Population` ,
                                      levels = rev(levels(merged_df$`Super Population`)))

eigen_val <- as.vector(data.table::fread(glue::glue("{out_dir}{prefix}.eigenval")))
var_explained <- eigen_val$V1 / sum(eigen_val$V1) * 100
pc_plot_thousand_genome_PC1_PC2 <- ggplot(data = merged_df) +
  aes(x=PC1,y=PC2,color=`Super Population`,shape=`Data Set`) +
  geom_point() + scale_shape_manual(values = c(20,4)) +
  xlab(paste0('PC1 (',signif(var_explained[1],2),'%)'))  +
  ylab(paste0('PC2 (',signif(var_explained[2],2),'%)'))
ggsave(pc_plot_thousand_genome_PC1_PC2,filename = glue::glue("{out_dir}PCA_Plot.png"))

if(sub_pop!='NA'){
  #Supopulation 1KG
  pc_df_Pop <- data.table::fread(glue::glue("{out_dir}{prefix}_{sub_pop}.eigenvec"))
  colnames(pc_df_Pop) <- c('FID','IID',paste0('PC',seq(1,ncol(pc_df_Pop) - 2)))
  merged_df_Pop<- pc_df_Pop %>%
    dplyr::left_join(tbl_1KG,c('IID'='Individual ID'))
  merged_df_Pop$`Data Set` <- as.factor(ifelse(is.na(merged_df_Pop$Population),'Study Cohort','1000 Genomes'))
  
  eigen_val_Pop <- as.vector(data.table::fread(glue::glue("{out_dir}{prefix}_{sub_pop}.eigenval")))
  var_explained_Pop <- eigen_val_Pop / sum(eigen_val_Pop) * 100
  
  pc_plot_Pop_pc1_pc2 <- ggplot(data = merged_df_Pop) +
    aes(x=PC1,y=PC2,color=`Population Description`,shape=`Data Set`) +
    geom_point() + scale_shape_manual(values = c(20,4)) +
    xlab(paste0('PC1 (',signif(var_explained_Pop[1],2),'%)'))  +
    ylab(paste0('PC2 (',signif(var_explained_Pop[2],2),'%)'))
  ggsave(pc_plot_Pop_pc1_pc2,filename = glue::glue("{out_dir}PCA_{sub_pop}_Plot.png")) 
}

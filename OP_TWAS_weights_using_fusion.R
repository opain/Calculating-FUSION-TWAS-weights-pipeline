#!/usr/bin/Rscript
library("optparse")

option_list = list(
  make_option("--PLINK_prefix", action="store", default=NA, type='character',
              help="Prefix of PLINK files for the target sample [required]"),
  make_option("--phenotype_file", action="store", default=NA, type='character',
              help="File name for normalised and adjusted expression data [required]"),
  make_option("--coordinate_file", action="store", default=NA, type='character',
              help="File name for coordinates of genes/transcripts [required]"),
  make_option("--gene_name", action="store", default=NA, type='character',
              help="Name of gene or transcript to be processed [required]"),
  make_option("--plink", action="store", default=NA, type='character',
              help="Path to PLINK [required]"),
  make_option("--gcta", action="store", default=NA, type='character',
              help="Path to gcta_nr_robust binary [required]"),
  make_option("--gemma", action="store", default=NA, type='character',
              help="Path to gemma binary [required]"),
  make_option("--ld_ref_dir", action="store", default=NA, type='character',
              help="FUSION LD reference directory [required]"),
  make_option("--fusion_software", action="store", default=NA, type='character',
              help="FUSION software directory [required]"),
  make_option("--output_dir", action="store", default=NA, type='character',
              help="Directory name for the output [required]")
)

opt = parse_args(OptionParser(option_list=option_list))
system(paste('mkdir ',opt$output_dir,'/Output',sep=''))

root_dir<-getwd()

library(data.table)

# Read in the gene expression data to extract the gene's chromosome and boundary coordinates
Gene_coordinates_file<-fread(opt$coordinate_file, header=T)

CHR<-Gene_coordinates_file$X.Chr[Gene_coordinates_file$ID == opt$gene_name]
Start<-Gene_coordinates_file$start[Gene_coordinates_file$ID == opt$gene_name]-0.5e6
Stop<-Gene_coordinates_file$end[Gene_coordinates_file$ID == opt$gene_name]+0.5e6

if (Start < 0) {
Start <- 0
}

# Extract gene from phenotype file
system(paste('mkdir ',opt$output_dir,'/temp',sep=''))
system(paste('awk -f ./t.awk c1=FID c2=IID c3=',opt$gene_name,' ',opt$phenotype_file,' > ',opt$output_dir,'/temp/temp_',opt$gene_name,'.pheno',sep=''))

# Using PLINK, extract variants within the specified gene +/- 500kb from start and stop coordinates
err_1<-system(paste(opt$plink,' --bfile ',opt$PLINK_prefix,' --make-bed --pheno ',opt$output_dir,'/temp/temp_',opt$gene_name,'.pheno', ' --out ', opt$output_dir,'/temp/temp_',opt$gene_name,' --chr ',CHR,' --from-bp ',Start,' --to-bp ',Stop,' --extract ',opt$ld_ref_dir,'/1000G.EUR.',CHR,'.bim', sep=''))
if (err_1 == 13) {
write.table('No SNPs within gene +/- 0.5Mb', paste(opt$output_dir,'/Output/',opt$gene_name,'.err',sep=''), col.names=F, row.names=F, quote=F)
} else {

# Using FUSION, calculate the weights for the genes expression using subset of genotypic data.
setwd(paste(opt$output_dir,'/temp', sep=''))
system(paste('ln -s ./ output', sep=''))
system(paste('Rscript ',opt$fusion_software,'/FUSION.compute_weights.R --bfile ', opt$output_dir,'/temp/temp_',opt$gene_name,' --tmp temp_',opt$gene_name,'.tmp --out ', opt$output_dir,'/Output/',opt$gene_name,' --verbose 0 --save_hsq --PATH_gcta ',opt$gcta,' --PATH_gemma ',opt$gemma,' --PATH_plink ',opt$plink,' --models top1,blup,lasso,enet', sep=''))
}

# Delete the temporary files
system(paste('rm ',opt$output_dir,'/temp/temp_',opt$gene_name,'*', sep=''))

#!/usr/bin/Rscript
library("optparse")

option_list = list(
  make_option("--ld_ref_dir", action="store", default=NA, type='character',
              help="Directory where FUSION LD reference files are stored [required]"),
  make_option("--PLINK_prefix", action="store", default=NA, type='character',
              help="Prefix of PLINK files for the target sample [required]"),
  make_option("--output", action="store", default=NA, type='character',
              help="File name for the output [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

library(data.table)

setwd(paste(opt$ld_ref_dir,"/",sep=''))
temp = list.files(pattern="*.bim")

Ref<-do.call(rbind, lapply(temp, function(x) data.frame(fread(x))))

# Read in the SNPs in CLOZUK
Target<-data.frame(fread(paste(opt$PLINK_prefix,'.bim',sep='')))

# Get intersect of the two based on RSID
Overlap<-intersect(Ref$V2, Target$V2)

sink(file = opt$output)
cat('Number of SNPs in FUSION LD Reference = ',length(Ref$V2),'
Number of SNPs in target PLINK files = ',length(Target$V2),'
Number of SNPs in both = ',length(Overlap),'
Percentage of SNPs in FUSION LD Reference that are in the target PLINK files = ',length(Overlap)/length(Ref$V2)*100,'%',sep='')
sink()


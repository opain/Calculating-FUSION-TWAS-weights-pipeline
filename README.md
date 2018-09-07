# Pipeline for creating SNP-weights for TWAS using FUSION

This pipeline is to create SNP-weights that can be used in TWAS using FUSION. I provide a series of R scripts to process the input and output files for the FUSION script called 'FUSION.compute_weights.R'.  I also provide all the software and data required for the pipeline to work. Download this repository to your ROCKS account, and follow the instructions below.



## Getting started

### Required data

1. Genetic data that has undergone standard quality control and imputation, and is in binary PLINK format (.bed/.bim/.fam).
2. A file containing phenotypic data (gene expression, methylation, chromatin modification etc.) in PLINK format (e.g. FID, IID, gene1, gene2, gene3...). This data should have been corrected for covariates and normalised. An example has been provided ('Example_phenotype_file.txt').
3. A file containing the coordinates of each feature (e.g. gene, methlyation site). For gene expression this could be start and stop coordinates. The file should be in the following order: Chromosome, start coordinate, end coordinate, and name of the feature (should match with phenotype columns names). The files should have a header but the column names will be ignored. An example has been provided ('Example_coordinate_file.txt').



### Required software

* R and the required packages:

```R
install.packages(c('data.table','optparse','foreach','doMC'))
```

* FUSION software:

```sh
git clone https://github.com/gusevlab/fusion_twas.git
```
* FUSION LD reference data ([download](https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2))

* [GEMMA 0.97 software](https://github.com/genetics-statistics/GEMMA)

  * Note: ROCKS users at Cardiff need to create a shell which specifies the correct compiler (see [example](http://gitlab.psycm.cf.ac.uk/mpmop/Calculating-FUSION-TWAS-weights-pipeline/blob/master/gemma.97.sh))

    

### Provided software

| Name                            | Description                                                  |
| ------------------------------- | ------------------------------------------------------------ |
| OP_fusion_ref_overlap_checker.R | A script that reads in the FUSION reference, the target sample .bim file, and produces a summary of the overlap. |
| OP_TWAS_weights_using_fusion.R  | Prepares and inputs the genotypic and phenotype data for the FUSION.compute_weights.R script. This script: 1) extracts the phenotype data for the specified feature and puts it into a PLINK phenotype file, 2) creates PLINK genotype files that only contains SNPs within the feature region +/-0.5Mb and are present in the FUSION 1KG reference, 3) inserts the feature specific phenotype data into the PLINK files, 4) inputs the PLINK files into the FUSION.compute_weights.R script, specifying the default settings (Note: BSLMM disabled), and 5) deletes the temporary files if the job completes. |
| OP_TWAS_weights_using_fusion.sh | A shell script for OP_TWAS_weights_using_fusion.R that allows all features to be run in an array. |
| OP_packaging_fusion_weights.R   | Creates folder containing SNP-weights in the same format as the FUSION released SNP-weights. |



## Instructions:

### 0 - Set working directory and create required variables

```sh
# Download this gitlab repository as .tar.gz
# Decompress 'Calculating-FUSION-TWAS-weights-pipeline.tar.gz'
tar -xzf Calculating-FUSION-TWAS-weights-pipeline.tar.gz

# Set the working directory to the folder 'Calculating-FUSION-TWAS-weights-pipeline'
cd Calculating-FUSION-TWAS-weights-pipeline

# Give permissions to fetal_weights_using_fusion-full_sample-gene.R.sh
chmod a+x fetal_weights_using_fusion-full_sample-gene.R.sh

# Create a directory for the output
mkdir <Desired name of output directory>

# Create the following variables
target_plink=<Path to prefix of PLINK files>
phenotype_file=<Path to file containing phenotypic data>
coordinate_file=<Path to file containing feature coordinate data>
fusion_ldref=<Path to folder containing fusion ld reference>
plink=<Path to plink binary>
gemma=<Path to gemma binary>
gcta=<Path to gcta_nr_robust binary>
output=<Name of the previously created output directory>

# For example, only account it could like this
cd /home/mpmop/Calculating-FUSION-TWAS-weights-pipeline
chmod a+x ./OP_TWAS_weights_using_fusion.sh
mkdir pipe_test
target_plink=/home/mpmop/Data/Heath/Fetal_all/genotypes
phenotype_file=/home/mpmop/Analyses/Fetal_GeneX_weights/Clean/GeneX/GeneX_norm_resid.pheno
coordinate_file=/home/mpmop/Analyses/Fetal_GeneX_weights/Clean/GeneX/Gene_locations.txt
fusion_ldref=/home/mpmop/Software/fusion_twas-master/LDREF
fusion_software=/home/mpmop/Software/fusion_twas-master
plink=/share/apps/plink2
gemma=/home/mpmop/Software/gemma.97.sh
gcta=/home/mpmop/Software/fusion_twas-master/gcta_nr_robust
output=/home/mpmop/Calculating-FUSION-TWAS-weights-pipeline/pipe_test
```



### 1. Check the percentage of SNPs in the FUSION LD reference available in the target data

```sh
Rscript ./OP_fusion_ref_overlap_checker.R \
    --ld_ref_dir ./fusion_twas-master/LDREF \
    --PLINK_prefix ${target_plink} \
    --output ${output}/overlap.txt
```



### 2. Run script that prepares and inputs the genotypic and phenotype data for the FUSION.compute_weights.R script.

##### - If you only want to create weights for one of the features in your phenotype file: 

```sh
# Assign the name of the feature to a variable called feature_name
feature_name=<Name of the feature>

# Create weights using 'OP_TWAS_weights_using_fusion.R' script
qsub -cwd -b y -l h_vmem=5G,mem_free=5G -e /dev/null -o /dev/null Rscript ./OP_TWAS_weights_using_fusion.R \
--PLINK_prefix ${target_plink} \
--phenotype_file ${phenotype_file} \
--coordinate_file ${coordinate_file} \
--gene_name ${feature_name} \
--plink ${plink} \
--gcta ${gcta} \
--gemma ${gemma} \
--ld_ref_dir ${fusion_ldref} \
--fusion_software ${fusion_software} \
--output_dir ${output}

```

##### - If you want to calculate weights for all features in your phenotype file:

```sh
# Use the shell script 'OP_TWAS_weights_using_fusion.sh' to submit each job in an array
# Or run the script for all features in an array.
qsub -t 2-$(wc -l ${coordinate_file} | cut -d' ' -f1) -N weights_calc -cwd -b y -l h_vmem=5G,mem_free=5G -e /dev/null -o /dev/null ./OP_TWAS_weights_using_fusion.sh \
--PLINK_prefix ${target_plink} \
--phenotype_file ${phenotype_file} \
--coordinate_file ${coordinate_file} \
--plink ${plink} \
--gcta ${gcta} \
--gemma ${gemma} \
--ld_ref_dir ${fusion_ldref} \
--fusion_software ${fusion_software} \
--output_dir ${output}

# Check that the jobs didn't crash. Could use qacct to idenitfy crashed jobs.
# Identify which job IDs failed and create a coordinates file only containing these features.
paste <(qacct -j weights_calc | grep 'taskid') <(qacct -j weights_calc | grep 'failed') | tr -s ' ' | cut -d ' ' -f 2,4 | awk -F" " '$2 != "0" { print $1 }' > failed
cat <(head -n 1 ${coordinate_file}) <(awk 'NR==FNR{ pos[$1]; next }FNR in pos' failed ${coordinate_file}) > rerun_coordinate_file

# Re-run scripts for features that crashed. Change the job name so they can be distinguished from the previous set of jobs.
qsub -t 2-$(wc -l rerun_coordinate_file | cut -d' ' -f1) -N weights_calc_2 -cwd -b y -l h_vmem=5G,mem_free=5G -e /dev/null -o /dev/null ./OP_TWAS_weights_using_fusion.sh \
--PLINK_prefix ${target_plink} \
--phenotype_file ${phenotype_file} \
--coordinate_file rerun_coordinate_file \
--plink ${plink} \
--gcta ${gcta} \
--gemma ${gemma} \
--ld_ref_dir ${fusion_ldref} \
--fusion_software ${fusion_software} \
--output_dir ${output}

# Check again until there are no crashed jobs
paste <(qacct -j weights_calc_2 | grep 'taskid') <(qacct -j weights_calc_2 | grep 'failed') | tr -s ' ' | cut -d ' ' -f 2,4 | awk -F" " '$2 != "0" { print $1 }' > failed

```



### 3. Create folder containing SNP-weights in the same format as the FUSION released SNP-weights

```sh
Rscript ./OP_packaging_fusion_weights.R \
--RDat_dir ${output}/Output \
--coordinate_file ${coordinate_file} \
--output_name TWAS_Weights \
--output_dir ${output}/TWAS_weights_package
```



## Help

This script was written by Dr Oliver Pain under the supervision of Dr Richard Anney whilst at the MRC Centre for Neuropsychiatric Genetics and Genomics, Cardiff University.

If you have any questions or comments please email Ollie (paino@cardiff.ac.uk).








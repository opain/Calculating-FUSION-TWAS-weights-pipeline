# Pipeline for creating SNP-weights for TWAS using FUSION

This pipeline is to create SNP-weights that can be used in TWAS using FUSION. I provide a series of R scripts to process the input and output files for the FUSION script called 'FUSION.compute_weights.R'.  I also provide all the software and data required for the pipeline to work. Download this repository to your ROCKS account, and follow the instructions below.



## Getting started

### Required data

1. 	Genetic data that has undergone standard quality control and imputation, and is in binary PLINK format (.bed/.bim/.fam).
2. 	A file containing phenotypic data (gene expression, methylation, chromatin modification etc.) in PLINK format (e.g. FID, IID, gene1, gene2, gene3...). This data should have been corrected for covariates and normalised. An example has been provided ('Example_phenotype_file.txt').

	3. 	A file containing the coordinates of each feature (e.g. gene, methlyation site). For gene expression this could be start and stop coordinates. The file should be in the following order: Chromosome, start coordinate, end coordinate, and name of the feature (should match with phenotype columns names). The files should have a header but the column names will be ignored. An example has been provided ('Example_coordinate_file.txt').
	



### Provided software and data

1. 	fusion_twas-master - This contains all the fusion scripts and the 1KG reference released by FUSION.

2. 	plink2 - PLINK v1.90b5.4

3. 	gcta_nr_robust - GCTA binary released by FUSION that enables robust non-linear optimization.

4. 	gemma.97.sh - GEMMA v0.97 binary in a wrapper that selects the correct compiler on ROCKS.

5. 	OP_fusion_ref_overlap_checker.R - A script that reads in the FUSION reference, the target sample bim file, and produces a summary of the overlap.

6. 	OP_TWAS_weights_using_fusion.R - Prepares and inputs the genotypic and phenotype data for the FUSION.compute_weights.R script. This script 1) extracts the phenotype data for the specified feature and puts it into a PLINK phenotype file, 2) creates PLINK genotype files that only contains SNPs within the feature region +/-0.5Mb and are present in the FUSION 1KG reference, 3) inserts the feature specific phenotype data into the PLINK files, 4) inputs the PLINK files into the FUSION.compute_weights.R script, specifying the default settings (Note: BSLMM disabled), and 5) deletes the temporary files if the job completes.

7. 	OP_TWAS_weights_using_fusion.sh - A shell script for OP_TWAS_weights_using_fusion.R that allows all features to be run in an array.

8. 	OP_packaging_fusion_weights.R- Creates folder containing SNP-weights in the same format as the FUSION released SNP-weights.



## Instructions:

### 0 - Set working directory and create required variables

```
# Decompress 'Calculating-FUSION-TWAS-weights-pipeline.tar.gz'
tar -xzf Calculating-FUSION-TWAS-weights-pipeline.tar.gz

# Set the working directory to the folder 'Calculating-FUSION-TWAS-weights-pipeline'
cd Calculating-FUSION-TWAS-weights-pipeline

# Create a directory for the output
mkdir <Desired name of output directory>

# Create 'target_plink' variable containing prefix to PLINK files
target_plink= <Path to prefix of PLINK files>

# Create 'phenotype_file' variable containing phenotypic data
phenotype_file= <Path to file containing phenotypic data>

# Create 'coordinate_file' variable containing coordinates of each feature data
coordinate_file= <Path to file containing coordinate data>

# Create 'output' variable containing the directory name for the various outputs
output= <Name of the previously created output directory>

# Give permissions to fetal_weights_using_fusion-full_sample-gene.R.sh
chmod a+x fetal_weights_using_fusion-full_sample-gene.R.sh

```



### 1. Check the percentage of SNPs in the FUSION LD reference available in the target data

```
Rscript ./OP_fusion_ref_overlap_checker.R \
--ld_ref_dir ./fusion_twas-master/LDREF \
--PLINK_prefix ${target_plink} \
--output ${output}/overlap.txt
```



### 2. Run script that prepares and inputs the genotypic and phenotype data for the FUSION.compute_weights.R script.

##### - If you only want to create weights for one of the features in your phenotype file: 

```
# Assign the name of the feature to a variable called feature_name
feature_name=<Name of the feature>

# Create weights using 'OP_TWAS_weights_using_fusion.R' script
qsub -cwd -b y -l h_vmem=5G,mem_free=5G -e /dev/null -o /dev/null Rscript ./OP_TWAS_weights_using_fusion.R \
--PLINK_prefix ${target_plink} \
--phenotype_file ${phenotype_file} \
--coordinate_file ${coordinate_file} \
--gene_name ${feature_name} \
--output_dir ${output}
```

##### - If you want to calculate weights for all features in your phenotype file:

```
# Use the shell script 'OP_TWAS_weights_using_fusion.sh' to submit each job in an array
qsub -t 2-$(wc -l ${coordinate_file} | cut -d' ' -f1) -cwd -b y -l h_vmem=5G,mem_free=5G -e /dev/null -o /dev/null ./OP_TWAS_weights_using_fusion.sh \
--PLINK_prefix ${target_plink} \
--phenotype_file ${phenotype_file} \
--coordinate_file ${coordinate_file} \
--output_dir ${output}

# Check that the jobs didn't crash. You can use qacct to idenitfy crashed jobs.
# e.g. qacct -j 405301 | grep 'taskid\|failed'
# Then re-run jobs that failed.
```



### 3. Create folder containing SNP-weights in the same format as the FUSION released SNP-weights

```
Rscript ./OP_packaging_fusion_weights.R \
--RDat_dir ${output}/Output \
--coordinate_file ${coordinate_file} \
--output_name TWAS_Weights \
--output_dir ${output}/TWAS_weights_package
```



## Help

This script was written by Dr Oliver Pain under the supervision of Dr Richard Anney whilst at the MRC Centre for Neuropsychiatric Genetics and Genomics, Cardiff University.

If you have any questions or comments please email Ollie (paino@cardiff.ac.uk).








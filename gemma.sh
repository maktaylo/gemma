#!/bin/bash
#SBATCH -D /home/maktaylo/fibr_gwas
#SBATCH -o /home/maktaylo/fibr_gwas/slurm-logs/mark-stdout-%j.txt
#SBATCH -e /home/maktaylo/fibr_gwas/slurm-logs/mark-stderr-%j.txt
#SBATCH -J gemma
#SBATCH --mail-type=END
#SBATCH --mail-user=maktaylor@ucdavis.edu

echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"
echo "SLURM_NODELIST = $SLURM_NODELIST"
echo "SLURM_NODE_ALIASES = $SLURM_NODE_ALIASES"
echo "SLURM_NNODES = $SLURM_NNODES"
echo "SLURM_TASKS_PER_NODE = $SLURM_TASKS_PER_NODE"
echo "SLURM_NTASKS = $SLURM_NTASKS"
echo "SLURM_JOB_ID = $SLURM_JOB_ID"

cd /home/maktaylo/fibr_gwas
traitname="dtb_ptu" #this phenotype is changed each time

mkdir $traitname #make a directory with the traitname that matches phenos[k] from the R code that creates the phenos[k].ped and phenos[k]_ecotypes1.txt files

path=$traitname
cd $traitname

#extract only the ecotypes I need from the master ecotypes file

file="/home/maktaylo/fibr_gwas/1001genomes/rawdata/K_2029/plink_format/k2029_plink.txt"
ecotype_file_direc="/home/maktaylo/fibr_gwas/phenotypes/"
ecotype_file_ext="_ecotypes1.txt" #this file has only 1 column of correctly ordered ecotype id's
ecotype_file=$ecotype_file_direc$traitname$ecotype_file_ext

pedext=".ped"
pedfile=$traitname$pedext 

grep -f $ecotype_file $file > $pedfile

#delete first column in pedfile because it is uneeded trait names by taking all columns after the first
cut -f2- testfile.txt $pedfile > $pedfile

#paste together ped file with phenotype data created in R and ped file with genotype data created above
pheno_file_direc="/home/maktaylo/fibr_gwas/phenotypes/"
pheno_file_ext=".ped" #this file has only 1 column of correctly ordered ecotype id's
pheno_file=$pheno_file_direc$traitname$pheno_file_ext

paste $pheno_file $pedfile > $pedfile

#make map file with correct trait name
mapext=".map"
mapfile=$traitname$mapext
cp /home/maktaylo/fibr_gwas/1001genomes/rawdata/snp-short-indel_only_ACGTN.map $mapfile

#make binary PLINK files to input into GEMMA
/home/maktaylo/plink/plink-1.07-x86_64/plink --file $traitname --noweb --make-bed --out $traitname

###########################################################################################
###########################################################################################
#gemma calls 
output="_output"
###########################################################################################
#run Bayesian sparse linear mixed model  (BSLMM)
out=$traitname$output
/home/maktaylo/GEMMA/gemma -bfile $traitname -bslmm 1 -maf 0.05 -o $out

################################################################################################
#use GEMMA to make relatedness matrix for use in LMM
/home/maktaylo/GEMMA/gemma -bfile $traitname -gk 1 -o $traitname
kmatext=".cXX.txt" #this is the file extension for the relatedness matrix created by gemma
kmat=$path$traitname$kmatext

cd output #must find the relatedness matrix in the output directory that GEMMA creates
cp *.cXX.txt $kmat

################################################################################################
#run EMMA (LMM in GEMMA)
cd $path
kfile=$traitname$kmatext
lmmout="_output_lmm"
lmout=$traitname$lmmout
/home/maktaylo/GEMMA/gemma -bfile $traitname -k $kfile -lmm 4 -maf 0.1 -o $lmout
################################################################################################
#clean out giant phenotype/genotype files that were made files that have been made
rm *ped
rm *map
rm *bim

################################################################################################



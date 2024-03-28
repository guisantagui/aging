#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=findMarksBlood
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 1
#SBATCH --time=01-00:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/aging/scripts/output_scRNAseq_integ_brain.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/aging/scripts/error_scRNAseq_integ_brain.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

module load devel/CMake/3.18.4-GCCcore-10.2.0
module load math/GMP/6.2.0-GCCcore-10.2.0
module load vis/cairo/1.16.0-GCCcore-10.2.0

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################

seuratFile="/home/users/gsantamaria/projects/aging/results/integration/blood_int_cca.rds"
resolution=2
nPCs=100

Rscript findMarkers.R --seuratFile $seuratFile --resolution $resolution --nPCs $nPCs


conda deactivate
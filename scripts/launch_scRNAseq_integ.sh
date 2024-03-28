#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=scRNAinteg
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=500GB
#SBATCH -c 64
#SBATCH --time=00-04:00:00
#Define sdout path
#SBATCH --output=/home/users/gsantamaria/projects/aging/scripts/output_scRNAseq_integ.txt
#Define sderr path
#SBATCH --error=/home/users/gsantamaria/projects/aging/scripts/error_scRNAseq_integ.txt
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p bigmem
#SBATCH --qos=normal

module load devel/CMake/3.18.4-GCCcore-10.2.0
module load math/GMP/6.2.0-GCCcore-10.2.0
module load vis/cairo/1.16.0-GCCcore-10.2.0

conda activate r-4.3.1

# Variables for the pipeline
########################################################################################################################

preProcDir="/home/users/gsantamaria/projects/aging/results/preprocessing/"
method="cca"

orgDirs=$(find "$preProcDir" -mindepth 1 -maxdepth 1 -type d)

for d in $orgDirs; do
  name=$(basename "$d")
  if [ -f "$d" ]; then
    echo "$name already has been integrated."
  else
    echo "Integrating $name samples" &
    Rscript integration_exec.R --preProcDir "$d" --method $method &
  fi
done

wait

conda deactivate
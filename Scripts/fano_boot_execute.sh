#!/bin/bash
#SBATCH --mail-user=user@mail.org
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --output=slurm%j.out
#SBATCH --error=slurm%j.err

module load R/4.1.2-foss-2021b

R CMD BATCH fano_bootstrapping.R

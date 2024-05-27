#!/bin/bash
#SBATCH --mail-user=user@site.org
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --output=slurm%j.out
#SBATCH --error=slurm%j.err


//path_to_software/cellranger count --id=Output_dir_name \
	--fastqs=//path_to_fastqs.dir/ \
	--sample=M15_M_GEX \
	--transcriptome=//path_to_transcriptome/PolyA_10X_genome

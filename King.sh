#!/bin/sh
#SBATCH --job-name=BWA
#SBATCH --output=BWA.out
#SBATCH --error=BWA.err
#SBATCH --time=40:00:00
#SBATCH --mem=10G
#SBATCH --ntasks=8
#SBATCH --partition=rome

# please change the number after --degree
singularity exec /usr/local/biotools/k/king\:2.2.7--hd03093a_0 king -b  PoritesPcyl2503_4_geno01_maf005_hwe_ed.bed --kinship --degree 1

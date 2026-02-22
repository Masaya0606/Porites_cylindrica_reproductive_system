#!/bin/bash
#SBATCH --job-name=refmap
#SBATCH --output=refmap.out
#SBATCH --error=refmap.err
#SBATCH --time=60:00:00
#SBATCH --mem=50G
#SBATCH --ntasks=3
#SBATCH --partition=rome
mkdir -p refmap
singularity exec /usr/local/biotools/s/stacks\:2.60--h9a82719_0 ref_map.pl -T 3 -o /path/refmap  --popmap popmap_2411.txt  --samples /path/bam/ -X "populations: --plink" -X "populations: --phylip-var" -X "populations: --vcf" -X "populations: --treemix"

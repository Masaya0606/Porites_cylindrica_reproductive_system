#!/bin/sh
#SBATCH --job-name=BWA
#SBATCH --output=BWA.out
#SBATCH --error=BWA.err
#SBATCH --time=40:00:00
#SBATCH --mem=10G
#SBATCH --ntasks=8
#SBATCH --partition=rome

module load singularity
singularity exec /usr/local/biotools/b/bwa-mem2\:2.2.1--h9a82719_1 bwa-mem2 index -p Pcylind /path/GCA_964035525.1_jaPorCyli1.1_genomic.fna

ls path | grep "^P*" | grep -v "Won"|grep -v "Nbus"|grep -v "Mup"|grep -v "Che"|grep -v "^A"|cut -f 1 -d "."|uniq > list.txt

mkdir -p sam

for i in `cat list.txt`
do

singularity exec /usr/local/biotools/b/bwa-mem2\:2.2.1--h9a82719_1 bwa-mem2 mem -t 8 Pcylind -R "@RG\tID:${i}\tSM:${i}\tPL:ILLUMINA\tLB:${i}" /path/"$i".1.fq.gz /path/"$i".2.fq.gz > sam/"$i".sam

singularity exec /usr/local/biotools/s/samtools\:1.9--h46bd0b3_0 samtools sort -O bam -o "$i".bam sam/"$i".sam
singularity exec /usr/local/biotools/s/samtools\:1.9--h46bd0b3_0 samtools index "$i".bam

done

mkdir -p bam
mv *.bam bam/
mv *.bai bam/

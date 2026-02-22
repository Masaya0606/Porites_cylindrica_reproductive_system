#!/bin/bash
##SBATCH --job-name=hetroIso
#SBATCH --time=60:00:00
#SBATCH --mem=20G
#SBATCH --ntasks=10
#SBATCH --partition=epyc

singularity exec /usr/local/biotools/b/bcftools\:1.9--ha228f0b_4 bcftools view -s Pc24_52,Pc24_78,Pc24_30,Pc24_34,Pc24_53,Pc24_72 PoritesPcyl2503_3_rm57_Pcyl_24_40_PcylYellow43_ed_geno005_maf005_hwe.vcf -Oz -o  Clone_group_2_M2_MH3_F1.vcf.gz
singularity exec /usr/local/biotools/b/bcftools\:1.9--ha228f0b_4 bcftools view -s Pc24_16,Pc24_19,Pc24_21,Pc24_29,Pc24_18,Pc24_22,Pc24_14,Pc24_55 PoritesPcyl2503_3_rm57_Pcyl_24_40_PcylYellow43_ed_geno005_maf005_hwe.vcf -Oz -o  Clone_group_7_M4_MH2_F2.vcf.gz
#index
singularity exec /usr/local/biotools/b/bcftools\:1.9--ha228f0b_4 bcftools index Clone_group_2_M2_MH3_F1.vcf.gz
singularity exec /usr/local/biotools/b/bcftools\:1.9--ha228f0b_4 bcftools index Clone_group_7_M4_MH2_F2.vcf.gz
#Isolation of the hetero GT
singularity exec /usr/local/biotools/b/bcftools\:1.9--ha228f0b_4 bcftools query -f '%CHROM\t%POS\t[%GT\t]\n'  Clone_group_2_M2_MH3_F1.vcf.gz >Clone_group_2_M2_MH3_F1_GT.txt
singularity exec /usr/local/biotools/b/bcftools\:1.9--ha228f0b_4 bcftools query -f '%CHROM\t%POS\t[%GT\t]\n'  Clone_group_7_M4_MH2_F2.vcf.gz >Clone_group_7_M4_MH2_F2_GT.txt

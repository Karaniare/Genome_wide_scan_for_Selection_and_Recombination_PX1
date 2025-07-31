#!/bin/bash
#SBATCH -J LDhat_batch1_run6
#SBATCH -t 120:00:00
#SBATCH -c 12
#SBATCH --mem-per-cpu=8g



#Importing ldhat package
cd $main_dir
PATH=$PATH:/nfs/jbailey5/baileyweb/kniare/Uganda_USCF_WGS/Uganda_WGS/Resolved_Ug_WGS_ART_resistance_study/Downstream_analysis/Recombination_analysis/LDhat


lkgen -nseq 192 -lk LDhat/lk_files/lk_n192_t0.001 -prefix tanzania_for_ldhat

for i in 7
	#8 9 10 11 12 13 14
do
	vcftools --gzvcf $main_dir/tanzania.snps.vcf.gz --out tanzania_chr"$i"_for_ldhat --ldhat-geno --chr Pf3D7_0"$i"_v3
        interval -seq tanzania_chr"$i"_for_ldhat.ldhat.sites -loc tanzania_chr"$i"_for_ldhat.ldhat.locs -its 10000000 -bpen 5 -samp 2000 -prefix tanzania_chr"$i" -lk tanzania_for_ldhatnew_lk.txt - -burn 100000
cd $main_dir/LDhat
      	./stat -input $main_dir/tanzania_chr"$i"rates.txt -loc $main_dir/tanzania_chr"$i"_for_ldhat.ldhat.locs -prefix $main_dir/tanzania_chr"$i"_summary
done



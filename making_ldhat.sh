#!/bin/bash
#SBATCH -J LDhat_batch1_run6
#SBATCH -t 120:00:00
#SBATCH -c 12
#SBATCH --mem-per-cpu=8g

### Specifying analysis directory
main_dir="/nfs/jbailey5/baileyweb/kniare/Uganda_USCF_WGS/Uganda_WGS/Resolved_Ug_WGS_ART_resistance_study/Downstream_analysis/Recombination_analysis"

cd $main_dir
 wget https://github.com/auton1/LDhat/archive/refs/heads/master.zip ./
 unzip master.zip

#########Importing ldhat package

PATH=$PATH:/nfs/jbailey5/baileyweb/kniare/Uganda_USCF_WGS/Uganda_WGS/Resolved_Ug_WGS_ART_resistance_study/Downstream_analysis/Recombination_analysis/LDhat

############# Generating likelihood lookup tables
lkgen -nseq 192 -lk LDhat/lk_files/lk_n192_t0.001 -prefix tanzania_for_ldhat

############ Computing per locus recombination rates
for i in 7
	#8 9 10 11 12 13 14
do
              # Converting VCF file into LDhat file formats
			  
	vcftools --gzvcf $main_dir/tanzania.snps.vcf.gz --out tanzania_chr"$i"_for_ldhat --ldhat-geno --chr Pf3D7_0"$i"_v3
           
		   # Now computing recombination rates 
		interval -seq tanzania_chr"$i"_for_ldhat.ldhat.sites -loc tanzania_chr"$i"_for_ldhat.ldhat.locs -its 10000000 -bpen 5 -samp 2000 -prefix tanzania_chr"$i" -lk tanzania_for_ldhatnew_lk.txt - -burn 100000
cd $main_dir/LDhat

######## Calculating recombinate rate summary statistics
      	./stat -input $main_dir/tanzania_chr"$i"rates.txt -loc $main_dir/tanzania_chr"$i"_for_ldhat.ldhat.locs -prefix $main_dir/tanzania_chr"$i"_summary
done



**Performing genome-wide scans for signatures of selection and recombnation**

***Features***
This code is condensed into two files. Each file contains step-by-step instructions on how to run analyses. The first file (PX1_Manuscript_github.R) contains R scripts that enable two independent analyses of selection signals in the genome of _Plasmodium falciparum_ using different measures.
The first analysis generates [_IsoRelate_]([url](https://github.com/bahlolab/isoRelate)) iR statistics which is based on identity-by-descent (IBD) sharing
The second analysis computes Integrated Haplotype Score (iHS), which is based on the extended haplotype homozygosity, using the [_rehh_]([url](https://cran.r-project.org/web/packages/rehh/index.html)) package.
The second file (making_ldhat.sh) is a batch script calculating the recombinate rates across chromosomes using the [LDhat]([url](https://github.com/auton1/LDhat/blob/master/lkgen.c)) package.
Each file contains step-by-step instructions on how to run each analysis.

***Input**
- The input is a filtered genome-wide VCF file (filtered for genotype missigness and minor allele frequency).

***How to install the package***

git clone git@github.com:Karaniare/Genome_wide_scan_for_Selection_and_Recombination_PX1.git
tar -xvf Genome_wide_scan_for_Selection_and_Recombination_PX1.git

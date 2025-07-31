#### Relevant packages to load

{
  library(isoRelate)
  library(radiator)
  library(GENESIS)
  library(randomForest)
  library(GWASdata)
  library(ggtext)
  library(scales)
  library(vcfR)
  library(McCOILR)
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  library(vroom)
  library(ape)
  library(hierfstat)
  library(adegenet)
  library(plyr)
  library(factoextra)
  library(FactoMineR)
  library(Rtsne)
  library(igraph)
  library(reshape2)
  library(rehh)
  library(dplyr)
  library(tidyverse)
  library(vroom)
  library(RColorBrewer)
  library(pROC)
  library(tidyverse)
  library(rstatix)
  library(circlize)
  library(viridis)
  library(ggtree)
  library(Biostrings)
  library(ggtext)
}

#### Defining some custom colors

  my_colors<-c("#0000FF", "#800080", "#006400","#FF69B4", "#7CFC00", "#000000",
               "#DAA520","#00BFFF","#FF4500", "#BC8F8F", "#8B0000", "#00FF7F",
               "#E0FFFF","#1E90FF")
  blind_safe7<-c("blue","skyblue","#b8227b","#088F8F","#E34234","orange","bluishgreen")
  blind_safe8<-c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  blind_safe15<-c("#000000","#004949","#009292","#ff6db6","#ffb6db",
                  "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
                  "#920000","red","#db6d00","#24ff24","#ffff6d")
  
  blind_safe15b<-c("#000000","#1AFF1A","#009292","#ff6db6","#ffb6db",
                   "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
                   "#920000","red","#db6d00","#24ff24","#ffff6d")
  
  



############ Processing VCF for ISORELATE analysis

    ### Reading in VCF and extracting important data
hib<-"AllChrs.pass.merged.annotated.updated.snps.miss10.vaf.norm.north.east.final.vcf.gz"
hib_header <- seqVCF_Header(hib)
hib_header$format$Number[hib_header$ID == "AD"] <- "."
info.import <- c("AC", "AF", "AN","RO","AO","MQM","MQMR","DPB", "DP", "DS",
                 "PAO", "PQA", "PQR", "MQ",
                 "PRO", "QD", "SAF","SAR", "SAP", "SRP", "DPRA")
format.import <- c("AD","AO","DP","GL","GQ","GT","MIN_DP","PL","QA","QR","RO")

###  Converting VCF to the GDS format

hib_gds<-seqVCF2GDS(hib,"north_east.gds",
                    header=hib_header, info.import=info.import,
                    fmt.import=format.import)
seqSummary(hib_gds)

mip_file <- seqOpen(hib_gds) 

seqSummary(mip_file)

sample.id <- seqGetData(mip_file, "sample.id")
coords <- moimix::getCoordinates(mip_file)

head(coords)

### Creating and saving map and bed files

ped_map_cor<-moimix::extractPED(mip_file, moi.estimates = NULL, use.hets = T,
                        outfile = NULL)  

saveRDS(ped_map_cor, file = "north_east_selected.rds")
full_maped<-readRDS("north_east_selected.rds")

 


######################################## Running ISORELATE Aanalyses #############################
#Making space for big data processing
library(doSNOW)
registerDoSNOW(makeCluster(8, type = "SOCK"))

stopCluster(makeCluster(8, type = "SOCK"))

## Extracting genotypes from the reformatted WGS data

full_maped<-readRDS("north_east_selected.rds")
full_genotypes <- isoRelate::getGenotypes(ped.map = full_maped,
                                          reference.ped.map = NULL,
                                          maf = 0.02,
                                          isolate.max.missing = 0.6,
                                          snp.max.missing = 0.2,
                                          chromosomes = NULL,
                                          input.map.distance = "cM",
                                          reference.map.distance = "cM")

## Running IBD analysis

st_parameters <- getIBDparameters(ped.genotypes = full_genotypes, 
                                  number.cores = 12)
write_tsv(x=st_parameters, file = "parameters_north_east.tsv")
st_parameters<-read.table("parameters_north_east.tsv", header = T)
head(st_parameters)

##Detecting IBD Segments
st_ibd<- getIBDsegments(ped.genotypes = full_genotypes,
                        parameters = st_parameters, 
                        number.cores = 12, 
                        minimum.snps = 20, 
                        minimum.length.bp = 50000,
                        error = 0.001)
write_tsv(x=st_ibd, file="ibd_north_east_selected.tsv")
st_ibd<-read.table("ibd_north_east_selected.tsv",header =T)
head(st_ibd)


###### Making population definition files

gp<-read.table("px1_genotypes_metadata_combined.tsv", header = T, sep = "\t")
gp$major_genotype[gp$major_genotype=="A675V"]<-"Mutant"
gp$major_genotype[gp$major_genotype=="C469Y"]<-"Mutant"

colnames(gp)[2]<-"fid"
my_groups <- full_genotypes[[1]][,1:3]

gp1<-merge(x=my_groups,y=gp, by = "fid", all = F)
my_groups<-data.frame(fid=gp1$fid,iid=gp1$iid,pid=gp1$region)

# generating a binary IBD matrix
st_matrix <- getIBDmatrix(ped.genotypes = full_genotypes, 
                          ibd.segments = st_ibd)

# calculating pairwise IBD fractions for each SNP

ibd_prop_raw <- getIBDproportion(ped.genotypes = full_genotypes, ### without stratification
                                 ibd.matrix = st_matrix, 
                                 groups = NULL)

st_proportion <- getIBDproportion(ped.genotypes = full_genotypes, ### with stratification by K13 mutations
                                  ibd.matrix = st_matrix, 
                                  groups = my_groups)
st_proportion0 <- getIBDproportion(ped.genotypes = full_genotypes, 
                                   ibd.matrix = st_matrix, 
                                   groups = my_groups)

####### Compoutings Isolate's iR metric of selection

st_iR <- getIBDiR(ped.genotypes = full_genotypes, ### without stratification
                  ibd.matrix = st_matrix, 
                  groups = NULL)
st_iR_grp <- getIBDiR(ped.genotypes = full_genotypes, ### with stratification by K13 mutations
                      ibd.matrix = st_matrix, 
                      groups = my_groups)

###### Reformatting IBD results
{
  ###
  st_iR$chr[st_iR$chr=="Pf3D7_01_v3"]<-1
  st_iR$chr[st_iR$chr=="Pf3D7_02_v3"]<-2
  st_iR$chr[st_iR$chr=="Pf3D7_03_v3"]<-3
  st_iR$chr[st_iR$chr=="Pf3D7_04_v3"]<-4
  st_iR$chr[st_iR$chr=="Pf3D7_05_v3"]<-5
  st_iR$chr[st_iR$chr=="Pf3D7_06_v3"]<-6
  st_iR$chr[st_iR$chr=="Pf3D7_07_v3"]<-7
  st_iR$chr[st_iR$chr=="Pf3D7_08_v3"]<-8
  st_iR$chr[st_iR$chr=="Pf3D7_09_v3"]<-9
  st_iR$chr[st_iR$chr=="Pf3D7_10_v3"]<-10
  st_iR$chr[st_iR$chr=="Pf3D7_11_v3"]<-11
  st_iR$chr[st_iR$chr=="Pf3D7_12_v3"]<-12
  st_iR$chr[st_iR$chr=="Pf3D7_13_v3"]<-13
  st_iR$chr[st_iR$chr=="Pf3D7_14_v3"]<-14
}

{
  ###
  st_iR_grp$chr[st_iR_grp$chr=="Pf3D7_01_v3"]<-1
  st_iR_grp$chr[st_iR_grp$chr=="Pf3D7_02_v3"]<-2
  st_iR_grp$chr[st_iR_grp$chr=="Pf3D7_03_v3"]<-3
  st_iR_grp$chr[st_iR_grp$chr=="Pf3D7_04_v3"]<-4
  st_iR_grp$chr[st_iR_grp$chr=="Pf3D7_05_v3"]<-5
  st_iR_grp$chr[st_iR_grp$chr=="Pf3D7_06_v3"]<-6
  st_iR_grp$chr[st_iR_grp$chr=="Pf3D7_07_v3"]<-7
  st_iR_grp$chr[st_iR_grp$chr=="Pf3D7_08_v3"]<-8
  st_iR_grp$chr[st_iR_grp$chr=="Pf3D7_09_v3"]<-9
  st_iR_grp$chr[st_iR_grp$chr=="Pf3D7_10_v3"]<-10
  st_iR_grp$chr[st_iR_grp$chr=="Pf3D7_11_v3"]<-11
  st_iR_grp$chr[st_iR_grp$chr=="Pf3D7_12_v3"]<-12
  st_iR_grp$chr[st_iR_grp$chr=="Pf3D7_13_v3"]<-13
  st_iR_grp$chr[st_iR_grp$chr=="Pf3D7_14_v3"]<-14
}

write_tsv(x=st_iR, file = "IR_allsnps_north_east_maf2.csv")
write_tsv(x=st_iR_grp, file = "IR_allsnps_north_east_stratification_maf2.csv")





{
  ###
  ibd_prop_raw$chr[ibd_prop_raw$chr=="Pf3D7_01_v3"]<-1
  ibd_prop_raw$chr[ibd_prop_raw$chr=="Pf3D7_02_v3"]<-2
  ibd_prop_raw$chr[ibd_prop_raw$chr=="Pf3D7_03_v3"]<-3
  ibd_prop_raw$chr[ibd_prop_raw$chr=="Pf3D7_04_v3"]<-4
  ibd_prop_raw$chr[ibd_prop_raw$chr=="Pf3D7_05_v3"]<-5
  ibd_prop_raw$chr[ibd_prop_raw$chr=="Pf3D7_06_v3"]<-6
  ibd_prop_raw$chr[ibd_prop_raw$chr=="Pf3D7_07_v3"]<-7
  ibd_prop_raw$chr[ibd_prop_raw$chr=="Pf3D7_08_v3"]<-8
  ibd_prop_raw$chr[ibd_prop_raw$chr=="Pf3D7_09_v3"]<-9
  ibd_prop_raw$chr[ibd_prop_raw$chr=="Pf3D7_10_v3"]<-10
  ibd_prop_raw$chr[ibd_prop_raw$chr=="Pf3D7_11_v3"]<-11
  ibd_prop_raw$chr[ibd_prop_raw$chr=="Pf3D7_12_v3"]<-12
  ibd_prop_raw$chr[ibd_prop_raw$chr=="Pf3D7_13_v3"]<-13
  ibd_prop_raw$chr[ibd_prop_raw$chr=="Pf3D7_14_v3"]<-14
}


write_tsv(x=st_iR, file = "IR_allsnps_north_east_maf2.csv")

############################################ Genome-wide analysis for iHS signals of selection ##################################################################
### loading data from VCF and scanning for EHH
ug_ihh_all<-data.frame()
chromo_list<-c("Pf3D7_01_v3","Pf3D7_02_v3","Pf3D7_03_v3","Pf3D7_04_v3","Pf3D7_05_v3","Pf3D7_06_v3",
               "Pf3D7_07_v3","Pf3D7_08_v3","Pf3D7_09_v3","Pf3D7_10_v3","Pf3D7_11_v3","Pf3D7_12_v3","Pf3D7_13_v3","Pf3D7_14_v3")
chr_list<-c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14")
for (i in chr_list) {
  ug_hh<-data2haplohh(hap_file =   "AllChrs.pass.merged.annotated.maf1.updated.major.cleaned.snps.mono.north.east.vcf.gz", min_perc_geno.mrk = 90,
                      vcf_reader = "data.table", verbose = T, polarize_vcf = F, min_maf=0.02,remove_multiple_markers=T,chr.name = i) 
  ug_ihh <- scan_hh(ug_hh,discard_integration_at_border = T,threads = 24,limehhs = 0.1, limhomohaplo = 4, polarized = F,maxgap = 2000000)
  ug_ihh<-data.frame(ug_ihh)
  ug_ihh_all <- dplyr::bind_rows(ug_ihh_all,ug_ihh)
 
}

## reformatting chromosome names in the iHS results

{ug_ihh_all$CHR[ug_ihh_all$CHR=="Pf3D7_01_v3"]<-1
  ug_ihh_all$CHR[ug_ihh_all$CHR=="Pf3D7_02_v3"]<-2
  ug_ihh_all$CHR[ug_ihh_all$CHR=="Pf3D7_03_v3"]<-3
  ug_ihh_all$CHR[ug_ihh_all$CHR=="Pf3D7_04_v3"]<-4
  ug_ihh_all$CHR[ug_ihh_all$CHR=="Pf3D7_05_v3"]<-5
  ug_ihh_all$CHR[ug_ihh_all$CHR=="Pf3D7_06_v3"]<-6
  ug_ihh_all$CHR[ug_ihh_all$CHR=="Pf3D7_07_v3"]<-7
  ug_ihh_all$CHR[ug_ihh_all$CHR=="Pf3D7_08_v3"]<-8
  ug_ihh_all$CHR[ug_ihh_all$CHR=="Pf3D7_09_v3"]<-9
  ug_ihh_all$CHR[ug_ihh_all$CHR=="Pf3D7_10_v3"]<-10
  ug_ihh_all$CHR[ug_ihh_all$CHR=="Pf3D7_11_v3"]<-11
  ug_ihh_all$CHR[ug_ihh_all$CHR=="Pf3D7_12_v3"]<-12
  ug_ihh_all$CHR[ug_ihh_all$CHR=="Pf3D7_13_v3"]<-13
  ug_ihh_all$CHR[ug_ihh_all$CHR=="Pf3D7_14_v3"]<-14
}

{ug_ihh_all$CHR[ug_ihh_all$CHR=="chr1"]<-1
  ug_ihh_all$CHR[ug_ihh_all$CHR=="chr2"]<-2
  ug_ihh_all$CHR[ug_ihh_all$CHR=="chr3"]<-3
  ug_ihh_all$CHR[ug_ihh_all$CHR=="chr4"]<-4
  ug_ihh_all$CHR[ug_ihh_all$CHR=="chr5"]<-5
  ug_ihh_all$CHR[ug_ihh_all$CHR=="chr6"]<-6
  ug_ihh_all$CHR[ug_ihh_all$CHR=="chr7"]<-7
  ug_ihh_all$CHR[ug_ihh_all$CHR=="chr8"]<-8
  ug_ihh_all$CHR[ug_ihh_all$CHR=="chr9"]<-9
  ug_ihh_all$CHR[ug_ihh_all$CHR=="chr10"]<-10
  ug_ihh_all$CHR[ug_ihh_all$CHR=="chr11"]<-11
  ug_ihh_all$CHR[ug_ihh_all$CHR=="chr12"]<-12
  ug_ihh_all$CHR[ug_ihh_all$CHR=="chr13"]<-13
  ug_ihh_all$CHR[ug_ihh_all$CHR=="chr14"]<-14
}

#### Calculating iHS and and saving final results
ug_ihs_all<-ihh2ihs(ug_ihh_all, freqbin = 50, verbose = F, standardize = T, min_nhaplo = 4)

write_tsv(x=ug_ihs_all$ihs, file = "ihs_Uganda_allsnps_north_east_raw.csv")

#####

######################################## Anlyzing EHH decay around drug resistance markers (here PX1)-- ################################################################################

### Reading data from VCFs (SNPs only and SNPs +indels)
px1<-data2haplohh(hap_file = "AllChrs.pass.merged.annotated.maf1.updated.major.cleaned.snps.vcf.gz", min_perc_geno.mrk = 90,
                  vcf_reader = "data.table", verbose = T, polarize_vcf = F, min_maf=0.05,remove_multiple_markers=T,chr.name = "chr7")

px1_indel<-data2haplohh(hap_file = "AllChrs.pass.merged.annotated.maf1.updated.major.cleaned.vcf.gz", min_perc_geno.mrk = 90,
                        vcf_reader = "data.table", verbose = T, polarize_vcf = F, min_maf=0.05,remove_multiple_markers=T,chr.name = "chr7")

#### only px1 gene signals
{   px1_res1<-calc_ehh(px1, mrk = "chr7:897322:897322:G:A", include_nhaplo = T,include_zero_values = T)
  anc<-data.frame(position=px1_res1$ehh$POSITION,ehh=px1_res1$ehh$EHH_A,Haplotype="WT",mutation="D1705N")
  derv1<-data.frame(position=px1_res1$ehh$POSITION,ehh=px1_res1$ehh$EHH_D,Haplotype="Mutant",mutation="D1705N")
  
  dat_ehh1<-rbind(anc,derv1)
  dat_ehh1<-dat_ehh1 %>%
    mutate(marker_position=897322)
  
  
  px1_res1<-calc_ehh(px1, mrk = "chr7:897312:897312:G:A", include_nhaplo = T,include_zero_values = T)
  anc<-data.frame(position=px1_res1$ehh$POSITION,ehh=px1_res1$ehh$EHH_A,Haplotype="WT",mutation="M1701I")
  derv2<-data.frame(position=px1_res1$ehh$POSITION,ehh=px1_res1$ehh$EHH_D,Haplotype="Mutant",mutation="M1701I")
  
  dat_ehh2<-rbind(anc,derv2)
  dat_ehh2<-dat_ehh2 %>%
    mutate(marker_position=897312)
  
  
  px1_res1<-calc_ehh(px1, mrk = "chr7:895874:895874:T:C", include_nhaplo = T,include_zero_values = T)
  anc<-data.frame(position=px1_res1$ehh$POSITION,ehh=px1_res1$ehh$EHH_A,Haplotype="WT",mutation="L1222P")
  derv3<-data.frame(position=px1_res1$ehh$POSITION,ehh=px1_res1$ehh$EHH_D,Haplotype="Mutant",mutation="L1222P")
  
  dat_ehh3<-rbind(anc,derv3)
  dat_ehh3<-dat_ehh3 %>%
    mutate(marker_position=895874)
  
  px1_res1<-calc_ehh(px1_indel, mrk = "chr7:897246:897246:CGTGGATAATATGTATAAT:C", include_nhaplo = T,include_zero_values = T)
  anc<-data.frame(position=px1_res1$ehh$POSITION,ehh=px1_res1$ehh$EHH_A,Haplotype="WT",mutation="del(V1680-1685)")
  derv3<-data.frame(position=px1_res1$ehh$POSITION,ehh=px1_res1$ehh$EHH_D,Haplotype="Mutant",mutation="del(V1680-1685)")
  
  dat_ehh4<-rbind(anc,derv3)
  dat_ehh4<-dat_ehh4 %>%
    mutate(marker_position=897246)
  
  px1_res1<-calc_ehh(px1_indel, mrk = "chr7:896071:896071:GTAAATAATGTGAACAATA:G", include_nhaplo = T,include_zero_values = T)
  anc<-data.frame(position=px1_res1$ehh$POSITION,ehh=px1_res1$ehh$EHH_A,Haplotype="WT",mutation="del(N1289-N1294)")
  derv3<-data.frame(position=px1_res1$ehh$POSITION,ehh=px1_res1$ehh$EHH_D,Haplotype="Mutant",mutation="del(N1289-N1294)")
  
  dat_ehh5<-rbind(anc,derv3)
  dat_ehh5<-dat_ehh5 %>%
    mutate(marker_position=896071)
  
  
  px1_res1<-calc_ehh(px1_indel, mrk = "chr7:894637:894637:GATAATGATGATAATAATTATAATGATGATAATAATT:G;chr7:894637:894637:G:GATAATGATGATAATAATT;chr7:894637:894637:GATAATGATGATAATAATT:G", include_nhaplo = T,include_zero_values = T)
  anc<-data.frame(position=px1_res1$ehh$POSITION,ehh=px1_res1$ehh$EHH_A,Haplotype="WT",mutation="del(N811-Y822)")
  derv3<-data.frame(position=px1_res1$ehh$POSITION,ehh=px1_res1$ehh$EHH_D1,Haplotype="Mutant",mutation="del(N811-Y822)")
  
  dat_ehh6<-rbind(anc,derv3)
  dat_ehh6<-dat_ehh6 %>%
    mutate(marker_position=894637)
  
  dat_ehh_px1<-rbind(dat_ehh1,dat_ehh2,dat_ehh3,dat_ehh4,dat_ehh6) 
  
  ## dat_ehh_px1 is final results
  
}




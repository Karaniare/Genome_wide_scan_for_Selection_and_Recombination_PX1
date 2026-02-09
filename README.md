# Genome-wide Scans for Selection and Recombination in *Plasmodium falciparum*

This repository provides scripts to perform **genome-wide scans for signatures of selection and recombination** in *Plasmodium falciparum* using population genomic data.

The workflows implemented here were developed in the context of the **PX1 manuscript**, but they are generalizable and can be applied to other datasets with minimal adaptation.

---

## Features

This repository contains **two main analysis pipelines**, implemented in two files.  
Each file includes **step-by-step instructions** to guide users through the analyses.

---

### 1. Genome-wide Selection Scans  
**`PX1_Manuscript_github.R`**

This R script performs **two independent and complementary tests for positive selection**:

#### IBD-based selection scan
- Computes *iR statistics* using  
  👉 [IsoRelate](https://github.com/bahlolab/isoRelate)
- Based on excess identity-by-descent (IBD) sharing across the genome

#### Haplotype-based selection scan
- Computes **Integrated Haplotype Score (iHS)**
- Based on extended haplotype homozygosity (EHH)
- Implemented using the  
  👉 [`rehh`](https://cran.r-project.org/web/packages/rehh/index.html) R package

The script covers the full workflow, from input processing to genome-wide visualization of selection signals.

---

### 2. Genome-wide Recombination Rate Estimation  
**`making_ldhat.sh`**

This bash script estimates **recombination rates across chromosomes** using **LDhat**:

- Chromosome-wise processing
- Estimation of population-scaled recombination rates (ρ)
- Designed for large-scale, genome-wide analyses

Inline comments provide clear guidance on how to run and adapt the pipeline.

---

## Input Data

- A **genome-wide VCF file**
- The VCF should be filtered for:
  - Genotype missingness
  - Minor allele frequency (MAF)

---

## Installation

Clone the repository:

```bash
git clone https://github.com/Karaniare/Genome_wide_scan_for_Selection_and_Recombination_PX1.git
cd Genome_wide_scan_for_Selection_and_Recombination_PX1```
```
## Resources

- **Example data for selection analysis**  
  *(link to be added](https://github.com/Karaniare/Genome_wide_scan_for_Selection_and_Recombination_PX1/tree/main/Examples)*

- **Example data for recombination rate estimation**  
  *(link to be added](https://github.com/Karaniare/Genome_wide_scan_for_Selection_and_Recombination_PX1/tree/main/Examples)*

- **LDhat documentation and resources**
  - GitHub repository: https://github.com/auton1/LDhat
  - Includes source code, likelihood tables, and usage instructions

---

## Dependencies

### R Packages

The following R packages are required for the selection analyses:

- `isoRelate`
- `moimix`
- `vcfR`
- `ggplot2`
- `ggsci`
- `ggpubr`
- `rehh`
- `dplyr`
- `SeqArray`
- `R.utils`

### External Software

- **LDhat**  
  Required for genome-wide recombination rate estimation  
  - https://github.com/auton1/LDhat

---

## Notes

- The pipelines were developed for *Plasmodium falciparum* but can be adapted to other organisms.
- Users are encouraged to inspect and adjust filtering thresholds depending on sample size and study design.

---

## Citation

If you use this code, please cite the **Niare K, Tafesse B, Treat M, Sadler J, Okitwi M, Orena S, Asua V, Kreutzfeld O, Legac J, Samuel NL, Yeka A, Rosenthal PJ, Juliano JJ, Bailey JA, Conrad MD. A novel locus associated with decreased susceptibility of Plasmodium falciparum to lumefantrine and dihydroartemisinin has emerged and spread in Uganda. bioRxiv [Preprint]. 2025 Aug 2:2025.07.30.667738. doi: 10.1101/2025.07.30.667738. PMID: 40766679; PMCID: PMC12324449** and the relevant software packages (**IsoRelate**, **rehh**, and **LDhat**).




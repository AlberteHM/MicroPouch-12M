# MicroPouch - 1-year follow-up
## Abstract

**Background**

Chronic pouchitis is a frequent complication following ileal pouch-anal anastomosis surgery, with microbial dysbiosis considered a key contributor to its pathogenesis. Fecal Microbiota Transplantation (FMT) has shown inconsistent therapeutic effects, and its long-term efficacy remains unclear. In this follow-up study of the MicroPouch trial, we assessed 1-year microbial changes and clinical outcomes of FMT in patients with chronic pouchitis.

**Methods**

Thirty patients were randomized 1:1 to receive FMT or placebo. Clinical outcomes including stool frequency and clinical Pouch Disease Activity Index (cPDAI) scores were assessed longitudinally during the 1-year follow-ups. Concurrent fecal samples were collected and subjected to metagenomic short read sequencing for gut microbiome profiling.  

**Results**

The results showed, no decrease in stool frequency in the FMT group at any follow-up, while a decrease was observed in the placebo group up to 3 months after treatment. Both groups demonstrated reductions in cPDAI scores; these improvements persisted up to 6 months in the FMT group and 12 months in the placebo group. However, no significant differences between groups were observed at the 1-year follow-up, except for just after the treatment period. Microbiome analyses revealed no sustained shift toward fecal donor composition after FMT, and microbial profiles at 1 year were comparable to baseline in both groups.

**Conclusions**

Overall, FMT treatment did not provide sustained clinical or microbial benefits in patients with chronic pouchitis up to one year after treatment. However, as the study was underpowered, the results should be considered hypothesis-generating. 

## Analysis
This repository contains the code used for clinical data plot generation and the microbiome analysis.

`data` contains: 
 - `data/MetaPhlAn_4.1.0.txt`: MetaPhlAn taxonomic profile generated from short-read metagenomic data. </li>          
 - `data/metadata.csv`: Contains the metadata for the included participants. </li>
 - `data/within_group_cpdai_stats.tsv`: Statistical data for within-group timepoint comparisons of cPDAI. </li>          
 - `data/within_group_stoolfre_stats.tsv`: Statistical data for within-group timepoint comparisons of stool frequency. </li>

`analysis` contains: 
 - [`plot_clinical_results.R`](analysis/plot_clinical_results.R): script used to generate Figure 2 and Figure 3. 
 - [`plot_microbiome_results.R`](analysis/plot_microbiome_results.R): script used to generate Figure 4, Figure 5, and Figure 6.

## Prerequisites
All data was analyzed using R (4.1.0) and RStudio.

## Data
All DNA sequences have been deposited at the European Nucleotide Archive under the accession numbers PRJEB66493 and PRJEB98069. 

## Cite 
Please cite our article:



# MicroPouch - 1-year follow-up
## Abstract
Chronic pouchitis is a frequent complication after ileal pouch–anal anastomosis surgery, and microbial dysbiosis is believed to play a central role in the pathogenesis of pouchitis. Fecal microbiota transplantation (FMT) has shown varying results in patients suffering from pouchitis, but long-term effects remain unclear. This study is a follow-up of the Micropouch trial, evaluating the clinical and microbial impact of FMT in patients with chronic pouchitis during a one-year period.

Thirty patients with a pouch were included in the original study and randomized 1:1 to receive FMT or placebo. Patients treated with FMT exhibited reduced stool frequency for up to 1 month, while the placebo group showed reduced stool frequency at the 3-month follow-up. Both groups demonstrated significant early improvements in cPDAI scores. The FMT group had significantly lower scores  up to the 6 month follow-up. In the placebo group, improvements persisted up to the 3 month follow-up. Microbiome analyses revealed no long-term shifts toward the donor profile following FMT. No significant microbial changes were observed in the placebo group.

Overall, FMT did not result in sustained clinical improvement, and the microbial composition did not differ from baseline after 1 year, in patients with chronic pouchitis.

## Analysis
This repository contains the code used for clinical data plot generation and the microbiome analysis.

`data` contains: 
 - `data/MetaPhlAn_4.1.0.txt`: MetaPhlAn taxonomic profile generated from short-read metagenomic data. </li>          
 - `data/metadata.csv`: Contains the metadata for the included participants. </li>

`analysis` contains: 
 - [`plot_clinical_results.R`](analysis/plot_clinical_results.R): script used to generate Figure 2 and Figure 3. 
 - [`plot_microbiome_results.R`](analysis/plot_microbiome_results.R): script used to generate Figure 4, Figure 5, and Figure 6.

## Prerequisites
All data was analyzed using R (4.1.0) and RStudio.

## Data
All DNA sequences have been deposited at the European Nucleotide Archive under the accession numbers PRJEB66493 and PRJEB98069. 

## Cite 
Please cite our article:



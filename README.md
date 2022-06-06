# ST372
Scripts to replicate analysis performed on ST372 *E. coli* associated with the manuscript "Identification of genes influencing the evolution of *Escherichia coli* ST372 in dogs and humans" by Elankumaran *et al*, 2022.

## Overview
The following repository allows conscientious readers of the manuscript to reproduce all the data processing, statistics and figures presented in the paper.

It comprises two directories: __`scripts`__ and __`data`__ (the contents of which should be self explanatory) and generates a __`data/outputs`__ folder with __`figures`__ and __`data`__ subdirectories.


## Installation
### Software requirements
These scripts are currently functional on mac OS Monterey 12.3.1 using RStudio 2022.02.2 (Build 485) and R version 4.1.3. We cannot guarantee they will work on other distributions of R or RStudio. Your OS should not be an issue provided you use these versions of R and RStudio though.

### Packages required
You will need to install the following packages and versions to work with the scripts:
- paletteer_1.4.0    
- plasmidmapR_0.1.0  
- abricateR_0.1.1    
- RColorBrewer_1.1-3 
- ggtree_3.5.0.900
- pheatmap_1.0.12
- reshape2_1.4.4     
- ggpubr_0.4.0       
- ggplot2_3.3.6      
- tibble_3.1.7
- purrr_0.3.4        
- readr_2.1.2        
- stringr_1.4.0      
- forcats_0.5.1      
- tidyr_1.2.0
- dplyr_1.0.9

### Issues with ggtree
There have recently been some issues with ggtree in the way it interacts with dplyr. The solution is to install the latest version of ggtree directly from github instead of via BiocManager. You can do this in the console on RStudio with the __`remotes`__ package like so:
```
install.packages("remotes")
remotes::install_github("YuLab-SMU/ggtree")
```

### plasmidmapR and abricateR
Please see [plasmidmapR](https://github.com/maxlcummins/plasmidmapR) and [abricateR](https://github.com/maxlcummins/abricateR) repositories to install these packages.

## Usage
Clone this repository
```
git clone https://github.com/CJREID/ST372.git
cd ST372
pwd
```
Open the ST372_ANALYSIS.R script in a text editor or RStudio and set the variable __`wrkdir`__ on line 23 to the output of __`pwd`__ above and save the script.

Run the ST372_ANALYSIS.R script and your __`outputs`__ folder will populate with figures and tables. Please note all figure legends were edited manually after export from R except for Figure 1.

## Outputs
### Figures
1. Figure 1. Genomic epidemiological features of the ST372 genome collection
2. Figure 2. Core gene maximum-likelihood phylogenetic tree with metadata
3. Figure 3. Low pairwise SNP distance heatmap
4. Figure 4. Map of cluster-associated genes aligned to the core gene phylogeny
5. Figure 5. Schematic representations of gene loci containing cluster-associated genes and predicted genomic islands

### Supplementary
#### Tables
1. Table S1. Metadata, accession numbers and gene screening results for 407 ST372 *E. coli* genomes used in the study
2. Table S2. Scoary results of genes significantly associated with phylogenetic clusters
3. Table S3. BLASTn results of MVC121 *pdu* operon to the *pdu* operon originally described in *Salmonella*


#### Figures
1. Fig S1. Alignment of all sequences to *pdu* island observed in MVC107 mapped to core gene phylogeny
2. Fig S2. Alignment of all sequences to *pdu* island observed in MVC121/MVC18 mapped to core gene phylogeny
3. Fig S3. Alignment of all sequences to *kps* island observed in MVC18 mapped to core gene phylogeny
4. Fig S4. Alignment of all sequences to *kps* island observed in MVAST6839 mapped to core gene phylogeny


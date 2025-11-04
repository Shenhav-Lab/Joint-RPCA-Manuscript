# Joint-RPCA Data & Code Repository
---
Joint-RPCA is a novel multi-omics dimensionality reduction tool designed to integrate multiple omic types from matched samples, uncovering patterns that drive ecosystem-wide differences across phenotypes. Joint-RPCA is part of the [Gemelli](https://github.com/cameronmartino/gemelli/tree/master/gemelli) tool box.

Here we provide the data and scripts used to conduct the case studies, simulations, and benchmarking presented in the Joint-RPCA publication. 

## Data
---
The following datasets were used in case studies and data-driven simulations: 

- **iHMP:** The metatranscriptomics, metagenomics, metabolomics, proteomics, and viromics data were downloaded from the integrated Human Microbiome Project [database](https://ibdmdb.org). See [Lloyd-Price, J. et al.](https://www.nature.com/articles/s41586-019-1237-9) for more details.
- **Additional IBD cohorts:** Instructions for accessing the independent inflammatory bowel disease (IBD) cohort data can be found [here](https://github.com/knightlab-analyses/uc-severity-multiomics). More details on the cohorts are provided by [Mills, R. H. et al](https://www.nature.com/articles/s41564-021-01050-3).
- **American Gut Project (AGP):** Metagenomics from controls with no IBD or irritable bowel syndrome (IBS). Instructions to access data can be found in the following [repository](https://github.com/knightlab-analyses/american-gut-analyses).
or in the original publication (McDonald, D. et al. mSystems, 2018).
- **Biocrust soils:** The true positive metabolite associations and evaluation code were obtained from [Morton, J. T. et al](https://github.com/knightlab-analyses/multiomic-cooccurrences/tree/rebuttal). Also see [Swenson, T. L. et al.](https://www.nature.com/articles/s41467-017-02356-9) for study details.
- **Mammalian Safari:** All data modalities were obtained and preprocessed following the original publication from [Gregor, R. et al.](https://academic.oup.com/ismej/article/16/5/1262/7474189?login=true).
- **Decomposer Network:** All data modalities were obtained from [Burcham, Z. M. et al.](https://github.com/Metcalf-Lab/2023-Universal-microbial-decomposer-network/tree/master/jrpca_network/network)

## Code
---
## Conda environment installation

Gemelli is most easily used inside of a [QIIME2](https://qiime2.org/) environment. The directions for creating a QIIME2 environment can be found [here](https://docs.qiime2.org/2024.10/install/native/#install-qiime-2-within-a-conda-environment). All of the versions are stored in the dockerfile in this capsule. 

To reduce compute time, we provide intermediate files under the data directory and include only the scripts required to reproduce the figures in `run.sh`. Scripts used for data preprocessing as well as other time-intensive tasks, such as machine learning and running the benchmarked tools, are provided under the *code* directory.

## Scripts

**Simulations:** Scripts used for the data-driven simulations are stored under `code/simulations-benchmarking/`. Scripts are numbered to indicate the order in which they need to be run, as well as the data they are associated with. For example, all **2._** scripts incorporate the biocrust wetting data from Swenson et al. (2018).

**Case Studies:** Similarly, scripts used for the case studies are stored under `code/case-studies/` and are numbered to indicate run order as well as the case study they are relevant to. 

## Specific scripts used for generating figures and tables
---
The following list documents which analysis scripts were used to generate each figure and table included in the main text and supplementary materials.

### **Main Figures**

#### **Figure 1**
- `simulations-benchmarking/0.0-introduction-figure.ipynb`

#### **Figure 2**
- `case-studies/1.1-iHMP-analysis-plotting.ipynb`: **Figure 2A**  
- `simulations-benchmarking/3.6-ihmp-benchmarks-plotting.ipynb`: **Figure 2B**

#### **Figure 3**
- `simulations-benchmarking/2.4-biocrust-benchmarks-plotting.ipynb`: **Figure 3A–B**  
- `simulations-benchmarking/1.0-finrisk-runtime-plotting.ipynb`: **Figure 3C–D**  
- `case-studies/1.5.5-runtime-plotting.ipynb`: **Figure 3E**

#### **Figure 4**
- `case-studies/1.1-iHMP-analysis-plotting.ipynb`: **Figure 4A, 4F**  
- `case-studies/1.3-UCSD-analysis-plotting.ipynb`: **Figure 4B–D**  
- `case-studies/1.4-AGP-analysis-plotting.ipynb`: **Figure 4E**

#### **Figure 5**
- `case-studies/3.2-decomposer-scatterplots.ipynb`: **Figure 5A**  
- `case-studies/3.7-decomposer-regression-plotting.ipynb`: **Figure 5B**  
- `case-studies/3.5-decomposer-venn.ipynb`: **Figure 5C**  
- `case-studies/2.2-mammalian-analysis-plotting.ipynb`: **Figure 5D–F**

### **Supplementary Figures**

#### **From `case-studies/`**
- `1.1-iHMP-analysis-plotting.ipynb`: **Supp. Fig. 1 and 6**  
- `3.4-decomposer-network.ipynb`: **Supp. Fig. 7**  
- `3.2-decomposer-scatterplots.ipynb`: **Supp. Fig. 8**  
- `2.2-mammalian-analysis-plotting.ipynb`: **Supp. Fig. 9**

#### **From `simulations-benchmarking/`**
- `3.6-ihmp-benchmarks-plotting.ipynb`: **Supp. Fig. 2A–B**  
- `3.8-ihmp-feature-overlap-plotting.ipynb`: **Supp. Fig. 2C**  
- `1.7-lowrank-sims-plotting.ipynb`: **Supp. Fig. 3**  
- `1.1-breakpoint-simulation.ipynb`: **Supp. Fig. 4**  
- `2.4-biocrust-benchmarks-plotting.ipynb`: **Supp. Fig. 5**

### **Supplementary Tables**

#### **From `simulations-benchmarking/`**
- `3.6-ihmp-benchmarks-plotting.ipynb`: **Supp. Table 1**  
- `1.7-lowrank-sims-plotting.ipynb`: **Supp. Table 2**  
- `2.4-biocrust-benchmarks-plotting.ipynb`: **Supp. Table 3–4**  
- `1.0-finrisk-runtime-plotting.ipynb`: **Supp. Table 5**

#### **From `case-studies/`**
- `1.5.5-runtime-plotting.ipynb`: **Supp. Table 6**  
- `3.1-decomposer-rpca.ipynb`: **Supp. Table 7**


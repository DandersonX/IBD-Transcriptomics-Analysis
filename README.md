# Identifying Transcriptomic Drivers of Active Inflammatory Bowel Disease (IBD)

An end-to-end computational biology pipeline using R and `limma` to analyze high-throughput transcriptomic data. This study screens **1,832 clinical biopsy samples** across **25,970 genes** to evaluate structural and immunological changes during active IBD flares.

---

## Project Overview & Objective
The goal of this analysis was to identify statistically robust differentially expressed genes (DEGs) that characterize the transition from inactive mucosal tissue to active inflammation in patients with IBD.

---

## Technical Challenges & Data Engineering Flex
Publicly available transcriptomic datasets are notoriously messy. During this analysis, a primary indexing failure occurred because the supplementary count matrix used short sample headers (`MSCCR_reGRID_1_Biopsy_1`), while the GEO metadata rows utilized complex trailing text descriptions (`MSCCR_reGRID_1_Biopsy_1, UC participants, Rectum...`).

**Solution:** I engineered a custom string-cleaning step in R using regular expressions (`sub()`) and name standardization (`make.names()`) to truncate trailing characters and realign all 1,832 samples seamlessly, preventing index corruption and recovering 100% of the clinical data.

---

## Methods & Tech Stack
* **Language:** R (v4.3+)
* **Core Packages:** `limma` (linear modeling and empirical Bayes moderation), `GEOquery` (data acquisition)
* **Statistical Thresholds:** Benjamini-Hochberg Adjusted p-value $< 0.05$

---

## Results & Key Findings

### 1. Differential Expression Profile
Our linear model successfully isolated highly significant transcriptomic shifts. Below is the Volcano Plot highlighting significantly upregulated genes (red) and downregulated genes (blue) during active flares.

### 2. Top Monitored Portfolio Genes
The analysis revealed a strict biological profile marking the shift to an active disease state:

### 2. Top Monitored Portfolio Genes

The analysis revealed a strict biological profile marking the shift to an active disease state:

| Ensembl ID | Gene Symbol | LogFC | Adj. P-Value | Biological Function |
| :--- | :--- | :--- | :--- | :--- |
| `ENSG00000136770` | **IL4R** | +0.17 | ~1.1e-16 | Interleukin-4 Receptor; drives mucosal inflammatory signaling. |
| `ENSG00000123352` | **CD14** | +0.21 | ~3.1e-15 | Monocyte marker; recognizes bacterial lipopolysaccharides. |
| `ENSG00000103257` | **LILRB2** | +0.60 | ~1.0e-14 | Regulates macrophage activation in inflamed tissue. |
| `ENSG00000133067` | **LIMS1** | -0.48 | ~5.8e-14 | Downregulated structural element; indicates epithelial barrier breakdown. |
---

## How to Run This Code
The complete, fully commented pipeline script can be reviewed directly in [analysis_pipeline.R](analysis_pipeline.R). 

To reproduce these results locally, clone this repository and run:
```R
source("analysis_pipeline.R")

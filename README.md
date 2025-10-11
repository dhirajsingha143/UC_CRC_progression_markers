# 🧬 Bioinformatics Exploration: UC–CRC Progression Markers

This repository is a part of the manuscript submitted to *British Journal of Cancer (BJC)* titled:  
**“From Inflammation to Malignancy: An Integrated Bioinformatics Exploration of Ulcerative Colitis – Associated Colorectal Cancer”**

It contains the complete workflow, analysis scripts, and processed data supporting the findings of this study.

---

<a name="readme-top"></a>

<div align="center">
  <h3>Repository Summary</h3>
  <p>Comprehensive bioinformatics pipeline integrating GEO transcriptomic datasets (GSE47908, GSE110224) to identify shared molecular markers and regulatory pathways linking Ulcerative Colitis (UC) with Colorectal Cancer (CRC) progression.</p>
</div>

---

## 📘 Table of Contents

- [About the Project](#about-the-project)
  - [Workflow Summary](#workflow-summary)
  - [Built With](#built-with)
- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
- [Repository Structure](#repository-structure)
- [Contributing](#contributing)
- [Contact](#contact)
- [Acknowledgements](#acknowledgements)
- [License](#license)

---

## About the Project

Inflammation due to Ulcerative Colitis (UC) increases the risk of developing Colorectal Cancer (CRC).  
This project performs an **integrated transcriptomic analysis** to identify **shared genes, pathways, and regulatory networks** involved in the progression from UC to CRC.

### Workflow Summary

1. **Data Retrieval**
2. **Pre-processing**
3. **DEG Analysis**
4. **Functional Enrichment**
5. **Network Construction**
6. **Validation**
7. **Diagnostic Evaluation**

---

### Built With

This analysis was conducted using:

<details>
  <summary>Tech Stack</summary>
  <ul>
    <li>R (v4.5.1)</li>
    <li>RStudio (v2025.09.0+387)</li>
    <li>Bioconductor (v3.21)</li>
    <li>R Packages: GEOquery, dplyr, tibble, limma, clusterProfiler, org.Hs.eg.db, STRINGdb, igraph, pROC, ggplot2, ..etc</li>
    <li>Data Sources: GEO, ENSEMBL, STRING, HPA</li>
  </ul>
</details>
---

## Getting Started

### Prerequisites

Ensure you have the following installed:

```r
R >= 4.5.1
Bioconductor >= 3.21
RStudio >= 2025.09.0+387
```

### Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/dhirajsingha143/UC_CRC_progression_markers.git
   cd UC_CRC_progression_markers
   ```
2. Reproduce the environment or Install required R packages:
   
   This project uses the [`renv`] (https://rstudio.github.io/renv/) package to ensure a fully reproducible R environment.
   To restore the same package versions used in the analysis:
   
   ```r
    # Step 1: Install renv (if not already installed) 
    install.packages("renv")

    # Step 2: Restore the project environment
    renv::restore()
   ```
   OR
   ```r
   install packages packages mentioned in main.R script
   ```

3. Run the main analysis script:

   ```r
   main.R
   ```
---

## Repository Structure
```text
UC_CRC_progression_markers/
│
├── LICENSE                           # MIT License © 2025 Dhiraj Singha
├── README.md                         # Project overview and usage guide
├── .gitignore                        # Ignored files for Git tracking
├── .Rhistory                         # R console history (auto-generated)
├── .Rprofile                         # R environment profile
├── UC_CRC_progression_markers.Rproj  # RStudio project file
├── renv.lock                         # R environment lock (for reproducibility)
├── main.R                            # Main pipeline script integrating all steps
│
├── Datasets/                         # Input data (GEO, meta, etc.)
│   ├── GSE47908_series_matrix.txt.gz
│   └── GSE110224_series_matrix.txt.gz
│
├── scripts/                          # R scripts for each analysis stage
│   ├── deg_functions.R
│   ├── enrichment.R
│   └── PPI_functions.R
│
├── results/                          # Output and analysis results
│   │
│   ├── DEG/                          # Differentially expressed gene results
│   │   ├── LSC_vs_HC   = LSC - HC_DEG_results.csv
│   │   ├── PC_vs_HC    = PC  - HC_DEG_results.csv
│   │   ├── CRC_vs_HC   = CRC - HC_DEG_results.csv
│   │   └── pdf files                 # Heat map, volcano plot, venn diagram pdf
│   │
│   ├── Enrichment/                   # GO/KEGG/GSEA functional analysis results
│   │   ├── CRC_vs_HC/
│   │   │   ├── UP/
│   │   │   └── DOWN/
│   │   ├── LSC_vs_HC/
│   │   │   ├── UP/
│   │   │   └── DOWN/
│   │   ├── PC_vs_HC/
│   │   │   ├── UP/
│   │   │   └── DOWN/
│   │   └── Input/
│   │       ├── UP/
│   │       └── DOWN/
│   │
│   ├── HPAdb/                        # HPA-based expression validation results
│   │
│   ├── PPI/                          # Network and hub gene analysis
│   │   ├── CRC_vs_HC/
│   │   │   ├── UP/
│   │   │   └── DOWN/
│   │   ├── LSC_vs_HC/
│   │   │   ├── UP/
│   │   │   └── DOWN/
│   │   └── PC_vs_HC/
│   │   │   ├── UP/
│   │   │   └── DOWN/
│   │   └── pdf files                 # Heat map, venn diagram pdf
│   │
│   ├── Pre_processing/               # Data normalization and QC outputs
│   │   ├── Raw/
│   │   ├── Normalized/
│   │   └── Batch_corrected/
│   │
│   ├── UC_to_CRC_HUB/                # Common hub gene analysis results
│   │   └── Boxplots/                 # Gene expression comparison plots
│   │
│   └── ROC/                          # ROC evaluation for diagnostic accuracy
│   │
│   └── sessionInfo.txt               # Session information
│
└── renv/                             # R environment management (auto-generated)
    ├── activate.R
    ├── library/
    └── settings.dcf
```
---
### Contributing

Contributions are welcome! Please fork the repository and create a pull request with your changes.

### Contact

<div align="left">
- Dhiraj Singha
- Department of Biotechnology, Amity University Uttar Pradesh, India
- Email: dhiraj.singha@s.amity.edu
- Github: [@dhirajsingha143](https://github.com/dhirajsingha143)
- LinkedIn: [Dhiraj Singha](https://www.linkedin.com/in/dhiraj-singha-b6871717a/)
</div>
<br>
<div align="left">
- Prof Kumud Bala
- Deputy Dean Student Welfare, Department of Biotechnology, Amity University Uttar Pradesh, India
- Email: kbala@amity.edu
- LinkedIn: [Prof Kumud Bala](https://www.linkedin.com/in/kumudbala/)
</div>

### Acknowledgements

- Department of Biotechnology, Amity University Uttar Pradesh, India
- Prof. Dr. Kumud Bala, Deputy Dean Student Welfare, Department of Biotechnology, Amity University Uttar Pradesh, India
- Original authors of the datasets
- R and Bioconductor communities for their invaluable tools and packages

### License
This repository is distributed under the MIT License copyright © 2025 DHIRAJ SINGHA. See `LICENSE` for more information.

<p align="right">(<a href="#readme-top">Back to top</a>)</p>

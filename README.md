
# Digital Twin Multi-organ Interactions!

This project is based on our study: Human Phenomics-based Digital Twin: Advancing Systematic Insights into Multi-organ Interactions
**The content of this repository is still being updated, and we will upload all the code after the article goes online.**

## Overview

This study developed an **organ-wide ‘digital human’ framework** to characterize multi-organ interaction patterns across both population and individual levels. Using organ-level phenotypes from the China Phenobank (2,361 phenotypes) and the UK Biobank (312 pheno-types) spanning seven major organs, our framework enables systematic observation and analysis of multi-organ interaction patterns. Unlike previous digital twin models that focused solely on individual organ functions, this organ-wide approach offers a comprehensive view of inter-organ dynamics and interactions.

## Highlights

1.	**A novel digital twin framework**: We provide the foundation to model and analyze multi-organ interac-tions, exploring key factors that influence organ interaction patterns and the underlying mechanisms of in-teraction perturbations.
2.	**Genomic and pathway discoveries**: Our analysis revealed dynamic inter-organ connections mediated by cellular, inflammatory, and metabolic pathways. We identified 44 genetic loci associated with multi-organ interactions, 18 of which are co-localized with known diseases.
3.	**Age- and sex-related variations**: We observed that variations in organ interaction patterns increase with age and are more pronounced in males than in females.
4.	**Predictive potential for disease**: Organ interaction perturbations show potential for predicting the long-term occurrence of major diseases.

## Supporting Packages

This repository contains supporting code for the study and utilizes the following major packages:

### R Packages:
- **ggplot2**: For data visualization and plotting.
- **dplyr**: For data manipulation and transformation.
- **tidyr**: For tidying data.
- **mediation**: For performing mediation analysis.
- **MendelianRandomization**: For Mendelian Randomization analysis.

### Python Packages:
- **hail**: For process gwas data
- **numpy**: For numerical computations.
- **pandas**: For data manipulation and analysis.
- **scipy**: For scientific computing and statistical analysis.
- **statsmodels**: For statistical modeling and regression analysis.
- **matplotlib**: For generating plots and visualizations.

---

## Installation and Setup

1. Clone the repository:

   ```bash
   git clone https://github.com/FD-ImageLab/OrganInteraction.git
   cd organ-interaction-atlas
   ```

2. Install required dependencies:

   For Python:
   for example
   ```bash
   pip install hail
   ```

   For R:
   for example:
   ```bash
   install.packages("forestploter")
   ```

3. Run the analysis:
   - For data preprocessing, correlation analysis, and Mendelian Randomization, use the respective scripts located under the `scripts` directory.


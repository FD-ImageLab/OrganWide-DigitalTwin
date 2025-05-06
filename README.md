# OrganInteraction
This project is based on our study: Atlas of the Human Organ Interaction Network.
## Overview

This study provides a **comprehensive atlas of human organ interactions**, exploring the dynamic relationships between 9 major organs and their implications in multi-organ diseases. By utilizing multi-omics and population-level phenotypic data, we reveal key genetic loci and molecular mediators influencing organ interactions and establish its predictive value for various diseases.

## Highlights

- **Comprehensive Atlas**: Maps out the interactions between 9 human organs.
- **Age and Gender Effects**: The perturbation of organ interactions increases with age, with more pronounced effects in males than females.
- **Genetic Insights**: Identifies 143 genetic loci regulating multi-organ interactions, with 24 loci co-localizing with common diseases.
- **Inflammatory and Metabolic Pathways**: Enrichment and mediation analyses suggest that inflammatory and metabolic factors are critical mediators of organ interactions.
- **Predictive Value**: Organ interaction perturbations serve as a new phenotype with predictive potential for major diseases.


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

---
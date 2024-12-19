# Install necessary packages if not already installed
# install.packages("remotes")
# remotes::install_github("MRCIEU/TwoSampleMR")
# install.packages("devtools")
# devtools::install_github("mrcieu/ieugwasr", force=TRUE)
# devtools::install_github("explodecomputer/plinkbinr")

# Load required libraries
library(tidyverse)      # For data manipulation
library(TwoSampleMR)    # For Two-sample Mendelian Randomization analysis
library("devtools")     # For installing GitHub packages
library(plinkbinr)      # For plink binary file manipulation
get_plink_exe()         # Get the path to the plink executable

# Load the 'ieugwasr' package for GWAS data manipulation
library("ieugwasr")

# Set the working directory to your project folder
setwd("/cpfs01/projects-HDD/cfff-6117e6302119_HDD/lm_21210880006/wangza/projects/Comorbidity")

# Source custom utility functions
source("scripts/utils/Func.R")

# List of exposure files (disease-related files) to analyze
exposure_list <- list.files("fuma/ukb/disease", pattern = "*.tsv.gz", full.names = TRUE)
# Extract exposure names from the filenames
exposure_names <- map_chr(str_split(basename(exposure_list), "_"), ~ .x[[1]][1])

# List of outcome files (phenotype files) to analyze
outcome_list <- list.files("fuma/ukb/new_pheno", pattern = "*.tsv.gz", full.names = TRUE)
# Extract outcome names from the filenames
outcome_names <- map_chr(str_split(basename(outcome_list), "_"), ~ .x[[1]][1])

# Initialize empty data frames to store results
dat_df <- data.frame()
mr_df <- data.frame()
mr_heterogeneity_df <- data.frame()
mr_pleiotropy_df <- data.frame()

# Loop through each exposure file for analysis
for (i in seq_along(exposure_list)) {
    try({
        # Load and process exposure data
        exposure_name <- exposure_names[i]
        exposure_dat <- read_exposure_data(exposure_list[i],
            sep = "\t",
            snp_col = "rsid",       # Column with SNP IDs
            beta_col = "beta",      # Column with SNP effect size (beta)
            se_col = "se",          # Column with standard error of the SNP effect
            effect_allele_col = "alt", # Column with effect allele
            other_allele_col = "ref", # Column with reference allele
            clump = FALSE
        )
        
        # Filter SNPs based on exposure p-value threshold
        exposure_dat <- exposure_dat %>% filter(pval.exposure < 5 * 10e-8)

        # Assign the exposure name to the data
        exposure_dat$id.exposure <- exposure_name

        # Create a temporary data frame for SNPs and p-values
        exposure_dat_temp <- dplyr::tibble(
            rsid = exposure_dat$SNP,
            pval = exposure_dat$pval.exposure,
            id = exposure_dat$id.exposure
        )

        # Perform LD clumping to remove correlated SNPs
        exposure_clump_temp <- ld_clump(exposure_dat_temp,
            plink_bin = get_plink_exe(),  # Path to PLINK executable
            bfile = "ld/EUR"              # Reference LD data
        )
        # Filter the exposure data to retain only the clumped SNPs
        exposure_dat_clump <- exposure_dat[exposure_dat$SNP %in% exposure_clump_temp$rsid, ]

        # Loop through each outcome file to perform MR analysis
        for (j in seq_along(outcome_list)) {
            # Load and process outcome data
            outcome_name <- outcome_names[j]
            outcome_dat <- read_outcome_data(outcome_list[j],
                snps = exposure_dat_clump$SNP,  # Use the SNPs from the exposure data
                sep = "\t",                    # Separator for the data file
                snp_col = "rsid",              # Column with SNP IDs
                beta_col = "beta",             # Column with SNP effect size (beta)
                se_col = "se",                 # Column with standard error of the SNP effect
                effect_allele_col = "alt",     # Column with effect allele
                other_allele_col = "ref",     # Column with reference allele
                pval_col = "P",                # Column with p-values
            )
            outcome_dat$id.outcome <- outcome_name
            outcome_dat$outcome <- outcome_name

            # Harmonize the exposure and outcome data to align the alleles
            dat <- harmonise_data(
                exposure_dat = exposure_dat_clump,
                outcome_dat = outcome_dat
            )
            # Append the harmonized data to the main dataframe
            dat_df <- rbind(dat_df, dat)

            ############### MR Analysis ###################
            # Perform Mendelian Randomization analysis
            single_mr <- mr(dat)
            mr_df <- rbind(mr_df, single_mr)

            ############## Heterogeneity Test ############
            # Perform heterogeneity test to check for inconsistency across SNPs
            single_mr_heterogeneity <- mr_heterogeneity(dat, method_list = c())  # Specify the methods for heterogeneity test
            mr_heterogeneity_df <- rbind(mr_heterogeneity_df, single_mr_heterogeneity)

            ############## Pleiotropy Test ##############
            # Perform pleiotropy test to check for violations of the MR assumption
            single_mr_pleiotropy <- mr_pleiotropy_test(dat)
            mr_pleiotropy_df <- rbind(mr_pleiotropy_df, single_mr_pleiotropy)

            # Save the MR results, heterogeneity results, and pleiotropy results as CSV files
            check_path(paste0("data/MR/", exposure_name, "_", outcome_name, "_", "mr.csv"))
            write.csv(single_mr, paste0("data/MR/", exposure_name, "_", outcome_name, "_", "mr.csv"))
            write.csv(single_mr_heterogeneity, paste0("data/MR/", exposure_name, "_", outcome_name, "_", "mr_heterogeneity.csv"))
            write.csv(single_mr_pleiotropy, paste0("data/MR/", exposure_name, "_", outcome_name, "_", "mr_pleiotropy.csv"))

            # Print progress
            print(paste0("i is ", i, ", j is ", j))
        }
    })
}

# Save the overall MR analysis results, heterogeneity, and pleiotropy test results
write.csv(mr_df, "data/MR/mr.csv")
write.csv(mr_heterogeneity_df, "data/MR/mr_heterogeneity.csv")
write.csv(mr_pleiotropy_df, "data/MR/mr_pleiotropy.csv")

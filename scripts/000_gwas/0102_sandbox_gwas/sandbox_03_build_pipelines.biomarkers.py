import sys
import hail as hl
import os

project_path = os.path.expanduser("/cpfs01/projects-HDD/cfff-6117e6302119_HDD/lm_21210880006/projects/indNet/")


############## phenotype ##################
ht_phenotypes = hl.read_table(project_path + 'data/processed/biomarker/biomarkers_gwas.ht')
phenotypes = list(ht_phenotypes.row_value)
chunk_size = 11
groups = [phenotypes[i:(i + chunk_size)] for i in range(0, len(phenotypes), chunk_size)]

i = 0
for group in groups:
    print(f"In pipeline {i}, the phenotypes include: {group}")
    ht = ht_phenotypes.select(*group)
    ht.write(project_path + f'data/processed/biomarker/biomarkers_gwas.pipeline_{i}.ht', overwrite=True)
    i += 1


######## disease #############
ht_phenotypes = hl.read_table(project_path + 'data/processed/biomarker/disease_gwas.ht')
phenotypes = list(ht_phenotypes.row_value)

chunk_size = 11
groups = [phenotypes[i:(i + chunk_size)] for i in range(0, len(phenotypes), chunk_size)]


i = 100
for group in groups:
    print(f"In pipeline {i}, the phenotypes include: {group}")
    ht = ht_phenotypes.select(*group)
    ht.write(project_path + f'data/processed/biomarker/biomarkers_gwas.pipeline_{i}.ht', overwrite=True)
    i += 1
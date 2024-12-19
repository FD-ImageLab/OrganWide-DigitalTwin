import hail as hl
import pandas as pd
from hail.utils import new_temp_file
import argparse
import os
import re


argparse = argparse.ArgumentParser(description="add arguments")
argparse.add_argument("--pipeline", type=int, help="pipeline number")
args = argparse.parse_args()
pipeline = args.pipeline

project_path = os.path.expanduser("/cpfs01/projects-HDD/cfff-6117e6302119_HDD/lm_21210880006/projects/indNet/")
pipeline_file_path = project_path + "scripts/pheno_dict.csv"


temp_pipeline_df = pd.read_csv(pipeline_file_path)

codes = dict(zip(temp_pipeline_df["description"], temp_pipeline_df["fid"]))

ht_autosomes = hl.read_table(project_path + f"data/processed/gwas/biomarker_gwas_results.autosomes.pipeline_{pipeline}.ht")

ht_results = ht_autosomes


ht_results = ht_results.annotate(variant=hl.delimit(hl.array([
    ht_results['locus'].contig,
    hl.str(ht_results['locus'].position),
    ht_results['alleles'][0],
    ht_results['alleles'][1]]), delimiter=':'))
ht_results = ht_results.key_by('locus','alleles')
ht_results = ht_results.repartition(116)
# ht_results = ht_results.cache()


GWAS_qc = project_path + 'data/raw/gene/variants.tsv.bgz'
ukb_qc = hl.import_table(GWAS_qc)
ukb_qc = ukb_qc.annotate(vstruct = hl.parse_variant(ukb_qc.variant))
ukb_qc = ukb_qc.annotate(locus = ukb_qc.vstruct.locus, alleles = ukb_qc.vstruct.alleles).key_by('locus','alleles')

ht_results = ht_results.join(ukb_qc, how="inner")

count = 1
phenotypes = ht_results['phenotypes'].collect()[0]


for i, phenotype in enumerate(phenotypes):
    pattern = re.compile(r'\(.*?\)')
    variable_type = re.sub(pattern, "", phenotype).replace(" ", "")
    code = codes[phenotype]
    print(f'Exporting LDSC sumstats for trait {code} ({count})...')
    ht_export = ht_results.annotate(
        A1 = ht_results.alleles[1],
        A2 = ht_results.alleles[0],
        N = ht_results['n'][i],
        Z = ht_results['t_stat'][i][0], 
        beta=ht_results['beta'][i][0],
        se=ht_results['standard_error'][i][0],
        P = ht_results["p_value"][i][0]
        )
    ht_export = ht_export.select('chr','pos','rsid','A1','A2','ref','alt','N','Z', 'beta', 'se', 'P')
    
    ht_export = ht_export.filter(hl.is_finite(ht_export["P"]))
    ht_export.export(project_path + f'data/processed/gwas/{code}_{variable_type}_imputed_v3_fuma.tsv.gz')
    count += 1


print('#######################')
print('## COMPLETED     ######')
print('## Pipeline number: {}'.format(pipeline))
print('#######################')

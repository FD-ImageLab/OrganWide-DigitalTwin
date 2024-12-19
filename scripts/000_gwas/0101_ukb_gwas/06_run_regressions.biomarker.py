import sys
import hail as hl
import argparse
import os

spark_conf = {
    'spark.master': 'local[30]',
    'spark.executor.memory': '128g',
    'spark.driver.memory': '128g',
}
hl.init(spark_conf=spark_conf)

project_path = os.path.expanduser("/cpfs01/projects-HDD/cfff-6117e6302119_HDD/lm_21210880006/projects/indNet/")

contig = "autosomes"
# 第几个block
argparse = argparse.ArgumentParser(description="add arguments")
argparse.add_argument("--pipeline", type=int, help="pipeline number")
args = argparse.parse_args()
pipeline = args.pipeline

if contig == 'autosomes':
    contig_expr = 'c{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}'
else:
    contig_expr = contig

if contig not in set(['autosomes', 'chrX', 'chrXY']):
    raise ValueError(f'Invalid contig argument "{contig}" - must be one of {{"autosomes", "chrX", "chrXY"}}.')


ht_phenotypes = hl.read_table(project_path + f'data/processed/biomarker/biomarkers_gwas.pipeline_{pipeline}.ht')
ht_covariates = hl.read_table(project_path + f"data/processed/biomarker/gwas_covariates.ht")
ht_variants = hl.read_table(project_path + f'data/processed/gene/gwas_variants.ht')

mt = hl.import_bgen(
    path=project_path+ f'data/raw/gene/bgen/ukb22828_{contig_expr}_b0_v3.bgen',
    sample_file=project_path+ f'data/raw/gene/sample/{contig}.sample',
    entry_fields=['dosage'],
    variants=ht_variants.variant)

mt = mt.annotate_cols(
    phenotypes=ht_phenotypes[mt.s],
    covariates=ht_covariates[mt.s])

# hail是否能够自动解决这个问题？暂时先加上这个过滤
mt = mt.filter_cols(hl.is_defined(mt.phenotypes) & hl.is_defined(mt.covariates))

phenotypes = list(mt['phenotypes'].keys())
print(phenotypes)

ht = hl.linear_regression_rows(
    y=[[mt['phenotypes'][y]] for y in phenotypes],
    x=mt.dosage,
    covariates=[1, *[mt['covariates'][x] for x in list(mt['covariates'].keys())]],
    pass_through=['varid', 'rsid'])

ht = ht.annotate_globals(phenotypes=phenotypes)

ht.write(project_path + f"data/processed/gwas/biomarker_gwas_results.{contig}.pipeline_{pipeline}.ht",
         overwrite=True)
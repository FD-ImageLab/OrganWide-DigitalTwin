import pandas as pd
import numpy as np
import csv
import hail as hl
import os
import datetime

from utils.tools import check_path
from utils.tools import write_pipeline_pheno
from utils.tools import pipeline_add_82_pheno
spark_conf = {
    'spark.master': 'local[30]',
    'spark.executor.memory': '128g',
    'spark.driver.memory': '128g',
}
hl.init(spark_conf=spark_conf)

project_path = os.path.expanduser("/cpfs01/projects-HDD/cfff-6117e6302119_HDD/lm_21210880006/projects/indNet/")



GWAS_bgen = [project_path + 'data/raw/gene/bgen/' + filename  for filename in os.listdir(project_path + 'data/raw/gene/bgen/') if filename.endswith('.bgen')]
GWAS_sample = project_path + 'data/raw/gene/sample/autosomes.sample'
GWAS_ht  = project_path + 'data/processed/gene/gwas_variants.ht'
GWAS_mfi = project_path + 'data/raw/gene/ukb_imp_mfi/gwas_variants.imputed_v3.mfi.ht'


ef = ['GT']
force_all = False



if (not check_path(GWAS_ht)) or force_all:
    ht_qc_variants = hl.import_table(project_path + f'data/raw/gene/variants.tsv', impute=True, missing="")
    ht_qc_variants = ht_qc_variants.annotate(variant = hl.parse_variant(ht_qc_variants.variant))
    ht_qc_variants = ht_qc_variants.annotate(locus = ht_qc_variants.variant.locus, alleles = ht_qc_variants.variant.alleles)
    ht_qc_variants = ht_qc_variants.filter((ht_qc_variants.info > 0.8) &
                        (ht_qc_variants.minor_AF >= 0.01) &
                        (ht_qc_variants.p_hwe > 1e-10))
    ht_qc_variants = ht_qc_variants.key_by('locus','alleles')
    ht_qc_variants.write(GWAS_ht, overwrite=True)


# 处理表型数据
if (not check_path(project_path + "data/processed/biomarker/pipeline_pheno.ht")) or force_all or True:
    type_dict = {
    "eid": hl.tstr,
    "y_res_Brain_Heart": hl.tfloat64,
    "y_res_Brain_Lung": hl.tfloat64,
    "y_res_Brain_Liver": hl.tfloat64,
    "y_res_Brain_Spleen": hl.tfloat64,
    "y_res_Brain_Pancreas": hl.tfloat64,
    "y_res_Brain_Kidney": hl.tfloat64,
    "y_res_Heart_Lung": hl.tfloat64,
    "y_res_Heart_Liver": hl.tfloat64,
    "y_res_Heart_Spleen": hl.tfloat64,
    "y_res_Heart_Pancreas": hl.tfloat64,
    "y_res_Heart_Kidney": hl.tfloat64,
    "y_res_Lung_Liver": hl.tfloat64,
    "y_res_Lung_Spleen": hl.tfloat64,
    "y_res_Lung_Pancreas": hl.tfloat64,
    "y_res_Lung_Kidney": hl.tfloat64,
    "y_res_Liver_Spleen": hl.tfloat64,
    "y_res_Liver_Pancreas": hl.tfloat64,
    "y_res_Liver_Kidney": hl.tfloat64,
    "y_res_Spleen_Pancreas": hl.tfloat64,
    "y_res_Spleen_Kidney": hl.tfloat64,
    "y_res_Pancreas_Kidney": hl.tfloat64
}
    pipeline_ht = hl.import_table(project_path + "data/ukb_new_pheno_without_regression/whole/ukb_new_pheno.csv", missing="NA", delimiter=",", types=type_dict)
    pipeline_ht = pipeline_ht.annotate(eid = hl.str(pipeline_ht.eid)) # convert column s to string
    pipeline_ht = pipeline_ht.key_by("eid")
    pipeline_ht.write(project_path + "data/processed/biomarker/pipeline_pheno.ht", overwrite=True)


# 处理疾病数据
if (not check_path(project_path + "data/processed/biomarker/pipeline_disease.ht")) or force_all:
    temp_df = pd.read_csv(project_path + "data/clean/disease.csv", sep=",")
    pipeline_ht = hl.import_table(project_path + "data/clean/disease.csv", impute=True, missing="", delimiter=",")
    pipeline_ht = pipeline_ht.annotate(eid = hl.str(pipeline_ht.eid)) # convert column s to string
    pipeline_ht = pipeline_ht.key_by("eid")
    pipeline_ht.write(project_path + "data/processed/biomarker/pipeline_disease.ht", overwrite=True)


##### covariates file generate #####
if not check_path(project_path + "data/processed/biomarker/gwas_covariates.ht") or force_all:
    pca = pd.read_csv(project_path + "data/ukb_sqc_v3.txt", sep=",")
    pca = pca.loc[:, ["eid", "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10",
                      "22009-0.11", "22009-0.12", "22009-0.13", "22009-0.14", "22009-0.15", "22009-0.16", "22009-0.17", "22009-0.18", "22009-0.19", "22009-0.20", 
                      "22006-0.0", "22019-0.0", "22020-0.0"]]
    cnames = ["eid", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10",
            "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", 
            "in.white.British.ancestry.subset", "putative.sex.chromosome.aneuploidy", "used.in.pca.calculation"]
    pca.columns = cnames
    pca["eid"] = pca["eid"].astype(str) # convert column eid to string
    qced_pca = pca[(pca["used.in.pca.calculation"] == 1) & (~(pca["in.white.British.ancestry.subset"].isna())) & (~(pca["putative.sex.chromosome.aneuploidy"] == 1))]
    qced_pca = qced_pca.drop(columns=["used.in.pca.calculation", "in.white.British.ancestry.subset", "putative.sex.chromosome.aneuploidy"])


    covar = hl.import_table(project_path + "data/clean/covariates.tsv", impute=True, missing="", delimiter="\t")
    covar = covar.annotate(eid = hl.str(covar.eid)) # convert column eid to string
    covar = covar.key_by("eid")
    covar = covar.select(covar.isMale, covar.age, covar.age_squared, covar.age_isMale, covar.age_squared_isMale)
    pca = hl.Table.from_pandas(qced_pca)
    pca = pca.annotate(eid = hl.str(pca.eid)) # convert column eid to string
    pca = pca.key_by("eid")
    
    newcovar = covar.join(pca, how="inner")
    newcovar = newcovar.key_by("eid")
    check_path(project_path + "data/processed/biomarker/gwas_covariates.ht")
    newcovar.write(project_path + "data/processed/biomarker/gwas_covariates.ht", overwrite=True)






# # 为bgen文件创建索引

# # 设置BGEN文件路径
# bgen_file_path = project_path + "data/raw/gene/bgen/"
# # 获取BGEN文件列表
# bgen_files = [f for f in os.listdir(bgen_file_path) if f.endswith(".bgen")]
# # 为每个BGEN文件创建索引文件
# for bgen_file in bgen_files:
#     # 构造BGEN文件路径和索引文件路径
#     bgen_path = os.path.join(bgen_file_path, bgen_file)
#     idx2_path = bgen_path + ".idx2" 
#     if not check_path(idx2_path) or force_all:
#         # 创建索引文件
#         hl.index_bgen(bgen_path, index_file_map={bgen_path: idx2_path}, 
#                     contig_recoding={"01": "1", "02": "2", "03": "3", "04": "4", "05": "5", "06": "6", "07": "7", "08": "8", "09": "9"}, 
#                     skip_invalid_loci=True)

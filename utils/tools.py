# utility to check if filepath exists
import os
import csv
import pandas as pd
import numpy as np
from typing import List

def check_path(filepath):
    # check for file
    stat = os.path.exists(filepath)
    return stat


project_path = os.path.expanduser("~/PheWAS/")
dataPath = project_path + "data/raw/biomarker/ukb670788_all.csv"
# utility to get FIA from FID, file from main dataset. e.g. get 31-0.0 from 31
def get_UDIs_from_fids(fids: List[int], dataPath=dataPath):
    with open(dataPath, 'r') as infile:
        UDIs = csv.DictReader(infile)
        UDIs = UDIs.fieldnames
    field_ids = []
    for i, UDI in enumerate(UDIs):
        field_id = UDI.split("-")[0]
        field_ids.append(field_id)
    field_ids = pd.Series(field_ids)
    a = field_ids.isin(list(map(str, fids)))
    b = ~field_ids.duplicated()
    idx = list(a & b)
    UDIs = np.array(UDIs)
    idx = np.array(idx)
    used_UDIs = list(map(str, UDIs[idx]))
    fids_in_UDIs = fids.astype(str).isin(list(map(str, field_ids)))
    return used_UDIs,fids_in_UDIs



raw_pheno_path = project_path + "data/raw/biomarker/ukb670788_all.csv"
pipeline_file_path = project_path + "scripts/pipeline_ecg.csv"
to_tsv_path = project_path + "data/processed/biomarker/pipeline_pheno_ecg.tsv"
# need to give main dataset path and pipeline file path
def write_pipeline_pheno(raw_pheno_path, pipeline_file_path, to_tsv_path=to_tsv_path):
    pipeline = pd.read_csv(pipeline_file_path)
    fids = pipeline["fid"]
    used_UDIs, fids_in_UDIs = get_UDIs_from_fids(fids, dataPath=raw_pheno_path)
    used_UDIs = ["eid"] + used_UDIs
    data = pd.read_csv(raw_pheno_path, usecols=used_UDIs)
    cnames = ["s"] + list(np.array(list(map(str, pipeline["names"])))[np.array(fids_in_UDIs)])
    data.columns = cnames
    data = qc(data)
    data.to_csv(to_tsv_path, index=False, sep="\t")
    print("Pipeline phenotype file is writting to: " + to_tsv_path)
    return None


def pipeline_add_82_pheno():
    bridge_path = project_path + "data/raw/biomarker/ukb96511bridge18545.txt"
    temp_df1 = pd.read_csv(project_path + "data/processed/biomarker/pipeline_pheno_ecg.tsv", sep="\t")
    temp_df2 = pd.read_excel(project_path + "data/raw/biomarker/bai_clinical_measures_39k.xlsx", index_col=0)
    temp_df2.reset_index(inplace=True)
    temp_df2.rename(columns={"index":"s"}, inplace=True)
    bridge = pd.read_csv(bridge_path, sep=" ", header=None, index_col=False, names=["AID96511", "AID18545"])
    bridge_dict = dict(bridge[["AID18545", "AID96511"]].values)
    temp_df2["s"] = temp_df2["s"].map(bridge_dict)
    temp_df3 = pd.merge(temp_df1, temp_df2, on="s", how="left")
    temp_df3 = temp_df3.set_index("s")
    temp_df3_idx = temp_df3.index
    temp_df3 = remove_outliers(temp_df3)
    temp_df3.index = temp_df3_idx
    temp_df3.reset_index(inplace=True)
    temp_df3.rename(columns={"index":"s"}, inplace=True)
    
    temp_df3.to_csv(project_path + "data/processed/biomarker/pipeline_pheno.tsv", sep="\t", index=False)
    print("Bai 82 phenotypes added to: " + project_path + "data/processed/biomarker/pipeline_pheno.tsv")
    return None

##################
### data qc
##################
def qc(df):
    if "12-lead ECG measuring method" in df.columns:
        df.loc[((df["12-lead ECG measuring method"] == 6) | (df["12-lead ECG measuring method"] == 7)), 
                ["Ventricular rate", "P duration", "QRS duration", "PQ interval", "QT interval", "QTC interval", 
                 "RR interval", "PP interval", "P axis", "R axis", "T axis"]] = np.nan
        df.drop(columns=["12-lead ECG measuring method"], inplace=True)
        return(df)
    else:
        return(df)
    
        
##################################
### outlier removal
##################################
def remove_outliers(df, id_col=True):
    if id_col:
        # 如果第一列是id列，则先将其保留
        id_series = df.iloc[:, 0]
        data_df = df.iloc[:, 1:]  # 不包含id列的数据
    else:
        data_df = df
    
    # 计算每列的标准差
    stds = data_df.std()
    
    # 遍历每列，替换大于3倍标准差的值为NaN
    for col in data_df.columns:
        data_df[col] = np.where(np.abs(data_df[col] - data_df[col].mean()) > 3 * stds[col], np.nan, data_df[col])
    
    # 合并id列和处理后的数据列
    if id_col:
        df_cleaned = pd.concat([id_series, data_df], axis=1)
    else:
        df_cleaned = data_df
    
    # 返回处理后的DataFrame
    return df_cleaned


if __name__ == "__main__":
    # used_UDIs = get_UDIs_from_fids([22423, 22420])
    # write_pipeline_pheno(raw_pheno_path, pipeline_file_path)
    
    # temp_df = pd.read_csv(project_path + "data/processed/biomarker/pipeline_pheno.tsv", sep="\t")
    # df = qc(temp_df)
    # print(df)
    # pipeline_add_82_pheno()

    # df = pd.DataFrame([[1,2], [2, 3], [3, 100], [100, 4], [2, 3], [2, 3], [2, 3], [2, 3], [2, 3], [2, 3], [2, 3], [2, 3]])
    # df = remove_outliers(df)
    # print(df)
    fids = pd.Series([25781])
    a, b = get_UDIs_from_fids(fids, dataPath=dataPath)
    print(a)
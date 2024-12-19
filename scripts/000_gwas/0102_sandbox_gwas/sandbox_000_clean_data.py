# 创建疾病/表型/协变量数据
import os
import pandas as pd
import numpy as np
import hail as hl
from functools import reduce
import datetime
import math
from utils.tools import check_path

project_path = os.path.expanduser("/cpfs01/projects-HDD/cfff-6117e6302119_HDD/lm_21210880006/projects/indNet/")
force_all = True



##### covariates file generate #####
dataPath = project_path + "data/clean/basic_info.csv"
if not check_path(project_path + "data/clean/covariates.csv") or force_all:
    covar_UDIs = ["eid", "31-0.0", "34-0.0", "52-0.0", "53-2.0", '21001-2.0', '21002-2.0']
    covar_UDIs_name = ["eid", "isMale", "birth_year", "birth_month", "attending_date", "bmi", "weight"]
    data = pd.read_csv(dataPath, usecols=covar_UDIs)
    data["eid"] = data["eid"].astype(str)
    raw_data = data.copy()
    data = data.dropna(axis=0)
    data["34-0.0"] = data["34-0.0"].astype(int)
    data["52-0.0"] = data["52-0.0"].astype(int)
    eid = data["eid"].values
    sex = data['31-0.0'].values
    age = np.zeros(len(data))
    for i in range(len(data)):
        # Calculate age
        d1 = datetime.date(data.iloc[i]['34-0.0'], data.iloc[i]['52-0.0'], 15)
        s = data.iloc[i]['53-2.0']
        d2 = datetime.date(int(s[:4]), int(s[5:7]), int(s[8:]))
        age[i] = np.round((d2 - d1).days / 365.25, 1)
    weight = data['21002-2.0'].values
    bmi = data['21001-2.0'].values
    height = np.round(np.sqrt(weight / bmi) * 100)
    # Confounding factors
    conf = np.stack((eid, sex, age, weight, height), axis=1)
    df_conf = pd.DataFrame(conf, index=data.index, columns=['eid', 'isMale', 'age', 'weight', 'height'])
    df_conf["age_squared"] = df_conf["age"] ** 2
    df_conf["age_isMale"] = df_conf["age"] * df_conf["isMale"]
    df_conf["age_squared_isMale"] = df_conf["age_squared"] * df_conf["isMale"]
    df_conf["eid"] = df_conf["eid"].astype(str)


    sqc = pd.read_csv(project_path + "data/ukb_sqc_v3.txt", sep=",")
    sqc = sqc.loc[:, ["eid", "22009-0.1", "22009-0.2", "22009-0.3", "22009-0.4", "22009-0.5", "22009-0.6", "22009-0.7", "22009-0.8", "22009-0.9", "22009-0.10",
                        "22009-0.11", "22009-0.12", "22009-0.13", "22009-0.14", "22009-0.15", "22009-0.16", "22009-0.17", "22009-0.18", "22009-0.19", "22009-0.20", 
                        "22006-0.0", "22019-0.0", "22020-0.0"]]
    cnames = ["eid", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10",
            "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20", 
            "in.white.British.ancestry.subset", "putative.sex.chromosome.aneuploidy", "used.in.pca.calculation"]
    sqc.columns = cnames
    sqc["eid"] = sqc["eid"].astype(str) # convert column s to string
    eur = sqc[(~(sqc["in.white.British.ancestry.subset"].isna()))]


    df_conf = df_conf[df_conf["eid"].isin(eur["eid"])]
    df_conf["eid"] = df_conf["eid"].astype(str)
    df_conf.sort_values(by="eid", inplace=True)
    df_conf.to_csv(project_path + "data/clean/covariates.tsv", index=False, sep="\t")
    df_conf.to_csv(project_path + "data/clean/covariates.csv", index=False)

    

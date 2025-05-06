import pandas as pd
import os

data = pd.read_csv("data/disease_days.csv", index_col="eid")
disease = pd.raed_csv("data/disease_with_pn.csv", index_col="eid")
n_nondisease = (data > 5000).sum()
n_disease_minus = (data < 0).sum()
n_disease_plus = (data > 0).sum() - (data > 5000).sum()
disease_num = pd.DataFrame({"n_nondisease": n_nondisease, "n_disease_minus": n_disease_minus, "n_disease_plus": n_disease_plus})
disease_num.rename_axis("disease", axis=0, inplace=True)
min_max_date = pd.read_csv("data/disease_min_max_date.csv", index_col=0)
min_max_date.rename_axis("disease", axis=0, inplace=True)
disease_meta_info = pd.merge(disease_num, min_max_date, on="disease", how="outer")
disease_meta_info.to_csv("data/disease_meta_info.csv", index=True)
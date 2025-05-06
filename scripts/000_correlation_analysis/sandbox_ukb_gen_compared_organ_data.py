import pandas as pd
import numpy as np
import os

heart_dict = pd.read_csv("Sandbox/ukb_category/Heart.csv")
heart_cols = heart_dict["UDI"][heart_dict["Category"] == "157"].tolist()
heart_cols = ["eid"] + heart_cols

brain_dict = pd.read_csv("Sandbox/ukb_category/Brain.csv")
brain_cols = brain_dict["UDI"][brain_dict["Category"] == 1101].tolist()
brain_cols = ["eid"] + brain_cols

kidney_dict = pd.read_csv("Sandbox/ukb_category/Kidney.csv")
kidney_cols = kidney_dict["UDI"][kidney_dict["Region"] == "Kidney"].tolist()
kidney_cols = ["eid"] + kidney_cols

liver_dict = pd.read_csv("Sandbox/ukb_category/Liver.csv")
liver_cols = liver_dict["UDI"][liver_dict["Region"] == "Liver"].tolist()
liver_cols = ["eid"] + liver_cols

lung_dict = pd.read_csv("Sandbox/ukb_category/Lung.csv")
lung_cols = lung_dict["UDI"][lung_dict["Region"] == "Lung"].tolist()
lung_cols = ["eid"] + lung_cols

pancreas_dict = pd.read_csv("Sandbox/ukb_category/Pancreas.csv")
pancreas_cols = pancreas_dict["UDI"][pancreas_dict["Region"] == "Pancreas"].tolist()
pancreas_cols = ["eid"] + pancreas_cols

spleen_dict = pd.read_csv("Sandbox/ukb_category/Spleen.csv")
spleen_cols = spleen_dict["UDI"][spleen_dict["Region"] == "Spleen"].tolist()
spleen_cols = ["eid"] + spleen_cols

heart = pd.read_csv("data/ukb674047.csv", index_col="eid", usecols=heart_cols)
brain = pd.read_csv("data/ukb674208.csv", index_col="eid", usecols=brain_cols)
kidney = pd.read_csv("data/ukb674208.csv", index_col="eid", usecols=kidney_cols)
liver = pd.read_csv("data/ukb674208.csv", index_col="eid", usecols=liver_cols)
lung = pd.read_csv("data/ukb674208.csv", index_col="eid", usecols=lung_cols)
pancreas = pd.read_csv("data/ukb674208.csv", index_col="eid", usecols=pancreas_cols)
spleen = pd.read_csv("data/ukb674208.acsv", index_col="eid", usecols=spleen_cols)

heart.to_csv("data/heart.csv", index=True)
brain.to_csv("data/brain.csv", index=True)
kidney.to_csv("data/kidney.csv", index=True)
liver.to_csv("data/liver.csv", index=True)
lung.to_csv("data/lung.csv", index=True)
pancreas.to_csv("data/pancreas.csv", index=True)
spleen.to_csv("data/spleen.csv", index=True)


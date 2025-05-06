import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
import re
import os
from scipy import stats
from scripts.utils.my_fdr import fdr_threshold
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt


df_new_pheno = pd.read_csv('data/sandbox_new_pheno/whole/sandbox_new_pheno.csv', index_col="eid")


df_disease = pd.read_csv('Sandbox/disease_sandbox.csv', index_col=0)



diseases = list(df_disease.columns)
phenotypes = list(df_new_pheno.columns)
phenotype_name = dict(zip(phenotypes, df_new_pheno.columns.str.removeprefix("y_res_")))



phenotype_unit = list(df_new_pheno.columns)

# For each disease of interest, perform logistic regression
#
# disease ~ Imaging_Phenotype
#
n_disease = len(diseases)
n_pheno = len(phenotypes)
df_beta = pd.DataFrame(index=diseases, columns=phenotypes, dtype=float)
df_conf = df_beta.copy()
df_p_val = df_beta.copy()
df_n = {}


odd_df = pd.DataFrame()
for pheno in phenotypes:
    for disease in diseases:

        # Dependent variables
        y = df_disease[disease]
        y = y.loc[((y == 0) | (y == -1))]
        y = y.replace(-1, 1)
        # Independent variables
        # X = pd.concat([df_info[['Sex', 'Age']], df[[pheno]]], axis=1)
        X = df_new_pheno[[pheno]]

        # Standard deviation of the phenotype for normalisation
        pheno_std = df_new_pheno[pheno].std()

        # Remove NaN values before regression
        valid_idx = pd.merge(X, y, on="buId", how="inner").dropna().index
        df_n[pheno] = len(valid_idx)
        
        # Logistic regression
        model = sm.Logit(y.loc[valid_idx], X.loc[valid_idx]).fit(disp=False)
        alpha = 0.05
        t = stats.t.isf(alpha / 2, model.df_resid)

        df_beta.loc[disease, pheno] = model.params[pheno] * pheno_std
        df_conf.loc[disease, pheno] = t * model.bse[pheno] * pheno_std
        df_p_val.loc[disease, pheno] = model.pvalues[pheno]

        odds = np.exp(df_beta.loc[disease, pheno])
        odds_lower = np.exp(df_beta.loc[disease, pheno] - df_conf.loc[disease, pheno])
        odds_upper = np.exp(df_beta.loc[disease, pheno] + df_conf.loc[disease, pheno])
        # print('{0}, per {1:.1f} {2}: Odds = {3:.2f} [{4:.2f}, {5:.2f}], p = {6:.3f}'.format(
        #     disease, pheno_std, pheno,
        #     odds, odds_lower, odds_upper, df_p_val.loc[disease, pheno]))
        temp = pd.DataFrame({"pheno": [pheno], "disease": [disease], "odds": [odds], "p": [df_p_val.loc[disease, pheno]]})
        odd_df = pd.concat([odd_df, temp])
odd_df.to_csv("data/sandbox_disease_odd.csv", index=False)




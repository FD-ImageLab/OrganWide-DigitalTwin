import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
import re
import os
from scipy import stats
from utils.my_fdr import fdr_threshold
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# heart = pd.read_csv("data/heart.csv", index_col="eid")
# heart = heart.dropna()

# brain = pd.read_csv("data/brain.csv", index_col="eid")
# brain = brain.dropna()

# kidney = pd.read_csv("data/kidney.csv", index_col="eid")
# kidney = kidney.dropna()

# liver = pd.read_csv("data/liver.csv", index_col="eid")
# liver = liver.dropna()

# lung = pd.read_csv("data/lung.csv", index_col="eid")
# lung = lung.dropna()

# pancreas = pd.read_csv("data/pancreas.csv", index_col="eid")
# pancreas = pancreas.dropna()

# spleen = pd.read_csv("data/spleen.csv", index_col="eid")
# spleen = spleen.dropna()


# pca = PCA(n_components=1)
# pca.fit(heart)
# heart_pca = pca.transform(heart)
# heart_pca = pd.DataFrame({"eid": heart.index, "heart_pca1": heart_pca[:, 0]})
# heart_pca.set_index("eid")

# pca = PCA(n_components=1)
# pca.fit(brain)
# brain_pca = pca.transform(brain)
# brain_pca = pd.DataFrame({"eid": brain.index, "brain_pca1": brain_pca[:, 0]})
# brain_pca.set_index("eid")

# pca = PCA(n_components=1)
# pca.fit(kidney)
# kidney_pca = pca.transform(kidney)
# kidney_pca = pd.DataFrame({"eid": kidney.index, "kidney_pca1": kidney_pca[:, 0]})
# kidney_pca.set_index("eid")

# pca = PCA(n_components=1)
# pca.fit(liver)
# liver_pca = pca.transform(liver)
# liver_pca = pd.DataFrame({"eid": liver.index, "liver_pca1": liver_pca[:, 0]})
# liver_pca.set_index("eid")

# pca = PCA(n_components=1)
# pca.fit(lung)
# lung_pca = pca.transform(lung)
# lung_pca = pd.DataFrame({"eid": lung.index, "lung_pca1": lung_pca[:, 0]})
# lung_pca.set_index("eid")

# pca = PCA(n_components=1)
# pca.fit(pancreas)
# pancreas_pca = pca.transform(pancreas)
# pancreas_pca = pd.DataFrame({"eid": pancreas.index, "pancreas_pca1": pancreas_pca[:, 0]})
# pancreas_pca.set_index("eid")

# pca = PCA(n_components=1)
# pca.fit(spleen)
# spleen_pca = pca.transform(spleen)
# spleen_pca = pd.DataFrame({"eid": spleen.index, "spleen_pca1": spleen_pca[:, 0]})
# spleen_pca.set_index("eid")

# pca_organs = [brain_pca, heart_pca, lung_pca, liver_pca, kidney_pca, pancreas_pca, spleen_pca]
# df_final = pd.DataFrame(columns=["eid"])
# for pca_organ in pca_organs:
#     df_final = df_final.merge(pca_organ, on="eid", how="outer")
# df_final = df_final.set_index("eid")
# df_final.to_csv("data/organ_pca.csv", index=True)

df_new_pheno = df_final
df_new_pheno.index.name = "s"
# df_new_pheno["y_res_heart_brain"] = abs(df_new_pheno["y_res_heart_brain"])

df_disease = pd.read_csv('data/disease_with_pn.csv', index_col=0)
df_disease = df_disease.drop('End stage renal', axis=1)
df_disease = df_disease.drop('Motor neurone', axis=1)
df_disease = df_disease.loc[:, ["Stroke", "Dementia", "Parkinson", "myocardial infarction", "hypertension", "nonrheumatic aortic valve disorders", "cardiomyopathy", 
                                "atrioventricular and left bundle-branch block", "Atrial_Fibrillation", "Heart_Failure", "other cardiac arrhythmias", 
                                "Asthma", "COPD", "Chronic kidney disease", "non-insulin-dependent diabetes mellitus"]]

df_covariate = pd.read_csv("data//covariates.tsv", index_col=0, sep="\t")


diseases = list(df_disease.columns)
phenotypes = list(df_new_pheno.columns)
phenotype_name = dict(zip(phenotypes, df_new_pheno.columns))
phenotype_unit = list(df_new_pheno.columns)


# For each disease of interest, perform logistic regression
#
# disease ~ Sex + Age + Weight + Height + Imaging_Phenotype
#
n_disease = len(diseases)
n_pheno = len(phenotypes)
df_beta = pd.DataFrame(index=diseases, columns=phenotypes, dtype=float)
df_conf = df_beta.copy()
df_p_val = df_beta.copy()
df_n = {}


for pheno in phenotypes:
    for disease in diseases:

        # Dependent variables
        y = df_disease[disease]

        # Independent variables
        # X = pd.concat([df_info[['Sex', 'Age']], df[[pheno]]], axis=1)
        X = pd.merge(df_covariate[['isMale', 'age', 'weight', 'height']], df_new_pheno[[pheno]], on="s", how="inner")

        # Standard deviation of the phenotype for normalisation
        pheno_std = df_new_pheno[pheno].std()

        # Remove NaN values before regression
        valid_idx = pd.merge(X, y, on="s", how="inner").dropna().index
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
        print('{0}, per {1:.1f} {2}: Odds = {3:.2f} [{4:.2f}, {5:.2f}], p = {6:.3f}'.format(
            disease, pheno_std, pheno,
            odds, odds_lower, odds_upper, df_p_val.loc[disease, pheno]))
        


# Bonferroni correction
p_bonf = 0.05 / (n_disease * n_pheno)

# FDR correction
p_fdr, _ = fdr_threshold(df_p_val.values.flatten(), 0.05)

# Number of phenotypes that is significantly associated with at least one of the IDPs
print('p_bonf = {0}'.format(p_bonf))
print('p_fdr = {0}'.format(p_fdr))

# Create the annotation text matrix
df_char = df_p_val.copy()
for disease in diseases:
    for pheno in phenotypes:
        if df_p_val.loc[disease, pheno] < p_bonf:
            df_char.loc[disease, pheno] = '**'
        elif df_p_val.loc[disease, pheno] < p_fdr:
            df_char.loc[disease, pheno] = '*'
        else:
            df_char.loc[disease, pheno] = ''








# Plot the heatmap for the odds ratios
plt.figure()
df_beta2 = df_beta.rename(columns=phenotype_name)
ax = sns.heatmap(np.exp(df_beta2), annot=df_char, annot_kws={'size': 13}, fmt='',
                 vmin=0.4, vmax=1.6, center=1, square=True,
                 cmap=sns.diverging_palette(240, 10, as_cmap=True),
                 cbar_kws={'shrink': .6, 'aspect': 15})
cbar_ax = ax.figure.get_axes()[-1]
cbar_ax.set_title('Odds\nratio', fontsize=13)
cbar_ax.set_yticklabels(cbar_ax.get_yticklabels(), fontsize=13)
ax.set_xticklabels(ax.get_xticklabels(), fontsize=13, rotation=30, ha="right")
ax.set_yticklabels(ax.get_yticklabels(), fontsize=13)
plt.gcf().set_size_inches(20, 8)
plt.tight_layout()
plt.savefig('plot/pca_plot_clinical_outcome_odds.png', bbox_inches='tight')

# Plot the heatmap for the logarithm of p-values
plt.figure()
df_p_val2 = df_p_val.rename(columns=phenotype_name)
ax = sns.heatmap(np.log(df_p_val2), annot=df_char, annot_kws={'size': 13}, fmt='',
                 vmin=-10, vmax=0, square=True,
                 cmap=sns.light_palette('#da3b46', reverse=True, as_cmap=True),
                 cbar_kws={"shrink": .6, 'aspect': 15})
cbar_ax = ax.figure.get_axes()[-1]
cbar_ax.set_title('log(p)', fontsize=13)
cbar_ax.set_yticklabels(cbar_ax.get_yticklabels(), fontsize=13)
ax.set_xticklabels(ax.get_xticklabels(), fontsize=13, rotation=30, ha="right")
ax.set_yticklabels(ax.get_yticklabels(), fontsize=13)
plt.gcf().set_size_inches(20, 8)
plt.tight_layout()
plt.savefig('plot/pca_plot_clinical_outcome_p_val.png', bbox_inches='tight')
# plt.show()

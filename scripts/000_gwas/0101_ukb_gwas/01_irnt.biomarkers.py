import hail as hl
import pandas as pd
import os

project_path = os.path.expanduser("/cpfs01/projects-HDD/cfff-6117e6302119_HDD/lm_21210880006/projects/indNet/")
hl.init(master='local[10]')

ht = hl.read_table(project_path + 'data/processed/biomarker/pipeline_pheno.ht')

# df = ht.to_pandas()

# df.index = df['s']
# df = df.drop('s', axis=1)

# dfp = df.rank()
# dfp = (dfp - 0.5) / (~dfp.isnull()).sum()
# dfp.columns = [x + '_prob' for x in dfp.columns]
# df.columns = [x + '_raw' for x in df.columns]

# df = pd.merge(df, dfp, how='inner', left_index=True, right_index=True)
# df.loc[:, 's'] = df.index

# ht = hl.Table.from_pandas(df, key='s')
# ht = ht.annotate(**{x.replace('_prob', '_irnt'): hl.qnorm(ht[x])
#                     for x in ht.row_value if x.endswith('_prob')})
# ht = ht.annotate(**{x: hl.or_missing(~hl.is_nan(ht[x]), ht[x]) for x in ht.row_value})
# ht = ht.drop(*[x for x in ht.row_value if x.endswith('_prob')])

ht.write(project_path + 'data/processed/biomarker/biomarkers_gwas.ht', overwrite=True)



ht = hl.read_table(project_path + 'data/processed/biomarker/pipeline_disease.ht')
ht.write(project_path + 'data/processed/biomarker/disease_gwas.ht', overwrite=True)


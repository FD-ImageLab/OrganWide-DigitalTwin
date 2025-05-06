import pandas as pd
from scripts.utils.func import draw_mht


data = pd.read_excel("data/mediation/mht.xlsx", sheet_name=None)

for fid, df in data.items():
    figname = "plot/med/mht/" + fid + ".png"
    draw_mht(df, figname=figname)
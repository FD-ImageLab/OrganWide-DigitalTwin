import pandas as pd
import re
import os

# main
# user-defined path
mainPath = "E:/Projects/indNet/"  # local


# feature_df = pd.read_csv("Sandbox/Category/feature_sly.csv")
# data = pd.read_csv("data/mediation/sandbox_med_filtered.csv")
# data = data[(data["cell_category"] != "Serum") & (data["cell_category"] != "IF&MF")]
# data["abs_prop.mediated"] = abs(data["prop.mediated"])
# data = data[data["x"].isin(feature_df[feature_df['Select']]['DispName'].tolist())]


# # feature_dict = feature_df.set_index("Short_name")["DispName"].to_dict()
# # data["x"] = data["x"].map(feature_dict)

# data = data[data["cell_category"] != "IF&MF"]

# data["x"] = data["y"].str.split("_").str[-1]
# data["y"] = data["y"].str.split("_").str[-2]

# rev_data = data.copy()
# rev_data["x"] = data["y"]
# rev_data["y"] = data["x"]
# data = pd.concat([data, rev_data], axis=0)


# a = data.groupby(["y", "cell_category"])
# df = pd.DataFrame()
# for x, y in a:
#     y_fix = x[0]
#     cell_cate = x[1]
#     category_num = len(y)
#     mean_prop = y["abs_prop.mediated"].mean()
#     temp_df = pd.DataFrame({"x_fix": [pd.NA], "y_fix": [y_fix], "cell_category": [cell_cate],
#                             "mean_prop.mediated": [mean_prop], "category_num": [category_num]})
#     df = pd.concat([df, temp_df], axis=0)

# b = data.groupby(["x", "cell_category"])
# for x, y in b:
#     x_fix = x[0]
#     cell_cate = x[1]
#     category_num = len(y)
#     mean_prop = y["abs_prop.mediated"].mean()
#     temp_df = pd.DataFrame({"x_fix": [x_fix], "y_fix": [pd.NA], "cell_category": [cell_cate],
#                             "mean_prop.mediated": [mean_prop], "category_num": [category_num]})
#     df = pd.concat([df, temp_df], axis=0)
# df.to_csv("data/sfig13/organ_cell_organ_mean.csv", index=False)











feature_df = pd.read_csv("Sandbox/Category/feature_sly.csv")
whole_data = pd.read_csv("data/mediation/sandbox_med_filtered.csv")

for feature in whole_data["x"].unique():
    print(feature)
    data = whole_data[whole_data["x"] == feature]
    data = data[(data["cell_category"] != "Serum") & (data["cell_category"] != "IF&MF")]
    data["abs_prop.mediated"] = abs(data["prop.mediated"])
    data = data[data["x"].isin(feature_df[feature_df['Select']]['DispName'].tolist())]


    # feature_dict = feature_df.set_index("Short_name")["DispName"].to_dict()
    # data["x"] = data["x"].map(feature_dict)

    data = data[data["cell_category"] != "IF&MF"]

    data["x"] = data["y"].str.split("_").str[-1]
    data["y"] = data["y"].str.split("_").str[-2]

    rev_data = data.copy()
    rev_data["x"] = data["y"]
    rev_data["y"] = data["x"]
    data = pd.concat([data, rev_data], axis=0)


    a = data.groupby(["y", "cell_category"])
    df = pd.DataFrame()
    for x, y in a:
        y_fix = x[0]
        cell_cate = x[1]
        category_num = len(y)
        mean_prop = y["abs_prop.mediated"].mean()
        temp_df = pd.DataFrame({"x_fix": [pd.NA], "y_fix": [y_fix], "cell_category": [cell_cate],
                                "mean_prop.mediated": [mean_prop], "category_num": [category_num]})
        df = pd.concat([df, temp_df], axis=0)

    b = data.groupby(["x", "cell_category"])
    for x, y in b:
        x_fix = x[0]
        cell_cate = x[1]
        category_num = len(y)
        mean_prop = y["abs_prop.mediated"].mean()
        temp_df = pd.DataFrame({"x_fix": [x_fix], "y_fix": [pd.NA], "cell_category": [cell_cate],
                                "mean_prop.mediated": [mean_prop], "category_num": [category_num]})
        df = pd.concat([df, temp_df], axis=0)
    if not os.path.exists("data/sfig13/feature/" + feature):
        os.makedirs("data/sfig13/feature/" + feature)
    df.to_csv("data/sfig13/feature/" + feature + "/organ_cell_organ_mean.csv", index=False)

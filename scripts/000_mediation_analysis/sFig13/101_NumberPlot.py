import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import math


# def plot_organ_number(data,color,figure_size,text_size,order):
#     #构建画布大小
#     plt.rcParams["figure.figsize"] = figure_size
#
#     #读入数据
#     organ = data[['individual']]
#     for i in organ.columns:
#         new_organ=organ[i].values.tolist()
#     organ=new_organ
#     organ.reverse()
#
#     number = data[['value']]
#     for i in number.columns:
#         new_number=number[i].values.tolist()
#     new_number.reverse()
#     number = new_number
#     number=list(number)
#     number=[round(x,2) for x in number]
#
#
#     color.reverse()
#     #画图
#     fig, ax = plt.subplots()
#     plt.tick_params(labelsize=text_size)
#     y=[i+1 for i in range(len(organ))]
#     ax.barh(y=y,width=number,color=color,zorder=1)
#
#     for i in range(len(organ)):
#         if number[i] != 0:
#             ax.text(x=number[i]/2,y=y[i],s=str(number[i]),verticalalignment='center', horizontalalignment="left",
#                     zorder=order,size=text_size, weight='bold')
#
#     ax.set_yticks(y)
#     ax.set_yticklabels(organ,verticalalignment='top',size=text_size,weight='bold',rotation=45)
#     plt.xticks(size=text_size, weight="bold")
#     ax.set_xticks([0, round(max(number) / (2 * 0.05)) * 0.05, math.ceil(max(number) / 0.05) * 0.05])
#     plt.tight_layout()
#main
if __name__ == '__main__':
    #main
    # user-defined path
    mainPath = "./"  # local
    
    # mainPath = "E:/projects/indNet/"  # local
    # mainPath = "/Users/apple/Seafile/Research/SandBox/"  # local
    # mainPath = "D:/Seafile/Seafile/SandBox/"
    # mainPath = "/public/sandbox/workdir/liumeng/SandBox"  # sandbox
    # mainPath = "/public/sandbox/workdir/wangchy/SharedFolder"  # sandbox
    codePath = mainPath+"/Code/R/MediationBarPlot/"
    resultPath = mainPath + "/Results/MediationBarPlot/"
    dataPath = mainPath+"/data/sfig13/100/Organ/"

    import sys
    sys.path.append(mainPath)
    from scripts.utils.func import plot_organ_number, bar_number

    color_map = pd.read_csv(mainPath + "/data/color_map/micro_cell_color_map.csv")

    file_names = os.listdir(dataPath)

    maximum = 0
    for file in file_names:
        data = pd.read_csv(os.path.join(dataPath, file), encoding="UTF-8")
        data["individual"] = data["individual"].str.replace("WB_Events_", "")
        data["individual"] = data["individual"].str.replace("WB_MFI_", "")
        tm = max([i for i in map(len, data['individual'])])
        maximum = max(maximum, tm)

    for file in file_names:
        data = pd.read_csv(os.path.join(dataPath, file), encoding="UTF-8")
        data["value"] = abs(data["value"])
        # delete WB_Events prefix
        data["individual"] = data["individual"].str.replace("WB_Events_", "")
        data["individual"] = data["individual"].str.replace("WB_MFI_", "")
        def f(string, length=20):
            if len(string) > length:
                string = string[:length]
            return string
        data["individual"] = [i for i in map(f, data["individual"])]
        # data["significant"] = data["significant"].replace(False, "")
        # data["significant"] = data["significant"].replace(True, "$ *$")
        # data['individual'] = data['individual'] + data['significant'].map(str)
        color = []
        for i, j in enumerate(data["group"]):
            if j in list(color_map["type"]):
                color.append(color_map["color"][color_map["type"] == j].tolist()[0])
            else:
                color.append("white")
        figure_size=[4,6]
        text_size=40
        order=2##whether to display number，2 for display，0 for not to display
        fig=plot_organ_number(data,color,figure_size,text_size,order, number_show=True, ylab_len=False, height=0.8, edgecolor="black", number_position=0.04)
        if not os.path.exists(os.path.join(mainPath + "plot/sfig13/100/")):
            os.makedirs(os.path.join(mainPath + "plot/sfig13/100/"))
        plt.savefig(os.path.join(os.path.join(mainPath + "plot/sfig13/100/"),
                                    os.path.splitext(file)[0] + ".png"),
                    dpi=600)





#     maximum = 0
#     for dataPath in os.listdir(mainPath + "/data/sfig13/100/Organ"):
#         data = pd.read_csv(mainPath + "/data/sfig13/100/Organ" + dataPath, encoding="UTF-8")
#         data["individual"] = data["individual"].str.replace("WB_Events_", "")
#         data["individual"] = data["individual"].str.replace("WB_MFI_", "")
#         tm = max([i for i in map(len, data['individual'])])
#         maximum = max(maximum, tm)

#     for dataPath in os.listdir(mainPath + "/data/sfig13/100/Organ"):
#         organ = os.path.splitext(dataPath)[0]
#         data = pd.read_csv(mainPath + "/data/sfig13/100/Organ" + dataPath, encoding="UTF-8")
#         data["individual"] = data["individual"].str.replace("WB_Events_", "")
#         data["individual"] = data["individual"].str.replace("WB_MFI_", "")
#         def f(string, length=40):
#             if len(string) > length:
#                 string = string[:length]
#             return string
#         data["individual"] = [i for i in map(f, data["individual"])]
        
#         color = []
#         for i, j in enumerate(data["group"]):
#             if j in list(color_map["type"]):
#                 color.append(color_map["color"][color_map["type"] == j].tolist()[0])
#             else:
#                 color.append("white")
#         figure_size = [9, 6]
#         text_size = 30
#         order = 2  ##数字是否显示，2为显示，0为不显示
#         fig = bar_number(data, color, figure_size, text_size, order, number_show=True, ylab_len=False, height=1, edgecolor="black", number_position=0.04)
#         if not os.path.exists(os.path.join(mainPath + "plot/sfig13/100/")):
#             os.makedirs(os.path.join(mainPath + "plot/sfig13/100/"))
#         plt.savefig(os.path.join(os.path.join(mainPath + "plot/sfig13/100/"), organ + ".png"),
#                     dpi=600)



# # import numpy as np
# # from matplotlib import pyplot as plt
# # x = np.linspace(1, 100, 1000)
# # y = np.log(x)
# # y1 = np.sin(x)
# # fig = plt.figure("Line plot")
# # legendFig = plt.figure("Legend plot")
# # ax = fig.add_subplot(111)
# # line1, = ax.plot(x, y, c="red", lw=10)
# # line2, = ax.plot(x, y1, c="green", lw=10)
# # legendFig.legend([line1, line2], ["y=log(x)", "y=sin(x)"], loc='center')
# # legendFig.savefig('legend.png')
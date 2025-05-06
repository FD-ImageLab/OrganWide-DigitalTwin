# -*- coding = utf-8 -*-
# @File_name = concatenate
import matplotlib.text
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os

def check_path(filepath):
    # check for file
    stat = os.path.exists(os.path.dirname(filepath))
    if not stat:
        os.makedirs(os.path.dirname(filepath))
        print("Path created: " + os.path.dirname(filepath))
    else:
        print("Path exists: " + os.path.dirname(filepath))


def concate_png(paths, cshape=[1, -1]):
    import os
    import cv2
    import numpy as np

    if cshape[1] == -1:
        cshape[1] = len(paths)
    concate_img = np.empty(0)
    for i in range(cshape[0]):
        if all([os.path.splitext(path)[1] == ".png" for path in paths]):
            img = [cv2.imread(path, -1) for path in paths[cshape[1] * i: cshape[1] * i + cshape[1]]]
            concate_row_img = np.concatenate(([i for i in img]), axis=1)
            if len(concate_img) == 0:
                concate_img = concate_row_img
            else:
                if concate_row_img.shape[1] != concate_img.shape[1]:
                    fill_img = np.full([concate_row_img.shape[0],
                                         concate_img.shape[1] - concate_row_img.shape[1],
                                         concate_row_img.shape[2]], 255)
                    concate_row_img = np.concatenate([concate_row_img, fill_img], axis=1)
                concate_img = np.concatenate([concate_img, concate_row_img], axis=0)
    return concate_img

# 文字放在y轴上
# def plot_organ_number(data, color, figure_size, text_size, order=0, number_show=False, ylab_len=16):
#     # 构建画布大小
#     import matplotlib.pyplot as plt
#     import math
#     plt.rcParams["figure.figsize"] = figure_size
#     plt.rc('font', family='Arial')
#
#     # 读入数据
#
#     # 调整y轴文字长度
#     if ylab_len is not False:
#         def f(s):
#             s = s.rjust(ylab_len, " ")
#             return s
#         data['individual'] = [i for i in map(f, data['individual'])]
#     organ = data[['individual']]
#     for i in organ.columns:
#         new_organ = organ[i].values.tolist()
#     organ = new_organ
#     organ.reverse()
#
#     number = data[['value']]
#     for i in number.columns:
#         new_number = number[i].values.tolist()
#     new_number.reverse()
#     number = new_number
#     number = list(number)
#     number = [round(x, 2) for x in number]
#
#     color.reverse()
#     # 画图
#     fig, ax = plt.subplots()
#     plt.tick_params(labelsize=text_size)
#     y = [i + 1 for i in range(len(organ))]
#     ax.barh(y=y, width=number, color=color, zorder=1)
#
#     if number_show:
#         for i in range(len(organ)):
#             if number[i] != 0:
#                 ax.text(x=number[i] / 2, y=y[i], s=str(number[i]), verticalalignment='center', horizontalalignment="center",
#                         zorder=order, size=text_size, weight='bold')
#
#     ax.set_yticks(y)
#     ax.set_yticklabels(organ, verticalalignment='top', size=text_size, weight='bold', rotation=35)
#     # for label in ax.get_yticklabels():
#     #     label.set_y(label.get_position()[1] + 2)
#     # [ax.get_yticklabels()[i].update({"y": ax.get_yticklabels()[i].get_position()[1] + 2}) for i in range(len(ax.get_yticklabels()))]
#     # label = ax.get_yticklabels()
#     # ax.set_yticklabels(label)
#     # matplotlib.text.Text
#     plt.xticks(size=text_size, weight="bold")
#     if round(max(number) / (2 * 0.05)) * 0.05 != 0:
#         ax.set_xticks([0, round(max(number) / (2 * 0.05)) * 0.05, math.ceil(max(number) / 0.05) * 0.05])
#         ax.set_xticklabels([0,
#                        "{:.2f}".format(round(max(number) / (2 * 0.05)) * 0.05),
#                        "{:.2f}".format(math.ceil(max(number) / 0.05) * 0.05)])
#     else:
#         ax.set_xticks([0, math.ceil(max(number) / 0.05) * 0.05])
#         ax.set_xticklabels([0,
#                             "{:.2f}".format(math.ceil(max(number) / 0.05) * 0.05)])
#     plt.tight_layout()

def plot_organ_number(data, color, figure_size, text_size, order=0, number_show=False, ylab_len=False, number_position=0.04, height=0.8, edgecolor=None):
    import matplotlib.pyplot as plt
    import math
    plt.rcParams["figure.figsize"] = figure_size
    plt.rc('font', family='Arial')

    if ylab_len is not False:
        def f(s):
            s = s.rjust(ylab_len, " ")
            return s
        data['individual'] = [i for i in map(f, data['individual'])]
    organ = data[['individual']]
    for i in organ.columns:
        new_organ = organ[i].values.tolist()
    organ = new_organ
    organ.reverse()

    number = data[['value']]
    for i in number.columns:
        new_number = number[i].values.tolist()
    new_number.reverse()
    number = new_number
    number = list(number)
    number = [round(x, 2) for x in number]

    color.reverse()
    fig, ax = plt.subplots()
    plt.tick_params(labelsize=text_size)
    y = [i + 1 for i in range(len(organ))]
    ax.barh(y=y, width=number, color=color, zorder=1, height=height, edgecolor=edgecolor)

    if number_show:
        for i in range(len(organ)):
            if number[i] != 0:
                ax.text(x=max(number) * number_position, y=y[i], s=str(organ[i]), verticalalignment='center', horizontalalignment="left",
                        zorder=order, size=text_size) # , weight='bold'

    ax.set_yticks([])
    # ax.set_yticklabels(organ, verticalalignment='top', size=text_size, weight='bold', rotation=35)
    # for label in ax.get_yticklabels():
    #     label.set_y(label.get_position()[1] + 2)
    # [ax.get_yticklabels()[i].update({"y": ax.get_yticklabels()[i].get_position()[1] + 2}) for i in range(len(ax.get_yticklabels()))]
    # label = ax.get_yticklabels()
    # ax.set_yticklabels(label)
    # matplotlib.text.Text
    plt.xticks(size=text_size) # , weight="bold"
    if round(max(number) / (2 * 0.05)) * 0.05 != 0:
        ax.set_xticks([0, round(max(number) / (2 * 0.05)) * 0.05, math.ceil(max(number) / 0.05) * 0.05])
        ax.set_xticklabels([0,
                       "{:.2f}".format(round(max(number) / (2 * 0.05)) * 0.05),
                       "{:.2f}".format(math.ceil(max(number) / 0.05) * 0.05)])
    else:
        ax.set_xticks([0, math.ceil(max(number) / 0.05) * 0.05])
        ax.set_xticklabels([0,
                            "{:.2f}".format(math.ceil(max(number) / 0.05) * 0.05)])
    plt.tight_layout()


# mediaiton_bar_plot and moderation_bar_plot num_barplot the text is on the yaxis
# def bar_number(data, color, figure_size, text_size, order=0, number_show=False, ylab_len=16):
#         # 构建画布大小
#         import matplotlib.pyplot as plt
#         import math
#         plt.rcParams["figure.figsize"] = figure_size
#         plt.rc('font', family='Arial')
#         # 读入数据
#         # 调整y轴文字长度
#         if ylab_len is not False:
#             def f(s):
#                 s = s.rjust(ylab_len, " ")
#                 return s
#             data['individual'] = [i for i in map(f, data['individual'])]
#         organ = data[['individual']]
#         for i in organ.columns:
#             new_organ = organ[i].values.tolist()
#         organ = new_organ
#         organ.reverse()
#
#         number = data[['value']]
#         for i in number.columns:
#             new_number = number[i].values.tolist()
#         new_number.reverse()
#         number = new_number
#         number = list(number)
#         number = [round(x, 2) for x in number]
#
#         color.reverse()
#         # 画图
#         fig, ax = plt.subplots()
#         plt.tick_params(labelsize=text_size)
#         y = [i + 1 for i in range(len(organ))]
#         ax.barh(y=y, width=number, color=color, zorder=1)
#
#         if number_show:
#             for i in range(len(organ)):
#                 if number[i] != 0:
#                     ax.text(x=number[i] / 2, y=y[i], s=str(number[i]), verticalalignment='center',
#                             horizontalalignment="center",
#                             zorder=order, size=text_size) # , weight='bold'
#
#         ax.set_yticks(y)
#         ax.set_yticklabels(organ, verticalalignment='top', size=text_size, weight='bold', rotation=35)
#         plt.xticks(size=text_size) # , weight="bold"
#         if round(max(number) / (2 * 0.05)) * 0.05 != 0:
#             ax.set_xticks([0, round(max(number) / 2), math.ceil(max(number))])
#             ax.set_xticklabels([0,
#                                 "{:.0f}".format(round(max(number) / 2)),
#                                 "{:.0f}".format(math.ceil(max(number)))])
#         else:
#             ax.set_xticks([0, math.ceil(max(number))])
#             ax.set_xticklabels([0,
#                                 "{:.0f}".format(math.ceil(max(number)))])
#         plt.tight_layout()


# mediaiton_bar_plot and moderation_bar_plot num_barplot the text is inside the bar
def bar_number(data, color, figure_size, text_size, order=0, number_show=False, ylab_len=16, number_position=0.04, height=0.8, edgecolor=None):
        # 构建画布大小
        import matplotlib.pyplot as plt
        import math
        plt.rcParams["figure.figsize"] = figure_size
        plt.rc('font', family='Arial')
        # 读入数据
        # 调整y轴文字长度
        if ylab_len is not False:
            def f(s):
                s = s.rjust(ylab_len, " ")
                return s
            data['individual'] = [i for i in map(f, data['individual'])]
        organ = data[['individual']]
        for i in organ.columns:
            new_organ = organ[i].values.tolist()
        organ = new_organ
        organ.reverse()

        number = data[['value']]
        for i in number.columns:
            new_number = number[i].values.tolist()
        new_number.reverse()
        number = new_number
        number = list(number)
        number = [round(x, 2) for x in number]

        color.reverse()
        # 画图
        fig, ax = plt.subplots()
        plt.tick_params(labelsize=text_size)
        y = [i + 1 for i in range(len(organ))]
        ax.barh(y=y, width=number, color=color, zorder=1, height=height, edgecolor=edgecolor)

        if number_show:
            for i in range(len(organ)):
                if number[i] != 0:
                    ax.text(x=max(number) * number_position, y=y[i], s=str(organ[i]), verticalalignment='center',
                            horizontalalignment="left",
                            zorder=order, size=text_size) # , weight='bold'

        ax.set_yticks([])
        plt.xticks(size=text_size) # , weight="bold"
        if round(max(number) / (2 * 0.05)) * 0.05 != 0:
            ax.set_xticks([0, round(max(number) / 2), math.ceil(max(number))])
            ax.set_xticklabels([0,
                                "{:.0f}".format(round(max(number) / 2)),
                                "{:.0f}".format(math.ceil(max(number)))])
        else:
            ax.set_xticks([0, math.ceil(max(number))])
            ax.set_xticklabels([0,
                                "{:.0f}".format(math.ceil(max(number)))])
        plt.tight_layout()


def grep(pattern, *args, return_info_flag=True, flags=0):
    import re
    return_list = []
    for string in args:
        if re.search(pattern, string, flags=flags):
            return_list += [string]
    if not return_list:
        if return_info_flag:
            print("This pattern can't match any item, please check your input!")
            return
    return return_list


def grepl(pattern, *args, return_info_flag=True, flags=0):
    import re
    return_list = []
    for string in args:
        if re.search(pattern, string, flags=flags):
            return_list += [True]
        else: 
            return_list += [False]
    if not return_list:
        if return_info_flag:
            print("This pattern can't match any item, please check your input!")
            return
    return return_list


def draw_mht(df, figname):
    mht = plt.figure()
    
    color_map = {"NMR metabolomics": "red", "Urine assays": "blue", "Blood chemistry": "orange"}
    sns.scatterplot(x="BP", y="prop.mediated", hue="molecular_category", data=df, alpha=1, legend="full", palette=color_map)
    
    first_occurrences = df.drop_duplicates(subset="y")
    x_before = np.array((first_occurrences["BP"] - 1))
    x = np.array((first_occurrences["BP"] - 1)[1:])
    x_after = np.array(list(x) + [len(df["BP"])])
    xticks = [i for i in (x_before + x_after) / 2]


    for i in x:
        plt.axvline(x=i, lw=0.5, c='black', ls='--', dashes=(5, 5))
    plt.axhline(y=0, lw=0.5, c='black', ls='--', dashes=(5, 5))
    plt.axvline(x=0, lw=0.5, c='black', ls='--', dashes=(5, 5))


    xticknames = df.y.unique()
    xticklabels = xticknames
    plt.xticks(xticks, xticklabels, fontsize=10, rotation=45, )

    plt.xlabel("")
    plt.ylabel("")
    plt.title("")
    plt.legend(loc="upper right")
    
    fig = plt.gcf()
    fig.set_size_inches(16, 8)
    plt.tight_layout()
    check_path(figname)
    plt.savefig(figname, bbox_inches='tight', dpi=300)


if __name__ == "__main__":
    import os
    import cv2
    
    check_path("plot/med/mht/1.png")
    # main
    # user-defined path
    mainPath = "E:/seafile/Seafile/SandBox/"  # local
    # mainPath = "/Users/apple/Seafile/Research/SandBox/"  # local
    # mainPath = "D:/Seafile/Seafile/SandBox/"
    # mainPath = "/public/sandbox/workdir/liumeng/SandBox"  # sandbox
    # mainPath = "/public/sandbox/workdir/wangchy/SharedFolder"  # sandbox
    codePath = mainPath + "/Code/R/GenderBarPlot/"
    resultPath = mainPath + "/Results/GenderBarPlot/"

    path = os.path.join(resultPath, "Organ")
    positions = [os.path.join(path, file_name) for file_name in os.listdir(path)]
    png = concate_png(positions, [2, 5])
    cv2.imwrite("temp.png", png)

    import re
    li = tuple(["A"+ str(i) for i in range(100)])
    print(grep("2", *li))

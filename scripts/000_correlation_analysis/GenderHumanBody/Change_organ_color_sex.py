import os
from PIL import Image
import numpy as np
import pandas as pd
import cairosvg
import xml.etree.ElementTree as ET


class FigureDrawer:
    def __init__(self, output_path: str):
        """Draw svg color in different parts
        """
        self.man_organ_path = 'SVG/man_organ.svg'
        self.woman_organ_path = 'SVG/woman_organ.svg'
        self.mixgender_organ_path = 'SVG/mixgender_organ.svg'
        self.colorbar_path = 'colorbar.png'
        self.colorbar_jet_path = 'jet.png'
        self.output_path = output_path

    def read_xml(self, in_path):
        """读取并解析xml文件

        Args:
            in_path (str): xml路径

        Returns:
            ElementTree
        """
        # 先要注册命名空间，不然会'ElementTree' object has no attribute 'register_namespace'
        ET.register_namespace("", "http://www.w3.org/2000/svg")
        ET.register_namespace("xlink", "http://www.w3.org/1999/xlink")
        ET.register_namespace("space", "preserve")
        tree = ET.ElementTree()
        tree.parse(in_path)
        return tree

    def write_xml(self, tree, out_path):
        """将xml文件写出

        Args:
            tree : xml树
            out_path (str): 写出路径
        """
        tree.write(out_path, encoding="utf-8", xml_declaration=True, method='xml')

    def RGB_to_Hex(self, rgb):
        """RGB格式颜色转换为16进制颜色格式颜色
        """
        color = '#'
        for i in rgb:
            num = int(i)
            # 将R、G、B分别转化为16进制拼接转换并大写  hex() 函数用于将10进制整数转换成16进制，以字符串形式表示
            color += str(hex(num))[-2:].replace('x', '0').upper()
        # print(color)
        return color

    def Hex_to_RGB(self, hex):
        """16进制颜色格式颜色转换为RGB格式
        """
        r = int(hex[1:3], 16)
        g = int(hex[3:5], 16)
        b = int(hex[5:7], 16)
        return (r, g, b)

    def rgb_show(self, rgb1, rgb2, ratio):
        rgb = [0, 0, 0]
        for i in range(len(rgb1)):
            if rgb1[i] < rgb2[i]:
                rgb[i] = int((rgb2[i] - rgb1[i]) * ratio) + rgb1[i]
            elif rgb1[i] > rgb2[i]:
                rgb[i] = int((rgb1[i] - rgb2[i]) * ratio) + rgb2[i]
            else:
                rgb[i] = rgb2[i]
        return (rgb[0], rgb[1], rgb[2])

    def get_ratio(self, num, max_num):
        return num / max_num

    def find_rgb(self, num, colorbar_name):
        if colorbar_name == 'jet':
            colorbar_path = self.colorbar_jet_path
        else:
            colorbar_path = self.colorbar_path
        colorbar = Image.open(colorbar_path)
        ratio = self.get_ratio(num, 1)
        if ratio > 1:
            ratio = 1
        shape = colorbar.size
        if ratio == 0:
            color = colorbar.getpixel((0, int(shape[1] / 2)))
        else:
            color = colorbar.getpixel((int(shape[0] * ratio) - 1, int(shape[1] / 2)))
        color = (color[0], color[1], color[2])
        return color

    def get_Hex(self, num, colorbar_name='1'):
        rgb = self.find_rgb(num, colorbar_name)
        hex = self.RGB_to_Hex(rgb)
        return hex

    def get_Hex_color_set(self, num):
        color_set = ['#d3311f', '#bea589', '#5f803f', '#b8802b', '#5360a0', '#d9712a', '#327487', '#387e24', '#dd7094']
        if num < 1 / 9:
            hex = color_set[0]
        elif num >= 1 / 9 and num < 2 / 9:
            hex = color_set[1]
        elif num >= 2 / 9 and num < 3 / 9:
            hex = color_set[2]
        elif num >= 3 / 9 and num < 4 / 9:
            hex = color_set[3]
        elif num >= 4 / 9 and num < 5 / 9:
            hex = color_set[4]
        elif num >= 5 / 9 and num < 6 / 9:
            hex = color_set[5]
        elif num >= 6 / 9 and num < 7 / 9:
            hex = color_set[6]
        elif num >= 7 / 9 and num < 8 / 9:
            hex = color_set[7]
        elif num >= 8 / 9:
            hex = color_set[8]
        return hex

    def get_norm_val(self, val):
        "input should be list"
        val_array = np.array(val)
        #        max_val = max(val_array)
        #        min_val = min(val_array)
        max_val = 1.2
        min_val = 0.4

        eps = np.finfo(min_val).eps
        val_norm_array = np.maximum((val_array - min_val) / np.maximum(max_val - min_val, eps), 0)
        val_norm = val_norm_array.tolist()
        return val_norm

    def change_man_color(self, svg_data: list, do_norm: bool = False, colorbar_name: str = '1'):
        """修改svg中颜色并保存
        svg_data: 数据list
        do_norm: {True, False} 是否对数据进行归一化
        colorbar_name: {'1', 'jet', 'color_set'}
                        '1': 使用默认colorbar
                        'jet': 使用jet作为colorbar
                        'color_set': 使用离散colorbar
        do_border: {True, False} 是否使用加边框的svg
        """
        # 计算各个部位的HEX颜色
        if do_norm:
            svg_data = self.get_norm_val(svg_data)
        if colorbar_name == 'color_set':
            hex_left_brain = self.get_Hex_color_set(svg_data[0])
            hex_right_brain = self.get_Hex_color_set(svg_data[1])
            hex_heart = self.get_Hex_color_set(svg_data[2])
            hex_left_lung = self.get_Hex_color_set(svg_data[3])
            hex_right_lung = self.get_Hex_color_set(svg_data[4])
            hex_liver = self.get_Hex_color_set(svg_data[5])
            hex_spleen = self.get_Hex_color_set(svg_data[6])
            hex_pancreas = self.get_Hex_color_set(svg_data[7])
            hex_left_kidney = self.get_Hex_color_set(svg_data[8])
            hex_right_kidney = self.get_Hex_color_set(svg_data[9])
            # hex_prostate = self.get_Hex_color_set(svg_data[10])
            hex_prostate = "#D3D3D3"
        else:
            hex_left_brain = self.get_Hex(svg_data[0], colorbar_name)
            hex_right_brain = self.get_Hex(svg_data[1], colorbar_name)
            hex_heart = self.get_Hex(svg_data[2], colorbar_name)
            hex_left_lung = self.get_Hex(svg_data[3], colorbar_name)
            hex_right_lung = self.get_Hex(svg_data[4], colorbar_name)
            hex_liver = self.get_Hex(svg_data[5], colorbar_name)
            hex_spleen = self.get_Hex(svg_data[6], colorbar_name)
            hex_pancreas = self.get_Hex(svg_data[7], colorbar_name)
            hex_left_kidney = self.get_Hex(svg_data[8], colorbar_name)
            hex_right_kidney = self.get_Hex(svg_data[9], colorbar_name)
            hex_prostate = self.get_Hex(svg_data[10], colorbar_name)
            hex_prostate = "#D3D3D3"

        # 修改svg代码来显示不同效果
        tree = self.read_xml(self.man_organ_path)
        new_path = self.output_path
        for element in tree.iter("{http://www.w3.org/2000/svg}path"):
            id_name = element.get("id")
            # 给不同id的部分加颜色
            if id_name == "right_lung":
                element.set("style", f"fill:{hex_right_lung}")
            elif id_name == "left_lung":
                element.set("style", f"fill:{hex_left_lung}")
            elif id_name == "heart":
                element.set("style", f"fill:{hex_heart}")
            elif id_name == "liver":
                element.set("style", f"fill:{hex_liver}")
            elif id_name == "spleen":
                element.set("style", f"fill:{hex_spleen}")
            elif id_name == "pancreas":
                element.set("style", f"fill:{hex_pancreas}")
            elif id_name == "right_kidney":
                element.set("style", f"fill:{hex_right_kidney}")
            elif id_name == "left_kidney":
                element.set("style", f"fill:{hex_left_kidney}")
            elif id_name == "man_genitals":
                element.set("style", f"fill:{hex_prostate}")
            elif id_name == "right_brain":
                element.set("style", f"fill:{hex_right_brain}")
            elif id_name == "left_brain":
                element.set("style", f"fill:{hex_left_brain}")
        self.write_xml(tree, new_path)
        print('man Success!')

    def change_woman_color(self, svg_data: list, do_norm: bool = False, colorbar_name: str = '1'):
        """修改svg中颜色并保存
        svg_data: 数据list
        do_norm: {True, False} 是否对数据进行归一化
        colorbar_name: {'1', 'jet', 'color_set'}
                        '1': 使用默认colorbar
                        'jet': 使用jet作为colorbar
                        'color_set': 使用离散colorbar
        do_border: {True, False} 是否使用加边框的svg
        """
        # 计算各个部位的HEX颜色
        if do_norm:
            svg_data = self.get_norm_val(svg_data)
        if colorbar_name == 'color_set':
            hex_left_brain = self.get_Hex_color_set(svg_data[0])
            hex_right_brain = self.get_Hex_color_set(svg_data[1])
            hex_heart = self.get_Hex_color_set(svg_data[2])
            hex_left_lung = self.get_Hex_color_set(svg_data[3])
            hex_right_lung = self.get_Hex_color_set(svg_data[4])
            hex_liver = self.get_Hex_color_set(svg_data[5])
            hex_spleen = self.get_Hex_color_set(svg_data[6])
            hex_pancreas = self.get_Hex_color_set(svg_data[7])
            hex_left_kidney = self.get_Hex_color_set(svg_data[8])
            hex_right_kidney = self.get_Hex_color_set(svg_data[9])
            # hex_uterus = self.get_Hex_color_set(svg_data[10])
            hex_uterus = "#D3D3D3"
        else:
            hex_left_brain = self.get_Hex(svg_data[0], colorbar_name)
            hex_right_brain = self.get_Hex(svg_data[1], colorbar_name)
            hex_heart = self.get_Hex(svg_data[2], colorbar_name)
            hex_left_lung = self.get_Hex(svg_data[3], colorbar_name)
            hex_right_lung = self.get_Hex(svg_data[4], colorbar_name)
            hex_liver = self.get_Hex(svg_data[5], colorbar_name)
            hex_spleen = self.get_Hex(svg_data[6], colorbar_name)
            hex_pancreas = self.get_Hex(svg_data[7], colorbar_name)
            hex_left_kidney = self.get_Hex(svg_data[8], colorbar_name)
            hex_right_kidney = self.get_Hex(svg_data[9], colorbar_name)
            # hex_uterus = self.get_Hex(svg_data[10], colorbar_name)
            hex_uterus = "#D3D3D3"

        # 修改svg代码来显示不同效果
        tree = self.read_xml(self.woman_organ_path)
        new_path = self.output_path

        for element in tree.iter("{http://www.w3.org/2000/svg}path"):
            id_name = element.get("id")
            # 给不同id的部分加颜色
            if id_name == "right_lung":
                element.set("style", f"fill:{hex_right_lung}")
            elif id_name == "left_lung":
                element.set("style", f"fill:{hex_left_lung}")
            elif id_name == "heart":
                element.set("style", f"fill:{hex_heart}")
            elif id_name == "liver":
                element.set("style", f"fill:{hex_liver}")
            elif id_name == "spleen":
                element.set("style", f"fill:{hex_spleen}")
            elif id_name == "pancreas":
                element.set("style", f"fill:{hex_pancreas}")
            elif id_name == "right_kidney":
                element.set("style", f"fill:{hex_right_kidney}")
            elif id_name == "left_kidney":
                element.set("style", f"fill:{hex_left_kidney}")
            elif id_name == "woman_genitals":
                element.set("style", f"fill:{hex_uterus}")
            elif id_name == "right_brain":
                element.set("style", f"fill:{hex_right_brain}")
            elif id_name == "left_brain":
                element.set("style", f"fill:{hex_left_brain}")
        self.write_xml(tree, new_path)
        print('woman Success!')


if __name__ == '__main__':

    # user-defined path
    mainPath = "/mnt/d/University/fdurop/remote_ws/indNet/"  # local

    dataPath = mainPath + "data/report_on_aging_and_gender/"
    os.chdir(mainPath + "scripts/000_correlation_analysis/GenderHumanBody/")

    do_norm = False
    do_norm2 = True
    do_border = True

    ##### 全身  ukb ######
     # # man (left_brain, right_brain, heart, left_lung, right_lung, liver, spleen, pancreas, left_kidney, right_kidney, prostate)
    output_path = mainPath + "plot/report_on_aging_and_gender/ukb/"
    data = pd.read_csv(dataPath + "ukb/male_humanbody.csv", index_col=0)
    
    # 调整性器官的值，使其不影响颜色
    data["value"][11] = min(data["value"])
    
     # save average data
    meta_data = data["value"].values
    tmpData2 = meta_data.tolist()
    svg_data = list(map(float, tmpData2))
    figure = FigureDrawer(output_path + 'man_MeanNorm.svg')
    figure.change_man_color(svg_data, do_norm=do_norm2)
    figure = FigureDrawer(output_path + '/' + 'man_Mean.svg')
    figure.change_man_color(svg_data)

    cairosvg.svg2png(url=output_path + 'man_MeanNorm.svg', write_to=output_path + 'man_MeanNorm.png', scale=10)
    cairosvg.svg2png(url=output_path + 'man_Mean.svg', write_to=output_path + 'man_Mean.png', scale=10)


     # #    # woman(left_brain, right_brain, heart, left_lung, right_lung, liver, spleen, pancreas, left_kidney, right_kidney, uterus)
    data = pd.read_csv(dataPath + "ukb/female_humanbody.csv")
    
    # 调整性器官的值，使其不影响颜色
    data["value"][11] = min(data["value"])
    
     # save average data
    meta_data = data["value"].values
    tmpData2 = meta_data.tolist()
    svg_data = list(map(float, tmpData2))
    figure = FigureDrawer(output_path + 'woman_MeanNorm.svg')
    figure.change_woman_color(svg_data, do_norm=do_norm2)
    figure = FigureDrawer(output_path + 'woman_Mean.svg')
    figure.change_woman_color(svg_data)
    
    
    cairosvg.svg2png(url=output_path + 'woman_MeanNorm.svg', write_to=output_path + 'woman_MeanNorm.png', scale=10)
    cairosvg.svg2png(url=output_path + 'woman_Mean.svg', write_to=output_path + 'woman_Mean.png', scale=10)
    
    
    
    
    
    
    
    ##### 全身  sandbox ######
     # # man (left_brain, right_brain, heart, left_lung, right_lung, liver, spleen, pancreas, left_kidney, right_kidney, prostate)
    output_path = mainPath + "plot/report_on_aging_and_gender/sandbox/"
    data = pd.read_csv(dataPath + "sandbox/male_humanbody.csv", index_col=0)
    
    # 调整性器官的值，使其不影响颜色
    data["value"][11] = min(data["value"])
    
     # save average data
    meta_data = data["value"].values
    tmpData2 = meta_data.tolist()
    svg_data = list(map(float, tmpData2))
    figure = FigureDrawer(output_path + 'man_MeanNorm.svg')
    figure.change_man_color(svg_data, do_norm=do_norm2)
    figure = FigureDrawer(output_path + '/' + 'man_Mean.svg')
    figure.change_man_color(svg_data)

    cairosvg.svg2png(url=output_path + 'man_MeanNorm.svg', write_to=output_path + 'man_MeanNorm.png', scale=10)
    cairosvg.svg2png(url=output_path + 'man_Mean.svg', write_to=output_path + 'man_Mean.png', scale=10)


     # #    # woman(left_brain, right_brain, heart, left_lung, right_lung, liver, spleen, pancreas, left_kidney, right_kidney, uterus)
    data = pd.read_csv(dataPath + "sandbox/female_humanbody.csv")
    
    # 调整性器官的值，使其不影响颜色
    data["value"][11] = min(data["value"])
     # save average data
    meta_data = data["value"].values
    tmpData2 = meta_data.tolist()
    svg_data = list(map(float, tmpData2))
    figure = FigureDrawer(output_path + 'woman_MeanNorm.svg')
    figure.change_woman_color(svg_data, do_norm=do_norm2)
    figure = FigureDrawer(output_path + 'woman_Mean.svg')
    figure.change_woman_color(svg_data)
    
    
    cairosvg.svg2png(url=output_path + 'woman_MeanNorm.svg', write_to=output_path + 'woman_MeanNorm.png', scale=10)
    cairosvg.svg2png(url=output_path + 'woman_Mean.svg', write_to=output_path + 'woman_Mean.png', scale=10)

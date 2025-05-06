# -*- coding = utf-8 -*-
# @File_name = concate_pic

# main
# user-defined path
mainPath = "E:/projects/indNet/"  # local
# mainPath = "/Users/apple/Seafile/Research/SandBox/"  # local
# mainPath = "D:/Seafile/Seafile/SandBox/"
# mainPath = "/public/sandbox/workdir/liumeng/SandBox"  # sandbox
# mainPath = "/public/sandbox/workdir/wangchy/SharedFolder"  # sandbox
resultPath = mainPath + "/plot/sfig7/function/"



import os
import cv2
import sys

from scripts.utils.func import concate_png
import pandas as pd
import numpy as np


# positions = [os.path.join(resultPath, file_name) for file_name in os.listdir(resultPath)]

# # concate_name = [i[0] for i in map(os.path.splitext, os.listdir(resultPath))]
# # sx = pd.read_csv(mainPath + "/data/color_map/organ_color_map.csv")
# # sx = {value:key for key, value in sx["organ"].items()}
# # concate_id = [sx.get(concate_name[i]) for i in range(len(concate_name))]
# # concate_order = np.argsort(concate_id)
# # positions = [positions[i] for i in concate_order]

# png = concate_png(positions, cshape=[len(positions), 1])
# cv2.imwrite(os.path.join(resultPath, "concate.png"), png)






# resultPath = mainPath + "/plot/sfig7/morphology/"
# positions = [os.path.join(resultPath, file_name) for file_name in os.listdir(resultPath)]
# png = concate_png(positions, cshape=[len(positions), 1])
# cv2.imwrite(os.path.join(resultPath, "concate.png"), png)





# resultPath = mainPath + "/plot/sfig7/texture/"
# positions = [os.path.join(resultPath, file_name) for file_name in os.listdir(resultPath)]
# png = concate_png(positions, cshape=[len(positions), 1])
# cv2.imwrite(os.path.join(resultPath, "concate.png"), png)




# resultPath = mainPath + "/plot/sfig7/vessel/"
# positions = [os.path.join(resultPath, file_name) for file_name in os.listdir(resultPath)]
# png = concate_png(positions, cshape=[len(positions), 1])
# cv2.imwrite(os.path.join(resultPath, "concate.png"), png)



resultPath = mainPath + "/plot/sfig7/modality/"
positions = [os.path.join(resultPath, file_name) for file_name in os.listdir(resultPath)]
png = concate_png(positions, cshape=[len(positions), 1])
cv2.imwrite(os.path.join(resultPath, "concate.png"), png)

import angles as angles
import matplotlib
import matplotlib.pyplot as plt
import math
import numpy as np
import os
import pandas as pd
matplotlib.rcParams['font.family']='Arial'
plt.rcParams["font.weight"] = "bold"

# user-defined path
mainPath = "E:/Projects/indNet/"
codePath = mainPath + "plot/RadarMap/"
dataPath = mainPath + "data/RadarMap/organ_change_by_age.csv"

if not os.path.exists(codePath):
    os.makedirs(codePath)   
# set path and load data
os.chdir(codePath)
data_= pd.read_csv(dataPath)

theta=data_[['organ']]
for i in theta.columns:
    new_theta = theta[i].values.tolist()
theta = new_theta

rader_labels=radar_labels=np.array(theta)

for i in range(1,len(data_.columns)):
    data=data_.iloc[:,i]
    angl=np.arange(0,2*math.pi,2*math.pi/7)
    fig=plt.figure(facecolor="white")
    ax=plt.subplot(111, polar=True)
    plt.plot(angl,data,'o-',linewidth=1.5, alpha=0.3,color="#FF6600")
    plt.fill(angl,data, alpha=0.25,color="#FF6600")
    ax.set_theta_zero_location('N')
    ax.set_ylim(0,0.8)
    ax.set_yticks([0,0.4,0.8])
    ax.spines['polar'].set_visible(False)
    plt.yticks(fontsize=15)
    plt.thetagrids(angl*180/math.pi,radar_labels,fontsize=17)
    plt.grid(True)
    plt.savefig(codePath+data_.columns[i]+".png",dpi=600)
    plt.clf()

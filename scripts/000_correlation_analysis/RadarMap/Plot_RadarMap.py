import angles as angles
import matplotlib
import matplotlib.pyplot as plt
import math
import numpy as np
import os
import pandas as pd

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

matplotlib.rcParams['font.family']='Arial'
plt.rcParams["font.weight"] = "bold"

# user-defined path
mainPath = '/mnt/d/University/fdurop/remote_ws/indNet'
# mainPath="E:/Projects/indNet"  # xumei
# mainPath = "/Users/apple/Seafile/Research/SandBox"  # local
# mainPath ="/public/sandbox/workdir/wangchy/SharedFolder"  # sandbox
# mainPath = "/public/sandbox/workdir/liumeng/SandBox"  # sandbox
codePath =mainPath+"/plot/RadarMap/sandbox/"
dataPath =mainPath+"/data/RadarMap/sandbox/organ_change_by_age.csv"

ensure_dir(codePath)

data_= pd.read_csv(dataPath)

theta=data_[['organ']]
for i in theta.columns:
    new_theta = theta[i].values.tolist()
theta = new_theta

rader_labels=radar_labels=np.array(theta)

for i in range(1,len(data_.columns)):
    data=data_.iloc[:,i]
    angl=np.arange(2*math.pi,0,-2*math.pi/7)
    fig=plt.figure(facecolor="white")
    ax=plt.subplot(111, polar=True)
    plt.plot(angl,data,'o-',linewidth=1.5, alpha=0.3,color="#1F78B4")
    plt.fill(angl,data, alpha=0.25,color="#1F78B4")
    ax.set_theta_zero_location('N')
    ax.set_ylim(0.5,0.9)
    ax.set_yticks([0.5, 0.6, 0.7, 0.8,0.9])
    ax.spines['polar'].set_visible(False)
    plt.yticks(fontsize=20)
    plt.thetagrids(angl*180/math.pi,radar_labels,fontsize=20)
    plt.grid(True)
    plt.savefig(codePath+data_.columns[i]+".png",dpi=600)
    plt.clf()




######ukb
codePath =mainPath+"/plot/RadarMap/ukb/"
dataPath =mainPath+"/data/RadarMap/ukb/organ_change_by_age.csv"

ensure_dir(codePath)
data_= pd.read_csv(dataPath)

theta=data_[['organ']]
for i in theta.columns:
    new_theta = theta[i].values.tolist()
theta = new_theta

rader_labels=radar_labels=np.array(theta)

for i in range(1,len(data_.columns)):
    data=data_.iloc[:,i]
    angl=np.arange(2*math.pi,0,-2*math.pi/7)
    fig=plt.figure(facecolor="white")
    ax=plt.subplot(111, polar=True)
    plt.plot(angl,data,'o-',linewidth=1.5, alpha=0.3,color="#1F78B4")
    plt.fill(angl,data, alpha=0.25,color="#1F78B4")
    ax.set_theta_zero_location('N')
    ax.set_ylim(0.5,0.8)
    ax.set_yticks([0.5,0.6,0.7,0.8])
    ax.spines['polar'].set_visible(False)
    plt.yticks(fontsize=20)
    plt.thetagrids(angl*180/math.pi,radar_labels,fontsize=20)
    plt.grid(True)
    plt.savefig(codePath+data_.columns[i]+".png",dpi=600)
    plt.clf()

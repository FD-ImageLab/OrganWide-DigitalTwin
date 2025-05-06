import matplotlib
import matplotlib.pyplot as plt
plt.rc('font',family='Arial')
import numpy as np
import pandas as pd
import os





# user-defined path
mainPath='/mnt/d/University/fdurop/remote_ws/indNet'  # xumei
organPath=mainPath+"/SandBox/Organ/"
savePath=mainPath+"/plot/sFig1/sandbox/"
sexPath=mainPath + "/正常队列扫描名单及异常病变记录-中山.csv"
# organ_name=["Brain", "Heart", "Lung", "Liver", "Spleen","Pancreas", "Kidney", "Prostate", "Uterus"]
organ_name=["Brain", "Heart", "Lung", "Liver", "Spleen","Pancreas", "Kidney"]
brain=[]
heart=[]
kidney=[]
liver=[]
lung=[]
pancreas=[]
prostate=[]
spleen=[]
uterus=[]
sexdata=pd.read_csv(sexPath)
sexdata = sexdata.rename(columns={'编号': 'buId'})
file=os.listdir(organPath)
for table in range(len(file)):
    if "brain" in file[table]:
        brain.append(pd.read_csv(organPath+file[table]).iloc[:,0])
    elif "heart" in file[table]:
        heart.append(pd.read_csv(organPath+file[table]).iloc[:,0].values.tolist())
    elif "kidney" in file[table]:
        kidney.append(pd.read_csv(organPath+file[table]).iloc[:,0])
    elif "liver" in file[table] and ".csv" in file[table]:
        liver.append(pd.read_csv(organPath+file[table]).iloc[:,0])
    elif "liver" in file[table] and ".xlsx" in file[table]:
        liver.append(pd.read_xlsx(organPath+file[table]).iloc[:,0])
    elif "lung" in file[table]:
        lung.append(pd.read_csv(organPath+file[table]).iloc[:,0])
    elif "pancreas" in file[table]:
        pancreas.append(pd.read_csv(organPath+file[table]).iloc[:,0])
    elif "prostate" in file[table]:
        prostate.append(pd.read_csv(organPath+file[table]).iloc[:,0]) 
    elif "spleen" in file[table] and ".csv" in file[table]:
        spleen.append(pd.read_csv(organPath+file[table]).iloc[:,0])
    elif "spleen" in file[table] and ".xlsx" in file[table]:
        spleen.append(pd.read_xlsx(organPath+file[table]).iloc[:,0])
    elif "uterus" in file[table]:
        uterus.append(pd.read_csv(organPath+file[table]).iloc[:,0]) 
brain=pd.DataFrame({'buId':list(set(brain[0])),})
heart=pd.DataFrame({'buId':list(set(heart[0])),})
kidney=pd.DataFrame({'buId':list(set(kidney[0])),})
liver=pd.DataFrame({'buId':list(set(liver[0])),})
lung=pd.DataFrame({'buId':list(set(lung[0])),})
pancreas=pd.DataFrame({'buId':list(set(pancreas[0])),})
prostate=pd.DataFrame({'buId':list(set(prostate[0])),})
spleen=pd.DataFrame({'buId':list(set(spleen[0])),})
uterus=pd.DataFrame({'buId':list(set(uterus[0])),})

a=list(set(brain['buId']).intersection(heart['buId'],kidney['buId'],liver['buId'],prostate['buId'],spleen['buId']))#         print(len(list(set(brain['buId']).intersection(heart['buId'],kidney['buId'],liver['buId'],spleen['buId'],uterus['buId']))))
b=list(set(brain['buId']).intersection(heart['buId'],kidney['buId'],liver['buId'],spleen['buId'],uterus['buId']))
number=list(set(a).union(set(b)))

# number=[len(brain),len(heart),len(lung),len(liver),len(spleen),len(pancreas),len(kidney),len(prostate),len(uterus)]
number=[len(brain),len(heart),len(lung),len(liver),len(spleen),len(pancreas),len(kidney)]


# 输出器官交集人数
set.intersection(set(brain["buId"]), set(heart["buId"]), set(lung["buId"]), set(liver["buId"]), set(kidney["buId"]), set(pancreas["buId"]), set(spleen["buId"]))

# figuresize=[9,7]

fig, ax = plt.subplots()
# plt.rcParams["figure.figsize"] = figuresize

ax.bar(organ_name,height=number,color=["#BDB0A5", "#EB8677", "#BDDD78", "#F2B670", "#7DBFA6", "#BDBBD7", "#EE924F"])
plt.xticks(rotation=30)
plt.tick_params(labelsize=16)
plt.subplots_adjust(bottom=0.2, left=0.15)
i=0
for x,y in zip(organ_name,number):
    plt.text(x,y/2,str(number[i]),ha='center',fontsize=16)
    i=i+1
ax.set_yticks([0,250,500,750,1000])
ax.set_ylabel("Number",size=20)
if not os.path.exists(savePath):
    os.makedirs(savePath)
plt.savefig(savePath+'/Organ_number.png',dpi=600) 
plt.clf()
        
        
        
        
        
# pheno_number = [453, 769, 136, 125, 88, 222, 162, 153, 210]

# fig, ax = plt.subplots()
# # plt.rcParams["figure.figsize"] = figuresize

# ax.bar(organ_name,height=pheno_number,color=["#BDB0A5", "#EB8677", "#BDDD78", "#F2B670", "#7DBFA6", "#BDBBD7", "#EE924F", "#7AADD2", "#DA8FC0"])
# plt.xticks(rotation=30)
# plt.tick_params(labelsize=16)
# plt.subplots_adjust(bottom=0.2, left=0.15)
# i=0
# for x,y in zip(organ_name,pheno_number):
#     plt.text(x,y/2,str(pheno_number[i]),ha='center',fontsize=16)
#     i=i+1
# ax.set_yticks([0,150,300,450,600,750])
# ax.set_ylabel("Pheno number",size=20)
# if not os.path.exists(savePath):
#     os.makedirs(savePath)
# plt.savefig(savePath+'/Organ_pheno_number.png',dpi=600) 
# plt.clf()




############## ukb ###################


savePath=mainPath+"/plot/sFig1/ukb/"
organ_name = ["Brain", "Heart", "Lung", "Liver", "Spleen","Pancreas", "Kidney"]
number = [46393, 32462, 23154, 27418, 28475, 29464, 36137]

fig, ax = plt.subplots()
# plt.rcParams["figure.figsize"] = figuresize

ax.bar(organ_name,height=number,color=["#BDB0A5", "#EB8677", "#BDDD78", "#F2B670", "#7DBFA6", "#BDBBD7", "#EE924F"])
plt.xticks(rotation=30)
plt.tick_params(labelsize=16)
plt.subplots_adjust(bottom=0.2, left=0.15)
i=0
for x,y in zip(organ_name,number):
    plt.text(x,y/2,str(number[i]),ha='center',fontsize=14)
    i=i+1
ax.set_yticks([0,10000,20000,30000,40000,50000])
ax.set_ylabel("Pheno number",size=20)
if not os.path.exists(savePath):
    os.makedirs(savePath)
plt.savefig(savePath+'/Organ_number.png',dpi=600, bbox_inches='tight')
plt.clf()




pheno_number = [139, 82, 11, 6, 3, 3, 6]

fig, ax = plt.subplots()
# plt.rcParams["figure.figsize"] = figuresize

ax.bar(organ_name,height=pheno_number,color=["#BDB0A5", "#EB8677", "#BDDD78", "#F2B670", "#7DBFA6", "#BDBBD7", "#EE924F"])
plt.xticks(rotation=30)
plt.tick_params(labelsize=24)
plt.subplots_adjust(bottom=0.2, left=0.15)
i=0
for x,y in zip(organ_name,pheno_number):
    plt.text(x,y/2,str(pheno_number[i]),ha='center',fontsize=24)
    i=i+1
ax.set_yticks([0,30,60,90,120,150])
ax.set_ylabel("Pheno number",size=20)
if not os.path.exists(savePath):
    os.makedirs(savePath)
plt.savefig(savePath+'/Organ_pheno_number.png',dpi=600, bbox_inches='tight')
plt.clf()
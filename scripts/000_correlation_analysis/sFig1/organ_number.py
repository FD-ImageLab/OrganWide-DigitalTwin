import matplotlib
import matplotlib.pyplot as plt
plt.rc('font',family='Arial')
import numpy as np
import pandas as pd
from numpy import *
import os

if __name__ == "__main__":
        # user-defined path
        # mainPath="E:/Projects/indNet/"  # xumei
        mainPath="/mnt/d/University/fdurop/remote_ws/indNet"      
        organPath=mainPath+"/Sandbox/Organ/"
        savePath=mainPath+"/plot/sfig1/"
        sexPath1=mainPath+"/SandBox/Demography/张江扫描名单.csv"
        sexPath2=mainPath+"/SandBox/Demography/中山扫描名单.csv"

        organ_name=["Brain", "Heart", "Lung", "Liver", "Spleen","Pancreas", "Kidney", "Prostate", "Uterus"]
        brain=[]
        heart=[]
        kidney=[]
        liver=[]
        lung=[]
        pancreas=[]
        prostate=[]
        spleen=[]
        uterus=[]
        sexdata1=pd.read_csv(sexPath1).iloc[:,[3,6,5]]
        sexdata1.columns=['ID','Age','Sex']
        sexdata2=pd.read_csv(sexPath2).iloc[:,[3,6,5]]
        sexdata2.columns=['ID','Age','Sex']
        sexdata=sexdata1._append(sexdata2)

        file=os.listdir(organPath)
        for table in range(len(file)):
            if "brain" in file[table] and ".csv" in file[table]:
                brain.append(pd.read_csv(organPath+file[table]).iloc[:,0])
            elif "heart" in file[table] and ".csv" in file[table]:
                heart.append(pd.read_csv(organPath+file[table]).iloc[:,0].values.tolist())
            elif "kidney" in file[table] and ".csv" in file[table]:
                kidney.append(pd.read_csv(organPath+file[table]).iloc[:,0])
            elif "liver" in file[table] and ".csv" in file[table]:
                liver.append(pd.read_csv(organPath+file[table]).iloc[:,0])
            elif "lung" in file[table] and ".csv" in file[table]:
                lung.append(pd.read_csv(organPath+file[table]).iloc[:,0])
            elif "pancreas" in file[table] and ".csv" in file[table]:
                pancreas.append(pd.read_csv(organPath+file[table]).iloc[:,0])
            elif "prostate" in file[table] and ".csv" in file[table]:
                prostate.append(pd.read_csv(organPath+file[table]).iloc[:,0]) 
            elif "spleen" in file[table] and ".csv" in file[table]:
                spleen.append(pd.read_csv(organPath+file[table]).iloc[:,0])
            elif "uterus" in file[table] and ".csv" in file[table]:
                uterus.append(pd.read_csv(organPath+file[table]).iloc[:,0]) 
        brain=pd.DataFrame({'ID':list(set(brain[0])),})
        heart=pd.DataFrame({'ID':list(set(heart[0])),})
        kidney=pd.DataFrame({'ID':list(set(kidney[0])),})
        liver=pd.DataFrame({'ID':list(set(liver[0])),})
        lung=pd.DataFrame({'ID':list(set(lung[0])),})
        pancreas=pd.DataFrame({'ID':list(set(pancreas[0])),})
        prostate=pd.DataFrame({'ID':list(set(prostate[0])),})
        spleen=pd.DataFrame({'ID':list(set(spleen[0])),})
        uterus=pd.DataFrame({'ID':list(set(uterus[0])),})
        a=list(set(brain['ID']).intersection(heart['ID'],kidney['ID'],liver['ID'],prostate['ID'],spleen['ID'],lung['ID']))#         print(len(list(set(brain['buId']).intersection(heart['buId'],kidney['buId'],liver['buId'],spleen['buId'],uterus['buId']))))
        b=list(set(brain['ID']).intersection(heart['ID'],kidney['ID'],liver['ID'],spleen['ID'],uterus['ID'],lung['ID']))
        number1=list(set(a).union(set(b)))
        union=prostate._append(uterus)
#         print(len(number1))
#         print(len(list(set(brain['ID'])))-number1)
        
#         print(len(list(set(heart['ID'])))-number1)
        
#         print(len(list(set(kidney['ID'])))-number1)
        
#         print(len(list(set(liver['ID'])))-number1)
        
#         print(len(list(set(lung['ID'])))-number1)
        
#         print(len(list(set(pancreas['ID'])))-number1)
        
#         print(len(list(set(prostate['ID'])))-number1)
        
#         print(len(list(set(spleen['ID'])))-number1)
        
#         print(len(list(set(uterus['ID'])))-number1)
        
    
        
        number=[len(brain),len(heart),len(lung),len(liver),len(spleen),len(pancreas),len(kidney),len(prostate),len(uterus)]
        
        brain_sex=pd.merge(brain,sexdata,on=["ID"])
        heart_sex=pd.merge(heart,sexdata,on=["ID"])
        kidney_sex=pd.merge(kidney,sexdata,on=["ID"])
        liver_sex=pd.merge(liver,sexdata,on=["ID"])
        lung_sex=pd.merge(lung,sexdata,on=["ID"])
        pancreas_sex=pd.merge(pancreas,sexdata,on=["ID"])
        prostate_sex=pd.merge(prostate,sexdata,on=["ID"])
        spleen_sex=pd.merge(spleen,sexdata,on=["ID"])
        uterus_sex=pd.merge(uterus,sexdata,on=["ID"])

        brain_male=brain_sex[brain_sex["Age"].isin(["男"])]
        brain_female=brain_sex[brain_sex["Age"].isin(["女"])]
        heart_male=heart_sex[heart_sex["Age"].isin(["男"])]
        heart_female=heart_sex[heart_sex["Age"].isin(["女"])]
        kidney_male=kidney_sex[kidney_sex["Age"].isin(["男"])]
        kidney_female=kidney_sex[kidney_sex["Age"].isin(["女"])]
        liver_male=liver_sex[liver_sex["Age"].isin(["男"])]
        liver_female=liver_sex[liver_sex["Age"].isin(["女"])]
        lung_male=lung_sex[lung_sex["Age"].isin(["男"])]
        lung_female=lung_sex[lung_sex["Age"].isin(["女"])]
        pancreas_male=pancreas_sex[pancreas_sex["Age"].isin(["男"])]
        pancreas_female=pancreas_sex[pancreas_sex["Age"].isin(["女"])]
        prostate_male=prostate_sex[prostate_sex["Age"].isin(["男"])]
        prostate_female=prostate_sex[prostate_sex["Age"].isin(["女"])]
        spleen_male=spleen_sex[spleen_sex["Age"].isin(["男"])]
        spleen_female=spleen_sex[spleen_sex["Age"].isin(["女"])]
        uterus_male=uterus_sex[uterus_sex["Age"].isin(["男"])]
        uterus_female=uterus_sex[uterus_sex["Age"].isin(["女"])]

#         uterus_female=uterus_female.append(prostate_female)
#         print(len(list(set(brain_female['ID']).intersection(heart_female['ID'],kidney_female['ID'],lung_female['ID'],pancreas_female['ID'],liver_female['ID'],spleen_female['ID'],uterus_female['ID']))))
# #         print(len(list(set(brain_female['buId']).intersection(heart_female['buId'],kidney_female['buId'],liver_female['buId'],spleen_female['buId'],uterus_female['buId']))))
# #         print(len(list(set(brain_female['buId']).intersection(heart_female['buId'],kidney_female['buId'],liver_female['buId'],spleen_female['buId']))))
        
#         print(len(list(set(brain_female['ID'])-set(heart_female['ID'])-set(kidney_female['ID'])-set(liver_female['ID'])-set(lung_female['ID'])-set(pancreas_female['ID'])-set(spleen_female['ID'])-set(uterus_female['ID']))))
        
#         print(len(list(set(heart_female['ID'])-set(brain_female['ID'])-set(kidney_female['ID'])-set(liver_female['ID'])-set(lung_female['ID'])-set(pancreas_female['ID'])-set(spleen_female['ID'])-set(uterus_female['ID']))))
        
#         print(len(list(set(kidney_female['ID'])-set(heart_female['ID'])-set(brain_female['ID'])-set(liver_female['ID'])-set(lung_female['ID'])-set(pancreas_female['ID'])-set(spleen_female['ID'])-set(uterus_female['ID']))))
        
#         print(len(list(set(liver_female['ID'])-set(heart_female['ID'])-set(kidney_female['ID'])-set(brain_female['ID'])-set(lung_female['ID'])-set(pancreas_female['ID'])-set(spleen_female['ID'])-set(uterus_female['ID']))))
        
#         print(len(list(set(lung_female['ID'])-set(heart_female['ID'])-set(kidney_female['ID'])-set(liver_female['ID'])-set(brain_female['ID'])-set(pancreas_female['ID'])-set(spleen_female['ID'])-set(uterus_female['ID']))))
        
#         print(len(list(set(pancreas_female['ID'])-set(heart_female['ID'])-set(kidney_female['ID'])-set(liver_female['ID'])-set(lung_female['ID'])-set(brain_female['ID'])-set(spleen_female['ID'])-set(uterus_female['ID']))))
        
# # #         print(len(list(set(prostate_female['buId'])-set(heart_female['buId'])-set(kidney_female['buId'])-set(liver_female['buId'])-set(lung_female['buId'])-set(pancreas_female['buId'])-set(brain_female['buId'])-set(spleen_female['buId'])-set(uterus_female['buId']))))
        
#         print(len(list(set(spleen_female['ID'])-set(heart_female['ID'])-set(kidney_female['ID'])-set(liver_female['ID'])-set(lung_female['ID'])-set(pancreas_female['ID'])-set(brain_female['ID'])-set(uterus_female['ID']))))
        
#         print(len(list(set(uterus_female['ID'])-set(heart_female['ID'])-set(kidney_female['ID'])-set(liver_female['ID'])-set(lung_female['ID'])-set(pancreas_female['ID'])-set(spleen_female['ID'])-set(brain_female['ID']))))
   
        print(len(list(set(brain_male['ID']).intersection(heart_male['ID'],kidney_male['ID'],liver_male['ID'],prostate_male['ID'],spleen_male['ID']))))
#         print(len(list(set(brain_male['buId']).intersection(heart_male['buId'],kidney_male['buId'],liver_male['buId'],spleen_male['buId'],uterus_male['buId']))))
#         print(len(list(set(brain_male['buId']).intersection(heart_male['buId'],kidney_male['buId'],liver_male['buId'],spleen_male['buId'],prostate_male['buId']))))
        
        print(len(list(set(brain_male['ID'])-set(heart_male['ID'])-set(kidney_male['ID'])-set(liver_male['ID'])-set(lung_male['ID'])-set(pancreas_male['ID'])-set(spleen_male['ID']))))
        
        print(len(list(set(heart_male['ID'])-set(brain_male['ID'])-set(kidney_male['ID'])-set(liver_male['ID'])-set(lung_male['ID'])-set(pancreas_male['ID'])-set(spleen_male['ID']))))
        
        print(len(list(set(kidney_male['ID'])-set(heart_male['ID'])-set(brain_male['ID'])-set(liver_male['ID'])-set(lung_male['ID'])-set(pancreas_male['ID'])-set(spleen_male['ID']))))
        
        print(len(list(set(liver_male['ID'])-set(heart_male['ID'])-set(kidney_male['ID'])-set(brain_male['ID'])-set(lung_male['ID'])-set(pancreas_male['ID'])-set(spleen_male['ID']))))
        
        print(len(list(set(lung_male['ID'])-set(heart_male['ID'])-set(kidney_male['ID'])-set(liver_male['ID'])-set(brain_male['ID'])-set(pancreas_male['ID'])-set(spleen_male['ID']))))
        
        print(len(list(set(pancreas_male['ID'])-set(heart_male['ID'])-set(kidney_male['ID'])-set(liver_male['ID'])-set(lung_male['ID'])-set(brain_male['ID'])-set(spleen_male['ID']))))
                
        print(len(list(set(spleen_male['ID'])-set(heart_male['ID'])-set(kidney_male['ID'])-set(liver_male['ID'])-set(lung_male['ID'])-set(pancreas_male['ID'])-set(brain_male['ID']))))
        
    
        brain_male_number=[]
        brain_female_number=[]
        heart_male_number=[]
        heart_female_number=[]
        kidney_male_number=[]
        kidney_female_number=[]
        liver_male_number=[]
        liver_female_number=[]
        lung_male_number=[]
        lung_female_number=[]
        pancreas_male_number=[]
        pancreas_female_number=[]
        prostate_male_number=[]
        prostate_female_number=[]
        spleen_male_number=[]
        spleen_female_number=[]
        uterus_male_number=[]
        uterus_female_number=[]
        for i in range(20,61,1):
            brain_male_number.append(list(brain_male["Age"]).count(i))
            brain_female_number.append(list(brain_female["Age"]).count(i))
            heart_male_number.append(list(heart_male["Age"]).count(i))
            heart_female_number.append(list(heart_female["Age"]).count(i))
            kidney_male_number.append(list(kidney_male["Age"]).count(i))
            kidney_female_number.append(list(kidney_female["Age"]).count(i))
            liver_male_number.append(list(liver_male["Age"]).count(i))
            liver_female_number.append(list(liver_female["Age"]).count(i))
            lung_male_number.append(list(lung_male["Age"]).count(i))
            lung_female_number.append(list(lung_female["Age"]).count(i))
            pancreas_male_number.append(list(pancreas_male["Age"]).count(i))
            pancreas_female_number.append(list(pancreas_female["Age"]).count(i))
            prostate_male_number.append(list(prostate_male["Age"]).count(i))
            prostate_female_number.append(list(prostate_female["Age"]).count(i))
            spleen_male_number.append(list(spleen_male["Age"]).count(i))
            spleen_female_number.append(list(spleen_female["Age"]).count(i))
            uterus_male_number.append(list(uterus_male["Age"]).count(i))
            uterus_female_number.append(list(uterus_female["Age"]).count(i))
        prostate_female_number=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        male=[brain_male_number,heart_male_number,lung_male_number,liver_male_number,spleen_male_number,pancreas_male_number,kidney_male_number,prostate_male_number,uterus_male_number]
        female=[brain_female_number,heart_female_number,lung_female_number,liver_female_number,spleen_female_number,pancreas_female_number,kidney_female_number,prostate_female_number,uterus_female_number]
        for i in range(len(male)):
            a=sum(male[i])
            print(a)
        for i in range(len(female)):
            b=sum(female[i])
            print(b)
        
        figuresize=[28,15]

        fig, ax = plt.subplots()
        plt.rcParams["figure.figsize"] = figuresize

        ax.bar(organ_name,height=number,color=["#BDB0A5", "#EB8677", "#BDDD78", "#F2B670", "#7DBFA6", "#BDBBD7", "#EE924F", "#7AADD2", "#DA8FC0"])
        plt.xticks(rotation=30)
        plt.tick_params(labelsize=14)
        i=0
        for x,y in zip(organ_name,number):
            plt.text(x,y/2,str(number[i]),ha='center',fontsize=10)
            i=i+1
        ax.set_yticks([0,175,350,525,700])
        ax.set_ylabel("Number",size=16)
        plt.savefig(savePath+'Organ_number.png',dpi=600) 
        plt.clf()
        
        for i in range(len(organ_name)):
            fig, ax = plt.subplots()
            plt.rcParams["figure.figsize"] = [14,8]
            age=list(range(20,61,1))
            ax.bar(age,height=male[i],color="#6BAED6",label="Male ($\mathit{N}$="+str(sum(male[i]))+")")
            ax.bar(age,bottom=male[i],height=female[i],color="#E78AC3",label="Female ($\mathit{N}$="+str(sum(female[i]))+")")
            plt.tick_params(labelsize=30)
            c=[male[i][a]+female[i][a] for a in range(min(len(male[i]),len(female[i])))]

            ax.set_yticks([0,int(1/4*max(c)),int(2/4*max(c)),int(3/4*max(c)),max(c)])
            ax.set_ylabel("Number",size=30)
            plt.title(organ_name[i],size=30,pad=12)
            plt.rcParams.update({'font.size':25})
            plt.legend(loc='upper right')
            plt.savefig(savePath+organ_name[i]+'.png',dpi=200)
            plt.clf()
  

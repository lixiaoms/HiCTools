import numpy as np
import pandas as pd
depth=[10,100,1000,5000]

##Juicer
juicer=[]
for line in open('Juicer/rao_chr10/raw.txt').readlines():
    juicer.append(float(line.strip()))
juicer_l=np.zeros([len(juicer),80])
for i in range(4):
    for j in range(20):
        sign=[]
        for line in open('Juicer/rao_chr10/'+str(depth[i])+'/eigen_'+str(j+1)+'.txt').readlines():
            sign.append(float(line.strip()))
        juicer_l[:,20*i+j]=sign

##CscoreTool
cscore=[]
for line in open('CscoreTool/rao_chr10/raw/Test_cscore.txt').readlines():
    cscore.append(float(line.strip().split('\t')[1]))
cscore_l=np.zeros([len(juicer),80])
for i in range(4):
    for j in range(20):
        sign=[]
        for line in open('CscoreTool/rao_chr10/'+str(depth[i])+'/Test_'+str(j+1)+'_cscore.txt').readlines():
            sign.append(float(line.strip().split('\t')[1]))
        cscore_l[:,20*i+j]=sign

##GeSICA
gesica=[]
for line in open('GeSICA/rao_chr10/raw/ratio.txt/ratio.txt_interaciton_ratio_res_100000_dis_400000.wig').readlines()[1:]:
    gesica.append(float(line.strip().split('\t')[1]))
gesica_l=np.zeros([len(juicer),80])
for i in range(4):
    for j in range(20):
        sign=[]
        for line in open('GeSICA/rao_chr10/'+str(depth[i])+'/ratio_'+str(j+1)+'.txt/ratio_'+str(j+1)+'.txt_interaciton_ratio_res_100000_dis_400000.wig').readlines()[1:]:
            sign.append(float(line.strip().split('\t')[1]))
        gesica_l[:,20*i+j]=sign
        
##concat the data
df=pd.DataFrame(juicer_l,columns=['J_seed_'+str(i) for i in range(80)])
df=pd.concat([df, pd.DataFrame(columns=['C_seed_'+str(i) for i in range(80)])],sort=False)
df[['C_seed_'+str(i) for i in range(80)]]=cscore_l
df=pd.concat([df, pd.DataFrame(columns=['G_seed_'+str(i) for i in range(80)])],sort=False)
df[['G_seed_'+str(i) for i in range(80)]]=gesica_l
df=df.fillna(0)
df1=pd.DataFrame(np.array([juicer,cscore,gesica]).T,columns=['Juicer','Cscore','GeSICA'])
df1=df1.fillna(0)

##caculate similarity between compartments identified in down sampled data and raw data
similarity=np.zeros([3,80])
for i in range(3):
    total=sum(df1[df1.columns[i]]!=0)
    for j in range(80):
        similarity[i,j]=float(max(sum(df1[df1.columns[i]]*df[df.columns[80*i+j]]>0),sum(df1[df1.columns[i]]*df[df.columns[80*i+j]]<0)))/total

##figure of example
import matplotlib.pyplot as plt
fig,axes=plt.subplots(figsize=(8,8))
for sf in range(3):
    a=similarity[sf]
    a.shape=(4,20)
    axes.boxplot(a.T,showmeans=True,showfliers=False)
    axes.plot(range(0,5),np.concatenate(([1],np.mean(a,axis=1))),label=df1.columns[sf])
plt.legend(loc=3)
plt.xlim((0,4.2))
plt.ylim(0.5,1)
plt.setp(axes, xticks=range(5),xticklabels=['raw depth','1_10 depth','1_100 depth','1_1000 depth','1_5000 depth'])
plt.title('Accuracy Index',fontsize=20)

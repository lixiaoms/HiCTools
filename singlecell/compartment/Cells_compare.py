import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import random

index=range(2,8)+range(9,19) ##cell#18 is pooled cells

juicer=[]
for line in open('../../review/Compartment/Juicer/rao_chr10/raw.txt').readlines():
    juicer.append(float(line.strip()))

juicer1=[]
for i in index:
    sign=[]
    for line in open('chr10_j/cell_'+str(i)+'_eigen.txt').readlines():
        sign.append(float(line.strip()))
    juicer1.append(sign)

df2=pd.DataFrame(np.array([juicer]+juicer1).T,columns=['raw']+['cell'+str(i) for i in index[:-1]]+['pool'])
df2=df2.fillna(0)

##caculate similarity between compartments identified in cells and ensemble data in Rao et al.
similarity=pd.DataFrame(np.zeros([16,2]),columns=['similarity','reads_num'],index=['cell'+str(i) for i in index[:-1]]+['pool'])
total=sum(df2[df2.columns[0]]!=0)
for i in range(16):
    similarity['similarity'][df2.columns[i+1]]=float(max(sum(df2[df2.columns[0]]*df2[df2.columns[i+1]]>0),sum(df2[df2.columns[0]]*df2[df2.columns[i+1]]<0)))/total
    similarity['reads_num'][df2.columns[i+1]]=pd.read_csv('../TAD/chr10/reads_num',header=None)[0][i]

##figure of example
plt.scatter(similarity['reads_num'],similarity['similarity'])
plt.xlabel('contacts num')
plt.ylabel('comparment accuracy')

import numpy as np
import pandas as pd
import re
import sklearn.metrics
import matplotlib.pyplot as plt
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.preprocessing import scale
#Mutual Information
def MI(windows,boundary1,boundary2):
    L1=[]
    L2=[]
    boundary1=sorted(list(set(boundary1[boundary1<windows]))+[0,windows])
    boundary2=sorted(list(set(boundary2[boundary2<windows]))+[0,windows])
    for x in range(len(boundary1)-1):
        for y in range(-boundary1[x]+boundary1[x+1]):
            L1.append(x)
    for x in range(len(boundary2)-1):
        for y in range(-boundary2[x]+boundary2[x+1]):
            L2.append(x)
    return adjusted_mutual_info_score(L1,L2)
#Weights Similarity
def WSlist(windows,L,l):
    L=sorted(list(set(L[L<windows]))+[0,windows])
    l=sorted(list(set(l[l<windows]))+[0,windows])
    return(min(WS(L[:-1],[i-1 for i in L[1:]],l[:-1],[i-1 for i in l[1:]]),WS(l[:-1],[i-1 for i in l[1:]],L[:-1],[i-1 for i in L[1:]])))

TAD={}
raw={}
index=range(2,8)+range(9,19)

raw['IS']=[]
file=open('chr10/raw.is520001.ids320001.insulation.boundaries').readlines()
for line in file[1:]:
    raw['IS'].append(int(line.split('\t')[3]))

TAD['IS']=[]
for i in index:
    TAD['IS'].append([])
    file=open('chr10/cell_'+str(i)+'.is520001.ids320001.insulation.boundaries').readlines()
    for line in file[1:]:
        TAD['IS'][-1].append(int(line.split('\t')[3]))

WSs=pd.DataFrame(np.zeros([16,3]),index=['cell_'+str(i) for i in index[:-1]]+['cell_pool'],columns=['DD','IS','reads_num'])
MIs=pd.DataFrame(np.zeros([16,3]),index=['cell_'+str(i) for i in index[:-1]]+['cell_pool'],columns=['DD','IS','reads_num'])

size=3389
for i in range(16):
    WSs['IS'][i]=WSlist(size,np.unique(TAD['IS'][i]),np.unique(raw['IS']))
    MIs['IS'][i]=MI(size,np.unique(TAD['IS'][i]),np.unique(raw['IS']))
for i in range(16):
    WSs['reads_num'][i]=pd.read_csv('chr10/reads_num',header=None)[0][i]
    MIs['reads_num'][i]=pd.read_csv('chr10/reads_num',header=None)[0][i]

##figure of example
plt.scatter(WSs['reads_num'],WSs['IS'])
plt.xlabel('contacts num')
plt.ylabel('weight similarity')

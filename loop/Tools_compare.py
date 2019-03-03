import numpy as np
import pandas as pd
import re
import sklearn.metrics
import matplotlib.pyplot as plt
import gzip
import seaborn as sns
import os

##nickname
name={}
name['HC']='HiCCUPS'
name['FT']='Fit-Hi-C'
name['FA']='FastHic'
name['HM']='HOMER'
name['DF']='diffHic'
name['CL']='cLoops'

software=['FT', 'CL', 'DF', 'FA', 'HM', 'HC']
loop=pd.DataFrame(np.zeros([51,6]),columns=software)

##HiCCUPS
for i in range(51):
    if i==0: 
        k='raw'
    else:
        k=str(i)
    try:
        loop['HC'][i]=len(open('hiccups/rao_chr10/'+k+'/merged_loops.bedpe').readlines())-1
    except:
        continue

##Fit-Hi-C
for i in range(51):
    if i==0: 
        k='raw'
    else:
        k=str(i)
    fithic=pd.read_csv('fithic/rao_chr10/'+k+'/FitHiC.spline_pass1.res10000.significances.txt',sep='\t')
    loop['FT'][i]=sum(fithic['q-value']<0.05)

##FastHiC
for i in range(25):
    for j in range(1,185):
        try:
            loop['FA'][2*i+1]+=sum(pd.read_csv('fasthic/rao_chr10/'+str(i+26)+'/TAD_'+str(j)+'_prob',sep='\t')['PeakProbability']>0.95)
        except:
            continue

##HOMER
for i in range(51):
    if i==0: 
        k='raw'
    else:
        k=str(i)
    try:
        loop['HM'][i]=len(open('homer/rao_chr10/'+k+'/0/0.loop.2D.bed').readlines())
    except:
        continue

##cLoops
for i in range(51):
    if i==0: 
        k='fast'
    else:
        k='_'+str(i)
    try:
        loop['CL'][i]=len(open('cloops/rao_chr10/hic'+k+'_loops_juicebox.txt').readlines())-1
    except:
        continue

##figure of example
fig = plt.figure(figsize=[8,8]) 
colors = ['b','g','r','orange','c','m']
for i in name.keys():
    plt.scatter(range(51),np.log10(loop[i]),label=name[i])
plt.legend(prop={'size':12})
plt.xlabel('sample rate',size=15)
plt.ylabel('Number of loops',size=15)
plt.ylim(-0.1,5.5)
plt.xticks( [i*np.log(4)/(np.log(10)-np.log(9)) for i in range(5)]+[50] , ('raw','1/4', '1/16', '1/64','1/256','') )
plt.plot( [50,50],[-0.1,np.log10(loop['FA'][49])] , color = 'blue', linewidth=2, linestyle=":" )
plt.text(-6.5,5.5,'xlog10')

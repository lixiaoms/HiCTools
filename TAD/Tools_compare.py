import numpy as np
import pandas as pd
import re
import sklearn.metrics
import matplotlib.pyplot as plt
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

#Weights Similarity for inconsecutive TADs
def WS(op1,ed1,op2,ed2):
    length=[]
    weight=[]
    for x in range(len(op1)):
        length.append(ed1[x]+1-op1[x])
        inter=[]
        for y in range(len(op2)):
            inter.append(len(list(set(range(op1[x],ed1[x]+1))&set(range(op2[y],ed2[y]+1))))*((ed1[x]+1-op1[x])*(ed2[y]+1-op2[y]))**(-0.5))
        weight.append(max(inter)*(ed1[x]+1-op1[x]))
    return(sum(weight)*(int(sum(length))**(-1)))

#Weights Similarity for consecutive TADs
def WSlist(windows,L,l):
    L=sorted(list(set(L[L<windows]))+[0,windows])
    l=sorted(list(set(l[l<windows]))+[0,windows])
    return(min(WS(L[:-1],[i-1 for i in L[1:]],l[:-1],[i-1 for i in l[1:]]),WS(l[:-1],[i-1 for i in l[1:]],L[:-1],[i-1 for i in L[1:]])))

##nickname
name={}
name['IS']='Insulation Score'
name['HS']='HiCSeg'
name['MT']='MrTADFinder'
name['AT']='Armatus'
name['DD']='deDoc(E)'
name['HM']='HOMER'
software=['IS','DD','MT','HM',"AT",'HS']

##example on down sampled chr18 data in Dixon et al.
n=18
size=1903

##start point of TADs
TAD={}
for sf in software:
    TAD[sf]=[]
raw={}
for sf in software:
    raw[sf]=[]

##end points of inconsecutive TADs
TAD1={}
software1=['MT','HM',"AT",'DD']
for sf in software1:
    TAD1[sf]=[]
raw1={}
for sf in software1:
    raw1[sf]=[]

##dataframe to save result
for i in [10,50,200]:
    for sf in software:
        TAD[sf].append([])
        for j in range(100,120):
            TAD[sf][-1].append([])
    for sf in software1:
        TAD1[sf].append([])
        for j in range(100,120):
            TAD1[sf][-1].append([])

##IS
file=open('IS/chr'+str(n)+'/raw.is520001.ids320001.insulation.boundaries').readlines()
for line in file[1:]:
    raw['IS'].append(int(line.split('\t')[3]))
for i in range(3):
    for j in range(20):
        file=open('IS/chr'+str(n)+'/1_'+str([10,50,200][i])+'/'+str(100+j)+'/matrix.is520001.ids320001.insulation.boundaries').readlines()
        for line in file[1:]:
            TAD['IS'][i][j].append(int(line.split('\t')[3]))

##HS
file=open('HiCseg/chr'+str(n)+'/raw').readlines()
for line in file:
    raw['HS'].append(int(line))
for i in range(3):
    for j in range(20):
        file=open('HiCseg/chr'+str(n)+'/chr'+str(n)+'_'+str([10,50,200][i])+'_'+str(100+j)).readlines()
        for line in file:
            TAD['HS'][i][j].append(int(line))

##DD
file=open('deDoc/raw/chr'+str(n)+'.graph.deDoc(E)').readlines()
for line in file:
    if len(line.strip().split(' '))>4:
        raw['DD'].append(int(line.strip().split(' ')[0]))
        raw1['DD'].append(int(line.strip().split(' ')[-1]))
for i in range(3):
    for j in range(20):
        file=open('deDoc/1_'+str([10,50,200][i])+'/chr'+str(n)+'_'+str(100+j)+'.graph.deDoc(E)').readlines()
        for line in file:
            if i==0:
                if len(line.strip().split(' '))>4:
                    TAD['DD'][i][j].append(int(line.strip().split(' ')[0]))
                    TAD1['DD'][i][j].append(int(line.strip().split(' ')[-1]))
            if i==1:
                if len(line.strip().split(' '))>4:
                    TAD['DD'][i][j].append(int(line.strip().split(' ')[0]))
                    TAD1['DD'][i][j].append(int(line.strip().split(' ')[-1]))
            if i==2:
                if len(line.strip().split(' '))>4:
                    TAD['DD'][i][j].append(int(line.strip().split(' ')[0]))
                    TAD1['DD'][i][j].append(int(line.strip().split(' ')[-1]))

##MT
file=open('MrTADFinder/raw/TAD_chr'+str(n)+'.bed').readlines()
for line in file:
    if int(line.strip().split('\t')[2])/40000-int(line.strip().split('\t')[1])/40000>0:
        raw['MT'].append(int(line.strip().split('\t')[1])/40000+1)
        raw1['MT'].append(int(line.strip().split('\t')[2])/40000)
for i in range(3):
    for j in range(20):
        file=open('MrTADFinder/1_'+str([10,50,200][i])+'/TAD_chr'+str(n)+'_'+str([10,50,200][i])+'_'+str(100+j)+'.bed').readlines()
        for line in file:
            if int(line.strip().split('\t')[2])/40000-int(line.strip().split('\t')[1])/40000>0:
                TAD['MT'][i][j].append(int(line.strip().split('\t')[1])/40000+1)
                TAD1['MT'][i][j].append(int(line.strip().split('\t')[2])/40000)

##HM
file=open('homer/chr'+str(n)+'/raw/raw.tad.2D.bed').readlines()
for line in file:
    raw['HM'].append(int(line.strip().split('\t')[1])/40000+1)
    raw1['HM'].append(int(line.strip().split('\t')[2])/40000+1)
for i in range(3):
    for j in range(20):
        file=open('homer/chr'+str(n)+'/1_'+str([10,50,200][i])+'_'+str(100+j)+'/1_'+str([10,50,200][i])+'_'+str(100+j)+'.tad.2D.bed').readlines()
        for line in file:
            TAD['HM'][i][j].append(int(line.strip().split('\t')[1])/40000+1)
            TAD1['HM'][i][j].append(int(line.strip().split('\t')[2])/40000+1)

##AT
file=open('armatus/chr'+str(n)+'/chr'+str(n)+'_raw.consensus.txt').readlines()
for line in file:
    if int(line.strip().split('\t')[2])/40000-int(line.strip().split('\t')[1])/40000>0:
        raw['AT'].append(int(line.strip().split('\t')[1])/40000)
        raw1['AT'].append(int(line.strip().split('\t')[2])/40000)
for i in range(3):
    for j in range(20):
        file=open('armatus/chr'+str(n)+'/chr'+str(n)+'_'+str([10,50,200][i])+'_'+str(100+j)+'.consensus.txt').readlines()
        for line in file:
            if int(line.strip().split('\t')[2])/40000-int(line.strip().split('\t')[1])/40000>0:
                TAD['AT'][i][j].append(int(line.strip().split('\t')[1])/40000)
                TAD1['AT'][i][j].append(int(line.strip().split('\t')[2])/40000)




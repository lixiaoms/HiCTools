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
    return(WS(L[:-1],[i-1 for i in L[1:]],l[:-1],[i-1 for i in l[1:]]))

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

##caculate the similarity between TADs identified from down sampled data and raw data
WSs={}
MIs={}
for sf in software:
    WSs[sf]=np.zeros([3,20])
    MIs[sf]=np.zeros([3,20])
for i in range(3):
    for j in range(20):
        for sf in software:
            if sf in software1:
                WSs[sf][i][j]=WS(TAD[sf][i][j],TAD1[sf][i][j],raw[sf],raw1[sf])
                MIs[sf][i][j]=MI(size,np.unique(TAD[sf][i][j]+[x+1 for x in TAD1[sf][i][j]]),np.unique(raw[sf]+[x+1 for x in raw1[sf]]))
            else:
                WSs[sf][i][j]=WSlist(size,np.unique(TAD[sf][i][j]),np.unique(raw[sf]))
                MIs[sf][i][j]=MI(size,np.unique(TAD[sf][i][j]),np.unique(raw[sf]))

##caculate the number of TADs
Num={}
for sf in software:
    Num[sf]=np.zeros([4,20])
for i in range(3):
    for j in range(20):
        for sf in software:
            Num[sf][i+1][j]=len(TAD[sf][i][j])
for j in range(20):
        for sf in software:
            Num[sf][0][j]=len(raw[sf])

##caculate the length of TADs
Length={}
for sf in software:
    Length[sf]=np.zeros([4,20])
for i in range(3):
    for j in range(20):
        for sf in software:
            if sf in software1:
                Length[sf][i+1][j]=np.mean(np.array(sorted(TAD1[sf][i][j]))-np.array(sorted(TAD[sf][i][j]))+1)
            else:
                bound=np.array(sorted(list(set(np.array(TAD[sf][i][j])[np.array(TAD[sf][i][j])<size]))+[0,size]))
                Length[sf][i+1][j]=np.mean(bound[1:]-bound[:-1])
for j in range(20):
        for sf in software:
            if sf in software1:
                Length[sf][0][j]=np.mean(np.array(sorted(raw1[sf]))-np.array(sorted(raw[sf]))+1)
            else:
                bound=np.array(sorted(list(set(np.array(raw[sf])[np.array(raw[sf])<size]))+[0,size]))
                Length[sf][0][j]=np.mean(bound[1:]-bound[:-1])

##caculate the enrichment of chip-seq peaks on TADs boundary
ctcf=[0 for i in range(size)]
for line in open('chip-seq/ctcf_IMR90_hg18.bed').readlines():
    if line.strip().split('\t')[0]=='chr'+str(n):
        pos=(int(line.strip().split('\t')[1])+int(line.strip().split('\t')[2]))/2
        ctcf[pos/40000]+=1
h3k4me3=[0 for i in range(size)]
for line in open('chip-seq/h3k4me3_IMR90_hg18.bed').readlines():
    if line.strip().split('\t')[0]=='chr'+str(n):
        pos=(int(line.strip().split('\t')[1])+int(line.strip().split('\t')[2]))/2
        h3k4me3[pos/40000]+=1
h3k36me3=[0 for i in range(size)]
for line in open('chip-seq/h3k36me3_IMR90_hg18.bed').readlines():
    if line.strip().split('\t')[0]=='chr'+str(n):
        pos=(int(line.strip().split('\t')[1])+int(line.strip().split('\t')[2]))/2
        h3k36me3[pos/40000]+=1

def enrich(bound):
    ratio=np.array([])
    for d in range(-10,11):
        ratio_sum=np.array([])
        for pos in bound:
            i=pos+d
            if i in range(1,size+1):
                try:
                    ratio_sum=np.row_stack((ratio_sum,np.array([ctcf[i-1],h3k4me3[i-1],h3k36me3[i-1]])))
                except:
                    ratio_sum=np.array([ctcf[i-1],h3k4me3[i-1],h3k36me3[i-1]])
        if np.mean(ratio_sum,axis=0).size==1:
            try:
                ratio=np.row_stack((ratio,ratio_sum))
            except:
                ratio=ratio_sum
        else:
            try:
                ratio=np.row_stack((ratio,np.mean(ratio_sum,axis=0)))
            except:
                ratio=np.mean(ratio_sum,axis=0)
    return ratio

def signal(sf):
    sample=[]
    if sf in software1:
        ratio=enrich(np.unique(raw[sf]+[x+1 for x in raw1[sf]]+[0]))
    else:
        ratio=enrich(raw[sf]+[0])
    sample.append(ratio)
    if sf == 'HS':
        return sample
    for i in range(3):
        if sf in software1:
            ratio=enrich(np.unique(TAD[sf][i][0]+[x+1 for x in TAD1[sf][i][0]]+[0]))
            for j in range(1,20):
                ratio=ratio+enrich(np.unique(TAD[sf][i][j]+[x+1 for x in TAD1[sf][i][j]]+[0]))
        else:
            ratio=enrich(TAD[sf][i][0]+[0])
            for j in range(1,20):
                ratio=ratio+enrich(TAD[sf][i][j]+[0])
        sample.append(ratio)
    return sample

signal_chipseq={}
for sf in software:
    signal_chipseq[sf]=signal(sf)

def zscore(signal):
    peak=np.array(signal[9:12])
    background=np.array(signal[[0,1,2,3,4,-5,-4,-3,-2,-1]])
    return (np.mean(peak)-np.mean(background))*(np.std(background))**(-1)

##figure of exsample
fig,axes=plt.subplots(figsize=(8,8))
for sf in software:
    axes.boxplot(MIs[sf].T,showmeans=True,showfliers=False)
    axes.plot(range(0,4),np.concatenate(([1],np.mean(MIs[sf],axis=1))),label=name[sf])
plt.legend(loc=3)
plt.xlim((0,3.2))
plt.ylim(0,1)
plt.setp(axes, xticks=range(4),xticklabels=['raw depth','1_10 depth','1_50 depth','1_200 depth'])
plt.title('Mutual Information',fontsize=20)

fig,axes=plt.subplots(figsize=(8,8))
for sf in software:
    axes.boxplot(WSs[sf].T,showmeans=True,showfliers=False)
    axes.plot(range(0,4),np.concatenate(([1],np.mean(WSs[sf],axis=1))),label=name[sf])
plt.legend(loc=1)
plt.xlim((0,3.2))
plt.setp(axes, xticks=range(4),xticklabels=['raw depth','1_10 depth','1_50 depth','1_200 depth'])
plt.title('Weight Similarity',fontsize=20)

fig,axes=plt.subplots(figsize=(8,8))
for sf in software:
    axes.plot(range(1,5),np.mean(Num[sf],axis=1),'x-',label=name[sf])
plt.legend(loc=1,handleheight=0)
plt.xlim((0.8,4.2))
plt.ylim(0,400)
plt.setp(axes, xticks=range(1,5),xticklabels=['raw depth','1_10 depth','1_50 depth','1_200 depth'])
plt.title('Number of TADs',fontsize=20)

fig,axes=plt.subplots(figsize=(8,8))
for sf in software:
    axes.plot(range(1,5),np.mean(Length[sf],axis=1),'x-',label=name[sf])
plt.legend(loc=1)
plt.xlim((0.8,4.2))
plt.ylim(0,200)
plt.setp(axes, xticks=range(1,5),xticklabels=['raw depth','1_10 depth','1_50 depth','1_200 depth'])
plt.setp(axes, yticks=range(25,200,25),yticklabels=['1M','','3M','','5M','','7M'])
plt.title('Length of TADs',fontsize=20)

plt.subplots(figsize=(15,6),nrows=1, ncols=3)
plt.subplot(131)
fc={'depth_10':[],'depth_200':[],'depth_50':[]}
for i in range(5):
    fc['depth_10'].append(zscore(signal_chipseq[software[i]][1][:,0]))
for i in range(5):
    fc['depth_50'].append(zscore(signal_chipseq[software[i]][2][:,0]))
for i in range(5):
    fc['depth_200'].append(zscore(signal_chipseq[software[i]][3][:,0]))
plt.bar(range(0,15,3),fc['depth_10'],label='1_10 depth')
plt.bar(range(1,16,3),fc['depth_50'],label='1_50 depth')
plt.bar(range(2,17,3),fc['depth_200'],label='1_200 depth')
plt.ylim(-2,10)
plt.legend()
plt.xticks([3*x+1 for x in range(5)],[name[sf] for sf in software],rotation=270)
plt.title('CTCF peaks z-score')
plt.subplot(132)
fc={'depth_10':[],'depth_200':[],'depth_50':[]}
for i in range(5):
    fc['depth_10'].append(zscore(signal_chipseq[software[i]][1][:,1]))
for i in range(5):
    fc['depth_50'].append(zscore(signal_chipseq[software[i]][2][:,1]))
for i in range(5):
    fc['depth_200'].append(zscore(signal_chipseq[software[i]][3][:,1]))
plt.bar(range(0,15,3),fc['depth_10'],label='1_10 depth')
plt.bar(range(1,16,3),fc['depth_50'],label='1_50 depth')
plt.bar(range(2,17,3),fc['depth_200'],label='1_200 depth')
plt.ylim(-2,10)
plt.legend()
plt.xticks([3*x+1 for x in range(5)],[name[sf] for sf in software],rotation=270)
plt.title('H3K4me3 peaks z-score')
plt.subplot(133)
fc={'depth_10':[],'depth_200':[],'depth_50':[]}
for i in range(5):
    fc['depth_10'].append(zscore(signal_chipseq[software[i]][1][:,2]))
for i in range(5):
    fc['depth_50'].append(zscore(signal_chipseq[software[i]][2][:,2]))
for i in range(5):
    fc['depth_200'].append(zscore(signal_chipseq[software[i]][3][:,2]))
plt.bar(range(0,15,3),fc['depth_10'],label='1_10 depth')
plt.bar(range(1,16,3),fc['depth_50'],label='1_50 depth')
plt.bar(range(2,17,3),fc['depth_200'],label='1_200 depth')
plt.ylim(-2,10)
plt.legend()
plt.xticks([3*x+1 for x in range(5)],[name[sf] for sf in software],rotation=270)
plt.title('H3K36me3 peaks z-score')

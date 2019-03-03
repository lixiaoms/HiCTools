import pandas as pd
import matplotlib.pyplot as plt  
import numpy as np

loop={}
peak=pd.DataFrame(columns=['Frag1','Frag2','rawp'])
for j in range(1,185):
    try:
        prob=pd.read_csv('../../review/fasthic/rao_chr10/0/TAD_'+str(j)+'_prob',sep='\t')
        prob['rawp']=prob['PeakProbability']
        peak=pd.merge(peak,prob[['Frag1','Frag2','rawp']],how='outer')
    except:
        continue
loop['raw']=peak

##cell#18 is pooled cells
for i in range(2,8)+range(9,18):
    peak=pd.DataFrame(columns=['Frag1','Frag2','pvalue','count'])
    for j in range(1,185):
        try:
            prob=pd.read_csv('rao_chr10/cell_'+str(i)+'/TAD_'+str(j)+'_prob',sep='\t')
            prob['pvalue']=prob['PeakProbability']
            prob['count']=prob['ObservedCount']
            peak=pd.merge(peak,prob[['Frag1','Frag2','pvalue','count']],how='outer')
        except:
            continue
    peak1=loop['raw']
    peak=pd.merge(peak,peak1,on=['Frag1','Frag2'],how='outer')
    peak['dis']=peak['Frag2']-peak['Frag1']
    loop[i]=peak

lap1=[]
num1=[]
lap2=[]
num2=[]

##caculate the recall and precision in each cells
for i in range(2,8)+range(9,18):
    peak=loop[i]
    peak['dis']=peak['Frag2']-peak['Frag1']
    num1.append(len(np.where(np.array(np.minimum(peak['pvalue'],peak['rawp'])>0.95)*np.array(peak['count']>0))[0])*len(np.where(np.array(peak['rawp']>0.95)*np.array(peak['count']>0))[0])**(-1))
    lap1.append(len(np.where(np.array(np.minimum(peak['pvalue'],peak['rawp'])>0.95)*np.array(peak['count']>0))[0])*len(np.where(np.array(peak['pvalue']>0.95)*np.array(peak['count']>0))[0])**(-1))
    num2.append(len(np.where(np.array(peak['count']>0)*np.array(peak['dis']>50)*np.array(peak['rawp']>0.95))[0])*len(np.where(np.array(peak['rawp']>0.95)*np.array(peak['count']>0))[0])**(-1))
    lap2.append(len(np.where(np.array(peak['count']>0)*np.array(peak['dis']>50)*np.array(peak['rawp']>0.95))[0])*len(np.where(np.array(peak['count']>0)*np.array(peak['dis']>50))[0])**(-1))

##figure of example
labels=['cell'+str(i) for i in [2, 3, 4, 5, 6, 7, 9, 10, 11, 12,13,14,15,16,17]]
angles = np.linspace(0, 2*np.pi, 15, endpoint=False)
angles = np.concatenate((angles, [angles[0]]))
fig = plt.figure(figsize=(10,6))
fig.subplots_adjust(hspace=0.48,wspace=0.3)
ax = fig.add_subplot(121, polar=True)
ax.plot(angles, lap1+[lap1[0]], 'ro-', linewidth=2,label='FastHiC')
ax.plot(angles, lap2+[lap2[0]], 'bo-', linewidth=2,label='naive caller')
ax.set_thetagrids(angles * 180/np.pi, labels)
ax.set_ylim(0,1)
ax.set_yticks([0,0.2,0.4,0.6,0.8])
ax.set_title("precision",pad=10,size=15)
ax.grid(True)
ax.legend(bbox_to_anchor=(1.4, 1.12))
ax = fig.add_subplot(122, polar=True)
ax.plot(angles, num1+[num1[0]], 'ro-', linewidth=2,label='FastHiC')
ax.plot(angles, num2+[num2[0]], 'bo-', linewidth=2,label='naive caller')
ax.set_ylim(0,0.4)
ax.set_yticks([0,0.1,0.2,0.3])
ax.set_thetagrids(angles * 180/np.pi, labels)
ax.set_title("recall",pad=10,size=15)
ax.grid(True)
plt.show()

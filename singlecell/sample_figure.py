## cell#11: chr10: 64-70Mbp
import os,sys
import platform
import matplotlib
%matplotlib inline
from math import sqrt, isnan, floor, ceil, pi
from numpy import log2, array, max
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker  import MultipleLocator
from matplotlib.patches import Polygon, Rectangle, Circle
from scipy.signal import argrelextrema
from scipy import ndimage
import scipy.sparse as sps
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import argparse
import bisect
import warnings
import logging

matrix1=np.genfromtxt('TAD/chr10/cell_11',skip_header=1,filling_values='0',skip_footer=0)
matrix2=np.genfromtxt('TAD/chr10/cell_18',skip_header=1,filling_values='0',skip_footer=0)
matrix3=np.genfromtxt('TAD/chr10/raw',skip_header=1,filling_values='0',skip_footer=0)
matrix4=np.tril(matrix1,-1)+np.triu(matrix2,-1)
matrix5=np.tril(matrix2,-1)*50+np.triu(matrix3,-1)

cell_11_pcc=[]
for line in open('compartment/chr10_j/cell_11_eigen.txt').readlines():
    cell_11_pcc.append(float(line.strip()))
cell_pool_pcc=[]
for line in open('compartment/chr10_j/cell_18_eigen.txt').readlines():
    cell_pool_pcc.append(float(line.strip()))
cell_rao_pcc=[]
for line in open('../review/Compartment/Juicer/rao_chr10/raw.txt').readlines():
    cell_rao_pcc.append(float(line.strip()))

TAD_11=[]
TAD_pool=[]
TAD_rao=[]
file=open('TAD/chr10/cell_11.is520001.ids320001.insulation.boundaries').readlines()
for line in file[1:]:
    TAD_11.append(int(line.split('\t')[3]))
file=open('TAD/chr10/cell_18.is520001.ids320001.insulation.boundaries').readlines()
for line in file[1:]:
    TAD_pool.append(int(line.split('\t')[3]))
file=open('TAD/chr10/raw.is520001.ids320001.insulation.boundaries').readlines()
for line in file[1:]:
    TAD_rao.append(int(line.split('\t')[3]))

numOfcols=2
numOfrows=6
cmaps = ['Greys','Reds','YlOrBr','YlOrRd','hot']
orientation='lower'
start1=1600
end1=1750
start2=1600
end2=1750
fig=plt.figure(figsize=(numOfcols*5+2.5, numOfrows+numOfrows/2+0.5))
fig.set_size_inches(numOfcols*5+2.5, numOfrows+numOfrows/2+0.5)
fig.subplots_adjust(hspace=0.48,wspace=0.5)
g=gridspec.GridSpec(numOfrows, 4*numOfcols)
ax1 = plt.subplot2grid((numOfrows, 4*numOfcols), (0, 0), rowspan=4,colspan=4)
length = end1-start1
img=ax1.imshow(matrix4[start1:end1,start1:end1],vmin=0,vmax=5,cmap=plt.get_cmap(cmaps[1]),origin=orientation,interpolation="nearest",extent=(int(start1 or 1) - 0.5,int(start1 or 1) + length - 0.5,int(start1 or 1) - 0.5,int(start1 or 1) + length - 0.5),aspect='auto')
divider = make_axes_locatable(ax1)
ax1.get_xaxis().set_label_coords(0.5,-0.06)
ax1.set_ylabel('cell #11')
ax1.set_xlabel('cell pool')
cax = divider.append_axes("bottom", size="2.5%", pad=0.42)
plt.colorbar(img,cax=cax,ticks=MultipleLocator(1),orientation='horizontal',extendfrac='auto',spacing='uniform')
ax1.set_xticks(range(start1,end1,25))
ax1.set_xticklabels([str(i)+'M' for i in range(start1/25,end1/25)])
ax1.set_yticks(range(start1,end1,25))
ax1.set_yticklabels([str(i)+'M' for i in range(start1/25,end1/25)])
a=1
for b in TAD_pool:
    if a in range(start1,end1) and b in range(start1,end1):
        ax1.plot(range(a,b),[a]*(b-a),linewidth=2,color='g')
        ax1.plot([b]*(b-a),range(a,b),linewidth=2,color='g')
    a=b
a=1
for b in TAD_11:
    if a in range(start1,end1) and b in range(start1,end1):
        ax1.plot(range(a,b),[b]*(b-a),linewidth=2,color='y')
        ax1.plot([a]*(b-a),range(a,b),linewidth=2,color='y')
    a=b
ax1.set_title('Hi-C matrix')
ax2 = plt.subplot2grid((numOfrows, 4*numOfcols), (4, 0), rowspan=1,colspan=4,sharex=ax1)
ax2.fill_between([x*2.5 for x in range(int(ceil(start1/2.5)),int(ceil(end1/2.5)))], array([cell_pool_pcc[x] for x in range(int(ceil(start1/2.5)),int(ceil(end1/2.5)))]),0, color='black', interpolate=True)
ax2.set_title('compartment score of cell pool',pad=2)
ax2.set_yticks([1,0,-1])
ax3 = plt.subplot2grid((numOfrows, 4*numOfcols), (5, 0), rowspan=1,colspan=4,sharex=ax1)
ax3.set_xlim(int(start1 or 1) - 0.5,int(start1 or 1) + length - 0.5)
ax3.fill_between([x*2.5 for x in range(int(ceil(start1/2.5)),int(ceil(end1/2.5)))], array([cell_11_pcc[x] for x in range(int(ceil(start1/2.5)),int(ceil(end1/2.5)))]),0, color='black', interpolate=True)
ax3.set_title('compartment score of cell #11',pad=2)
ax3.set_yticks([1,0,-1])
ax4 = plt.subplot2grid((numOfrows, 4*numOfcols), (0, 4), rowspan=4,colspan=4,sharex=ax1)
length = end2-start2
img=ax4.imshow(matrix5.T[start2:end2,start2:end2],vmin=0,vmax=250,cmap=plt.get_cmap(cmaps[1]),origin=orientation,interpolation="nearest",extent=(int(start2 or 1) - 0.5,int(start2 or 1) + length - 0.5,int(start2 or 1) - 0.5,int(start2 or 1) + length - 0.5),aspect='auto')
divider = make_axes_locatable(ax4)
ax4.get_xaxis().set_label_coords(0.5,-0.06)
ax4.set_ylabel('ensemble data in Rao et al.')
ax4.set_xlabel('cell pool')
cax = divider.append_axes("bottom", size="2.5%", pad=0.42)
plt.colorbar(img,cax=cax,ticks=MultipleLocator(50),orientation='horizontal',extendfrac='auto',spacing='uniform')
ax4.set_xticks(range(start2,end2,25))
ax4.set_xticklabels([str(i)+'M' for i in range(start2/25,end2/25)])
#ax4.set_yticks(range(start2,end2,25))
#ax4.set_yticklabels([str(i)+'M' for i in range(start2/25,end2/25)])
ax4.set_yticks([])
a=1
for b in TAD_pool:
    if a in range(start2,end2) and b in range(start2,end2):
        ax4.plot(range(a,b),[a]*(b-a),linewidth=2,color='g')
        ax4.plot([b]*(b-a),range(a,b),linewidth=2,color='g')
    a=b
a=1
for b in TAD_rao:
    if a in range(start2,end2) and b in range(start2,end2):
        ax4.plot(range(a,b),[b]*(b-a),linewidth=2,color='y')
        ax4.plot([a]*(b-a),range(a,b),linewidth=2,color='y')
    a=b
ax4.set_xlim(int(start2 or 1) - 0.5,int(start2 or 1) + length - 0.5)
ax4.set_ylim(int(start2 or 1) - 0.5,int(start2 or 1) + length - 0.5)
ax4.set_title('Hi-C matrix')
ax5 = plt.subplot2grid((numOfrows, 4*numOfcols), (4, 4), rowspan=1,colspan=4,sharex=ax4)
ax5.fill_between([x*2.5 for x in range(int(ceil(start2/2.5)),int(ceil(end2/2.5)))], array([cell_pool_pcc[x] for x in range(int(ceil(start2/2.5)),int(ceil(end2/2.5)))]),0, color='black', interpolate=True)
ax5.set_title('compartment score of cell pool',pad=2)
ax5.set_yticks([1,0,-1])
ax6 = plt.subplot2grid((numOfrows, 4*numOfcols), (5, 4), rowspan=1,colspan=4,sharex=ax4)
ax6.set_xlim(int(start2 or 1) - 0.5,int(start2 or 1) + length - 0.5)
ax6.fill_between([x*2.5 for x in range(int(ceil(start2/2.5)),int(ceil(end2/2.5)))], array([cell_rao_pcc[x] for x in range(int(ceil(start2/2.5)),int(ceil(end2/2.5)))]),0, color='black', interpolate=True)
ax6.set_title('compartment score of ensemble data',pad=2)
ax2.set_ylim(-0.05,0.05)
ax3.set_ylim(-0.04,0.04)
ax5.set_ylim(-0.05,0.05)
ax6.set_ylim(-0.05,0.05)
ax6.set_yticks([0])
plt.show()

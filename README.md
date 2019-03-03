Methods for the Hi-C data normalization

1 ICE

Iterative correction was performed by our own R codes. We removed from consideration contacts within a bin and with adjacent bins. We then removed the bins having no contacts and 2% of bins having the fewest number of contacts, which is suggested by the original paper. 

2 HiCNorm

The codes are available at http://www.people.fas.harvard.edu/~junliu/HiCNorm/. Required genomic feature files were downloaded from http://enhancer.sdsc.edu/yunjiang/resources/genomic_features/. Since both methods are executed by R, we can compare its running time.

Methods for the analysis of compartments

1 Juicer

Eigenvector is part of Juicer Tools (https://github.com/aidenlab/juicer/wiki/Download) for being used to delineate compartments. We calculated the eigenvector with None normalization. Bins with NA eigenvalue are not assigned with compartment A or B.

2 CscoreTool

CscoreTool is available at https://github.com/scoutzxb/CscoreTool. Parameter minDis was set to 1M, as suggested by author. Bins with score 0 are not assigned with compartment A or B.

3 GeSICA

GeSICA is available at https://zhanglab.tongji.edu.cn/Softwares.htm. All parameters were left as default. Bins with score 0 are not assigned with compartment A or B.

Methods for the identification of TAD

1 InsulationScore

InsulationScore is available at https://github.com/dekkerlab/crane-nature-2015. All parameters were left as default.

2 HOMER

FindTADsAndLoops is part of HOMER (http://homer.ucsd.edu/homer/) for finding TADs and loops. Parameters res and window were both set to 40kb.

3 MrTADFinder

MrTADFinder is available at https://github.com/gersteinlab/MrTADFinder. Parameter res was set to 2.875, as suggested by the original paper. 

4 HiCseg

HiCseg (https://cran.r-project.org/web/packages/HiCseg/) is an R package. Arguments distrib and model were set to ‘B’ and ‘D’, as used by the original paper.

5 Armatus

Armatus (https://github.com/kingsfordgroup/armatus) was running in ubuntu14.04 while other methods’ code can be running in newer version. Parameters gammaMax and stepSize were set to 0.5 and 0.05, as suggested by the original paper. 

6 deDoc(E)

deDoc is available at https://github.com/yinxc/structural-information-minimisation. TAD have less than 5 bins were removed, as suggested by the author.

Methods for the identification of loop

1 HiCCUPS

HiCCUPS is part of Juicer Tools for identifying loop. We used modified version of HiCCUPS which can run on CPUs with the addition of the –-ignore_sparsity flag.

2 Fit-Hi-C

Fit-Hi-C is available at https://github.com/ay-lab/fithic. Bias file generated by HiCKRy.py was incorporated in the command by argument –t. Only interactions with a FDR<0.05 were considered as loop, as suggested by the author.

3 FastHiC

FastHiC is available at http://www.unc.edu/~yunmli/FastHiC/. Expected counts under random collision events were estimated ICE. We used InsulationScore to identify TAD and analyzed each TAD separately to detect intra-TAD chromatin interactions. Only interactions with a PeakProbability>0.95 were considered as loop.

4 diffHiC

diffHic24 (v.1.8.1) is an R package available in Bioconductor. We set the size of the neighborhood (flank parameter) at 50kb. Function filterPeaks‘s arguments min.enrich, min.count and min.diag were set to 0.5, 5 and 2L respectively, as suggested in the user guide.

5 cLoops

cLoops is available at https://github.com/YaqiangCao/cLoops. Key parameters were set to ‘-eps 5000,10000 -minPts 50,100 –hic’ that is suggested by the author. We considered loops_juicebox.txt as final loops.

6 HOMER

Parameters res and window were both set to 10kb.

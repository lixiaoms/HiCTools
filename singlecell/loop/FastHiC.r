library(parallel)
library(doParallel)

chrlen=135534747
index=c(2:7,9:17)

decay=function(G){
decay=c()
for ( x in 1:250){
    decay=c(decay,mean(G[col(G)-row(G)==x]))
}
return(decay)}

fastin=function(F,E,start){
A=data.frame(frag1=c(),frag2=c(),Oij=c(),Eij=c())
if (ncol(F)>3){
for ( i in 1:(ncol(F)-2)){
    for (j in (i+2):ncol(F)){
        A=rbind(A,c(i+start,j+start,F[i,j],E[j-i]))
    }
}
    }
colnames(A)=c('frag1','frag2','Oij','Eij')
return(A)}

insulation=read.table('../../review/fasthic/rao_chr10/raw.is500001.ids260001.insulation.boundaries')
TAD=as.numeric(as.vector(insulation[2:nrow(insulation),4]))
TAD=c(1,TAD,ncol(F))
TADindex=which(TAD[-1]-TAD[-length(TAD)]>1)
TADindex=TADindex[-47]

##multicore-caculate the expected number of reads and employ FastHiC tool in each cells 
for (i in index){
    system(paste0('mkdir rao_chr10/cell_',i))
    hic=read.table(paste0('../sunnyxie/clean_',i))
    F=matrix(0,13554,13554)
    for (j in 1:nrow(hic)){
    a=strsplit(as.character(hic[j,1]),',')[[1]]
    b=strsplit(as.character(hic[j,2]),',')[[1]]
    if (a[1]=='10' & b[1]=='10'){
        pos1=ceiling(as.numeric(a[2])/10000)
        pos2=ceiling(as.numeric(b[2])/10000)
        F[pos1,pos2]=F[pos1,pos2]+1
    }
    }
    F[lower.tri(F)]=t(F)[lower.tri(t(F))]
    E=decay(F)
    mc <- getOption("mc.cores", 6)
    func=function(x){
        A=fastin(F[TAD[x]:TAD[x+1],TAD[x]:TAD[x+1]],E,TAD[x]-1)
        write.table(A,file=paste0('rao_chr10/cell_',i,'/TAD_',x),quote = FALSE,row.names=FALSE)
        system(paste0('java -jar ../../software/callpeak_fasthic/FastHiC.jar --interaction rao_chr10/cell_',i,'/TAD_',x,' --outprob rao_chr10/cell_',i,'/TAD_',x,'_prob --outestimator rao_chr10/cell_',i,'/TAD_',x,'_esti'))
    }
    mclapply(TADindex,func,mc.cores=mc)
}
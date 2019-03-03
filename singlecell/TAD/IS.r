chrlen=135534747

for (j in c(2:7,9:17)){
hic=read.table(paste0('../sunnyxie/clean_',j))
F=matrix(0,ceiling(chrlen/40000),ceiling(chrlen/40000))
for (i in 1:nrow(hic)){
    a=strsplit(as.character(hic[i,1]),',')[[1]]
    b=strsplit(as.character(hic[i,2]),',')[[1]]
    if (a[1]=='19' & b[1]=='19'){
        pos1=ceiling(as.numeric(a[2])/40000)
        pos2=ceiling(as.numeric(b[2])/40000)
        F[pos1,pos2]=F[pos1,pos2]+1
    }
}
F[lower.tri(F)]=t(F)[lower.tri(t(F))]
assign(paste0('F',j),F)}

##pooled cell index 18
index=c(2:7,9:18)
F18=matrix(0,ceiling(chrlen/40000),ceiling(chrlen/40000))
for (i in index[-length(index)]){
    F18=F18+get(paste0('F',i))
}

##employ InsulationScore tool
for (j in index){
    G=get(paste0('F',j))
    rownames(G)=coor
    colnames(G)=coor
    colnames(G)[1]=paste0('	',coor[1])
    write.table(G,paste0('chr19/cell_',j),quote=FALSE,sep='\t')
    system(paste0('perl ../../software/crane-nature-2015/scripts/matrix2insulation.pl -i chr19/cell_',j))
}

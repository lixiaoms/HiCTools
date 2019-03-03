index=c(2:7,9:17) ##index of cells
for (j in index){
hic=read.table(paste0('../sunnyxie/clean_',j))
cat('',file=paste0('chr10_j/cell_',j,'.txt'))
for (i in 1:nrow(hic)){
    a=strsplit(as.character(hic[i,1]),',')[[1]]
    b=strsplit(as.character(hic[i,2]),',')[[1]]
    if (a[1]=='10' & b[1]=='10'){
        pos1=as.numeric(a[2])
        pos2=as.numeric(b[2])
        data=paste0(0,'\t','chr10','\t',pos1,'\t',0,'\t',0,'\t','chr10','\t',pos2,'\t',1,'\n')
        cat(data,file=paste0('chr10_j/cell_',j,'.txt'),append=TRUE)
    }
}
system(paste0('gzip chr10_j/cell_',j,'.txt'))
system(paste0('java -Xmx30g -jar /home/lixiao/software/juicer_tools.jar pre chr10_j/cell_',j,'.txt.gz chr10_j/cell_',j,'.hic hg19'))
system(paste0('gzip -d chr10_j/cell_',j,'.txt.gz'))
system(paste0('java -Xmx30g -jar /home/lixiao/software/juicer_tools.jar eigenvector -p  NONE chr10_j/cell_',j,'.hic 10 BP 100000 chr10_j/cell_',j,'_eigen.txt'))
}

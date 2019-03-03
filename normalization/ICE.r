ICEnorm=function(A, max_iter=200, eps=1e-5, sparse.filter=0.02){
    #change diag and nearby bins to 0
    A[abs(col(A)-row(A))<=1]=0
    #remove col and row with sum equal 0 or few reads
    zeros = unique(which(colSums(A) ==0), which(rowSums(A) ==0))
    sparse= which(rank(colSums(A))<(ncol(A)-length(zeros))*sparse.filter+length(zeros))
    print(length(zeros))
    print(length(sparse))
  if (length(sparse) > 0) {
    N = A[-sparse, -sparse]
    message(paste0('Cols/Rows removed: '))
    message(paste(" ", sparse, sep = " "))
  }
    bias=rep(1,ncol(A))
    #iteration
    for (i in 1:max_iter){
    S=colSums(N)
    S=S/mean(S)
    N= matrix(as.vector(N)*kronecker(S**(-1), S**(-1)),nrow(N),ncol(N))
    bias[-sparse]=bias[-sparse]*S**(-1)
    if (max(abs(log(S)))<eps){break}
    }  
    #recover col and row deleted
    A[-sparse,-sparse]=N
    A[sparse,]=0
    A[,sparse]=0
    return(c(i,bias))
}


##example on down sampled chr19 data in Rao et al.
downsample=function(matrix,nx){
    reads=matrix[!lower.tri(matrix)]
    sample=rep(0,length(reads))
    index=sort(sample(1:sum(reads),sum(reads)/nx))
    i=0
    sum=reads[1]
    for (j in 1:length(index)){
        for (k in i:length(reads)){
            if (sum>=index[j]){
                i=k
                break
            }
            sum=sum+reads[k+1]}
        sample[k]=sample[k]+1
    }
    S=matrix
    S[!lower.tri(S)]=sample
    S[lower.tri(S)]=t(S)[lower.tri(t(S))]
    return(S)
    }

F=read.table(paste('../GM12878/hic/GM12878_combined/10kb_resolution_intrachromosomal/chr19/MAPQGE30/chr19_10kb.RAWobserved',sep=''))
G=matrix(0,1479,1479)
for (x in 1:nrow(F)){
    a=F[x,1]%/%40000
    b=F[x,2]%/%40000
    if (a!=b){
    G[a,b]=G[a,b]+F[x,3]
    G[b,a]=G[b,a]+F[x,3]}
    if (a==b){
    G[a,a]=G[a,a]+F[x,3]}
}

ICEtest=matrix(0,50,4)
for (i in 1:50){
    sample=matrix(0,11,3)
    icenormtime=c()
    for (j in 1:11){
    set.seed(j)
    G1=downsample(as.matrix(G),(1.2)**i)
    ptm <- proc.time()
    res=ICEnorm(G1)
    icenormtime=c(icenormtime,proc.time() - ptm)
    sample[j,]=c(cor(log(C[-1]),log(res[-1])),var(log(res[-1])),res[1])}
    ##record median result of 11 seeds
    ICEtest[i,]=c(median(sample[,1]),median(sample[,2]),median(sample[,3]),median(icenormtime[seq(3,53,5)]))
}

##figure of example
rate <- 1:50
t <- 1.2
par(mar=c(5, 4, 4, 6) + 0.1)
plot(rate, ICEtest[,1], pch=16, axes=FALSE, ylim=c(0,1), xlab="", ylab="", 
   type="o",col="black", main="biases property")
axis(2, ylim=c(0,0.3),col="black",las=1)
mtext("biases correlation",side=2,line=2.5)
box()
par(new=TRUE)
plot(rate, ICEtest[,2], pch=15,  xlab="", ylab="",  
    axes=FALSE, type="o", col="black")
mtext("biases variance",side=4,col="black",line=2.5) 
axis(4, ylim=c(0,10), col="black",col.axis="black",las=1)
axis(1,c(log(4)/log(t), log(16)/log(t), log(64)/log(t), log(256)/log(t), log(1024)/log(t), log(4096)/log(t)),c('1_4','1_16','1_64','1_256','1_1024','1_4096'),las = 1)
mtext('reads propotion',side=1,col="black",line=2.5)
legend("left",legend=c("biases correlation","biases variance",'minimal depth to converge'),
  text.col=c("black","black",'red'),pch=c(16,15,15),col=c("black","black",'red'))
abline(v=47, lty=3, lwd=2, col="red")

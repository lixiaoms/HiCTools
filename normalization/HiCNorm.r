##modified version of the code from http://www.people.fas.harvard.edu/~junliu/HiCNorm/
v<-read.table("../software/HiCNorm/data/chr19_hg19_40kb") ##genomic feature file
for (i in range(ncol(v))){v[,i]=as.numeric(v[,i])}

Hicnorm=function(G){
u<-G

#change matrix into vector
u_vec<-u[upper.tri(u)]

#get cov matrix
if (length(which(v[,6]==0))>0){
len_m<-as.matrix(log(v[-which(v[,6]==0),4]%o%v[-which(v[,6]==0),4]))
gcc_m<-as.matrix(log(v[-which(v[,6]==0),5]%o%v[-which(v[,6]==0),5]))
map_m<-as.matrix(log(v[-which(v[,6]==0),6]%o%v[-which(v[,6]==0),6]))
}
if (length(which(v[,6]==0))==0){
len_m<-as.matrix(log(v[,4]%o%v[,4]))
gcc_m<-as.matrix(log(v[,5]%o%v[,5]))
map_m<-as.matrix(log(v[,6]%o%v[,6]))
}

#centralize cov matrix of enz, gcc
len_m<-(len_m-mean(c(len_m)))/sd(c(len_m))
gcc_m<-(gcc_m-mean(c(gcc_m)))/sd(c(gcc_m))

#change matrix into vector
len_vec<-len_m[upper.tri(len_m)]
gcc_vec<-gcc_m[upper.tri(gcc_m)]
map_vec<-map_m[upper.tri(map_m)]

#fit Poisson regression: u~len+gcc+offset(map)
fit<-glm(u_vec~len_vec+gcc_vec+offset(map_vec),family="poisson")

#user can use the following two lines to fit negative binomial regression: u~len+gcc+offset(map).
#library("MASS")
#fit<-glm.nb(u_vec~len_vec+gcc_vec+offset(map_vec))

#summary(fit)
coeff<-round(fit$coeff,4)
res<- round(u/exp(coeff[1]+coeff[2]*len_m+coeff[3]*gcc_m+map_m), 4)

return(fit$coefficient)
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

G1=G
normtest=matrix(0,101,2)
hicnormtime=c()
ptm <- proc.time()
coe=Hicnorm(G1)
hicnormtime=c(hicnormtime,proc.time() - ptm)
normtest[1,]=c(coe[[2]],coe[[3]])
k=2
for (i in c(10,50,200,1000,5000)){
    for (j in 100:119){
        set.seed(j)
        G2=downsample(G1,i)
        ptm <- proc.time()
        coe=Hicnorm(G2)
        hicnormtime=c(hicnormtime,proc.time() - ptm)
        normtest[k,]=c(coe[[2]],coe[[3]])
        k=k+1
    }
}

##figure of example
library(ggplot2)
test=matrix(0,101,3)
test[,1:2]=normtest
test[,3]=hicnormtime[seq(3,503,5)]
data <-data.frame(length=normtest[,1], gc=normtest[,2],reads_depth=factor(c('raw',rep(c('1_10','1_50','1_200','1_1000','1_5000'),each=20)),levels=c('raw','1_10','1_50','1_200','1_1000','1_5000')))
ggplot(data, aes(x = length, y = gc, colour = reads_depth)) +geom_point()+ theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+xlab('Effective length')+ylab('GC content')+scale_fill_discrete(name="reads depth")+xlim(0.07,0.27)+ylim(-0.1,-0.28)+theme_bw()+theme(legend.position=c(0.1,0.15))


#+eval=FALSE
#Library
library(fpc)
library(R.matlab)
library('cluster')
library(teigen)
library(mvtnorm)
library(MixGHD)

#Functions(these are self written)
#-------------------------------------------------------------------------------------------------------------------------------------
#K Medoids benchmark test using silhouette index
kmed.si.bench=function(X,g,G){
  if(is.list(X)&g%%1==0&G%%1==0&g>=2&G>g){
    SI=vector(length=length(X))
    nclust=vector(length=length(X))
    cluster.idx<-list()
    for(i in 1:length(X)){
      temp.si<-vector(length=G-g+1)
      for(j in g:G){
        pam.out<-pam(X[[i]],k=j)
        temp.si[j-g+1]<-pam.out$silinfo$avg.width
      }
      pam.out<-pam(X[[i]],k=which.max(temp.si)+g-1)
      cluster.idx[[i]]<-pam.out$clustering
      nclust[i]<-which.max(temp.si)+g-1
      SI[i]<-pam.out$silinfo$avg.width
    }
    names(cluster.idx)<-as.character(1:length(X))
    bench.list<-list()
    bench.list[[1]]<-SI
    bench.list[[2]]<-nclust
    bench.list[[3]]<-cluster.idx
    names(bench.list)<-c('SI','G','Cluster Index')
    return(bench.list)
  }else{
    return('Error: Double check arguments')
  }
}

#K medoids benchmark test using ARI
kmed.ari.bench=function(X,Y,g,G){
  if(is.list(X)&is.list(Y)&g%%1==0&G%%1==0&g>=2&G>g){
    adjusted.rand.index=vector(length=length(X))
    nclust=vector(length=length(X))
    cluster.idx<-list()
    for(i in 1:length(X)){
      temp.ari<-vector(length=G-g+1)
      for(j in g:G){
        pam.out<-pam(X[[i]],k=j)
        temp.ari[j-g+1]<-ARI(pam.out$clustering[1:length(Y[[i]])],Y[[i]])
      }
      pam.out<-pam(X[[i]],k=which.max(temp.ari)+g-1)
      cluster.idx[[i]]<-pam.out$clustering
      nclust[i]<-which.max(temp.ari)+g-1
      adjusted.rand.index[i]<-max(temp.ari)
    }
    names(cluster.idx)<-as.character(1:length(X))
    bench.list<-list()
    bench.list[[1]]<-adjusted.rand.index
    bench.list[[2]]<-nclust
    bench.list[[3]]<-cluster.idx
    names(bench.list)<-c('ARI','G','Cluster Index')
    return(bench.list)
  }else{
    return('Error: Double check arguments')
  }
}

#Student T distribution benchmark tests(SI and ARI in one since optimization is chose via BIC)
t.bench=function(X,Y,g,G){
  if(is.list(X)&is.list(Y)&g%%1==0&G%%1==0&g>=2&G>g){
    SI<-vector(length=length(X))
    adjusted.rand.index<-vector(length=length(X))
    nclust<-vector(length=length(X))
    cluster.idx<-list()
    count<-0
    for(i in 1:length(X)){
      t.out<-teigen(x=X[[i]],Gs=g:G,models="UUUU")
      if(t.out$G==1){
        SI[i]<-0
        adjusted.rand.index[i]<-0
        nclust[i]<-1
        cluster.idx[[i]]<-NA
      }else{
      d<-dist(X[[i]])
      temp.sval<-silhouette(t.out$classification,dist=dist(X[[i]]))
      SI[[i]]<-mean(temp.sval[,3])
      adjusted.rand.index[i]<-ARI(t.out$classification[1:length(Y[[i]])],Y[[i]])
      nclust[i]<-t.out$G
      cluster.idx[[i]]<-t.out$classification
      }
    }
    names(cluster.idx)<-as.character(1:length(X))
    bench.list<-list()
    bench.list[[1]]<-SI
    bench.list[[2]]<-adjusted.rand.index
    bench.list[[3]]<-nclust
    bench.list[[4]]<-cluster.idx
    names(bench.list)<-c('SI','ARI','G','Cluster Index')
    return(bench.list)
  }else{
    return('Error:Double check arguments')
  }
}

#DB Scan benchmark test(uses knn to pick epsilon, MinPts is set to 10) 
db.bench=function(X,Y,k,Pts){
  if(is.list(X)&is.list(Y)){
  SI<-vector(length=length(X))
  adjusted.rand.index<-vector(length=length(X))
  G<-vector(length=length(X))
  cluster.idx<-list()
  for(i in 1:length(X)){
    d<-dist(X[[i]])
    knn<-sort(d[1:(nrow(X[[i]])-1)])
    db.out<-dbscan(X[[i]],eps=knn[k],MinPts = Pts)
    if(max(db.out$cluster)==1|max(db.out$cluster)==0){
      SI[i]<-0
      adjusted.rand.index[i]<-0
      G[i]<-1
      cluster.idx[[i]]<-db.out$cluster
    }else{
    cluster.idx[[i]]<-db.out$cluster
    index<-db.out$cluster!=0
    temp.data<-X[[i]][index,]
    temp.cluster<-db.out$cluster[index]
    temp.sval<-silhouette(temp.cluster,dist=dist(temp.data))
    SI[[i]]<-mean(temp.sval[,3])
    adjusted.rand.index[i]<-ARI(temp.cluster[1:length(Y[[i]])],Y[[i]])
    G[i]<-max(db.out$cluster)
    }
  }
  names(cluster.idx)<-as.character(1:length(X))
  bench.list<-list()
  bench.list[[1]]<-SI
  bench.list[[2]]<-adjusted.rand.index
  bench.list[[3]]<-G
  bench.list[[4]]<-cluster.idx
  names(bench.list)<-c('SI','ARI','G','Cluster Index')
  return(bench.list)
  }else{
  return('Error: Double check arguments')
  }
}

#Contaiminated Normal Simulation Function
sim.func<-function(n,p,g,mu,sigma1,sigma2,alpha){
  if(alpha>=0&alpha<=0.5&is.matrix(mu)&ncol(mu)==p&nrow(mu)==g){
    n.data<-ceiling(alpha*n)+floor((1-alpha)*n)
    temp.mat<-matrix(0,nrow=5*n.data,ncol=p)
    temp.list<-list()
    for(i in 1:g){
      temp.mat.1<-rmvnorm(floor((1-alpha)*n),mu[i,],sigma1)
      temp.mat.2<-rmvnorm(ceiling(alpha*n),mu[i,],sigma2)
      temp.mat.3<-rbind(temp.mat.1,temp.mat.2)
      temp.list[[i]]<-temp.mat.3
    }
    simulation<-matrix(ncol=p)
    for(i in 1:g){
      simulation<-rbind(simulation,temp.list[[i]]) 
    }
    return(simulation[-1,])
  }else{
    return('Error: Double check arguments')
  }
}

#UPSP MNIST Handwritten Digits Tests 
#---------------------------------------------------------------------------------------------------------------------------------------

#Data set up (change number of pca tests here)
load('pcaMNIST.rdata')
dimensions<-c(24,31,40,53,77,130,276)
pcaMNIST.y<-vector()
for(i in 1:10){
  pcaMNIST.y<-c(pcaMNIST.y,rep(i,100))
}
pcaMNIST.y.list<-list()
for(i in 1:length(dimensions)){
  pcaMNIST.y.list[[i]]<-pcaMNIST.y
}

#Tests
pcaMNIST.kmed.si<-kmed.si.bench(pcaMNIST,2,10)
pcaMNIST.kmed.ari<-kmed.ari.bench(pcaMNIST,pcaMNIST.y.list,2,10)
pcaMNIST.t<-t.bench(pcaMNIST[1:5],pcaMNIST.y.list[1:5],2,10)
pcaMNIST.db<-db.bench(pcaMNIST,pcaMNIST.y.list,10,10)

#Plots
plot(dimensions,pcaMNIST.kmed.si$SI,ylim=c(0,max(pcaMNIST.db$SI)+0.1),ylab='Silhouette Index', main='Silhouette Index Comparison on MNIST Handwritten Digits')
lines(dimensions,pcaMNIST.kmed.si$SI,col='red')
points(dimensions[1:5],pcaMNIST.t$SI)
lines(dimensions[1:5],pcaMNIST.t$SI,col='blue')
points(dimensions,pcaMNIST.db$SI)
lines(dimensions,pcaMNIST.db$SI,col='green')
legend(180,0.75,legend=c('K medoids(G=2)','Student T model(G=5)','DBSCAN(G=2)'),col=c('red','blue','green'),lty=1)

plot(dimensions,pcaMNIST.kmed.ari$ARI,ylim=c(0,max(pcaMNIST.kmed.ari$ARI)+0.1),ylab='ARI',main='ARI Comparison on MNIST Handwritten Digits')
lines(dimensions,pcaMNIST.kmed.ari$ARI,col='red')
points(dimensions[1:5],pcaMNIST.t$ARI)
lines(dimensions[1:5],pcaMNIST.t$ARI,col='blue')
legend(180,0.4,legend=c('K medoids(G=10)','Student T model(G=5)','DBSCAN(G=NA)'),col=c('red','blue','green'),lty=1)


#Simulations outliers
#---------------------------------------------------------------------------------------------------------------------------------------
#Set up
sim.out.G<-5
sim.out.p<-10
sim.out.mu<-matrix(0,nrow=sim.out.G,ncol=sim.out.p)
for(i in 1:sim.out.G){
  for(j in 1:sim.out.p){
    sim.out.mu[i,j]<-runif(1,min=100*(i-1),max=100*i)
  }
}
sim.out.sig.1<-diag(100,sim.out.p,sim.out.p)
sim.out.sig.2<-diag(1000,sim.out.p,sim.out.p)

sim.out.data<-list()
sim.out.alpha<-seq(0.05,0.3,by=0.05)
for(i in 1:length(sim.out.alpha)){
  sim.out.data[[i]]<-sim.func(100,sim.out.p,sim.out.G,sim.out.mu,sim.out.sig.1,sim.out.sig.2,sim.out.alpha[i])
}
sim.out.y<-vector()
for(i in 1:sim.out.G){
  sim.out.y<-c(sim.out.y,rep(i,100))
}
sim.out.y.list<-list()
for(i in 1:7){
  sim.out.y.list[[i]]<-sim.out.y
}

#Tests
sim.out.kmed.si<-kmed.si.bench(sim.out.data,2,7)
sim.out.kmed.ari<-kmed.ari.bench(sim.out.data,sim.out.y.list,2,7)
sim.out.t<-t.bench(sim.out.data,sim.out.y.list,2,7)
sim.out.db<-db.bench(sim.out.data,sim.out.y.list,10,10)

#Plots
plot(sim.out.alpha,sim.out.kmed.si$SI,ylim=c(0,max(sim.out.kmed.si$SI)+0.1),xlab='Alpha',ylab='SI', main='Silhouette Index Comparison on Outlier Simulation')
lines(sim.out.alpha,sim.out.kmed.si$SI,col='red')
points(sim.out.alpha,sim.out.t$SI)
lines(sim.out.alpha,sim.out.t$SI,col='blue')
points(sim.out.alpha,sim.out.db$SI)
lines(sim.out.alpha,sim.out.db$SI,col='green')
legend(0.05,0.2,legend=c('K medoids(G=5)','Student T model(G=5)','DBSCAN(G=5)'),col=c('red','blue','green'),lty=1)

plot(sim.out.alpha,sim.out.kmed.ari$ARI,ylim=c(0,max(sim.out.kmed.ari$ARI)+0.1),xlab='Alpha',ylab='ARI',main='ARI Comparison on Outlier Simulation')
lines(sim.out.alpha,sim.out.kmed.ari$ARI,col='red')
points(sim.out.alpha,sim.out.t$ARI)
lines(sim.out.alpha,sim.out.t$ARI,col='blue')
points(sim.out.alpha,sim.out.db$ARI)
lines(sim.out.alpha,sim.out.db$ARI,col='green')
legend(0.05,0.3,legend=c('K medoids(G=5)','Student T model(G=5)','DBSCAN(G=5)'),col=c('red','blue','green'),lty=1)

#Simulation Number of Clusters
#------------------------------------------------------------------------------------------------------------------------------------
#Set up
sim.nclust.n<-50
sim.nclust.p<-5
sim.nclust.alpha<-0.05

sim.nclust.sig.1<-diag(100,sim.nclust.p,sim.nclust.p)
sim.nclust.sig.2<-diag(1000,sim.nclust.p,sim.nclust.p)

sim.nclust.data<-list()
sim.nclust.y.list<-list()
sim.nclust.G<-seq(5,30,by=5)
for(i in 1:length(sim.nclust.G)){
  sim.nclust.mu<-matrix(0,nrow=sim.nclust.G[i],ncol=sim.nclust.p)
  for(j in 1:sim.nclust.G[i]){
      sim.nclust.mu[j,]<-runif(sim.nclust.p,min=100*(i-1),max=100*i)
    }
  sim.nclust.data[[i]]<-sim.func(sim.nclust.n,sim.nclust.p,sim.nclust.G[i],sim.nclust.mu,sim.nclust.sig.1,sim.nclust.sig.2,sim.nclust.alpha)
  temp.vec<-vector()
  for(j in 1:sim.nclust.G[i]){
    temp.vec<-c(temp.vec,rep(j,sim.nclust.n))
  }
  sim.nclust.y.list[[i]]<-temp.vec
}

#Tests
sim.nclust.kmed.si<-kmed.si.bench(sim.nclust.data,2,30)
sim.nclust.kmed.ari<-kmed.ari.bench(sim.nclust.data,sim.nclust.y.list,5,30)
sim.nclust.t<-t.bench(sim.nclust.data,sim.nclust.y.list,2,30)
sim.nclust.db<-db.bench(sim.nclust.data,sim.nclust.y.list,10,10)

#Plots
plot(sim.nclust.G,sim.nclust.kmed.si$SI,ylim=c(0,max(sim.nclust.kmed.si$SI)+0.1),ylab='SI',xlab='Number of Clusters',main='Silhouette Index Comparison on Number of Clusters Simulation')
lines(sim.nclust.G,sim.nclust.kmed.si$SI,col='red')
points(sim.nclust.G,sim.nclust.t$SI)
lines(sim.nclust.G,sim.nclust.t$SI,col='blue')
points(sim.nclust.G,sim.nclust.db$SI)
lines(sim.nclust.G,sim.nclust.db$SI,col='green')
legend(5,0.15,legend=c('K medoids','Student T model','DBSCAN'),col=c('red','blue','green'),lty=1)

plot(sim.nclust.G,sim.nclust.kmed.ari$ARI,ylim=c(0,max(sim.nclust.kmed.ari$ARI)+0.1),ylab='ARI',xlab='Number of Clusters',main='ARI Comparison on Number of Clusters Simulation')
lines(sim.nclust.G,sim.nclust.kmed.ari$ARI,col='red')
points(sim.nclust.G,sim.nclust.t$ARI)
lines(sim.nclust.G,sim.nclust.t$ARI,col='blue')
points(sim.nclust.G,sim.nclust.db$ARI)
lines(sim.nclust.G,sim.nclust.db$ARI,col='green')
legend(8,0.23,legend=c('K medoids','Student T model','DBSCAN'),col=c('red','blue','green'),lty=1)

#Simulations Correlations
#------------------------------------------------------------------------------------------------------------------------------------
#Notes: Using var=1 to create correlation matrix. 2 sigmas is redundant here, but didn't want to write new function.
sim.cor.G<-5
sim.cor.p<-10
sim.cor.mu<-matrix(0,nrow=sim.cor.G,ncol=sim.cor.p)
for(i in 1:sim.cor.G){
  sim.cor.mu[i,]<-runif(1,min=100*(i-1),max=100*i)
}

sim.cor.vec<-seq(0,1,by=0.2)
sim.cor.sig<-list()
for(i in 1:length(sim.cor.vec)){
  sim.cor.sig[[i]]<-matrix(sim.cor.vec[i],nrow=sim.cor.p,ncol=sim.cor.p)
  for(j in 1:sim.cor.p){
    sim.cor.sig[[i]][j,j]<-1
  }
}

sim.cor.data<-list()
for(i in 1:length(sim.cor.sig)){
  sim.cor.data[[i]]<-sim.func(100,sim.cor.p,sim.cor.G,sim.cor.mu,sim.cor.sig[[i]],sim.cor.sig[[i]],0.5)
}
sim.cor.y<-vector()
for(i in 1:sim.cor.G){
  sim.cor.y<-c(sim.cor.y,rep(i,100))
}
sim.cor.y.list<-list()
for(i in 1:7){
  sim.cor.y.list[[i]]<-sim.cor.y
}

#Tests
sim.cor.kmed.si<-kmed.si.bench(sim.cor.data,2,8)
sim.cor.kmed.ari<-kmed.ari.bench(sim.cor.data,sim.cor.y.list,2,8)
sim.cor.t<-t.bench(sim.cor.data[-6],sim.cor.y.list[-6],2,8)
sim.cor.db<-db.bench(sim.cor.data,sim.cor.y.list,10,10)

#Plots
plot(sim.cor.vec,sim.cor.kmed.si$SI,ylim=c(0.5,max(sim.cor.kmed.si$SI)+0.1),ylab='SI',xlab='Correlation',main='Silhouette Index Comparison on Correlation Simulation')
lines(sim.cor.vec,sim.cor.kmed.si$SI,col='red')
points(sim.cor.vec[-6],sim.cor.t$SI)
lines(sim.cor.vec[-6],sim.cor.t$SI,col='blue')
points(sim.cor.vec,sim.cor.db$SI)
lines(sim.cor.vec,sim.cor.db$SI,col='green')
legend(0,0.7,legend=c('K medoids(G=5)','Student T model(G=5)','DBSCAN(G=5)'),col=c('red','blue','green'),lty=1)

plot(sim.cor.vec,sim.cor.kmed.ari$ARI,ylim=c(0,max(sim.cor.kmed.ari$ARI)+0.1),ylab='ARI',xlab='Correlation',main='ARI Comparison on Correlation Simulation')
lines(sim.cor.vec,sim.cor.kmed.ari$ARI,col='red')
points(sim.cor.vec[-6],sim.cor.t$ARI)
lines(sim.cor.vec[-6],sim.cor.t$ARI,col='blue')
points(sim.cor.vec,sim.cor.db$ARI)
lines(sim.cor.vec,sim.cor.db$ARI,col='green')
legend(0,0.3,legend=c('K medoids(G=5)','Student T model(G=5)','DBSCAN(G=5)'),col=c('red','blue','green'),lty=1)

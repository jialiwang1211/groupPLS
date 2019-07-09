m=100
MSE<-rep(NA,m)

for(i in 1:m){
  ##setup parameters###
  n=100 ##number of observations
  p=50  ##number of X variables in each cluster
  np=10 ##number of cluster
  pc=2  ##number of latent comp
  set.seed(i+100)
  ##no. of clusters##
  nop=rep(p,np)
  nopc=rep(pc,np)
  eigenval=NULL
  Xclu=NULL
  Zclu=NULL
  
  #set.seed(1)
  # corM=NULL
  # for(c in 1:length(nopc)){
  #   corM[[c]]<-matrix(0,nopc[c],nopc[c])
  #   for(i in 1:nopc[c]) {
  #     for(j in 1:nopc[c]) {
  #        corM[[c]][i,j]=0.1^(abs((i-j)))
  #       #corM[[c]][i,j]<-0.8^(1-I(i==j))
  #     }}
  # 
  gpindex=rep(1:np,each=pc)
  ##c is no. of clusters##
  for(c in 1:length(nopc)){
    eigenval[[c]]=diag(c(abs(rnorm(nopc[c],10,1))))
    #eigenval[[c]]=abs(rnorm(nopc[c],10,1))
    eigenvec=randortho(nop[c],type ="orthonormal")
    Z=rbind(mvrnorm(n/2, rep(0,nopc[c]), diag(nopc[c]),empirical =TRUE),mvrnorm(n/2, rep(0,nopc[c]), diag(nopc[c]),empirical =TRUE))%*%sqrt(eigenval[[c]])
    noise=matrix(rnorm(prod(n,nop[c]),0,0.2), nrow = n)
    
    Xnew=Z%*%t(eigenvec[,1:nopc[c]])+noise
    
    Xclu=cbind(Xclu,Xnew)
    Zclu=cbind(Zclu,Z)
  }
  
  ##number of responses
  q=3
  
  
  beta=cbind(rep(c(0,2,3,0,3,0,0,0,rep(0,12)),np/10),rep(c(3,0,0,3,0,0,3,0,rep(0,12)),np/10),rep(c(0,2,0,0,0,3,0,3,rep(0,12)),np/10))
  
  
  #beta<-mvrnorm(sum(nopc), rep(0,q), 4*diag(q))
  #beta=normalize.vector(beta)
  #Zclu=scale(Zclu,scale=TRUE, center=TRUE)
  
  y=Zclu%*%beta+mvrnorm(n,rep(0,q),1*diag(q))
xtrain=Xclu[1:(n/2),]
xtest=Xclu[(n/2+1):n,]
ytrain=as.matrix(y[1:(n/2),])
ytest=as.matrix(y[(n/2+1):n,])

xmean=apply(xtrain,2,mean)
ymean=apply(ytrain,2,mean)
xctrain=sweep(xtrain,2,xmean, "-")
yctrain=sweep(ytrain,2,ymean, "-")
xctest=sweep(xtest,2,xmean, "-")
yctest=sweep(ytest,2,ymean, "-")

trueindex=rep(1:np,each=p)
ncomp=30

rs<-groupPLS(yctrain~xctrain,scale.Y=FALSE,scale.X=FALSE,group.index=trueindex,k=FALSE,
             seed=100+i,ncomp=ncomp)

pred<-groupPLS.predict(rs,newdata =xctest)
MSE[i]<-sqrt(sum((pred-yctest)^2)/length(yctest))
}

vip<-vipscore(groupPLS.obj=rs,y=yctrain)

boxplot(MSE)
mean(MSE)

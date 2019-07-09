#' Two-step group penalized partial least squares regression
#'
#' Performs group penalized partial least squares regression.
#'
#' More details.
#'
#' @param formula model formula, \code{y~x}, where \code{y} are the response variables and x are the predictors.
#'
#' @param scale.X scale \code{x} by the sample standard deviation of each variable. TRUE or FALSE
#'
#' @param scale.Y scale \code{y} by the sample standard deviation of each variable. TRUE or FALSE
#'
#' @param k number of groups. If \code{group.index} is not provided, then perform the hierarchical clustering based on a similarity matrix to define \code{k} clusters.
#'
#' @param group.index pre-defined group index for each predictor.
#'
#' @param group.pen group penalization on latent components. TRUE or FALSE.
#'
#' @param lambda.min lower bound of group penalty paramter.
#'
#' @param lambda.max upper bound of group penalty paramter.
#'
#' @param ncomp number of latent components in partial least squares.
#'
#' @param seed an integer.
#'
#' @return Return \code{groupPLS} object.
#'
#' \item{ncomp.group}{indices of latent components from PLS in each group.}
#'
#' \item{latent.group}{combined latent components from PLS in each group.}
#'
#' \item{coef}{regression coefficients of latent components on responses.}
#'
#' \item{scale.X}{if scaling of X is requested.}
#'
#' \item{scale.Y}{if scaling of Y is requested.}
#'
#' \item{Xmean}{mean of each variable in X.}
#'
#' \item{Xsd}{standard deviation of each variable in X.}
#'
#' \item{Ymean}{mean of each variable in Y.}
#'
#' \item{Ysd}{standard deviation of each variable in Y.}
#'
#' \item{k}{number of groups.}
#'
#' \item{group.index}{group index for each predictor.}
#'
#' \item{group.lambda}{group penalty paramter selected by cross validation.}
#'
#' \item{call}{function call.}
#'
#' @examples data<-read.csv("Spectra_eps.csv")
#'data<-data[complete.cases(data),]
#'rs<-groupPLS(data[,c(2,4:6)] ~ data[,7:dim(data)[2]], k=3)
#'
#' @author Jiali Wang (\email{jiali.wang@@data61.csiro.au}); Le Chang (\email{le.chang@@anu.edu.au})
#'
#' @export
#' @import pls
#' @import MSGLasso
groupPLS<-function(formula,scale.X=FALSE,scale.Y=TRUE,k,group.index=NA,group.pen=TRUE,lambda.min=0.01,lambda.max=0.1,ncomp=20,seed=1){

set.seed(seed)

x<-eval(formula[[3]])
y<-eval(formula[[2]])

xmean<-apply(x,2,mean)
xsd<-apply(x,2,sd)
ymean<-apply(y,2,mean)
ysd<-apply(y,2,sd)

xc<-as.matrix(sweep(x,2,xmean, "-"))
if(scale.X==TRUE){
xc<-sweep(xc,2,xsd, "/")
}

yc<-as.matrix(sweep(y,2,ymean, "-"))
if(scale.Y==TRUE){
yc<-sweep(yc,2,ysd, "/")
}

dismat<-as.dist(1-abs(cor(xc)))
iniclu<-hclust(dismat, method = "average")

if(k){
group.index<-cutree(iniclu,k=k)
}

q<-dim(y)[2]

compall=NULL
compno=NULL
W<-NULL
for (k in unique(group.index)) {
  ink=which(group.index==k)
  plsink=plsr(yc~xc[,ink], ncomp=min(ncomp/length(unique(group.index)),length(ink)), method="simpls", validation = "LOO")
  scoreink<-apply(plsink$validation$PRESS/plsink$validation$PRESS0,2,mean)
  nopk<-(1:length(scoreink))[scoreink==min(scoreink)]

  plsinf=plsr(yc~xc[,ink], nopk, method="simpls")
  compall=cbind(compall,xc[,ink]%*%plsinf$projection)
  compno=c(compno,rep(k,nopk))
  W[[k]]<-plsinf$projection
}

#W<-bdiag(W)

lamG.v<-NA
if(group.pen==FALSE){
  coef<-coef(lm(yc~compall-1))
}else{

##
r<-dim(yc)[2]
PQ.grps<-as.matrix(cbind(c(compno-1,0:(r-1)),999))
GR.grps<-matrix(999,nrow=length(unique(compno))+r,ncol=max(table(compno),r)+1)
for (h in 1:length(unique(compno))){
  GR.grps[h,1:length(which(compno==h))]=which(compno==h)-1
}
GR.grps[(length(unique(compno))+1):(length(unique(compno))+r),1]<-0:(r-1)

PenL<-matrix(0,length(compno),r)
PenG<-matrix(1,length(unique(compno)),r)
grp_Norm0 <- matrix(0, length(unique(compno)),r)
grpWTs=matrix(1,length(unique(compno)), r)

lamG.v <- seq(lambda.min,lambda.max,length=20)
try.cv<-MSGLasso.cv(compall,yc,grpWTs,PenL,PenG,PQ.grps,GR.grps,0,lamG.v, fold=5,seed=seed)

MSGLassolamG<-try.cv$lams.c[which.min(as.vector(try.cv$rss.cv))][[1]]$lam3
MSGLassolamG.m<-matrix(MSGLassolamG,length(unique(compno)),r)

mglassofitcom<-MSGLasso(compall,yc,grpWTs,PenL,PenG,PQ.grps,GR.grps, grp_Norm0,0,MSGLassolamG.m,Beta0=NULL)

coef<-mglassofitcom$Beta
}
return(list(ncomp.group=compno,latent.group=compall,coef=coef,W=W,scale.X=scale.X,scale.Y=scale.Y,Xmean=xmean,Xsd=xsd,Ymean=ymean,Ysd=ysd,k=k,group.index=group.index,group.lambda=lamG.v,call=formula))
}


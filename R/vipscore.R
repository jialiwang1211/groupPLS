#' Variable importance in projection
#'
#' Compute variable importance in projection scores (VIP) for response variables.
#'
#' @param groupPLS.obj the object from fitting \code{groupPLS}.
#'
#' @param y data frame contains response variables.
#'
#' @return Variable importance in projection scores for response variables.
#'
#' @examples vip<-vipscore(groupPLS.obj=rs,y=ytrain)
#' par(mfrow=c(3,1))
#' for (i in 1:3){
#' plot(vip[,i],type="h")
#' }
#'
#' @export
#' @import Matrix
vipscore<-function(groupPLS.obj,y){
  p<-length(groupPLS.obj$Xmean)
  q<-length(groupPLS.obj$Ymean)
  nocomp<-length(groupPLS.obj$ncomp.group)
  comp<-groupPLS.obj$latent.group
  W<-bdiag(groupPLS.obj$W)
  Q<-groupPLS.obj$coef

  y<-as.matrix(sweep(y,2,groupPLS.obj$Ymean, "-"))
  if(groupPLS.obj$scale.Y==TRUE){
    y<-sweep(y,2,groupPLS.obj$Ysd, "/")
  }

  Rsq<-matrix(NA,nocomp,q)
  vip<-matrix(0,p,q)

  for (j in 1:q){
    select<-Q[,j]!=0
    for (i in 1:nocomp){
      if(i==1){
        Rsq[i,j]<-summary(lm(y~comp[,1]-1))[[j]]$r.squared
      }else{
        Rsq[i,j]<-summary(lm(y~comp[,(1:i)[select[1:i]]]-1))[[j]]$r.squared
      }
    }
  }

  for (j in 1:q){
    rsq<-c(Rsq[1,j],diff(Rsq[,j]))
    rsq<-rsq*(Q[,j]!=0)
    w_normsq<-(apply(W,2,function(x) x/sqrt(sum(x^2))))^2
    vip[,j]<-sqrt(w_normsq%*%rsq/sum(rsq)*p)
  }
  return(vip)
}

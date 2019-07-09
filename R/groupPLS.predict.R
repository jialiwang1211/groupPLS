#' Predict for new observations
#'
#' Prediction for new observations from the fit of two-step group penalized partial least squares regression.
#'
#' @param groupPLS.obj the object from fitting \code{groupPLS}.
#'
#' @param newdata predictors of new observations.
#'
#' @return Matrix of predicted values.
#'
#'
#' @export
groupPLS.predict<-function(groupPLS.obj,newdata=newdata){

  x<-as.matrix(newdata)

  xmean<-groupPLS.obj$Xmean
  xsd<-groupPLS.obj$Xsd
  ymean<-groupPLS.obj$Ymean
  ysd<-groupPLS.obj$Ysd

  xc<-as.matrix(sweep(x,2,xmean, "-"))
  if(groupPLS.obj$scale.X==TRUE){
    xc<-sweep(xc,2,xsd, "/")
  }

  ##
  compalltest=NULL
  for (k in unique(groupPLS.obj$group.index)) {
    ink<-which(groupPLS.obj$group.index==k)
    compalltest<-cbind(compalltest,xc[,ink]%*%groupPLS.obj$W[[k]])
  }
  ##

  yhat=compalltest%*%groupPLS.obj$coef

  if(groupPLS.obj$scale.Y==TRUE){
    yhat<-sweep(yhat,2,ysd, "*")
  }
  yhat<-sweep(yhat,2,ymean, "+")

  return(pred=yhat)
}

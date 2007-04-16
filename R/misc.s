scoref<-function(r,tau){
  tt<-signr<-sign(r)
  tt[signr>0]<-tau
  tt[signr<0]<-tau-1
  return(tt)
}

plotDim<-function(p){
  temp<-round(sqrt(p))
  temp1<-(temp+1)
  ifelse(temp*temp1<p, dims<-c(temp1,temp1), dims<-c(temp,temp1))
  if(temp^2==p || p<temp*temp){dims<-c(temp,temp)}
  return(list(dims=dims))
}


#meansd<-function(obj){
#  means <- apply(obj$theta, 2, mean)
#  sds <- sqrt(apply(obj$theta, 2, var))
#  res<-matrix(cbind(means, sds), ncol=2, byrow=FALSE)
#  dimnames(res)<-list(NULL, c("Mean","SD"))
#  return(res)
#}

rqmcmb.ci<-function(obj, alpha=0.1){
  sds <- sqrt(apply(obj$theta, 2, var))
  ciL<-obj$coef+qnorm(alpha/2)* sds
  ciU<-obj$coef-qnorm(alpha/2)* sds
  res<-matrix(cbind(obj$coef, sds, ciL,ciU ), ncol=4, byrow=FALSE)
  dimnames(res)<-list(NULL, c("Coef","SD", "L", "U"))
  return(res)
}


rqmcmb.plot<-function(obj, alpha=.10){
  p <- dim(obj$theta)[2]
  K <- length(obj$theta[,1])
  ci <- rqmcmb.ci(obj)
  lci <- ci[,3]
  uci <- ci[,4]
  #theta <- obj$theta
  thetaOrig <- obj$theta[1,]
  z<-qnorm(alpha/2)
  dev<-dev.cur()
  if(dev==1) {motif()}
  par(mfrow=plotDim(p)$dims)
  plotr <- apply(obj$theta,2,range)
  plotlim <- rbind(pmin(plotr[1,],lci), pmax(plotr[2,],uci))
  for (i in 1:p){
    plot(obj$theta[,i], ylab=paste("theta[",i,"]"),
      ylim=plotlim[,i])
    lines(seq(0, K+1,,10), rep(thetaOrig[i],10),type="b",pch=1,lty=2)
    lines(c(0, K+1), rep(lci[i],2))
    lines(c(0, K+1), rep(uci[i],2))
    title(paste("theta[",i-1,"]"))
  }
  par(mfrow=c(1,1))

}

.First.lib <-
function(lib, pkg) {
   library.dynam("rqmcmb2", pkg, lib)
   print("rqmcmb2 library loaded")}

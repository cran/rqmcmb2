rqmcmb <- function(x=x, y=y, tau=0.5, K=100, int=TRUE, plotTheta=FALSE) {

if(exists(".Random.seed") == FALSE){ 
  rnorm(1)
  seed <- .Random.seed[2]  
}
else{
  seed <- .Random.seed[2]
}

n<-length(y)
x<-as.matrix(x)
if (nrow(x)!=n) stop("Warning: Dimensions of x and y do not match") 

if(int==TRUE) x <- cbind(rep(1,n), x)

if (n>=200000){
  print("Sample size is too large. Maximum allowed is 200000")}
else {
  if (ncol(x)>=100){
    print("Number of parameters is too large. Maximum allowed is 100")}


#main code


else { 

  if ( n <= 5000) {fit<-rq.fit(x,y, tau=tau, ci=FALSE)} else {
    if ( n <=50000) {fit<-rq.fit.fn(x,y, tau=tau) } else {
      fit<-rq.fit.pfn(x,y, tau=tau) } }


  thetaOrig<-fit$coef
  p<-length(thetaOrig)
  thetaTilda<-thetaOrig


  if ((n*tau < 5*p+1) || (n*(1-tau) < 5*p+1)) {
    print("Warning: May not have enough data to estimate the")
    print("requested quantile reliably")}

  Z <- t(x)%*%x
  cov.svd<-svd(Z)
  A<-cov.svd$u%*%(diag(sqrt(cov.svd$d)^-1))%*%t(cov.svd$v)
  Ainv<-cov.svd$u%*%(diag(sqrt(cov.svd$d)))%*%t(cov.svd$v)
  thetaTilda<-matrix(Ainv%*%thetaOrig, nrow=length(thetaOrig))

  cn<-sqrt(max(cov.svd$d)/min(cov.svd$d))
  if (cn>100){
    print("Warning: Nearly singular design detected;")
    print("the results from rqmcmb may be unreliable")
  }

  psi<-scoref(fit$resid,tau)
  psimat<-matrix(psi,nrow=n,ncol=p,byrow=FALSE)
  ZTilda<-(x%*%A)*psimat

  # re-define x on the theta_tilda scale
  x<-x%*%A
  sumxij<-apply(x,2,sum)
  sumabsxij<-apply(abs(x),2,sum)
  zstar<-.C("rqmcmb",
        array(as.double(t(x))),
        array(as.double(y)),
        array(as.double(tau)),
        array(as.double(thetaTilda)),
        array(as.double(t(A))),
	array(as.double(ZTilda)),
        array(as.double(sumxij)),
        array(as.double(sumabsxij)),
        as.integer(n),
        as.integer(p),
	success=as.integer(1),
        theta=array(as.double(rep(0,K*p+p)),c(p,K+1)),
        as.integer(K),
        as.integer(seed),
        PACKAGE="rqmcmb2"
        )


  if(zstar$success==0){return(list(success=0))}

  thetaTildaMCMB<-t(zstar$theta)
  thetaOrigMCMB<-t(A%*%zstar$theta)

  if(plotTheta==TRUE){
    dev<-dev.cur()
    if(dev==1) {x11()}
    par(mfrow=plotDim(p)$dims)

  for (i in 1:p){
    plot(thetaOrigMCMB[,i], ylab=paste("theta[",i-1,"]"))
    lines(c(0, K+1), rep(thetaOrig[i],2))
    title(paste("theta[",i-1,"]"))
  }

  par(mfrow=c(1,1))
}

return(list(coef=thetaOrig, cov=var(thetaOrigMCMB), theta=thetaOrigMCMB, 
  success=zstar$success, cn=cn))
}

}
}


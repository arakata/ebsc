library(SemiPar)
library(splines)
library(forecast)
library(nlme)

bspline <- function(x, xl, xr, ndx, bdeg,L){
    dx<-(xr-xl)/ndx
    knots<-seq(xl-bdeg*dx,xr+bdeg*dx,by=dx)
    #B <- spline.des(knots, x, bdeg + 1, 0 * x,outer.ok=TRUE)$design
    B <- spline.des(knots, x, bdeg + 1, 0 * x,outer.ok=TRUE)$design
    solve(L)%*%B
}

mycorMatrix <-function(noise){
    noise.model<-auto.arima(noise)
    mycorcoefs<-noise.model$coef
    DATA <- cbind(index=seq(1,length.out=length(noise)),noise=as.numeric(noise))
    myDATA <- data.frame(DATA)
    index <- myDATA$index
    string <- names(mycorcoefs)
    substring <- substring(string,1,2)
    ar.ord <- length(substring[substring=="ar"])
    ma.ord <- length(substring[substring=="ma"])
    csARMA <- corARMA(mycorcoefs, form = ~index, p = ar.ord,q=ma.ord);
    csARMA <- Initialize(csARMA, data = myDATA);
    if((ar.ord==0)&(ma.ord==0)){
        print("white noise");
        R=diag(length(noise));
        sigma2=sd(noise)^2}else{
        R=corMatrix(csARMA)
        sigma2=noise.model$sigma2    
    }
    list(R=R,sigma2=sigma2);
}

fitting<-function(num.knots,bdeg,pord,L,y){
    n=length(y)
    x=1:n
    xl=1
    xr=n
    ndx=num.knots

    B <- bspline(x, xl, xr, ndx, bdeg,L)
    D <- diag(ncol(B))
    for (k in 1:pord) D <- diff(D)

    log10lambdas=seq(-10,-0.1,length.out=1000)
    m=length(log10lambdas)
    gcv=rep(NA,m)
    edf=rep(NA,m)
    count=0
    for(log10lambda in log10lambdas){
       count=count+1
       S <-B%*%solve(t(B) %*% B + 10^log10lambda * t(D) %*% D)%*%t(B)
       edf[count]<-sum(diag(S))
       a <- solve(t(B) %*% B + 10^log10lambda * t(D) %*% D,t(B) %*% y)
       yhat<- B%*%a
       s <- sum((y - yhat)^2)
       Q <- solve(t(B) %*% B + 10^log10lambda * t(D) %*% D)
                                        # matrix inversion
       t <- sum(diag(Q %*% (t(B) %*% B)))
       gcv[count]<-s/(nrow(B)-t)^2
   }
    index=which.min(log(gcv))
    lambda.star=10^log10lambdas[index]
    gcv.star=gcv[index]
    a <- solve(t(B) %*% B + 10^lambda.star * t(D) %*% D,t(B) %*% y)
    yhat<-B%*%a
    error<-y-yhat
    
    list(yhat=yhat,edf=edf,gcv=gcv,index=index,error=error)
}

####################################
#step1: working independence
####################################
fit0<-fitting(num.knots=5,bdeg=2,pord=2,L=diag(length(y)),y=y)

####################################
#step2: covariance function
####################################    
error<-fit0$error
covnoise=mycorMatrix(error)
R.hat=covnoise$R
sigma2.hat=covnoise$sigma2
    
L<-t(chol(R))
ytilde<-solve(L)%*%y

####################################
#step3: update estimation
####################################    

fitu<-fitting(num.knots=5,bdeg=2,pord=2,L=L,y=ytilde)

error<-fitu$error
covnoise=mycorMatrix(error)
R.hat=covnoise$R
sigma2.hat=covnoise$sigma2
    
L<-t(chol(R))
ytilde<-solve(L)%*%y

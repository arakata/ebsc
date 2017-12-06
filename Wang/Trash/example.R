library(forecast)
library(nlme)
library(SemiPar)
library(splines)


bspline <- function(x, xl, xr, ndx, bdeg){
    dx<-(xr-xl)/ndx
    knots<-seq(xl-bdeg*dx,xr+bdeg*dx,by=dx)
    B <- spline.des(knots, x, bdeg + 1, 0 * x,outer.ok=TRUE)$design
    B
}

fit.noise <-function(noise,p,q){
    n<-length(noise)
    noise.model<-arima(noise,order=c(p,0,q),include.mean=FALSE)
    coefs<-noise.model$coef
    data <- data.frame(index=1:n,noise=as.numeric(noise))
    index <- data$index
    string <- names(coefs)
    substring <- substring(string,1,2)
    ar.ord <- length(substring[substring=="ar"])
    ma.ord <- length(substring[substring=="ma"])
    csARMA <- corARMA(coefs, form = ~index, p = ar.ord,q=ma.ord);
    csARMA <- Initialize(csARMA, data = data);
    if((ar.ord==0)&(ma.ord==0)){
        print("white noise");
        R=diag(length(noise));
        sigma2=sd(noise)^2}else{
        R=corMatrix(csARMA)
        sigma2=noise.model$sigma2    
    }
    list(R=R,sigma2=sigma2);
}

getAIC<-function(num.knots,bdeg,pord,Btilde,ytilde,p){
    D <- diag(ncol(Btilde))
    for (k in 1:pord) D <- diff(D)

    log10lambdas=seq(-10,-0.1,length.out=1000)
    m=length(log10lambdas)
    AICmat=rep(NA,m)

        count=0
        for(log10lambda in log10lambdas){
            count=count+1
            S <-Btilde%*%solve(t(Btilde) %*% Btilde + 10^log10lambda * t(D) %*% D)%*%t(Btilde)
            edf<-sum(diag(S))
            a <- solve(t(Btilde) %*% Btilde + 10^log10lambda * t(D) %*% D,t(Btilde) %*% ytilde)
            yhat<- Btilde%*%a
            s <- sum((ytilde - yhat)^2)
            Q <- solve(t(Btilde) %*% Btilde + 10^log10lambda * t(D) %*% D)
            t <- sum(diag(Q %*% (t(Btilde) %*% Btilde)))
            
            a <- solve(t(Btilde) %*% Btilde + 10^log10lambda * t(D) %*% D,t(Btilde) %*% ytilde)
            yhat<-Btilde%*%a
            error<-ytilde-yhat
       
            sigma2.hat<-sd(error)^2
            AICmat[count]<-n*log(sigma2.hat)+2*(p+edf)
        }
    index=which.min(AICmat)
    lambda.star=10^log10lambdas[index]
    AIC.star=AICmat[index]

    a <- solve(t(Btilde) %*% Btilde + 10^lambda.star * t(D) %*% D,t(Btilde) %*% ytilde)
    yhat<-Btilde%*%a
    error<-ytilde-yhat

    list(AIC.star=AIC.star,lambda.star=lambda.star,error=error)
}

fit.mean<-function(num.knots,bdeg,pord,Btilde,ytilde,lambda.star){
    D <- diag(ncol(Btilde))
    for (k in 1:pord) D <- diff(D)
    a <- solve(t(Btilde) %*% Btilde + 10^lambda.star * t(D) %*% D,t(Btilde) %*% ytilde)
    yhat<-Btilde%*%a
    error<-ytilde-yhat

    list(yhat=yhat,error=error)
}

############################################################

x=seq(0,1,length.out=200)
mu=(exp(1)^(10*x-5))/(1+exp(1)^(10*x-5))
noise=arima.sim(n=length(mu),list(ar=c(0.5)),sd=0.1)
y=mu+noise

num.knots=5
bdeg=2
pord=2

n=length(y)
B=bspline(x=1:n, xl=1, xr=n, ndx=num.knots, bdeg=bdeg)

ytilde=y



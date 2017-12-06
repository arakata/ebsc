library(SemiPar)
library(splines)
library(forecast)
library(nlme)

bspline <- function(x, xl, xr, ndx, bdeg){
    dx<-(xr-xl)/ndx
    knots<-seq(xl-bdeg*dx,xr+bdeg*dx,by=dx)
    #B <- spline.des(knots, x, bdeg + 1, 0 * x)$design
    B <- spline.des(knots, x, bdeg + 1, 0 * x,outer.ok=TRUE)$design
    B
}

fit.noise <-function(noise,ar.ord){
    noise.model<-arima(noise,order=c(ar.ord,0,0),include.mean=FALSE)
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

fit.mean<-function(num.knots,bdeg,pord,y,ar.ord){
    n=length(y)
    x=1:n
    xl=1
    xr=n
    ndx=num.knots
    if(L==NULL){L=diag(n)}
    
    B <- bspline(x, xl, xr, ndx, bdeg)
    D <- diag(ncol(B))
    for (k in 1:pord) D <- diff(D)

    log10lambdas=seq(-10,-0.1,length.out=1000)
    m=length(log10lambdas)

    for(p in 1:4){
        count=0
        for(log10lambda in log10lambdas){
            count=count+1
            S <-B%*%solve(t(B) %*% B + 10^log10lambda * t(D) %*% D)%*%t(B)
            #edf[count]<-sum(diag(S))
            a <- solve(t(B) %*% B + 10^log10lambda * t(D) %*% D,t(B) %*% y)
            yhat<- B%*%a
            s <- sum((y - yhat)^2)
            Q <- solve(t(B) %*% B + 10^log10lambda * t(D) %*% D)
            t <- sum(diag(Q %*% (t(B) %*% B)))
            #gcv[count]<-s/(nrow(B)-t)^2

            a <- solve(t(B) %*% B + 10^lambda * t(D) %*% D,t(B) %*% y)
            yhat<-B%*%a
            error<-y-yhat
       
            noise.model=fit.noise(noise=error,ar.ord=ar.ord)
            R.hat=noise.model$R
            sigma2.hat=noise.model$sigma2
    
            L<-t(chol(R.hat))
            ytilde<-solve(L)%*%y
            Btilde<-solve(L)%*%B
       
            AIC[count,p]<-n*log(sigma2.hat)+2(ar.ord+edf)
        }
    }
    index=which(AIC == min(AIC), arr.ind = TRUE)
    p.star=index[2]
    lambda.star=10^log10lambdas[index[1]]

    a <- solve(t(B) %*% B + 10^lambda.star * t(D) %*% D,t(B) %*% y)
    yhat<-B%*%a
    error<-y-yhat


    #noise

    #AIC
    
    
    list(yhat=yhat,edf=edf,gcv=gcv,index=index,error=error,L=L,ytilde=ytilde,Btilde=Btilde)
}

CochranOrcutt<-function(num.knots,bdeg,pord,y){

    #step 1: working independence
    
    fit<-fit.mean(num.knots=num.knots,bdeg=bdeg,pord=pord,y=y,ar.ord=ar.ord)

    #step2: covariance function
    error<-fit0$error
    covnoise=mycorMatrix(error)
    R.hat=covnoise$R
    sigma2.hat=covnoise$sigma2
    
    L<-t(chol(R))
    ytilde<-solve(L)%*%y
}


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

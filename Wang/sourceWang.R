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

Sigma.eps <-function(eps,p,q){
    index<-1:length(eps)
    data <- data.frame(index=index,eps=as.numeric(eps))

    #white noise
    if((p==0)&(q==0)){R=diag(n);sigma2=sd(eps)^2}else{
        
        #not white noise
        eps.model<-arima(eps,order=c(p,0,q),include.mean=FALSE)
        csARMA <- corARMA(eps.model$coef, form = ~index, p = p,q=q);
        csARMA <- Initialize(csARMA, data = data);
        R<-corMatrix(csARMA);
        sigma2<-eps.model$sigma2;
    }
    list(R=R,sigma2=sigma2);
}

arma.eps<-function(eps){
    eps.model<- auto.arima(eps,d=0)
    string <- names(eps.model$coef)
    substring <- substring(string,1,2)
    p <- length(substring[substring=="ar"])
    q <- length(substring[substring=="ma"])
    list(p=p,q=q)
}

search.opt<-function(B,y,p){
    pord<-2
    D <- diag(ncol(B))
    for (k in 1:pord) D <- diff(D)

    log10lambdas=seq(-10,-0.1,length.out=1000)
    m=length(log10lambdas)
    AICs=rep(NA,m)

        count=0
        for(log10lambda in log10lambdas){
            count<-count+1
            S <- B%*%solve(t(B) %*% B + 10^log10lambda * t(D) %*% D)%*%t(B)
            edf.c<- sum(diag(S))            
            coefs <- solve(t(B) %*% B + 10^log10lambda * t(D) %*% D,t(B) %*% y)       
            sigma2.c<-sd(y-B%*%coefs)^2
            AICs[count]<-n*log(sigma2.c)+2*(p+edf.c)
        }
    
    lambda.star=10^log10lambdas[which.min(AICs)]
    AIC.star=AICs[which.min(AICs)]

    list(AIC.star=AIC.star,lambda.star=lambda.star)
}

stable.opt<-function(B,y,p){
    n<-length(y)
    R.storage <- array(data=NA, dim = c(5,n,n))
    R.storage[1,,]<-diag(n)
    yhat<-spline.mu(B=B,y=y,lambda=search.opt(B=B,y=y,p=0)$lambda.star)$yhat
    error<-y-yhat

    for(i in 2:5){
        R.storage[i,,]<-Sigma.eps(eps=error,p=p,q=0)$R
        L<-t(chol(R.storage[i,,]))
        ytilde<-solve(L)%*%y
        Btilde<-solve(L)%*%B
        param<-search.opt(B=Btilde,y=ytilde,p=p)
        mu<-spline.mu(Btilde,ytilde,param$lambda.star)$yhat
        l2norm<-norm(matrix(eigen(R.storage[i,,])$values-eigen(R.storage[i-1,,])$values,1),type="2")
        if(l2norm<10){break}else{error<-ytilde-mu}
    }
    list(Btilde=Btilde,ytilde=ytilde)
}

spline.mu<-function(B,y,lambda){
    pord<-2
    D <- diag(ncol(B))
    for (k in 1:pord) D <- diff(D)
    coefs <- solve(t(B) %*% B + lambda* t(D) %*% D,t(B) %*% y)
    yhat=B%*%coefs
    
    list(yhat=yhat,error=y-yhat)
}

CochraneOrcut<-function(num.knots,bdeg,y){    

    #working independence
    n<<-length(y)
    B<-bspline(x=1:n, xl=1, xr=n, ndx=num.knots, bdeg=bdeg)

    yhat<-spline.mu(B=B,y=y,lambda=search.opt(B=B,y=y,p=0)$lambda.star)$yhat
    error<-y-yhat
    print("ok")

    #try correlation models
    AICs<-rep(NA,4)
    lambdas<-rep(NA,4)
    for(p in 1:4){
        out<-stable.opt(B=B,y=y,p=p)
        param<-search.opt(B=out$Btilde,y=out$ytilde,p=p)
        AICs[p]<-param$AIC.star
        lambdas[p]<-param$lambda.star
    }
    p.star=which.min(AICs)
    lambda.star=lambdas[p.star]
    
    #re-estimate noise
    eps.model=Sigma.eps(eps=error,p=p.star,q=0)
    R.hat=eps.model$R
    L<-t(chol(R.hat))
    ytilde<-solve(L)%*%y
    Btilde<-solve(L)%*%B
    
    mu.star<-spline.mu(B=Btilde,y=ytilde,lambda=lambda.star)

    yhat<-L%*%mu.star$yhat
    
    list(yhat=yhat,Rhat=R.hat,p.star=p.star,lambda.star=lambda.star)
}





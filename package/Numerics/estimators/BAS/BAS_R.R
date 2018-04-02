###############################
## estimating equations and auxiliary functions

Elambda	= function(log10lambda,coefs,eigens,n,q)
{ # estimating equation for lambda
  lambda<-10^log10lambda
  r.mat1 <- 1/(1 + eigens*lambda)
  r.mat <- eigens*lambda/(1 + eigens*lambda)
  sum(coefs^2*r.mat*r.mat1) - (sum(r.mat1)-q)*(sum(coefs^2*r.mat)+1)/(n+1)
}

Eq=function(log10lambda,coefs,eigens,n,q,rho)
{#estimating equation for q
  lambda<-10^log10lambda
  r.mat1 <- 2*q*log((eigens/rho)^(1/(2*q)))/(1 + eigens*lambda)
  r.mat <- eigens*lambda/(1 + eigens*lambda)
  sum(coefs^2*r.mat*r.mat1) - sum(r.mat1[-c(1:q)])*(sum(coefs^2*r.mat)+1)/(n+1)
}

smooth_rhos = function(rhos,q=q)
{ #smoothing rhos
  n=length(rhos)
  x=seq(1,n,length=n)
  Basis = drbasis(nn=n,qq=q)
  PP = Basis$eigenvectorsQR
  ee = as.vector(n*Basis$eigenvalues)
  coefs = t(PP) %*% rhos
  temp = try(uniroot(Elambda,zero_range,coefs=coefs,eigens=ee,n=n,q=q,tol = eps,maxiter=max_iter,extendInt="yes")$root, silent = TRUE)
  if(class(temp)=="try-error"){temp = 5}
  rho.fit = PP %*% (coefs/(1+(10^temp)*ee)) 
  pmax( 10^-5, rho.fit)
}


############################
## generate data and set values

library(eBsc)

eps = .Machine$double.eps
max_iter = 10^5
zero_range = c(-25,1) 
q.max=6

n=500
x=seq(0,1,length=n)


Basis=list()
for (i in 1:q.max)
  Basis[[i]]=drbasis(nn=n,qq=i)

sigma=0.07
pp=1
qq=0
coef=c(0.6) 
set.seed(1121)

## function 1
beta <- function(x,p,q){beta <- gamma(p+q)/(gamma(p)*gamma(q))*x^(p-1)*(1-x)^(q-1) ; beta}
x <- seq(0,1,length.out=n)
mu <- (6*beta(x,30,17)+4*beta(x,3,11))/10;
mu <- (mu-min(mu))/(max(mu)-min(mu))
y=mu+sigma*arima.sim(n = n, list(ar = coef), sd = 1)

## function 2
A5 <- Basis[[5]]$eigenvectors
ev5 <- as.vector(Basis[[5]]$eigenvalues)
coef5 <- c(rep(0,5),(pi*(3:(n-3)))^(-5.1))*(cos(2*(1:n)))
mu5 <- A5%*%coef5
mu=mu5/(max(mu5)-min(mu5))
y=mu+sigma*arima.sim(n = n, list(ar = coef), sd = 1)

## function 3
A3 <- Basis[[3]]$eigenvectors
ev3 <- as.vector(Basis[[3]]$eigenvalues)
coef3 <- c(rep(0,3),(pi*(2:(n-2)))^(-3.1))*(cos(2*(1:n)))
mu3 <- A3%*%coef3
mu=mu3/(max(mu3)-min(mu3))
y=mu+sigma*arima.sim(n = n, list(ar = coef), sd = 1)


R.true=corMatrix(Initialize(corARMA(coef,p=pp,q=qq),data=diag(1:n)))
Tq=rep(NA,q.max-1)
Fit=matrix(NA,n,q.max-1)
R.fit=list()
Phi1=Basis[[1]]$eigenvectors
V.t=eigen(R.true)$vectors

#####
## Iterate for each q

for (q in 2:q.max)
{
  
  Phi=Basis[[q]]$eigenvectorsQR
  eta=as.vector(Basis[[q]]$eigenvalues)
  
  #set R to some initial value (I or something else) and get initial estimates 
  R.hat=diag(n)
  R.ev=eigen(R.hat)
  V=R.ev$vectors
  r=R.ev$values
  
  
  R.q=V%*%diag((r)^(1/(4*q)))%*%t(V)
  SS=eigen(R.q%*%Phi%*%diag((eta*n)^(1/(2*q)))%*%t(Phi)%*%t(R.q))
  U=(SS$vectors)[,n:1]
  eigens=((SS$values)^(2*q))[n:1]
  coefs=as.vector(t(U)%*%V%*%diag(1/sqrt(r))%*%t(V)%*%y)
  temp = try(uniroot(Elambda,zero_range,coefs=coefs,eigens=eigens,n=n,q=q,tol=eps,maxiter=max_iter,extendInt="yes")$root, silent=TRUE) # try to find zero of estimating eaquation
  if(class(temp)=="try-error") log10lambda = 5 else log10lambda = temp # set log-lambda to 5 (lambda is 10^5) if no zero is found
  Lambda=10^log10lambda
  R.hat.inv=V%*%diag(1/r)%*%t(V)
  ridge=chol(t(Phi)%*%R.hat.inv%*%Phi+Lambda*diag(n*eta))
  ridge.inv=backsolve(ridge,diag(rep(1,nrow(ridge))))
  Ridge.inv=ridge.inv%*%t(ridge.inv)
  S=Phi%*%Ridge.inv%*%t(Phi)%*%R.hat.inv
  fit=S%*%y
  
  if (isTRUE(all.equal(R.hat,diag(n)))) V=Phi1
 
  plot(x,y)
  lines(x,mu,lwd=2)
  lines(x,fit,col=2,lwd=2)
   count=1
   

repeat
  {
    
    #get raw rhos
    d=0
    resid=(r^d)*t(V)%*%(diag(n)-S)%*%y
    den=(r^d)*diag(t(V)%*%(diag(n)-S)%*%V)
    nom=diag(resid%*%t(resid))
    rhos=(nom/den)^(1/(d+1))
    rhos=rhos/mean(rhos)

    
    #smooth rhos
    rho.fit=smooth_rhos(rhos,q=2)
    plot(rhos)
    lines(rho.fit,lwd=2)
    lines(eigen(R.true)$values,lwd=2,col=2)
    
    #estimate R
    R.hat1=(Phi1)%*%diag(rho.fit)%*%t(Phi1)
    R.hat1=toeplitz(R.hat1[1,]/R.hat1[1,1])
    
  #  R.hat1=toeplitz(sapply(1:n, function(i) mean(cos(pi*x*(i-1))*rho.fit))) ## works a bit worse
   
    #get new eigenvectors and eigenvalues
    R.hat1.ev=eigen(R.hat1)
    V1=R.hat1.ev$vectors
    r1=R.hat1.ev$values
    
    plot(R.hat1[1:100,40],type="l")
    lines(R.true[1:100,40],col=2)
    
    #get the fit with the new R
    R.q=V1%*%diag((r1)^(1/(4*q)))%*%t(V1)
    SS=eigen(R.q%*%Phi%*%diag((eta*n)^(1/(2*q)))%*%t(Phi)%*%t(R.q))
    U1=(SS$vectors)[,n:1]
    eigens1=((SS$values)^(2*q))[n:1]
    coefs1=as.vector(t(U1)%*%V1%*%diag(1/sqrt(r1))%*%t(V1)%*%y)
    temp = try(uniroot(Elambda,zero_range,coefs=coefs1,eigens=eigens1,n=n,q=q,tol=eps,maxiter=max_iter,extendInt="yes")$root, silent=TRUE) # try to find zero of estimating eaquation
    if(class(temp)=="try-error") log10lambda = 5 else log10lambda = temp # set log-lambda to 5 (lambda is 10^5) if no zero is found
    
    Lambda1=10^log10lambda
    R.hat1.inv=V1%*%diag(1/r1)%*%t(V1)
    ridge=chol(t(Phi)%*%R.hat1.inv%*%Phi+Lambda1*diag(n*eta))
    ridge.inv=backsolve(ridge,diag(rep(1,nrow(ridge))))
    Ridge.inv=ridge.inv%*%t(ridge.inv)
    S1=Phi%*%Ridge.inv%*%t(Phi)%*%R.hat1.inv
    fit1=S1%*%y
   # fit1=V1%*%diag(sqrt(r1))%*%t(V1)%*%U1%*%(coefs1/(1+Lambda1*eigens1))
   
    plot(x,y)
    lines(x,mu,lwd=2)
    lines(x,fit1,col=2,lwd=2)
    lines(x,fit,col=4,lwd=2)
    
    Lambda.diff=abs(log10(Lambda1)-log10(Lambda))
    R.diff=sqrt(sum((r-r1)^2))
    print(q)
    print(Lambda.diff)
    print(R.diff)
  
    if ((Lambda.diff<0.1)&&(R.diff<1)) break
    
    count=count+1
    if (count==15) break
    fit=fit1
    V=V1
    S=S1
    r=r1
    Lambda=Lambda1
    R.hat.inv=R.hat1.inv
    R.hat=R.hat1
  
  }
  
  
  Fit[,q-1]=fit1
  Tq[q-1]=Eq(log10(Lambda1),coefs1,eigens1,n,q,r1)
  R.fit[[q-1]]=R.hat1

  }


plot(2:q.max,Tq,pch=18,ylim=c(-0.5,max(Tq)))
lines(2:q.max,Tq)
lines(2:q.max,rep(0,(q.max-1)),col="grey",lwd=3)

par(mfrow=c(3,2))
for (i in 1:5)
{
  plot(x,y)
  lines(x,mu,lwd=3)
  lines(x,Fit[,i],lwd=3,col=2)
}
par(mfrow=c(1,1))

par(mfrow=c(3,2))
for (i in 1:5)
{
plot((R.fit[[i]])[1:100,40],type="l",ylim=c(-0.5,1),lwd=3)
lines(R.true[1:100,40],col=2)
}
par(mfrow=c(1,1))

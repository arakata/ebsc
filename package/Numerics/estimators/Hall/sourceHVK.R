library(funtimes)
library(wavethresh)
library(nlme)
library(funtimes)

###############
#HVK estimator
##############
#based on:
#Hall and Van Keilegom 2003. 

Hall<-function(xvec,m1,m2,Y){
  
  # only consider AR 1 errors (can be extended)
  # xvec is [0,1] equidistant
  # Y is the signal  

  # => returns h.opt, R.eigenvalues, fhat

  #requires (to compute the Gassler et al. kernel estimators)  
  library(funtimes) #to get the phi's (already implemented by peps. at the University of maryland)
  library(nlme) #to get the covariance given phi's. Known and implemented for ARMA's in nlme 
    
  n <- length(xvec)
  if(is.null(m2)) m2 <- n^0.5
  if(is.null(m1)) m1 <- n^0.4
  h.set<-seq(0.06,0.94,length.out = 45)

#hat estimator#############################################
K<-function(x){
  x[abs(x)>1]<-NA
  o<-((15/16)*(1-x^2)^2)
  o[is.na(o)]<-0
  o
}

s1 <- function(x,xvec,h){
  n<-length(xvec)
  ((n*h)^(-1))*sum(K((x-xvec)/h)*(x-xvec)^1)
}

s2 <- function(x,xvec,h){
  n<-length(xvec)
  ((n*h)^(-1))*sum(K((x-xvec)/h)*((x-xvec)^2))
}

w <- function(x,xvec,k,h){
  K((x-xvec[k])/h)*(s2(x,xvec,h)-(x-xvec[k])*s1(x,xvec,h))
}

ghat <- function(xvec,h,Y){
  n<-length(xvec)
  ones<-rep(1,n)
  W<-matrix(NA,n,n)
  for(i in 1:n){W[i,]=w(x=xvec[i],xvec=xvec,k=1:n,h)}
  
  W%*%Y/W%*%ones
} 

#ghatloo is ghat with leave one out (for cross-validation)
ghatloo <- function(xvec,h,Y){
  n<-length(xvec)
  ones<-rep(1,n)
  W<-matrix(NA,n,n)
  o<-rep(NA,n)
  
  for(i in 1:n){W[i,]=w(x=xvec[i],xvec=xvec[-i],k=1:n,h)}
  W<-W[,-n]
  
  for(i in 1:n){
    o[i]<-(W[i,]%*%Y[-i])/(W[i,]%*%ones[-i])
  }
  o
}

Hall.gamma0 <- function(m2,m1,Y){
    n=length(Y);
    gamma0=0;
    for(i in m1:m2){
        temp1 = Y[1:(n-i)]
        temp2 = Y[(i+1):n]
        temp = mean( (temp2-temp1)^2 )/2
        gamma0 = gamma0 + temp
    }
    return( gamma0/(m2-m1+1) )
}
  
#end of hat estimator######################################

#tilde estimator###########################################
L<-function(x){
  x[abs(x)>1]<-NA
  x[x<0]<-NA
  o<-((15/16)*(1-x^2)^2)
  o[is.na(o)]<-0
  o
}

r1 <- function(x,xvec,h){
  n<-length(xvec)
  ((n*h)^(-1))*sum(L((x-xvec)/h)*(x-xvec)^1)
}

r2 <- function(x,xvec,h){
  n<-length(xvec)
  ((n*h)^(-1))*sum(L((x-xvec)/h)*((x-xvec)^2))
}

v <- function(x,xvec,k,h){
  L((x-xvec[k])/h)*(r2(x,xvec,h)-(x-xvec[k])*r1(x,xvec,h))
}
#end of tilde estimator####################################

phi.hat <- pmax(-.999, pmin(.999, HVK(X=as.vector(Y),m1=m1,m2=m2,ar.order=1)))
myDATA <- data.frame(X=Y,index=1:n)
csARMA <- corARMA(phi.hat, form = ~index, p = 1,q=0);
csARMA <- Initialize(csARMA, data = myDATA)
R <- corMatrix(csARMA)
R.eigenvalues<-Re(eigen(R)$values)
  
CV<-rep(NA,length(h.set))  
for(i in 1:length(h.set)){
  fit<-ghatloo(xvec=xvec,h=h.set[i],Y=Y)
  error<-Y-fit
  difference<-c(na.omit(error-phi.hat*c(NA,error[-n])))  
  CV[i]<-sum((difference)^2)/length(difference)
}

h.hat<-h.set[which.min(CV)]
f.hat <- ghat(xvec=xvec,h=h.hat,Y=Y)
sigma2.hat<-Hall.gamma0(m2=m2,m1=m1,Y=Y)
  
#list(h.hat=h.hat, R.eigenvalues=R.eigenvalues,f.hat=f.hat,sigma2.hat=sigma2.hat)
list(Rhat=R,yhat=f.hat)
}

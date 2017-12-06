library(funtimes)
###############
#HERR estimator
##############
#based on:
#Herrmann et al. 1992. 

Her<-function(xvec,m,Y){

  # considers m-dependent errors
  # xvec is [0,1] equidistant
  # Y is the signal  
  # m is the dependence degree. 

  # => returns h.opt, R.eigenvalues, fhat

  #requires (to compute the Gassler et al. kernel estimators)  
  library(lokern)
  
  n<-length(xvec)
  xgrid<-seq(-1,1,length.out = n)
  dx<-diff(xvec)[1]
  dxgrid<-diff(xgrid)[1]
  
  #functions needed
 
  aux.poly <- function(d){
    return( 1 + d + d^2 )
  }
  
  gap.estimator <- function(m,d,v){ # Difference-based estimate of length 2(m+1) and length 3 calculated on vector v
    n<- length(v)
    sum((v[1:(n-2*(m+1))] - (1+d)*v[(m+2):(n-m-1)] + d*v[(2*m+3):n])^2)/( 2*aux.poly(d)*(n-2*(m+1)) )
  }
  
  gkH<-function(m,k,v){
    n=length(v)
    
    aux<-gap.estimator(m,1,v)
    
    alphak<-(m+k+1)/(2*(m+1))
    
    betak<--alphak*( (k/(m+k+1))^2 + ((m+1)/(m+k+1))^2 + 1)
    
    temp<-alphak*sum( ( (k+m+1)*v[(k+1):(n-(m+1))] - k*v[(m+k+2):n] - (m+1)*v[1:(n-k-m-1)] )^2 )/((k+m+1)^2*(n-m-k-1))
    
    -(aux*betak+temp)
  }
  
  gkHn<-function(m,k,v){
    n<-length(v)
    aux<-gap.estimator(m,1,v)
    
    alphak<-(m+k+1)/(2*(m+1))
    
    betak<--alphak*( (k/(m+k+1))^2 + ((m+1)/(m+k+1))^2 + 1)
    
    temp<- (alphak/((m+k+1)^2*(n-m-1-k)))*sum( ((m+k+1)*v[(m+2):(n-k)] - (m+1)*v[(m+2+k):n] - k*v[1:(n-k-m-1)])^2 )
    
    -(aux*betak+temp)
  }
  
  covarianceH<-function(m,v){
    # l=length(v)
    cov<- numeric(m)
    for(i in 1:m){
      cov[i] <- gkH(m,i,v)
    }
    # cov[ cov < 0 ] <- 0
    c(gap.estimator(m,1,v),cov)
  }
  
  covarianceHn<-function(m,v){
    # l=length(v)
    cov<- numeric(m)
    for(i in 1:m){
      cov[i] <- gkHn(m,i,v)
    }
    # cov[ cov < 0 ] <- 0
    c(gap.estimator(m,1,v),cov)
  }
  

  # non-negative function with support in x and non-vanishing derivatives at the borders
  v<-function(x){
      mu=0.5
      sigma=0.2
      (1/(sqrt(2*pi)*sigma))*exp(-((x-mu)^2)/(2*sigma^2))
  }

  # a fourth degree kernel (need to use the second derivative)
  W<-function(x){
      x[abs(x)>1]<-NA
      o<-((15/16)*(1-x^2)^2)
      o[is.na(o)]<-0
      o
  }

  ghat<-function(xvec,b,sigma,Y){
      glkerns(x=xvec, y=Y, deriv = 0, x.out=xvec,
        korder= 4, hetero=FALSE, is.rand=FALSE,
        inputb= FALSE,
        m1 = 400, xl=NULL, xu=NULL,
        s=NULL, sig=sigma, bandwidth=NULL, trace.lev = 0)$est
  }

  g2hat<-function(xvec,b,sigma,Y){
      glkerns(x=xvec, y=Y, deriv = 2, x.out=xvec,
        korder= 4, hetero=FALSE, is.rand=FALSE,
        inputb= FALSE,
        m1 = 400, xl=NULL, xu=NULL,
        s=NULL, sig=sigma, bandwidth=NULL, trace.lev = 0)$est
  }

  #####################
  #bandwidth selection 
  #####################

  #input
  C1 <- sum(v(x=xvec)*dx)*sum((W(x=xgrid)^2)*dxgrid)
  C2 <- sum(((xgrid)^2)*W(x=xgrid)*dxgrid)
  
  covvec <- covarianceH(m=m,v=y)
  sigma<-sqrt(covvec[1])
  gammavec<-covvec[-1]
  gammavec[gammavec<10^(-2)] <- 0
  
  fr<-rep(0,n)
  fr[1:length(covarianceH(m,Y))]<-covarianceH(m,Y)/covarianceH(m,Y)[1]
  R<-toeplitz(fr)
  R.eigenvalues<-Re(eigen(R)$values)
  
  S <- sigma^2+2*sum(gammavec)

  sigma2.hat<-sigma^2
  
  # loop
  bvec<-rep(NA,12)
  bvec[1]<-1/n
  for(i in 2:12){
      g2vec<-g2hat(xvec=xvec,b=bvec[i-1]*(n^0.1),sigma=sqrt(sigma2.hat),Y=y)
      bvec[i]<-((C1*S)/(n*(C2^2)*(sum((v(xvec)*(g2vec^2))*dx))))^(1/5)
  }
  h.hat<-(bvec)[12]

  f.hat<-ghat(xvec=xvec,b=b.hat,sigma=sqrt(sigma2.hat),Y=y)

#list(h.hat=h.hat, R.eigenvalues=R.eigenvalues,f.hat=f.hat,sigma2.hat=sigma2.hat)
list(Rhat=R,yhat=f.hat)

}

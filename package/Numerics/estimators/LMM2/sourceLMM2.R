library(Brobdingnag)
library(nlme)
library(wavethresh)

lmm2<-function(y,ARpMAq, Basis){
  
  mycorMatrix <-function(mycorcoefs,myDATA){
    myDATA <- data.frame(myDATA)
    index <- myDATA$index
    string <- names(mycorcoefs)
    substring <- substring(string,1,3)
    ar.ord <- length(substring[substring=="Phi"])
    ma.ord <- length(substring[substring=="The"])
    csARMA <- corARMA(mycorcoefs, form = ~index, p = ar.ord,q=ma.ord);
    csARMA <- Initialize(csARMA, data = myDATA);
    corMatrix(csARMA);}
  
  lmm <-function (y, q, parameters){
    n <- length(y)
    basis <- Basis[[q]]
    col.ones <- rep(1,n)
    DATA <- cbind(index=seq(1,length.out=n),y=as.numeric(y))
    X <- basis$eigenvectorsQR[,1:q]
    Z <- basis$eigenvectorsQR[,(q+1):n]%*%diag(1/sqrt(basis$eigenvalues[(q+1):n]))
    N <- basis$eigenvectorsQR
    D <- diag(basis$eigenvalues)
    if(is.null(parameters)){
      lme.fit <- lme(y ~ -1+X, random = list(col.ones=(pdIdent(~-1+Z))),correlation=NULL, control = lmeControl(singular.ok=TRUE, returnObject=TRUE, opt = "optim"))
      sigma2.e<- lme.fit$sigma^2
      sigma2.b<-(lme.fit$sigma^2)*exp(2*unlist(lme.fit$modelStruct))
      lambda<- ((sigma2.e)/(sigma2.b))[1] 
      V <- diag(n)}else{
        lme.fit <- lme(y ~ -1+X, random = list(col.ones=(pdIdent(~-1+Z))),correlation=eval(parse(text=parameters)), control = lmeControl(singular.ok=TRUE, returnObject=TRUE, opt = "optim"))
        sigma2.e <- lme.fit$sigma^2
        sigma2.b<-(lme.fit$sigma^2)*exp(2*unlist(lme.fit$modelStruct))
        lambda <- ((sigma2.e)/(sigma2.b))[1]
        corcoefs <- coef(lme.fit$modelStruct$corStruct,unconstrained=FALSE)
        V <- mycorMatrix(corcoefs,DATA)}
    list(sigma2.hat=sigma2.e,lambda.hat=log(lambda/n,10),R.hat=V,X=X,Z=Z,N=N,D=D,yhat=lme.fit$fitted[,2])}
  
  q=2
  n=length(y)
  E=exp(1)
  x=seq(0,1,length=n)
  
  if(!(exists("Basis")==TRUE&&nrow(Basis[[1]]$eigenvectorsQR)==n)){
    Basis <- list()
    for(i in 1:6) Basis[[i]] <- drbasis(nn=n,qq=i)
  }
  
  if((ARpMAq[1]==0)&&(ARpMAq[2]==0)){
    correlation <-NULL
    mm <- lmm(y=y,q=q,parameters=correlation)
  }else{
    correlation <-parse(text = paste("corARMA(","p=",ARpMAq[1],",q=",ARpMAq[2],")",sep=""))
    mm <- lmm(y=y,q=q,parameters=correlation)
  }
  
  list(Rhat=mm$R.hat,yhat=mm$yhat)
  
}
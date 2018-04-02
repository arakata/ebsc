bspline <- function(x, xl, xr, ndx, bdeg){
    dx<-(xr-xl)/ndx
    knots<-seq(xl-bdeg*dx,xr+bdeg*dx,by=dx)
    B <- spline.des(knots, x, bdeg + 1, 0 * x)$design
    B
}
    
eilers<-function(num.knots,bdeg,pord,y){
    n=length(y)
    x=1:n
    xl=1
    xr=n
    ndx=num.knots
    #L<-t(chol(R))

    B <- bspline(x, xl, xr, ndx, bdeg)
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

    list(yhat=yhat,edf=edf,gcv=gcv,index=index)
}



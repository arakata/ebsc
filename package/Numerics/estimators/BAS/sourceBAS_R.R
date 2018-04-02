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

smooth_rhos = function(rhos, q)
{ #smoothing rhos
  n=length(rhos)
  x=seq(1,n,length=n)
  Basis = drbasis(nn=n,qq=q)
  PP = Basis$eigenvectorsQR
  ee = as.vector(n*Basis$eigenvalues)
  coefs = t(PP) %*% rhos
  temp = try(uniroot(Elambda,zero_range,coefs=coefs,eigens=ee,n=n,q=q,tol = eps,maxiter=max_iter,extendInt="yes")$root, silent = TRUE)
  if(class(temp)=="try-error"  | temp >zero_range[2]){temp = -0} ## ???
  # if(class(temp)=="try-error"){temp = 5}
  rho.fit = PP %*% (coefs/(1+(10^temp)*ee)) 
  
  return( pmax( 10^-5, rho.fit) )
}




EBSC = function(data, DRbases, q = NULL, iterations = 2)
{ ## applies our method; DR basis for q=1,...,6 must be passed as a list argument for the relevant n
  
  n = length(data)
  
  if(is.null(q)) # if no q is specified pick it suing Tq
  {
    
    
    
    
    
  }else{ # if Tq is specified then use that value of q
    
    result = NULL
    result$lambdas = rep(1, iterations)
    result$rhos = matrix(1, nrow = n, ncol = iterations)
    result$sigmasqs = rep(1, iterations)
    result$qs = rep(q, iterations)
    result$estimates = matrix(1, nrow = n, ncol = iterations)
    
    Phi=DRbases[[q]]$eigenvectorsQR
    Phi1=DRbases[[1]]$eigenvectors
    eta=as.vector(DRbases[[q]]$eigenvalues)
    Phi.q=Phi[,1:q]
    Phi.nq=Phi[,(q+1):n]
    eta.q=eta[(q+1):n]
    
    #set R to some initial value (I or something else) and get initial estimates 
    R.hat=diag(n)
    R.ev=eigen(R.hat)
    V=R.ev$vectors
    r=R.ev$values

    R.q=V%*%diag((r)^(1/(4*q)))%*%t(V)
    SS=eigen(R.q%*%Phi%*%diag((eta*n)^(1/(2*q)))%*%t(Phi)%*%t(R.q))
    U=(SS$vectors)[,n:1]
    eigens=((SS$values)^(2*q))[n:1]
    coefs=as.vector(t(U)%*%V%*%diag(1/sqrt(r))%*%t(V)%*%data)
    temp = try(uniroot(Elambda,zero_range,coefs=coefs,eigens=eigens,n=n,q=q,tol=eps,maxiter=max_iter,extendInt="yes")$root, silent=TRUE) # find zero
    if(class(temp)=="try-error" | temp >zero_range[2]) log10lambda = 5 else log10lambda = temp # set log-lambda to 5 if no zero is found ???
    Lambda=10^log10lambda; result$lambdas[1] = Lambda

    W=solve(R.hat+Phi.nq%*%diag(1/(Lambda*n*eta.q))%*%t(Phi.nq))
    IS=R.hat%*%W%*%(diag(n)-Phi.q%*%solve(t(Phi.q)%*%W%*%Phi.q)%*%t(Phi.q)%*%W)
    S=diag(n)-IS
    
    # result$sigmasqs[1] = (sum(coefs^2*(eigens*Lambda/(1 + eigens*Lambda)))+1)/(n+1)
    result$sigmasqs[1] = (t(data)%*%W%*%(diag(n)-Phi.q%*%solve(t(Phi.q)%*%W%*%Phi.q)%*%t(Phi.q)%*%W)%*%data+1)/(n+1)
    
    # R.hat.inv=V%*%diag(1/r)%*%t(V)
    # ridge=chol(t(Phi)%*%R.hat.inv%*%Phi+Lambda*diag(n*eta))
    # ridge.inv=backsolve(ridge,diag(rep(1,nrow(ridge))))
    # Ridge.inv=ridge.inv%*%t(ridge.inv)
    # S=Phi%*%Ridge.inv%*%t(Phi)%*%R.hat.inv
    
    fit=S%*%data; result$estimates[, 1] = fit
    
    if (isTRUE(all.equal(R.hat,diag(n)))) V=Phi1
    
    count=1
    
    repeat
    {
      
      #get raw rhos
      # d=0
      # resid=(r^d)*t(V)%*%(diag(n)-S)%*%data
      # den=(r^d)*diag(t(V)%*%(diag(n)-S)%*%V)
      # nom=diag(resid%*%t(resid))
      # rhos=(nom/den)^(1/(d+1))
      # rhos=rhos/mean(rhos)
      # plot(rhos)
      resid=t(V)%*%IS%*%data
      den=diag(t(V)%*%IS%*%V)
      nom=diag(resid%*%t(resid))
      rhos=nom/den
      rhos=rhos/mean(rhos)
      
      rhos[1:q] = rhos[q+1] # ???
      
      #smooth rhos
      rho.fit=smooth_rhos(rhos,q=2); result$rhos[,count+1] = rho.fit
      # lines(rho.fit)
      
      #estimate R
      R.hat1=(Phi1)%*%diag(rho.fit)%*%t(Phi1)
      R.hat1=toeplitz(R.hat1[1,]/R.hat1[1,1])
      
      #get new eigenvectors and eigenvalues
      # R.hat1.ev=eigen(R.hat1)
      # V1=R.hat1.ev$vectors
      # r1=R.hat1.ev$values
      R.hat1.ev=eigen(R.hat1)
      V1=R.hat1.ev$vectors
      r1=R.hat1.ev$values
      if (rho.fit[n]-rho.fit[1]>0)
      {
        r1=r1[n:1]
        V1=V1[,n:1]
      }
      
      #get the fit with the new R
      R.q=V1%*%diag((r1)^(1/(4*q)))%*%t(V1)
      SS=eigen(R.q%*%Phi%*%diag((eta*n)^(1/(2*q)))%*%t(Phi)%*%t(R.q))
      U1=(SS$vectors)[,n:1]
      eigens1=((SS$values)^(2*q))[n:1]
      coefs1=as.vector(t(U1)%*%V1%*%diag(1/sqrt(r1))%*%t(V1)%*%data)
      
      temp = try(uniroot(Elambda,zero_range,coefs=coefs1,eigens=eigens1,n=n,q=q,tol=eps,maxiter=max_iter,extendInt="yes")$root, silent=TRUE) # find zero
      if(class(temp)=="try-error" | temp >zero_range[2]) log10lambda = 5 else log10lambda = temp # set log-lambda to 5 (lambda is 10^5) if no zero is found ???
      Lambda1=10^log10lambda; result$lambdas[count+1] = Lambda1
        
      # R.hat1.inv=V1%*%diag(1/r1)%*%t(V1)
      # ridge=chol(t(Phi)%*%R.hat1.inv%*%Phi+Lambda1*diag(n*eta))
      # ridge.inv=backsolve(ridge,diag(rep(1,nrow(ridge))))
      # Ridge.inv=ridge.inv%*%t(ridge.inv)
      # S1=Phi%*%Ridge.inv%*%t(Phi)%*%R.hat1.inv
      W1=solve(R.hat1+Phi.nq%*%diag(1/(Lambda1*n*eta.q))%*%t(Phi.nq))
      IS1=R.hat1%*%W1%*%(diag(n)-Phi.q%*%solve(t(Phi.q)%*%W1%*%Phi.q)%*%t(Phi.q)%*%W1)
      S1=diag(n)-IS1
      
      # result$sigmasqs[count+1] = (sum(coefs1^2*(eigens1*Lambda1/(1 + eigens1*Lambda1)/rho.fit))+1)/(n+1)
      result$sigmasqs[count+1] = (t(data)%*%W1%*%(diag(n)-Phi.q%*%solve(t(Phi.q)%*%W1%*%Phi.q)%*%t(Phi.q)%*%W1)%*%data+1)/(n+1)
      
      fit1=S1%*%data; result$estimates[, count+1] = fit1
      # fit1=V1%*%diag(sqrt(r1))%*%t(V1)%*%U1%*%(coefs1/(1+Lambda1*eigens1))

      Lambda.diff=abs(log10(Lambda1)-log10(Lambda))
      R.diff=sqrt(mean((r-r1)^2))
      
      if ((Lambda.diff<0.1)&&(R.diff<0.1)) break
      
      count=count+1
      if (count==iterations) break
      fit=fit1
      V=V1
      IS=IS1
      S=S1
      r=r1
      Lambda=Lambda1
      R.hat=R.hat1
    }
    
    result$eigens = eigens1
    result$iterations = count
    result$estimate = fit1
    result$lambda = Lambda1
    result$sigmasq = result$sigmasqs[iterations]
    result$rho = rho.fit
    result$R = R.hat1
    result$S = S1
    
    ## see if the following can be computed from the new matrices
    R.hat1.inv=V1%*%diag(1/r1)%*%t(V1)
    ridge=chol(t(Phi)%*%R.hat1.inv%*%Phi+Lambda1*diag(n*eta))
    ridge.inv=backsolve(ridge,diag(rep(1,nrow(ridge))))
    Ridge.inv=ridge.inv%*%t(ridge.inv)
    result$SR = Phi%*%Ridge.inv%*%t(Phi)
    # result$SR = R.hat-IS%*%R.hat ## this does not work
    
    result$q = q
    result$Tq = rep(NA, q.max)
  }
  
  return( result )
}

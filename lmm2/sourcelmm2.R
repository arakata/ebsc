library(Brobdingnag)
library(nlme)
library(wavethresh)

lmm2<-function(y,ARpMAq){

    drbasis <-function(nn,qq){
        options(expressions=10000)

        #inputs
        n<<-nn
        x<<- seq(0,1,length.out=n)
        kset<<-seq(1,n);
        tset<<-seq(0,1,length.out=n);
        q<-qq
        q.global<<-qq
    
        #constants
        E<-exp(1)
        I<-1.i
        Pi<-pi

        const0<-1
        const1<- -sqrt(3) 
        const2<- 2*sqrt(3)
        const3<- -sqrt(5)
        const4<- 6*sqrt(5)
        const5<- - 6*sqrt(5)
        const6<- -sqrt(7)
        const7<- 12*sqrt(7)
        const8<- -30*sqrt(7)
        const9<- 20*sqrt(7)
        const10<- 3
        const11<- -60
        const12<- 270
        const13<- -420
        const14<- 210
        const15<- -sqrt(11)
        const16<- 30*sqrt(11)
        const17<- -210*sqrt(11)
        const18<- 560*sqrt(11)
        const19<- 630*sqrt(11)
        const20<- 252*sqrt(11)    

        #funcs
        Power<-function(x,y){
            if(length(x)<2){
                if(x==E){
                    ee<-as.brob(E);
                    yRe<-Re(y);
                    yIm<-complex(real=0,imaginary=Im(y));
                    ans<-(ee^yRe)*(ee^yIm)}else{ans<-x^y}
            }else{
                ans<-x^y
            }
            ans
        }
        Complex<-function(x,y) complex(real=x,imaginary=y)
    
        turn.Mo<-function(M,Mo){
            norm<-function(x){(sum(x^2))^(1/2)}
            sign.change<-function(x) {if(x[1]<=x[2]){-1}else{1}}
            D<-cbind(apply(M+Mo,2,norm),apply(M-Mo,2,norm))
            d<-apply(D,1,sign.change)
            cMo<-Mo%*%diag(d)
            cMo
        }

        #############
        #Eigenvectors
        #############
    
        if(q==1){
        #NULL solution
            nullspace<-rep(1,n)
        
        #ODE solution
            kset<-(q.global+1):n
            phi<-function(t,k){sqrt(2)*cos((-1 + k)*Pi*t)}        

        #output
            M<-cbind(nullspace,outer(tset,kset,FUN=phi))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-turn.Mo(M,Mo)
    }
    
        if(q==2){
        #NULL solution
            const0<-1
            const1<--sqrt(3) 
            const2<-2*sqrt(3) 
            nullspace<-cbind(const0,const1 + const2*tset)
        
        #ODE solution
            kset<-(q.global+1):n
            phi<-function(t,k){
                ans<-as.complex(
                    Power(as.complex(-1),1 + k)/Power(E,(-1.5 + k)*Pi*(1 - t)) + 
                        Power(E,-((-1.5 + k)*Pi*t)) + sqrt(2)*cos(Pi/4. + (-1.5 + k)*Pi*t)
                )
                Re(ans)
            }
        
        #output
            M<-cbind(nullspace,outer(tset,kset,FUN=phi))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-turn.Mo(M,Mo)
        }
    
        if(q==3){
        #NULL solution
            nullspace<-cbind(const0,const1 + const2*tset,const3 + const4*tset+const5*tset^2)

        #ODE solution
            kset<-(q.global+1):n
            phi<-function(t,k){
                Re(as.complex(
                    (sqrt(1.5) - Complex(0,1)/sqrt(2))*
                        (Power(as.complex(-1),1 + k)/
                             Power(E,Power(as.complex(-1),1/6)*(-2 + k)*Pi*(1 - t)) + 
                                 Power(E,-(Power(as.complex(-1),1/6)*(-2 + k)*Pi*t))) + 
                                     (sqrt(1.5) + Complex(0,1)/sqrt(2))*
                                         (Power(as.complex(-1),1 + k)/
                                              Power(E,(Complex(0,-0.5) + sqrt(3)/2.)*(-2 + k)*Pi*(1 - t)) + 
                                                  Power(E,-((Complex(0,-0.5) + sqrt(3)/2.)*(-2 + k)*Pi*t))) - 
                                                      sqrt(2)*sin((-2 + k)*Pi*t)
                ))
            }        

        #output
            M<-cbind(nullspace,outer(tset,kset,FUN=phi))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-turn.Mo(M,Mo)
        }
        
        if(q==4){
        #NULL solution
            nullspace<-cbind(const0,const1 + const2*tset,const3 + const4*tset+const5*tset^2,const6+const7*tset+const8*tset^2+const9*tset^3)

        #ODE solution
            kset<-(q.global+1):n
            phi<-function(t,k){
                Re(as.complex(
                    (1 + sqrt(2))*(Power(as.complex(-1),1 + k)/Power(E,(-2.5 + k)*Pi*(1 - t)) + 
                                       Power(E,-((-2.5 + k)*Pi*t))) + 
                                           (1/sqrt(2) - Complex(0,1)*(1 + 1/sqrt(2)))*
                                               (Power(as.complex(-1),1 + k)/Power(E,Power(as.complex(-1),0.25)*(-2.5 + k)*Pi*(1 - t)) + 
                                                    Power(E,-(Power(as.complex(-1),0.25)*(-2.5 + k)*Pi*t))) + 
                                                        (1/sqrt(2) + Complex(0,1)*(1 + 1/sqrt(2)))*
                                                            (Power(as.complex(-1),1 + k)/
                                                                 Power(E,(Complex(1,-1)*(-2.5 + k)*Pi*(1 - t))/sqrt(2)) + 
                                                                     Power(E,(Complex(-1,1)*(-2.5 + k)*Pi*t)/sqrt(2))) - 
                                                                         sqrt(2)*cos(Pi/4. - (-2.5 + k)*Pi*t)
                ))
            }
            
        #output
            M<-cbind(nullspace,outer(tset,kset,FUN=phi))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-turn.Mo(M,Mo)
        }
    
        if(q==5){
        #NULL solution
            nullspace<-cbind(const0,const1 + const2*tset,const3 + const4*tset+const5*tset^2,
                             const6+const7*tset+const8*tset^2+const9*tset^3,
                             const10+const11*tset+const12*tset^2+const13*tset^3+const14*tset^4)
            
        #ODE solution
            kset<-(q.global+1):n
            phi<-function(t,k){
                Re(as.complex(
                    (sqrt(2)*(1 + sqrt(5)/2.) - (Complex(0,0.5)*
                                                     (sqrt(10 - 2*sqrt(5)) + sqrt(2*(5 + sqrt(5)))))/sqrt(2))*
                                                         (Power(as.complex(-1),1 + k)/Power(E,Power(as.complex(-1),0.1)*(-3 + k)*Pi*(1 - t)) + 
                                                              Power(E,-(Power(as.complex(-1),0.1)*(-3 + k)*Pi*t))) + 
                                                                  (-(1/sqrt(2)) - (Complex(0,0.5)*
                                                                                       (sqrt(10 - 2*sqrt(5)) + sqrt(2*(5 + sqrt(5)))))/sqrt(2))*
                                                                                           (Power(as.complex(-1),1 + k)/Power(E,Power(as.complex(-1),0.3)*(-3 + k)*Pi*(1 - t)) + 
                                                                                                Power(E,-(Power(as.complex(-1),0.3)*(-3 + k)*Pi*t))) + 
                                                                                                    (sqrt(2)*(1 + sqrt(5)/2.) + 
                                                                                                         (Complex(0,0.5)*(sqrt(10 - 2*sqrt(5)) + sqrt(2*(5 + sqrt(5)))))/sqrt(2)
                                                                                                     )*(Power(as.complex(-1),1 + k)/
                                                                                                            Power(E,(sqrt(0.625 + sqrt(5)/8.) - Complex(0,0.25)*(-1 + sqrt(5)))*
                                                                                                                      (-3 + k)*Pi*(1 - t)) + 
                                                                                                                          Power(E,-((sqrt(0.625 + sqrt(5)/8.) - Complex(0,0.25)*(-1 + sqrt(5)))*
                                                                                                                                        (-3 + k)*Pi*t))) + (-(1/sqrt(2)) + 
                                                                                                                                                                (Complex(0,0.5)*(sqrt(10 - 2*sqrt(5)) + sqrt(2*(5 + sqrt(5)))))/sqrt(2)
                                                                                                                                                            )*(Power(as.complex(-1),1 + k)/
                                                                                                                                                                   Power(E,(sqrt(0.625 - sqrt(5)/8.) - Complex(0,0.25)*(1 + sqrt(5)))*
                                                                                                                                                                             (-3 + k)*Pi*(1 - t)) + 
                                                                                                                                                                                 Power(E,-((sqrt(0.625 - sqrt(5)/8.) - Complex(0,0.25)*(1 + sqrt(5)))*
                                                                                                                                                                                               (-3 + k)*Pi*t))) - sqrt(2)*cos((-3 + k)*Pi*t)
                ))
            }
                
        #output
            M<-cbind(nullspace,outer(tset,kset,FUN=phi))/sqrt(n)
            Mo<-qr.Q(qr(M));Mo<-turn.Mo(M,Mo)
        }
    
    if(q==6){
        #NULL solution
        nullspace<-cbind(const0,const1 + const2*tset,const3 + const4*tset+const5*tset^2,
                         const6+const7*tset+const8*tset^2+const9*tset^3,
                         const10+const11*tset+const12*tset^2+const13*tset^3+const14*tset^4,
                         const15+const16*tset+const17*tset^2+const18*tset^3+const19*tset^4+const20*tset^5)
        
        #ODE solution
        kset<-(q.global+1):n
        phi<-function(t,k){
            Re(as.complex(
            (3 + 2*sqrt(3))*(Power(as.complex(-1),1 + k)/Power(E,(-3.5 + k)*Pi*(1 - t)) + 
                                 Power(E,-((-3.5 + k)*Pi*t))) + 
                                     ((1 + sqrt(3))/2. - Complex(0,0.5)*(5 + 3*sqrt(3)))*
                                         (Power(as.complex(-1),1 + k)/
                                              Power(E,Power(as.complex(-1),1/6)*(-3.5 + k)*Pi*(1 - t)) + 
                                                  Power(E,-(Power(as.complex(-1),1/6)*(-3.5 + k)*Pi*t))) + 
                                                      ((-3 - sqrt(3))/2. - Complex(0,0.5)*(1 + sqrt(3)))*
                                                          (Power(as.complex(-1),1 + k)/
                                                               Power(E,Power(as.complex(-1),1/3)*(-3.5 + k)*Pi*(1 - t)) + 
                                                                   Power(E,-(Power(as.complex(-1),1/3)*(-3.5 + k)*Pi*t))) + 
                                                                       ((-3 - sqrt(3))/2. + Complex(0,0.5)*(1 + sqrt(3)))*
                                                                           (Power(as.complex(-1),1 + k)/
                                                                                Power(E,(0.5 - Complex(0,0.5)*sqrt(3))*(-3.5 + k)*Pi*(1 - t)) + 
                                                                                    Power(E,-((0.5 - Complex(0,0.5)*sqrt(3))*(-3.5 + k)*Pi*t))) + 
                                                                                        ((1 + sqrt(3))/2. + Complex(0,0.5)*(5 + 3*sqrt(3)))*
                                                                                            (Power(as.complex(-1),1 + k)/
                                                                                                 Power(E,(Complex(0,-0.5) + sqrt(3)/2.)*(-3.5 + k)*Pi*(1 - t)) + 
                                                                                                     Power(E,-((Complex(0,-0.5) + sqrt(3)/2.)*(-3.5 + k)*Pi*t))) - 
                                                                                                         sqrt(2)*cos(Pi/4. + (-3.5 + k)*Pi*t)
            ))
        }
        
        #output
        M<-cbind(nullspace,outer(tset,kset,FUN=phi))/sqrt(n)
        Mo<-qr.Q(qr(M));Mo<-turn.Mo(M,Mo)
    }

#############
#Eigenvalues
#############
    
        eigenvalues <- function(n,q){
        Q <- function(x,q){
            if (q==1) Q=1
            if (q==2) Q=1/3+2*cos(pi*x)^2/3
            if (q==3) Q=2/15+11*cos(pi*x)^2/15+2*cos(pi*x)^4/15
            if (q==4) Q=(17+180*cos(pi*x)^2+114*cos(pi*x)^4+4*cos(pi*x)^6)/315
            if (q==5) Q=(62+1072*cos(pi*x)^2+1452*cos(pi*x)^4+247*cos(pi*x)^6+2*cos(pi*x)^8)/2835
            if (q==6) Q=(1382+35396*cos(pi*x)^2+83021*cos(pi*x)^4+34096*cos(pi*x)^6+2026*cos(pi*x)^8+4*cos(pi*x)^10)/155295
            Q
        }
            
        sinc <- function(x) {sin(pi*x)/(pi*x)}
        
        s.vec <- (pi*(1:(n-q)+1/2*( (q+1)%%2 )+floor((q-1)/2)))^(2*q)/n
        j=1:(2*n)
        atten <- ((sinc(j/(2*n)))^(2*q)/Q(j/(2*n),q))[(q+1):n]
        s.vec <- s.vec*atten
            s.vec
    }
        ev<-c(rep(0,q),eigenvalues(n,q))[1:n]
    
        eigenvalues=ev
        eigenvectorsQR=Mo
        eigenvectors=M
        
    #final output
        list(eigenvectors=eigenvectors,eigenvectorsQR=eigenvectorsQR,eigenvalues=eigenvalues,x=x)
    }

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
        lme.fit <- lme(y ~ -1+X, random = list(col.ones=(pdIdent(~-1+Z))),correlation=NULL)
        sigma2.e<- lme.fit$sigma^2
        sigma2.b<-(lme.fit$sigma^2)*exp(2*unlist(lme.fit$modelStruct))
        lambda<- ((sigma2.e)/(sigma2.b))[1] 
        V <- diag(n)}else{
            lme.fit <- lme(y ~ -1+X, random = list(col.ones=(pdIdent(~-1+Z))),correlation=eval(parse(text=parameters)))
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

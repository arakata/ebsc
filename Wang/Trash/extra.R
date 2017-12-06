library(SemiPar)
library(splines)

COiterations<-function(num.knots,degree,y){
    n=length(y)
    x=seq(1,n,length.out=n)
    ndx=num.knots+1
    dx=(max(x)-min(x))/ndx
    knots=seq(min(x)-degree*dx,max(x)+degree*dx,by=dx)
    B=spline.des(knots,x,degree+1,0*y,outer.ok=T)$design
    y.hat=B%*%solve(t(B)%*%B)%*%t(B)%*%y
    list(y.hat=y.hat)
}

y=unique(log.income)
y.hat=COiterations(num.knots=5,degree=3,y=y)$y.hat


function B = bspline(x, xl, xr, ndx, bdeg)
dx=(xr-xl)/ndx;
t = xl + dx * [-bdeg:ndx-1];
T=(0*x+1)*t;
X=x*(0*t+1);
P=(X-T)/dx;
B=(T<=X)&(X<(T+dx));
r = [2:length(t) 1];
for k = 1:bdeg
B=(P.*B+(k+1-P).*B(:,r))/k;
end;
end;

B = bspline(x, xl, xr, ndx, bdeg);
[m n] = size(B);
D = diff(eye(n), pord);
a=(B’*B+lambda*D’*D)\(B’*y);
yhat=B*a;
Q=inv(B’*B+lambda*D’*D);
s=sum((y-yhat).^2)
t = sum(diag(Q * (B’ * B)));
gcv=s/(m-t)^2;

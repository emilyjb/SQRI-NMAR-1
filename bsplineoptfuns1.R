 
 ## optimize the penalized quantile regression under B-splines
 ## initial values are coefficient of unpenalized quantile regression
 ploss=function(b,tau,y,B,lambda,D)
 {
   resid=y-B%*%b
   U=tau-as.numeric(resid<0)
   loss=sum(resid*U)+lambda*b%*%t(D)%*%D%*%b
   return(loss)
 }
 
 coef.lam <- function(j, ploss, b.lam.init, delta.s, y.s, tau.j, lam.select, D, B){
   nlm(ploss,p=as.vector(b.lam.init),y=y.s[delta.s==1],tau=tau.j[j],lambda=lam.select,D=D,B=B[delta.s==1,])$est
 }	
 
 






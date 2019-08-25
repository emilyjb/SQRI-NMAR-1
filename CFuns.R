


Cjfun <- function(j, yimp1, phi22, wyimp1, gimp1, dgyimp1, E2J, BS){
	cyij <- exp(phi22*yimp1[,j])*dgyimp1[,j] + gimp1[,j]*exp(phi22*yimp1[,j])*phi22
	ctildeyij <- cyij/apply(exp(phi22*yimp1),1,sum) - E2J*wyimp1[,j]*phi22
	apply(BS*ctildeyij, 2, sum)/sum(groupind[,1])
}

CovFun <- function(wyimp1, yimp1, gimp1,n){
	sum(apply(wyimp1*yimp1*gimp1, 1, sum) - apply(wyimp1*yimp1, 1, sum) *apply(wyimp1*gimp1, 1,sum))/n
}


 

 


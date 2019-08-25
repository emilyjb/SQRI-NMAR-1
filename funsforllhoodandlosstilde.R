

p2fun <- function(y, phivec, x){
	num2 <- exp(phivec[1] + phivec[2]*x + phivec[3]*y)
	num3 <- exp(phivec[4] + phivec[5]*x + phivec[6]*y)
	den <-  1 + num2 + num3
	num2/den
}

p2funsimp2 <- function(y, phivec, x){
	num2 <- exp(phivec[1] + phivec[2]*x + phivec[3]*y)
	#num3 <- exp(phivec[4] + phivec[5]*x + phivec[6]*y)
	den <-  1 + num2  
	num2/den
}

p3fun <- function(x, phivec, y){
	num2 <- exp(phivec[1] + phivec[2]*x + phivec[3]*y)
	num3 <- exp(phivec[4] + phivec[5]*x + phivec[6]*y)
	den <-  1 +   num3
	num3/den
}

p3funsimp3 <- function(x, phivec, y){
	#num2 <- exp(phivec[1] + phivec[2]*x + phivec[3]*y)
	num3 <- exp(phivec[1] + phivec[2]*x + phivec[3]*y)
	den <-  1 +   num3
	num3/den
}

p1fun <-  function(x, phivec, y){
	num2 <- exp(phivec[1] + phivec[2]*x + phivec[3]*y)
	num3 <- exp(phivec[4] + phivec[5]*x + phivec[6]*y)
	den <-  1 + num2 + num3
	1/den
}

p1funsimp2 <-  function(x, phivec, y){
	num2 <- exp(phivec[1] + phivec[2]*x + phivec[3]*y)
	num3 <- exp(phivec[4] + phivec[5]*x + phivec[6]*y)
	den <-  1 + num2  
	1/den
}


lphinew <- function(phivec,  phicur, x, y, groupind, ximp1, yimp1){
	term1 <-p1fun( x[groupind[,1]==1], phivec, y[groupind[,1]==1])
	term2num <- apply(yimp1, 2, p2fun, phivec, x[groupind[,2]==1])
	term2nummod <- exp(phicur[3]*yimp1)*log(term2num)
	term2 <- apply(term2nummod, 1, mean)/apply(exp(phicur[3]*yimp1), 1, mean)
	term3num <- apply(ximp1, 2, p3fun, phivec, y[groupind[,3]==1])
	term3nummod <- exp(phicur[5]*ximp1)*log(term3num)
	term3 <- apply(term3nummod, 1, mean)/apply(exp(phicur[5]*ximp1), 1, mean)
	-(sum(log(term1)) + sum(term2) + sum(term3))
}

lphinewsimp2 <- function(phivec,  phicur, x, y, groupind, ximp1, yimp1){
	term1 <- p1funsimp2( x[groupind[,1]==1], phivec, y[groupind[,1]==1])
	term2num <- apply(yimp1, 2, p2funsimp2, phivec, x[groupind[,2]==1])
	term2nummod <- exp(phicur[3]*yimp1)*log(term2num)
	term2 <- apply(term2nummod, 1, mean)/apply(exp(phicur[3]*yimp1), 1, mean)
	 
	-(sum(log(term1)) + sum(term2)  )
}


lphinewsimp3 <- function(phivec,  phicur, x, y, groupind, ximp1, yimp1){
	term1 <- p1funsimp2( x[groupind[,1]==1], phivec, y[groupind[,1]==1])
	term2num <- apply(ximp1, 2, p3funsimp3, phivec, y[groupind[,3]==1])
	term2nummod <- exp(phicur[2]*ximp1)*log(term2num)
	term2 <- apply(term2nummod, 1, mean)/apply(exp(phicur[2]*ximp1), 1, mean)
	 
	-(sum(log(term1)) + sum(term2)  )
}



pitildefunx <- function(phi, x2,   beta, sigma2){
	num <- exp(phi[1] + phi[2]*x2 + phi[3]*beta[1] + beta[2]*phi[3]*x2 + beta[3]*phi[3]*x2^3 - phi[3]^2*sigma2/2)
	den <- 1+ num
	num/den
}


losstilde <- function(phi, groupind, x, ys1, beta, sigma2){
	#term1 <- pifun1(phi, x[groupind==1], ys1)
	term1 <- pitildefunx(phi, x[groupind==1], beta, sigma2)
	term2 <- pitildefunx(phi, x[groupind==0], beta, sigma2)
	-(sum(log(term1)) + sum( log(1-term2)))
}

losstildequant2 <- function(phi, groupind, x, yimpall){
	num <- exp(phi[1] + phi[2]*x - log( apply( exp(-yimpall*phi[3]), 1, mean)))
	den <- 1 + num
	pitilde <- num/den
	-(sum(log(pitilde[groupind==0])) + sum(log(1-pitilde[groupind==1])))
}



losstildequant3 <- function(phi, groupind, y, ximpall){
	num <- exp(phi[1] + phi[3]*y - log( apply( exp(-ximpall*phi[2]), 1, mean)))
	den <- 1 + num
	pitilde <- num/den
	-(sum(log(pitilde[groupind==0])) + sum(log(1-pitilde[groupind==1])))
}


pitildeq12 <- function(phi, groupind, x, yimpall){
  num <- exp(phi[1] + phi[2]*x - log( apply( exp(-yimpall*phi[3]), 1, mean)))
  den <- 1 + num
  num/den
}

pitildeq13 <- function(phi, groupind, y, ximpall){
  num <- exp(phi[1] + phi[3]*y - log( apply( exp(-ximpall*phi[2]), 1, mean)))
  den <- 1 + num
  num/den
}








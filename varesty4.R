
###  New Si:
pi12A12 <- pitildeq12(-phicur[1:3],  groupind = groupind[,2][groupind[,3]==0], x = x[groupind[,3]==0], yimpall = yimp1all) 
wyimp1all <- exp(-yimp1all*phivec2$par[3])/apply( exp(-yimp1all*phivec2$par[3]), 1, sum)
E0yall <- apply(wyimp1all*yimp1all, 1, sum)
Z2 <- -cbind(1, x[groupind[,3]==0], E0yall)
Siall <- (groupind[,1][groupind[,3]==0] - pi12A12)*Z2

bdx <- bw.ucv(x[groupind[,1]==1])
bdy <- bw.ucv(y[groupind[,1]==1])

Y1 <- kronecker(rep(1,sum(groupind[,1])), matrix(y[groupind[,1]==1], nrow = 1))
X1 <- kronecker(rep(1,sum(groupind[,1])), matrix(x[groupind[,1]==1], nrow = 1))
devxj <- X1 - x[groupind[,1]==1]
fhatdenj <- apply( exp(-devxj^2/2/bdx^2), 1, sum)
e1J <- y[groupind[,1]==1] - BS.X[groupind[,1] == 1,]%*%coefygivenx
#  fhat1ygx <- diff(tauvecforvar)/t(apply(qy1forvar, 1, diff))
Hjlist <- vector("list", length(tauvec))
Llist <- vector("list", length(tauvec))
fhatJmat <- c()
DSQJ <- 0
DBar <- 0
for(j in 1:length(tauvec)){


 qj <- (BS.X[groupind[,1] == 1,]%*%coefygivenx)[,j]
  devyj <- Y1 - qj
  qjall <- (BS.X[groupind[,3] == 0,]%*%coefygivenx)[,j]

  fhatnumj <- apply(exp(-devyj^2/2/bdy^2)*exp(-devxj^2/2/bdx^2), 1, sum)/bdy/sqrt(2*pi)
  
  fhat1ygxj <- fhatnumj/fhatdenj
  fhatJmat <- cbind(fhatJmat, fhat1ygxj)
  Hjlist[[j]] <-   1/sum(groupind[,1])*t(BS.X[groupind[,1]==1,])%*%diag( fhat1ygxj)%*%BS.X[groupind[,1]==1,] + lam.select/sum(groupind[,1])*t(D)%*%D
  
  psie1J <- tauvec[j] - ifelse( e1J[,j] < 0, 1, 0) 
  Llist[[j]] <-    1/sum(groupind[,1])*t(solve(Hjlist[[j]])%*%t(BS.X[groupind[,1]==1,]) )*psie1J 

  dsqjj <- wyimp1all[,j]*( 1 +  qjall*phicur[3] - phicur[3]*E0yall)  
  dsj1 <-  -pi12A12*(1-pi12A12)*wyimp1all[,j]*phicur[3]
  dsj2 <- dsj1*x[groupind[,3]==0]
  dsj3 <- dsqjj - pi12A12*dsqjj + E0yall*(-pi12A12*(1-pi12A12))*wyimp1all[,j]*phicur[3]
  
  dsj3star <- E0yall*(-pi12A12*(1-pi12A12))*wyimp1all[,j]*phicur[3]

  dsj <- cbind(dsj1, dsj2, dsj3star)
  DSQJ <- DSQJ +  dsj*as.vector(BS.X[groupind[,3]==0,]%*%apply(Llist[[j]], 2, sum)  )
  
  d1Bar <- apply(dsj1*BS.X[groupind[,3]==0,] , 2, sum)/nrow(groupind)
  d2Bar <- apply(dsj2*BS.X[groupind[,3]==0,] , 2, sum)/nrow(groupind)
  d3Bar <- apply(dsj3star*BS.X[groupind[,3]==0,] , 2, sum)/nrow(groupind)

  DBarj <-  rbind(as.vector(Llist[[j]]%*%d1Bar), as.vector(Llist[[j]]%*%d2Bar), as.vector(Llist[[j]]%*%d3Bar)) 

  DBar  <- DBar + DBarj

}

DBar <- DBar*nrow(groupind)

#Iphi <- phivec2$hessian
Iphi <- t(Z2)%*%diag(pi12A12*(1-pi12A12))%*%Z2# + diag(c(0,0,sum(pi12A12*apply(wyimp1all*yimp1all^2,1,sum) - apply(wyimp1all*yimp1all,1,sum)^2)))

 #solve(Iphi)%*%(t(Siall  + DSQJ)%*%(Siall  + DSQJ))%*%solve(Iphi)



Vy2 <- sum((apply(yimp1^2*wyimp1, 1, sum) - apply(yimp1 *wyimp1, 1, sum)^2))/n

#WB20 <- matrix(0, sum(groupind[,2]), length(tauvec))

WB2  <-  wyimp1*( 1 +  yimp1*phicur[3] - phicur[3]*apply(wyimp1*yimp1,1,sum)) 
WBSum <- t(t(BS.X[groupind[,2]==1,])%*%WB2)

T32 <-apply(WB2[,1]*BS.X[groupind[,2]==1,],2,sum)%*%t(Llist[[1]])
for(j in 2:length(tauvec)){
  WDj <- apply(WB2[,j]*BS.X[groupind[,2]==1,],2,sum)/n
  T32 <- T32 + WDj%*%t(Llist[[j]]*n)
}
#T32[1,]


#r <- rep(0,  sum(groupind[,1] + groupind[,2]))
groupind1A12 <- groupind[,1][groupind[,3]==0]
#r[groupind1A12==1] <- y[groupind[,3]==0][groupind[,1]==1]  -  meanyhatimp + as.vector(Vy2*matrix(c(0,0,1),nrow=1)%*%solve(Iphi)%*%(Siall  + DSQJ) + T32[1,]

phidel <- t(solve(Iphi/n)%*%t(Siall  + DSQJ))

rphi1 <-  solve(Iphi/n)%*%(t(Siall[groupind1A12==1,]) + DBar)
rphi2 <-  solve(Iphi/n)%*%(t(Siall[groupind1A12==0,])  )
phidel <- cbind(rphi1, rphi2)

r1 <-  y[ groupind[,1]==1] - meanyhatimp1  +  as.vector(Vy2*matrix(c(0,0,1),nrow=1)%*%(solve(Iphi/n))%*%(t(Siall[groupind1A12 ==1,] ) + DBar))  + T32[1,]  
r2 <- apply(wyimp1*yimp1, 1,sum) - meanyhatimp1  +  as.vector(Vy2*matrix(c(0,0,1),nrow=1)%*%(solve(Iphi/n))%*%t(Siall[groupind1A12 ==0,] ))  
r3 <- y[ groupind[,3]==1] - meanyhatimp1 



r <- c(r1, r2, r3)

vhaty1 <- var(r)/dim(groupind)[1]
vhaty1s <-c(vhaty1s, var(r)/dim(groupind)[1])


vhatphi <- diag( phidel%*%t(phidel))/5000^2

vhatphis <- rbind(vhatphis, vhatphi)



vhatbetahattaucheck <- diag(t(Llist[[10]])%*%Llist[[10]])
vhatbetahattauchecks <- rbind(vhatbetahattauchecks, vhatbetahattaucheck)

betahatygxs <- rbind(betahatygxs, coefygivenx[,10])

xydat <- data.frame(x,y)
xydat[,2][groupind[,2]==1] <- NA
xydat[,1][groupind[,3]==1] <- NA






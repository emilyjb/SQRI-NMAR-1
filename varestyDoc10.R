
#####  Note: phivec2 and phivec3 contain probabilities for delta1, while phicur contains probabilities for delta2 and delta3:

###  Compute S_i and Ti:
pi12A12 <- pitildeq12(-phicur[1:3],  groupind = groupind[,2][groupind[,3]==0], x = x[groupind[,3]==0], yimpall = yimp1all) 
pi13A13 <- pitildeq13(-phicur[4:6],  groupind = groupind[,3][groupind[,2]==0], y = y[groupind[,2]==0], ximpall = ximp1all) 

wyimp1all <- exp(-yimp1all*phivec2$par[3])/apply( exp(-yimp1all*phivec2$par[3]), 1, sum)
E0yall <- apply(wyimp1all*yimp1all, 1, sum)
Z2 <- -cbind(1, x[groupind[,3]==0], E0yall)
Siall <- (groupind[,1][groupind[,3]==0] - pi12A12)*Z2

wximp1all <- exp(-ximp1all*phivec3$par[3])/apply( exp(-ximp1all*phivec3$par[3]), 1, sum)
E0xall <- apply(wximp1all*ximp1all, 1, sum)
Z3 <- -cbind(1,  E0xall, y[groupind[,2]==0])
Tiall <- (groupind[,1][groupind[,2]==0] - pi13A13)*Z3

###  Compute beta-hat component:

bdx <- bw.ucv(x[groupind[,1]==1])
bdy <- bw.ucv(y[groupind[,1]==1])

Y1 <- kronecker(rep(1,sum(groupind[,1])), matrix(y[groupind[,1]==1], nrow = 1))
X1 <- kronecker(rep(1,sum(groupind[,1])), matrix(x[groupind[,1]==1], nrow = 1))
devx1j <- X1 - x[groupind[,1]==1]
fhatdenjygx <- apply( exp(-devx1j^2/2/bdx^2), 1, sum)

devy1j <- Y1 - y[groupind[,1]==1]
fhatdenjxgy <- apply( exp(-devy1j^2/2/bdy^2), 1, sum)

Llist <- vector("list", length(tauvec)); Mlist <- vector("list", length(tauvec))
dqLSum <- 0; fqMSum <- 0

CyEy <- CxEx <- CyEy2 <- CxEx2 <- CyExy <- CxExy  <- CyEydel <- CxExdel <- 0

for(j in 1:length(tauvec)){

      qjy <- (BS.X[groupind[,1] == 1,]%*%coefygivenx)[,j]
      devyj <- Y1 - qjy
      qjx <- (BS.Y[groupind[,1] == 1,]%*%coefxgiveny)[,j]
      devxj <- X1 - qjx

	fhatnumjy <- apply(exp(-devyj^2/2/bdy^2)*exp(-devx1j^2/2/bdx^2), 1, sum)/bdy/sqrt(2*pi)
      fhat1ygxj <- fhatnumjy/fhatdenjygx	
      
	Hyj <- 1/sum(groupind[,1])*t(BS.X[groupind[,1]==1,])%*%diag( fhat1ygxj)%*%BS.X[groupind[,1]==1,] + lam.select/sum(groupind[,1])*t(D)%*%D

	fhatnumjx <- apply(exp(-devxj^2/2/bdx^2)*exp(-devy1j^2/2/bdy^2), 1, sum)/bdx/sqrt(2*pi)
      fhat1xgyj <- fhatnumjx/fhatdenjxgy
	
	Hxj <- 1/sum(groupind[,1])*t(BS.Y[groupind[,1]==1,])%*%diag( fhat1xgyj)%*%BS.Y[groupind[,1]==1,] + lam.select/sum(groupind[,1])*t(D)%*%D

      psitau1ygx <- tauvec[j] - ifelse( y[groupind[,1]==1] - BS.X[groupind[,1]==1,]%*%coefygivenx[,j] < 0, 1, 0) 
	psitau1xgy <-  tauvec[j] - ifelse( x[groupind[,1]==1] - BS.Y[groupind[,1]==1,]%*%coefxgiveny[,j]< 0, 1, 0) 
	
	Llist[[j]] <- (BS.X[groupind[,1]==1,]%*% solve(Hyj))*as.vector(psitau1ygx)
	Mlist[[j]] <- (BS.Y[groupind[,1]==1,]%*%t( solve(Hxj)))*as.vector(psitau1xgy)

	dqj <- pi12A12*(1-pi12A12)*wyimp1all[,j]*phicur[3]*Z2
	fqj <- pi13A13*(1-pi13A13)*wximp1all[,j]*phicur[6]*Z3

	dqjB1 <- dqj[,1]*BS.X[groupind[,1] == 1 | groupind[,2] ==1,]
	dqjB2 <- dqj[,2]*BS.X[groupind[,1] == 1 | groupind[,2] ==1,]
	dqjB3 <- dqj[,3]*BS.X[groupind[,1] == 1 | groupind[,2] ==1,]

	fqjB1 <- fqj[,1]*BS.Y[groupind[,1] == 1 | groupind[,3] ==1,]
	fqjB2 <- fqj[,2]*BS.Y[groupind[,1] == 1 | groupind[,3] ==1,]
	fqjB3 <- fqj[,3]*BS.Y[groupind[,1] == 1 | groupind[,3] ==1,]
		
	dqBBar <- rbind(apply(dqjB1,2,sum),apply(dqjB2,2,sum), apply(dqjB3,2,sum))/sum(groupind[,1])
	fqBBar <- rbind(apply(fqjB1,2,sum),apply(fqjB2,2,sum), apply(fqjB3,2,sum))/sum(groupind[,1])
	
	dqLSum <- dqLSum + dqBBar%*%t(Llist[[j]])
	fqMSum <- fqMSum + fqBBar%*%t(Mlist[[j]])

######  Input the C-terms for estimating the variance of the mean:

# Cyj for EY:
CyjEy <- Cjfun(j, yimp1, phicur[3], wyimp1, yimp1, matrix(1,nrow = nrow(yimp1), ncol = ncol(yimp1)), apply(wyimp1*yimp1, 1, sum), BS.X[groupind[,2]==1,])
CyEy <- CyEy +  Llist[[j]]%*%CyjEy
# Cxj for EX:
CxjEx <- Cjfun(j, ximp1, phicur[6], wximp1, ximp1, matrix(1,nrow = nrow(ximp1), ncol = ncol(ximp1)), apply(wximp1*ximp1, 1, sum), BS.Y[groupind[,3]==1,])
CxEx <-  CxEx + Mlist[[j]]%*%CxjEx
# Cyj for EY2:
CyjEy2 <- Cjfun(j, yimp1, phicur[3], wyimp1,  yimp1^2, 2*yimp1, apply(wyimp1*yimp1^2, 1, sum), BS.X[groupind[,2]==1,])
CyEy2 <-CyEy2 +  Llist[[j]]%*%CyjEy2
# Cxj for Ex2:
CxjEx2 <- Cjfun(j, ximp1, phicur[6], wximp1,  ximp1^2, 2*ximp1, apply(wximp1*ximp1^2, 1, sum), BS.Y[groupind[,3]==1,])
CxEx2 <- CxEx2 + Mlist[[j]]%*%CxjEx2

# Cyj for Exy:
CyjExy <- Cjfun(j, yimp1, phicur[3], wyimp1,  yimp1*x[groupind[,2]==1] ,  kronecker(x[groupind[,2]==1], matrix(1, ncol = length(tauvec))) , apply(wyimp1*x[groupind[,2]==1]*yimp1, 1, sum), BS.X[groupind[,2]==1,])
CyExy <- CyExy + Llist[[j]]%*%CyjExy
# Cxj for Exy:
CxjExy <- Cjfun(j, ximp1, phicur[6], wximp1,  ximp1*y[groupind[,3]==1] ,  kronecker(y[groupind[,3]==1], matrix(1, ncol = length(tauvec))) , apply(wximp1*y[groupind[,3]==1]*ximp1, 1, sum), BS.Y[groupind[,3]==1,])
CxExy <- CxExy + Mlist[[j]]%*%CxjExy

# Cyj for Eydel:
CyjEydel <- Cjfun(j, yimp1, phicur[3], wyimp1,  yimp1*domainInd[groupind[,2]==1] ,  kronecker(domainInd[groupind[,2]==1], matrix(1, ncol = length(tauvec))) , apply(wyimp1*domainInd[groupind[,2]==1]*yimp1, 1, sum), BS.X[groupind[,2]==1,])
CyEydel <- CyEydel +  Llist[[j]]%*%CyjEydel
# Cxj for Exdel:
CxjExdel <- Cjfun(j, ximp1, phicur[6], wximp1,  ximp1*domainInd[groupind[,3]==1] ,  kronecker(domainInd[groupind[,3]==1], matrix(1, ncol = length(tauvec))) , apply(wximp1*domainInd[groupind[,3]==1]*ximp1, 1, sum), BS.Y[groupind[,3]==1,])
CxExdel <- CxExdel + Mlist[[j]]%*%CxjExdel


}

Inphi2 <- t(Z2)%*%diag(pi12A12*(1-pi12A12))%*%Z2/n
Inphi3 <- t(Z3)%*%diag(pi13A13*(1-pi13A13))%*%Z3/n

Uphi2g2 <-  Siall[groupind[,2][groupind[,3]==0]==1,] 
Uphi2g1 <-  Siall[groupind[,1][groupind[,3]==0]==1,] + t(dqLSum)
Uphi2 <- matrix(0, nrow = n, ncol = 3)
Uphi2[groupind[,1] ==1,] <- Uphi2g1	
Uphi2[groupind[,2] == 1,] <- Uphi2g2

Tphi3g3 <-  Tiall[groupind[,3][groupind[,2]==0]==1,] 
Tphi3g1 <-  Tiall[groupind[,1][groupind[,2]==0]==1,] + t(fqMSum)
Tphi3 <- matrix(0, nrow = n, ncol = 3)
Tphi3[groupind[,1] ==1,] <- Tphi3g1	
Tphi3[groupind[,3] == 1,] <- Tphi3g3

Vhatphi2 <- solve(Inphi2)%*%t(Uphi2)%*%Uphi2%*%solve(Inphi2)/n^2
Vhatphi3 <- solve(Inphi3)%*%t(Tphi3)%*%Tphi3%*%solve(Inphi3)/n^2

rEy <- rep(0,n)
rEy1 <- y[groupind[,1]==1] - meanyhatimp1  + CovFun(wyimp1, yimp1, yimp1,n)*(Uphi2g1%*%solve(Inphi2))[,3]   +  CyEy 
rEy2 <- apply(wyimp1*yimp1,1,sum) - meanyhatimp1  +  CovFun(wyimp1, yimp1, yimp1,n)*(Uphi2g2%*%solve(Inphi2))[,3]
rEy3 <- y[groupind[,3] ==1] - meanyhatimp1 
rEy[groupind[,1] ==1] <- rEy1
rEy[groupind[,2] ==1] <- rEy2
rEy[groupind[,3] ==1] <- rEy3
vhatEy <- var(rEy)/n

rEy2 <- rep(0,n)
rEy21 <- (y^2)[groupind[,1]==1] - Ey2  + CovFun(wyimp1, yimp1, yimp1^2,n)*(Uphi2g1%*%solve(Inphi2))[,3]   +  CyEy2 
rEy22 <- apply(wyimp1*yimp1^2,1,sum) - Ey2  +  CovFun(wyimp1, yimp1, yimp1^2,n)*(Uphi2g2%*%solve(Inphi2))[,3]
rEy23 <- (y^2)[groupind[,3] ==1] - Ey2 
rEy2[groupind[,1] ==1] <- rEy21
rEy2[groupind[,2] ==1] <- rEy22
rEy2[groupind[,3] ==1] <- rEy23
vhatEy2 <- var(rEy2)/n

rEx <- rep(0,n)
rEx1 <- x[groupind[,1]==1] - meanxhatimp1  + CovFun(wximp1, ximp1, ximp1,n)*(Tphi3g1%*%solve(Inphi3))[,3]   +  CxEx 
rEx2 <-  x[groupind[,2]==1] - meanxhatimp1  
rEx3 <-apply(wximp1*ximp1,1,sum) - meanxhatimp1 +  CovFun(wximp1, ximp1, ximp1,n)*(Tphi3g3%*%solve(Inphi3))[,3]
rEx[groupind[,1] ==1] <- rEx1
rEx[groupind[,2] ==1] <- rEx2
rEx[groupind[,3] ==1] <- rEx3
vhatEx <- var(rEx)/n

rEx2 <- rep(0,n)
rEx21 <- (x^2)[groupind[,1]==1] - Ex2  + CovFun(wximp1, ximp1, ximp1^2,n)*(Tphi3g1%*%solve(Inphi3))[,3]   +  CxEx2 
rEx22 <-  (x^2)[groupind[,2]==1] - Ex2
rEx23 <-apply(wximp1*(ximp1^2),1,sum) -  Ex2 +  CovFun(wximp1, ximp1, ximp1^2,n)*(Tphi3g3%*%solve(Inphi3))[,3]
rEx2[groupind[,1] ==1] <- rEx21
rEx2[groupind[,2] ==1] <- rEx22
rEx2[groupind[,3] ==1] <- rEx23
vhatEx2 <- var(rEx2)/n

rExy <- rep(0,n)
rExy1 <-  (x*y)[groupind[,1]==1] - Exy  +  CovFun(wyimp1, yimp1, yimp1*(x[groupind[,2]==1]),n)*(Uphi2g1%*%solve(Inphi2))[,3]   +  CovFun(wximp1, ximp1, ximp1*y[groupind[,3]==1],n)*(Tphi3g1%*%solve(Inphi3))[,3]   +  CxExy + CyExy
rExy2 <-  apply(wyimp1*yimp1*x[groupind[,2]==1], 1, sum) - Exy + CovFun(wyimp1, yimp1, yimp1*(x[groupind[,2]==1]),n)*(Uphi2g2%*%solve(Inphi2))[,3]
rExy3 <-  apply(wximp1*ximp1*y[groupind[,3]==1],1,sum) -  Exy +  CovFun(wximp1, ximp1, ximp1*y[groupind[,3]==1],n)*(Tphi3g3%*%solve(Inphi3))[,3]
rExy[groupind[,1] ==1] <- rExy1
rExy[groupind[,2] ==1] <- rExy2
rExy[groupind[,3] ==1] <- rExy3
vhatExy <- var(rExy)/n

rydel <- rep(0,n)
rydel1 <- (y*domainInd)[groupind[,1]==1] - DomMeanY +  CovFun(wyimp1, yimp1, yimp1*(domainInd[groupind[,2]==1]),n)*(Uphi2g1%*%solve(Inphi2))[,3] + CyEydel
rydel2 <- apply(wyimp1*yimp1*x[groupind[,2]==1], 1, sum) - DomMeanY +  CovFun(wyimp1, yimp1, yimp1*(domainInd[groupind[,2]==1]),n)*(Uphi2g2%*%solve(Inphi2))[,3] 
rydel3 <- (y*domainInd)[groupind[,3]==1] - DomMeanY
rydel[groupind[,1]==1] <- rydel1
rydel[groupind[,2]==1] <- rydel2
rydel[groupind[,3]==1] <- rydel3

rydel2 <- 1/mean(domainInd)*(rydel - DomMeanY*domainInd)
vhatEydel <- var(rydel2)/n

rxdel <- rep(0,n)
rxdel1 <- (x*domainInd)[groupind[,1]==1] - DomMeanX +  CovFun(wximp1, ximp1, ximp1*(domainInd[groupind[,3]==1]),n)*(Tphi3g1%*%solve(Inphi3))[,3] + CxExdel
rxdel2 <- (x*domainInd)[groupind[,2]==1] - DomMeanX  
rxdel3 <- apply(wximp1*ximp1*y[groupind[,3]==1], 1, sum) - DomMeanX +  CovFun(wximp1, ximp1, ximp1*(domainInd[groupind[,3]==1]),n)*(Tphi3g3%*%solve(Inphi3))[,3] 
rxdel[groupind[,1]==1] <- rxdel1
rxdel[groupind[,2]==1] <- rxdel2
rxdel[groupind[,3]==1] <- rxdel3

rxdel3 <- 1/mean(domainInd)*(rxdel - DomMeanX*domainInd)
vhatExdel <- var(rxdel3)/n
 
if(vhatEx > 1 | vhatEy  > 1){ save.image("BadSample.Rdata")}

vhatimps <- rbind(vhatimps, c(vhatEy, vhatEx, vhatEy2, vhatEx2, vhatExy, vhatEydel, vhatExdel))

vhatphis <- rbind(vhatphis, c(diag(Vhatphi2), diag(Vhatphi3)))

















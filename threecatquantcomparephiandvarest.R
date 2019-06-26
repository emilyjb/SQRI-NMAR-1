############  Consider estimation of phi for full response:
rm(list = ls(all = TRUE))
library("quantreg")
library("splines")

source("bsplineoptfuns1.R")
source("computescoretermsmissingcovs1.R")
source("funsforllhoodandlosstilde.R")

phi20 <- -0.8
phi21 <- 0.15
phi22 <- 0.15
phi30 <- -0.7
phi31 <- 0.1
phi32 <- 0.1

phitrue <- c(phi20, phi21, phi22, phi30, phi31, phi32)

n <- 5000

phihats <- c()
itersim <- 0

p <- 3
K.n <- 5
####knots.select <- ((-p + 1):(K.n + p))/K.n
knots.select <- (1:(K.n-1))/K.n
J <- 100

meanys <- c()
meanxs <- c()

meanyhatimps <- c()
meanxhatimps <- c()

meanyhatimp1s <- c()
meanxhatimp1s <- c()

meanyccs <-c()
meanxccs <- c()

phihats <- c()
R <- 3
vhaty1s <- c()
vhatbetahattauchecks <- c()
vhatphis <- c()

vhatbetahattauchecks <- c()

betahatygxs <- c()

repeat{
  
  itersim <- itersim + 1
  
    x <- runif(n, -1,1)
  # y <- -1 + x + 2*(x - mean(x))^2+ -3*(x-mean(x))^3 +  rnorm(n, 0, 0.9)
#    y <-   x + 2*x^3 + 3*x^5 +  rnorm(n, 0, sd = 0.75*(1+2*abs(x)))

   y <- exp(2*x) +  rnorm(n, 0, sd = 0.5)


  quantx <- quantile(x, probs = seq(0, 1, by = 0.2))
  quanty <- quantile(y, probs = seq(0, 1, by = 0.2))

  BS.X <- model.matrix(x~bs(x,   degree =15, knots = quantx))
  BS.Y <- model.matrix(y~bs(y,   degree = 15, knots = quanty))
  
  denpi <-  1 + exp(phi20 + phi21*x + phi22*y ) + exp(phi30 + phi31*x + phi32*y)
  numpi1 <- 1
  numpi2 <- exp(phi20 + phi21*x + phi22*y ) 
  numpi3 <- exp(phi30 + phi31*x + phi32*y)
  
  pi1 <- numpi1/denpi
  pi2 <- numpi2/denpi
  pi3 <- numpi3/denpi
  
  pimat <- cbind(pi1, pi2, pi3)
  
  groupind <- t(apply(pimat, 1, function(x){rmultinom(1, 1, prob = x)}))
  #### head(t(groupind))
  
  #### Initial phi2
  groupind13 <- groupind[,1] + groupind[,3]
  groupind12 <- groupind[,1] + groupind[,2]
  groupind2 <- groupind[,2]
  groupind3 <- groupind[,3]
  
  
  groupind13 <- groupind[,1] + groupind[,3]
  groupind2 <- groupind[,2]
  
  ##############################  Initial phi value:
  
  D <- diag(ncol(BS.X))
  for (k in 1:2) D <- diff(D)
  Omega <- t(D)%*%D
  lams.seq <- seq(0,0.5, length = 500)
 # lam.select <- smooth.spline(x[groupind[,1]==1], y[groupind[,1]==1],nknots =length(quantx), spar = NULL)$lambda
 # gcvlams <- sapply(lams.seq, gcv.lam, B, Omega/2, y, delta.s)
 # lam.select <- lams.seq[which.min(gcvlams)]
 # lam.select <- min(sapply(tau.j, cobs.tau.nk, 26, 2))
  
  lam.select <-  2
  
  bx.lam.init <- solve(t(BS.X[groupind[,1]==1,])%*%BS.X[groupind[,1]==1,] + lam.select*Omega)%*%t(BS.X[groupind[,1]==1,])%*%(y[groupind[,1]==1])
  by.lam.init <- solve(t(BS.Y[groupind[,1]==1,])%*%BS.Y[groupind[,1]==1,] + lam.select*Omega)%*%t(BS.Y[groupind[,1]==1,])%*%(x[groupind[,1]==1])
  
  
  
  #tauvecforvar <- sort(runif(101))
 # tauvec <- (tauvecforvar[2:101] + tauvecforvar[1:100])/2
#   tauvec <- sort(runif(100))
  tauvec1 <- runif(1, 0, 0.01)
 
  tauvec <- tauvec1 + seq(0.01, 0.99, by = 0.01)
  coefygivenx <- sapply(1:length(tauvec), coef.lam,  ploss, bx.lam.init, groupind[,1], y, tauvec, lam.select, D, BS.X)
  coefxgiveny <- sapply(1:length(tauvec), coef.lam,  ploss, by.lam.init, groupind[,1], x, tauvec, lam.select, D, BS.Y)

  yimp1all <-    BS.X[groupind[,2] == 1 | groupind[,1]==1 ,]%*%coefygivenx
  ximp1all <-    BS.Y[groupind[,3] == 1 | groupind[,1]==1, ]%*%coefxgiveny

  yimp1 <- yimp1all[groupind[groupind[,2] == 1 | groupind[,1]==1 ,][,2]==1,]
  ximp1 <- ximp1all[groupind[groupind[,3] == 1 | groupind[,1]==1 ,][,3]==1,]

 
  yimp <- kronecker(rep(1, sum(groupind[,2]==1)), t(y[groupind[,1]==1]))
  ximp <- kronecker(rep(1, sum(groupind[,3]==1)), t(x[groupind[,1]==1]))
 
  yimphat1 <- apply(yimp1all, 1, mean)
  ximphat1 <- apply(ximp1all, 1, mean)
 
  phihat2cur <- glm(cbind(groupind[,2], 1- groupind[,2])[groupind[,1]==1 | groupind[,2]==1,]~  x[groupind[,3]==0] + yimphat1, family = binomial(link = "logit"))$coef
  phihat3cur <- glm(cbind(groupind[,3], 1- groupind[,3])[groupind[,1]==1 | groupind[,3]==1,]~  ximphat1 + y[groupind[,2]==0], family = binomial(link = "logit"))$coef
  
 
#  plot(x, y)
 # lines(x[order(x)], (BS.X%*%bx.lam.init)[order(x)],col = 2, lwd = 2)
# yimp1allforplot <- BS.X%*%coefygivenx
#for(i in 1:length(tauvec)){
# lines(x[order(x)], yimp1allforplot[order(x),i])
#}

 # plot(y, x)
 # lines(y[order(y)], (BS.Y%*%by.lam.init)[order(y)],col = 2, lwd = 2)
 #ximp1allforplot <- BS.Y%*%coefxgiveny
#for(i in 1:length(tauvec)){
 #lines(y[order(y)], ximp1allforplot[order(y),i])
#}


iterphi <- 0
phicur <- phivec <- c(phihat2cur, phihat3cur)

 

iterphi <- 0
phicur <- phivec <- c(phihat2cur, phihat3cur)

repeat{

	iterphi <- iterphi + 1

	phiupdateopt2 <- optim(phivec[1:3] ,lphinewsimp2   , phicur=phicur[1:3]  , x=x, y=y, groupind=groupind, ximp1=ximp1, yimp1 = yimp1)

################
	phiupdateopt3 <- optim(phivec[4:6] ,lphinewsimp3   , phicur=phicur[4:6]  , x=x, y=y, groupind=groupind, ximp1=ximp1, yimp1 = yimp1)
################

	
	phivec <- phicur <-c( phiupdateopt2$par,  phiupdateopt3$par)


	if(iterphi == 3){break}

}

#####  Start iteration with old initial value:

  phivec2 <- optim(- phihat2cur, losstildequant2, groupind = groupind[,2][groupind[,3]==0], x = x[groupind[,3]==0], yimpall = yimp1all, hessian = TRUE) 
  
  phivec3 <- optim(- phihat3cur, losstildequant3, groupind = groupind[,3][groupind[,2]==0], y = y[groupind[,2]==0], ximpall = ximp1all) 
  
  phihats <- rbind(phihats, c(phivec, -phivec2$par, -phivec3$par))
    
  print(phivec)  
  print(paste(itersim))


    
  wyimp <- exp(yimp1*phicur[3])/apply( exp(yimp1*phicur[3]), 1, sum)
  wximp <- exp(ximp1*phicur[5])/apply( exp(ximp1*phicur[5]), 1, sum)
 

  wyimp1 <- exp(-yimp1*phivec2$par[3])/apply( exp(-yimp1*phivec2$par[3]), 1, sum)
  wximp1 <- exp(-ximp1*phivec3$par[3])/apply( exp(-ximp1*phivec3$par[3]), 1, sum)
  meanyhatimp <- mean(c(y[groupind[,2]==0], apply(yimp1*wyimp, 1, sum)  )  )
  meanxhatimp <- mean(c(x[groupind[,3]==0], apply(ximp1*wximp, 1, sum)  )  )
 
    meanyhatimps <- c(meanyhatimps, meanyhatimp)
  meanxhatimps <- c(meanxhatimps, meanxhatimp)
 
  meanyhatimp1 <- mean(c(y[groupind[,2]==0], apply(yimp1*wyimp1, 1, sum)  )  )
  meanxhatimp1 <- mean(c(x[groupind[,3]==0], apply(ximp1*wximp1, 1, sum)  )  )
  
  meanys <- c(meanys, mean(y))
  meanxs <- c(meanxs, mean(x))
  
   
  meanyhatimp1s <- c(meanyhatimp1s, meanyhatimp1)
  meanxhatimp1s <- c(meanxhatimp1s, meanxhatimp1)
  
  meanyccs <- c(meanyccs, mean(y[groupind[,2]==0]))
  meanxccs <- c(meanxccs, mean(x[groupind[,3]==0]))
  
  phicur <- c(-phivec2$par, -phivec3$par)
#############     source("varesty3.R")
  
  
  if(itersim == 100){break}
  
}


xtable(cbind(c(apply(phihats[1:100,],2, var)[c(7,8,9)], var(meanyhatimp1s[1:100])), c(apply(vhatphis[1:100,] , 2, mean),  mean(vhaty1s[1:100]))), digits = 4)


phit <- (apply(phihats, 2, mean)[1:6] - c(phi20, phi21, phi22, phi30, phi31, phi32))/apply(phihats, 2, sd)[1:6]*sqrt(itersim)
phimat1 <- cbind(c(phi20, phi21, phi22, phi30, phi31, phi32), apply(phihats, 2, mean)[1:6], apply(phihats, 2, sd)[1:6], phit)
my1 <- c( mean(meanys), mean(meanyhatimps), sd(meanyhatimps), mean(meanys - meanyhatimps)/sd(meanys - meanyhatimps)*sqrt(length(meanys)))
mx1 <- c( mean(meanxs), mean(meanxhatimps), sd(meanxhatimps), mean(meanxs - meanxhatimps)/sd(meanxs - meanxhatimps)*sqrt(length(meanxs)))

rbind(phimat1, my1, mx1)

phit <- (apply(phihats, 2, mean)[1:6] - c(phi20, phi21, phi22, phi30, phi31, phi32))/apply(phihats, 2, sd)[7:12]*sqrt(itersim)
phimat2 <- cbind(c(phi20, phi21, phi22, phi30, phi31, phi32), apply(phihats, 2, mean)[7:12], apply(phihats, 2, sd)[7:12], phit)
my2 <- c( mean(meanys), mean(meanyhatimp1s), sd(meanyhatimp1s), mean(meanys - meanyhatimp1s)/sd(meanys - meanyhatimp1s)*sqrt(length(meanys)))
mx2 <- c( mean(meanxs), mean(meanxhatimp1s), sd(meanxhatimp1s), mean(meanxs - meanxhatimp1s)/sd(meanxs - meanxhatimp1s)*sqrt(length(meanxs)))

model1out <- rbind(rbind(phimat1, my1, mx1),
rbind(phimat2, my2, mx2))

msecompare <- matrix(sqrt((model1out[,2]-model1out[,1])^2 + model1out[,3]^2), nrow = 8, byrow = FALSE)
biascompare <- matrix( model1out[,2]-model1out[,1] , nrow = 8, byrow = FALSE)
sdcompare <- matrix(  model1out[,3], nrow = 8, byrow = FALSE)

xtable(cbind(model1out[,1][1:8], biascompare, sdcompare, msecompare))



msecompare[,2]/msecompare[,1]


mean(meanys - meanyhatimp1s)/sd(meanys - meanyhatimps)*sqrt(length(meanys))
mean(meanxs - meanxhatimp1s)/sd(meanxs - meanxhatimps)*sqrt(length(meanxs))

mean(meanys - meanyccs)/sd(meanys - meanyccs)*sqrt(length(meanys))
mean(meanxs - meanxccs)/sd(meanxs - meanxccs)*sqrt(length(meanxs))



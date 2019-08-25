
#x2 <- x[groupind[,2]==1]


####  compute score equation term:

computeScoreterm <- function(phivec, y1, y3, x1, x2, yimp, ximp, wyimp, wximp){
 # wyimp <- exp(yimp*phivec[2])/apply(exp(yimp*phivec[2]), 1, sum)
 # wximp <- exp(ximp*phivec[4])/apply(exp(ximp*phivec[4]), 1, sum)
  sG2J <- sapply(1:ncol(yimp), computescoretermGroup2a , phivec, y1, y3, x1, x2, yimp, wyimp)
  sG2 <- apply(sG2J, 1, sum)
  sG3J <- sapply(1:ncol(ximp), computescoretermGroup3a , phivec, y1, y3, x1, x2, ximp, wximp)
  sG3 <- apply(sG3J, 1, sum)
  sG1 <- computescoretermGroup1a( phivec, y1, y3, x1, x2)
  
  Svec <- as.vector(sG1 + sG2 + sG3)
  Svec
  
  #sum(Svec^2)
  
}

computeQuadterm <- function(phivec, y1, y3, x1, x2, yimp, ximp, wyimp, wximp){
 # wyimp <- exp(yimp*phivec[2])/apply(exp(yimp*phivec[2]), 1, sum)
 # wximp <- exp(ximp*phivec[4])/apply(exp(ximp*phivec[4]), 1, sum)
  sG2J <- sapply(1:ncol(yimp), computescoretermGroup2a , phivec, y1, y3, x1, x2, yimp, wyimp)
  sG2 <- apply(sG2J, 1, sum)
  sG3J <- sapply(1:ncol(ximp), computescoretermGroup3a , phivec, y1, y3, x1, x2, ximp, wximp)
  sG3 <- apply(sG3J, 1, sum)
  sG1 <- computescoretermGroup1a( phivec, y1, y3, x1, x2)
  
  Svec <- as.vector(sG1 + sG2 + sG3)
  #Svec
  
  sum(Svec^2)
  
}

computescoretermGroup2bnorm <- function(phivec, y1, y3, x1, x2, yimp, ximp, wyimp, wximp){

  Z1 <- kronecker(cbind(1, yimp[,k]), c(0, 1, 0))
  Z2 <- kronecker(cbind(1, x2), c(0, 0, 1))
  Z <- cbind(Z1, Z2)
  fitZ <- as.vector(Z%*%phivec)
  nummat <- matrix(exp(fitZ), nrow = 3, byrow = FALSE)
  phatmatphi <- nummat/apply(nummat, 2, sum)
  ZDPG2 <- apply(sapply(1:sum(groupind[,2]==1), function(j){  t( Z[((j-1)*R + 1):(j*R),][-1,]  )%*%(wyimp[j,k]*(groupind[groupind[,2]==1,-1][j,] - phatmatphi[-1,j])) }), 1, sum)
  as.vector(ZDPG2)
  
}




computeQuadterm2 <- function(phivec, y1, y3, x1, x2, yimp, ximp, wyimp, wximp){
  # wyimp <- exp(yimp*phivec[2])/apply(exp(yimp*phivec[2]), 1, sum)
  # wximp <- exp(ximp*phivec[4])/apply(exp(ximp*phivec[4]), 1, sum)
  sG2J <- sapply(1:length(tauvec), computescoretermGroup2 , phivec, y1, y3, x1, x2, yimp, wyimp)
  sG2 <- sum(sG2J)
  sG3J <- sapply(1:length(tauvec), computescoretermGroup3 , phivec, y1, y3, x1, x2, ximp, wximp)
  sG3 <- sum(sG3J)
  sG1 <- computescoretermGroup1( phivec, y1, y3, x1, x2)
  
  Svec <- as.vector(sG1 + sG2 + sG3)
  Svec
  
  # sum(Svec^2)
  
}



computescoretermGroup2a <- function(k, phivec, y1, y3, x1, x2, yimp, wyimp){
 
  Z1 <- kronecker(cbind(1, yimp[,k]), c(0, 1, 0))
  Z2 <- kronecker(cbind(1, x2), c(0, 0, 1))
  Z <- cbind(Z1, Z2)
  fitZ <- as.vector(Z%*%phivec)
  nummat <- matrix(exp(fitZ), nrow = 3, byrow = FALSE)
  phatmatphi <- nummat/apply(nummat, 2, sum)
  ZDPG2 <- apply(sapply(1:sum(groupind[,2]==1), function(j){  t( Z[((j-1)*R + 1):(j*R),][-1,]  )%*%(wyimp[j,k]*(groupind[groupind[,2]==1,-1][j,] - phatmatphi[-1,j])) }), 1, sum)
  as.vector(ZDPG2)
  
}




computescoretermGroup3a <- function(k, phivec, y1, y3, x1, x2, ximp, wximp){
  
  Z1 <- kronecker(cbind(1, y3), c(0, 1, 0))
  Z2 <- kronecker(cbind(1, ximp[,k]), c(0, 0, 1))
  Z <- cbind(Z1, Z2)
  fitZ <- as.vector(Z%*%phivec)
  nummat <- matrix(exp(fitZ), nrow = 3, byrow = FALSE)
  phatmatphi <- nummat/apply(nummat, 2, sum)
  ZDPG3 <- apply(sapply(1:sum(groupind[,3]==1), function(j){  t( Z[((j-1)*R + 1):(j*R),][-1,]  )%*%(wximp[j,k]*(groupind[groupind[,3]==1,-1][j,] - phatmatphi[-1,j])) }), 1, sum)
  as.vector(ZDPG3)
  
}


computescoretermGroup1a <- function( phivec, y1, y3, x1, x2 ){
  
  Z1 <- kronecker(cbind(1, y1), c(0, 1, 0))
  Z2 <- kronecker(cbind(1, x1), c(0, 0, 1))
  Z <- cbind(Z1, Z2)
  fitZ <- as.vector(Z%*%phivec)
  nummat <- matrix(exp(fitZ), nrow = 3, byrow = FALSE)
  phatmatphi <- nummat/apply(nummat, 2, sum)
  ZDPG1 <- apply(sapply(1:sum(groupind[,1]==1), function(j){  t( Z[((j-1)*R + 1):(j*R),][-1,]  )%*%((groupind[groupind[,1]==1,-1][j,] - phatmatphi[-1,j])) }), 1, sum)
  as.vector(ZDPG1)
  
}



computescoretermGroup2 <- function(k, phivec, y1, y3, x1, x2, yimp, wyimp){
 
  Z1 <- kronecker(cbind(1, yimp[,k]), c(0, 1, 0))
  Z2 <- kronecker(cbind(1, x2), c(0, 0, 1))
  Z <- cbind(Z1, Z2)
  fitZ <- as.vector(Z%*%phivec)
  nummat <- matrix(exp(fitZ), nrow = 3, byrow = FALSE)
  phatmatphi <- nummat/apply(nummat, 2, sum)
  ZDPG2 <- apply(sapply(1:sum(groupind[,2]==1), function(j){  t( Z[((j-1)*R + 1):(j*R),][-1,]  )%*%(wyimp[j,k]*(groupind[groupind[,2]==1,-1][j,] - phatmatphi[-1,j])) }), 1, sum)
  sum(as.vector(ZDPG2)^2)
  
}




computescoretermGroup3 <- function(k, phivec, y1, y3, x1, x2, ximp, wximp){
  
  Z1 <- kronecker(cbind(1, y3), c(0, 1, 0))
  Z2 <- kronecker(cbind(1, ximp[,k]), c(0, 0, 1))
  Z <- cbind(Z1, Z2)
  fitZ <- as.vector(Z%*%phivec)
  nummat <- matrix(exp(fitZ), nrow = 3, byrow = FALSE)
  phatmatphi <- nummat/apply(nummat, 2, sum)
  ZDPG3 <- apply(sapply(1:sum(groupind[,3]==1), function(j){  t( Z[((j-1)*R + 1):(j*R),][-1,]  )%*%(wximp[j,k]*(groupind[groupind[,3]==1,-1][j,] - phatmatphi[-1,j])) }), 1, sum)
  sum(as.vector(ZDPG3)^2)
  
}


computescoretermGroup1 <- function( phivec, y1, y3, x1, x2 ){
  
  Z1 <- kronecker(cbind(1, y1), c(0, 1, 0))
  Z2 <- kronecker(cbind(1, x1), c(0, 0, 1))
  Z <- cbind(Z1, Z2)
  fitZ <- as.vector(Z%*%phivec)
  nummat <- matrix(exp(fitZ), nrow = 3, byrow = FALSE)
  phatmatphi <- nummat/apply(nummat, 2, sum)
  ZDPG1 <- apply(sapply(1:sum(groupind[,1]==1), function(j){  t( Z[((j-1)*R + 1):(j*R),][-1,]  )%*%((groupind[groupind[,1]==1,-1][j,] - phatmatphi[-1,j])) }), 1, sum)
  sum(as.vector(ZDPG1)^2)
  
}


computePK22 <- function(phivecy, phiy1, y1, y3, x1, x2){
  #Z1 <- kronecker(cbind(1, y1), c(0, 1, 0))
  denpi <- 1/(1 + exp(phiy1 + phivecy*y1))
  sum( x1/denpi  ) - sum(c(x1, x2))
  
}

computePK20 <- function(phiy1,phivecy,  y1, y3, x1, x2){
  #Z1 <- kronecker(cbind(1, y1), c(0, 1, 0))
  denpi <- 1/(1 + exp(phiy1 + phivecy*y1))
  sum( 1/denpi  ) - length(c(x1, x2))
  
}




computePK31 <- function(phivecy, phiy1, y1, y3, x1, x2){
 # Z1 <- kronecker(cbind(1, x1), c(0, 1, 0))
  denpi <- 1/(1 + exp(phiy1 + phivecy*x1))
  sum( y1/denpi  ) - sum(c(y1, y3))
  
}

computePK30 <- function( phiy1,phivecy, y1, y3, x1, x2){
 # Z1 <- kronecker(cbind(1, x1), c(0, 1, 0))
  denpi <- 1/(1 + exp(phiy1 + phivecy*x1))
  sum( 1/denpi  ) - length(c(y1, y3))
  
}










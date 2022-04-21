rep.row <- function(x, n){
  matrix(rep(x, each = n), nrow = n)
}

covCor <- function(Y, k = -1) {
  dim.Y <- dim(Y)
  N <- dim.Y[1]
  p <- dim.Y[2]
  if (k < 0) {    # demean the data and set k = 1
    Y <- scale(Y, scale = F)
    k <- 1
  }
  n <- N - k    # effective sample size
  c <- p / n    # concentration ratio
  sample <- (t(Y) %*% Y) / n   
  
  # compute shrinkage target
  samplevar <- diag(sample)
  sqrtvar <- sqrt(samplevar)
  rBar <- (sum(sample / outer(sqrtvar, sqrtvar)) - p) / (p * (p - 1))
  target <- rBar * outer(sqrtvar, sqrtvar)
  diag(target) <- samplevar
  
  # estimate the parameter that we call pi in Ledoit and Wolf (2003, JEF)
  Y2 <- Y^2
  sample2 <- (t(Y2) %*% Y2) / n   
  piMat <- sample2 - sample^2
  pihat <- sum(piMat)
  
  # estimate the parameter that we call gamma in Ledoit and Wolf (2003, JEF)
  gammahat <- norm(c(sample - target), type = "2")^2
  
  # diagonal part of the parameter that we call rho 
  rho_diag <- sum(diag(piMat))
  
  # off-diagonal part of the parameter that we call rho 
  term1 <- (t(Y^3) %*% Y) / n;
  term2 <- rep.row(samplevar, p) * sample;
  term2 <- t(term2)
  thetaMat <- term1 - term2
  diag(thetaMat) <- 0
  rho_off <- rBar * sum(outer(1/sqrtvar, sqrtvar) * thetaMat)
  
  # compute shrinkage intensity
  rhohat <- rho_diag + rho_off
  kappahat <- (pihat - rhohat) / gammahat
  shrinkage <- max(0, min(1, kappahat / n))
  
  # compute shrinkage estimator
  sigmahat <- shrinkage * target + (1 - shrinkage) * sample
}
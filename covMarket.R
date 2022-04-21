rep.col <- function(x, n){
  matrix(rep(x, times = n), ncol = n, byrow = F)
}

covMarket <- function(Y, k = -1) {
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
  Ymkt <- matrix(apply(Y, 1, mean),ncol = 1)
  covmkt <- as.vector((t(Ymkt) %*% Y) / n)
  varmkt <- c((t(Ymkt) %*% Ymkt) / n)
  target <- outer(covmkt, covmkt) / varmkt
  diag(target) <- diag(sample)
  
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
  temp <- Y * rep.col(Ymkt, p)
  v1 <- (1/n) * t(Y2) %*% temp - rep.col(covmkt, p) * sample
  roff1 <- sum(v1 * t(rep.col(covmkt, p))) / varmkt - sum(diag(v1) * covmkt) / varmkt
  v3 <- (1/n) * t(temp) %*% temp - varmkt * sample
  roff3 <- sum(v3 * (covmkt %*% t(covmkt))) / varmkt^2 - sum(diag(v3) * covmkt^2) / varmkt^2
  rho_off <- 2 * roff1 - roff3
  
  # compute shrinkage intensity
  rhohat <- rho_diag + rho_off
  kappahat <- (pihat - rhohat) / gammahat
  shrinkage <- max(0, min(1, kappahat / n))
  
  # compute shrinkage estimator
  sigmahat <- shrinkage * target + (1 - shrinkage) * sample
}
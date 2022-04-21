cov1Para <- function(Y, k = -1) {
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
  meanvar <- mean(diag(sample))
  target <- meanvar * diag(p)
  
  # estimate the parameter that we call pi in Ledoit and Wolf (2003, JEF)
  Y2 <- Y^2
  sample2 <- (t(Y2) %*% Y2) / n   
  piMat <- sample2 - sample^2
  pihat <- sum(piMat)
  
  # estimate the parameter that we call gamma in Ledoit and Wolf (2003, JEF)
  gammahat <- norm(c(sample - target), type = "2")^2
  
  # diagonal part of the parameter that we call rho 
  rho_diag <- 0
  
  # off-diagonal part of the parameter that we call rho 
  rho_off <- 0
  
  # compute shrinkage intensity
  rhohat <- rho_diag + rho_off
  kappahat <- (pihat - rhohat) / gammahat
  shrinkage <- max(0, min(1, kappahat / n))
  
  # compute shrinkage estimator
  sigmahat <- shrinkage * target + (1 - shrinkage) * sample
}
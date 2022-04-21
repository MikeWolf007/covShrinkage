cov2Para <- function(Y, k = -1) {
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
  id_p <- diag(p)
  one_p <- matrix(rep(1, p^2), ncol = p)
  
  # compute shrinkage target
  meanvar <- mean(diag(sample))
  meancovar <- sum(sample - diag(diag(sample))) / (p * (p- 1))
  target <- meanvar * id_p + meancovar * (one_p - id_p)
  
  # estimate the parameter that we call pi in Ledoit and Wolf (2003, JEF)
  Y2 <- Y^2
  sample2 <- (t(Y2) %*% Y2) / n   
  piMat <- sample2 - sample^2
  pihat <- sum(piMat)
  
  # estimate the parameter that we call gamma in Ledoit and Wolf (2003, JEF)
  gammahat <- norm(c(sample - target), type = "2")^2
  
  # diagonal part of the parameter that we call rho 
  rho_diag <- sum(sample2) / p - (sum(diag((sample))))^2 / p
  
  # off-diagonal part of the parameter that we call rho 
  sum1 <- apply(Y, 1, sum)
  sum2 <- apply(Y2, 1, sum)
  rho_off1 <- sum((sum1^2 - sum2)^2)/ p / n
  rho_off2 <- (sum(sample) - sum(diag(sample)))^2 / p
  rho_off <- (rho_off1 - rho_off2) / (p - 1)
  
  # compute shrinkage intensity
  rhohat <- rho_diag + rho_off
  kappahat <- (pihat - rhohat) / gammahat
  shrinkage <- max(0, min(1, kappahat / n))
  
  # compute shrinkage estimator
  sigmahat <- shrinkage * target + (1 - shrinkage) * sample
}
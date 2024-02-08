rep.row <- function(x, n){
  matrix(rep(x, each = n), nrow = n )
}

lis <- function(Y, k = -1) {
  dim.Y <- dim(Y)
  N <- dim.Y[1]
  p <- dim.Y[2]
  if (k < 0) {    # demean the data and set k = 1
    Y <- scale(Y, scale = F)
    k <- 1
  }
  n <- N - k    # effective sample size
  c <- p / n    # concentration ratio
  sample <- (t(Y) %*% Y) / n    # sample covariance matrix    
  sample <- (t(sample) + sample) / 2   # enforce symmetry (even more)
  spectral <- eigen(sample, symmetric = T)    # spectral decompositon
  lambda <- spectral$values[p:1]    # sort eigenvalues in ascending order
  u <- spectral$vectors[,p:1]    # eigenvectors follow their eigenvalues
  h <- min(c^2, 1/c^2)^0.35 / p^0.35    # smoothing parameter
  invlambda <- 1 / lambda[max(1, p-n+1):p]    # inverse of non-null eigenvalues   
  Lj <- rep.row(invlambda, min(p, n))    # like 1 / lambda_j
  Lj.i <- Lj - t(Lj)    # like (1 / lambda_j) - (1 / lambda_i)
  theta <- rowMeans(Lj * Lj.i / (Lj.i^2 + h^2 * Lj^2))    # smoothed Stein shrinker
  if (p <= n) {   # case where sample covariance matrix is not singular
    delta <- (1 - c) * invlambda + 2 * c * invlambda * theta # shrunk inverse eigenvalues
    deltaLIS <- pmax(delta, min(invlambda))
  }
  else {    # case where sample covariance matrix is singular
    stop("p must be <= n for Stein''s loss")
  }
  sigmahat <- u %*% diag(1 / deltaLIS) %*% t(u)    #reconstruct covariance matrix
}

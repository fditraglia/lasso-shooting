#First attempt to code the shooting algorithm for LASSO
#Frank DiTraglia

#softmax function
softmax_R <- function(x, y){
  sign(x) * pmax(abs(x) - y, 0)
}

#Lasso via coordinate descent: the "Shooting Algorithm" of Fu (1998). Adapted from pseudocode algorithm 13.1 of Murphy (2012) and matlab code LassoShooting.m by Mark Schmidt.
lasso_shoot_R <- function(X, y, lambda, tol = 1e-5, 
                          max_iter = 10000){  
  p <- ncol(X)
  XX <- crossprod(X, X)
  XX2 <- 2 * XX
  Xy <- crossprod(X, y)
  Xy2 <- 2 * Xy
  
  beta <- solve(XX + diag(lambda, p, p), Xy)
  
  converged <- FALSE
  iteration <- 0
  
  while (!converged & (iteration < max_iter)){
    
    beta_prev <- beta
    
    for (j in 1:p){
      aj <- XX2[j,j]
      cj <- Xy2[j] - sum(XX2[j,] %*% beta) + beta[j] * XX2[j,j]
      beta[j] <- softmax_R(cj / aj, lambda / aj)
    }
    
    iteration <- iteration + 1
    converged <- sum(abs(beta - beta_prev)) < tol
  }
  out <- list(beta = beta, n_iter = iteration, 
              converged = converged)
  return(out)
}
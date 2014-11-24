library(MASS)

design1 <- function(Rsq.1, Rsq.2, n = 100, p = 200, 
                    a.0 = 0.5, rho = 0.5, s.nu = 1, s.xi = 1){
#---------------------------------------------------
# Generate one simulation draw from Design 1 of
# Belloni, Chernozhukov and Hansen (2014)
#---------------------------------------------------
# Arguments:
#     Rsq.1   first-stage R-squared
#     Rsq.2   second-stage R-squared
#     n       sample size
#     p       number of "control" regressors (x)
#     a.0     treatment effect
#     rho     controls covariance matrix of x
#     s.nu    std dev of first-stage errors
#     s.xi    std dev of second-stage errors
#
# Returns:    list (y, d, x) where y is a vector of
#             outcomes, d is a vector of treatments
#             and x is a matrix of control regressors
#
# Notes:      Requires the library MASS
#---------------------------------------------------
  b.0 <- ( 1 / (1:p) )^2 #decay profile of coefficients
  Sigma <- toeplitz(rho^(0:(p-1))) #cov matrix of x
  
  #calculate first-stage coefficients
  bSb <- crossprod(b.0, Sigma) %*% b.0
  c.d <- sqrt(Rsq.1 / ( (1 - Rsq.1) * bSb))
  theta.m <- c.d * b.0 #first-stage coefficients
  
  #calculate second-stage coefficients
  z.2 <- (1 - Rsq.2) * bSb
  z.1 <- 2 * (1 - Rsq.2) * a.0 * c.d * bSb
  z.0 <- (1 - Rsq.2) * (a.0 * c.d)^2 * bSb - Rsq.2 * (a.0^2 * s.nu^2 + s.xi^2) 
  discriminant <- z.1^2 - 4 * z.2 * z.0
  discriminant <- ifelse(abs(discriminant) > 1e-12, 
                         discriminant, 0)
  c.y <- (-z.1 + sqrt(discriminant)) / (2 * z.2)
  theta.g <- c.y * b.0 #second-stage coefficients
  
  #draw simulations
  xi <- rnorm(n, mean = 0, sd = s.xi)
  nu <- rnorm(n, mean = 0, sd = s.nu)
  x <- mvrnorm(n, rep(0, p), Sigma)
  d <- x %*% theta.m + nu
  y <- a.0 * d + x %*% theta.g + xi
  
  return(list(y = y, d = d, x = x))
}


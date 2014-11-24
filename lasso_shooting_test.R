library(Rcpp)
library(RcppArmadillo)
source("lasso_shooting.R")
sourceCpp("lasso_shooting.cpp")


#Test softmax_R and softmax_cpp
w <- seq(from = -3, 3, 0.01)
lam <- rep(1, length(w)) 
wsoft_R <- softmax_R(w, lam)
wsoft_cpp <- softmax_cpp(w, lam)
all.equal(wsoft_R, as.vector(wsoft_cpp))
plot(w, wsoft_R, type = 'l')
plot(w, wsoft_cpp, type = 'l')


#Test lasso_shoot_R and lasso_shoot_cpp
source("simulation_designs.R")
set.seed(3827)
sim_data <- design1(Rsq.1 = 0.2, Rsq.2 = 0.2)
x <- sim_data$x
d <- sim_data$d
lambda <- 1
system.time(results_R <- lasso_shoot_R(x, d, lambda))
system.time(results_cpp <- lasso_shoot_cpp(x, d, lambda))
all.equal(results_R, results_cpp)

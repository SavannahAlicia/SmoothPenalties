#demonstrate simple smooth model of three dimensional data with Duchon spline
# with separate penalty for non-null basis (second order derivatives+) and
# null space 
library(mgcv)
library(Rcpp)
sourceCpp("src/functions.cpp")

set.seed(42)
n <- 400

x1 <- runif(n, -1, 1)
x2 <- runif(n, -1, 1)
x3 <- runif(n, -1, 1)

f_true <- sin(pi * x1) + x2^2 - x1 * x3
y <- f_true + rnorm(n, sd = 0.3)

dat <- data.frame(x1, x2, x3, y)

fit_duc <- gam(y ~ s(x1, x2, x3, bs = "ds", m = c(2, 0)),
               family = gaussian(link = "identity"),
               data = dat)

sm <- smoothCon(s(x1, x2, x3, bs = "ds", m = c(2, 0)),
                data = dat)[[1]]

X <- sm$X
S <- sm$S[[1]]

X_full <- cbind(1, X)  # add functional intercept
S_full <- matrix(0, ncol(X_full), ncol(X_full))
S_full[-1,-1] <- S # I think this means I have two penalty matrix rows for intercepts?

beta_init <- rep(1, ncol(X_full))
lambda = 10

opt <- optim(beta_init, 
             negllk_norm,
             X = X, 
             y = y,
             S = S,
             lambda = lambda, method = "BFGS")


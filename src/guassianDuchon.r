#demonstrate simple smooth model of three dimensional data with Duchon spline
# with separate penalty for non-null basis (second order derivatives+) and
# null space 
library(mgcv)
library(tidyverse)
library(Rcpp)
sourceCpp("src/functions.cpp")

set.seed(42)
n <- 400

x1 <- runif(n, -1, 1)
x2 <- runif(n, -1, 1)
#x3 <- runif(n, -1, 1)

f_true <- sin(pi * x1) + x2^2 - x1 #* x3
y <- f_true + rnorm(n, sd = 0.3)

dat <- data.frame(x1, x2, #x3,
                  y)

fit_duc <- gam(y ~ s(x1, x2, #x3, 
                     bs = "ds", m = c(2, 0)),
               family = gaussian(link = "identity"),
               data = dat)
lambda = fit_duc$sp

sm <- smoothCon(s(x1, x2, #x3, 
                  bs = "ds", m = c(2, 0)),
                data = dat)[[1]]

X <- sm$X
S <- sm$S[[1]]
E <- eigen(S, symmetric = TRUE)
U <- E$vectors
D <- E$values

tol <- 1e-10
null_idx <- which(D < tol)
pen_idx  <- which(D >= tol)
# Rotate design matrix
X_rot <- X %*% U
# Split
X_null <- X_rot[, null_idx, drop = FALSE]
X_pen  <- X_rot[, pen_idx,  drop = FALSE]
# Rescale penalized part so penalty = identity
D_pen <- D[pen_idx]
X_pen_scaled <- X_pen %*% diag(1 / sqrt(D_pen))
X_new <- cbind(X_null, X_pen_scaled)
p_null <- ncol(X_null)
p_pen  <- ncol(X_pen_scaled)

S_new <- diag(c(rep(0, p_null), rep(1, p_pen)))

beta_init <- rep(1, ncol(X_new))


fitmy <- optim(beta_init, 
             negllk_norm,
             X = X_new, 
             y = y,
             Snewdiag = diag(S_new),
             lambda = lambda,
             method = "BFGS", control = list(maxit = 100))

beta_new <- fitmy$par
alpha_null <- beta_new[1:p_null]
theta      <- beta_new[(p_null + 1):(p_null + p_pen)]
alpha_pen <- theta / sqrt(D_pen)
U_null <- U[, null_idx, drop = FALSE]
U_pen  <- U[, pen_idx,  drop = FALSE]
beta_orig <- U_null %*% alpha_null + U_pen %*% alpha_pen
pred_my <- as.vector(X %*% beta_orig)
comparedat <- data.frame(y = y,
                         x1 = x1,
                         x2 = x2, 
                        # x3 = x3,
                         my = pred_my,
                         mgcv = predict.gam(fit_duc, dat)
)
plotdat <- comparedat |>
  pivot_longer(
    cols = c(my, mgcv, y),
    names_to = "source",
    values_to = "value"
  )                         
ggplot() +
  geom_point(plotdat, mapping = aes(x = x1, y = x2, color = value)) +
  scale_color_viridis_c(option = "magma")+
  facet_wrap(~ source, nrow = 1) +
  theme_minimal()

#with regular grid
grid <- expand.grid(
  x1 = seq(min(x1), max(x1), length.out = 100),
  x2 = seq(min(x2), max(x2), length.out = 100)
)
X_pred <- PredictMat(sm, data = grid)
grid$my <- X_pred %*% beta_orig
grid$mgcv <- predict(fit_duc, grid)
grid$y <- sin(pi * grid$x1) + grid$x2^2 - grid$x1 

plotdatg <- grid |>
  pivot_longer(cols = c(my, mgcv, y),
               names_to = "model",
               values_to = "value")

ggplot() +
  geom_raster(plotdatg, mapping = aes(x = x1, y = x2, fill = value)) +
  facet_wrap(~ model) +
  scale_fill_viridis_c() +
  coord_equal() +
  theme_minimal()

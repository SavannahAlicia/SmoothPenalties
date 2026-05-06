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
               data = dat,
               optimizer = c("outer", "bfgs"))
lambda = fit_duc$sp

sm <- smoothCon(s(x1, x2, #x3, 
                  bs = "ds", m = c(2, 0)),
                data = dat)[[1]]

X <- sm$X
S <- sm$S[[1]]
beta_init <- rep(1, ncol(X))

fitmy <- optim(beta_init, 
               negllk_norm,
               X = X, 
               y = y,
               S = S,
               lambda = lambda,
               method = "BFGS", control = list(maxit = 100))

fitmy_nopen <- optim(beta_init, 
               negllk_norm,
               X = X, 
               y = y,
               S = (S),
               lambda = 0,
               method = "BFGS", control = list(maxit = 100))
#check if this gets same issue when i fit basis with lm
Xdf <- as.data.frame(X)
Xdf$y <- y
fit_lm <- lm(y ~ V1+ V2+ V3+ V4+ V5+ V6+ V7+ V8+ V9+ V10+ V11+ V12+ V13+ 
               V14+ V15+ V16+ V17+ V18+ V19+ V20+ V21+ V22+ V23+ V24+ V25+
               V26+ V27+ V28+ V29+ V30+ V31+ V32+ V33 -1, 
             data = Xdf, )

#with rotation like mgcv does maybe?
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
fitmyrot <- optim(beta_init, 
               negllk_norm,
               X = X_new, 
               y = y,
               S = S_new,
               lambda = lambda,
               method = "BFGS", control = list(maxit = 100))




#data points from simulation
beta_my <- fitmy$par
beta_nopen <- fitmy_nopen$par
beta_rot <- fitmyrot$par
pred_rot <- as.vector(X_new %*% beta_rot)
pred_my <- as.vector(X %*% beta_my)
pred_nopen <- as.vector(X %*% beta_nopen)
comparedat <- data.frame(y = y,
                         x1 = x1,
                         x2 = x2, 
                        # x3 = x3,
                         my = pred_my,
                        rot = pred_rot,
                        nopen = pred_nopen,
                         mgcv = predict.gam(fit_duc, dat)
)
plotdat <- comparedat |>
  pivot_longer(
    cols = c(my, mgcv, nopen,
             y),
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
grid$my <- as.vector(X_pred %*% beta_my)
grid$mgcv <- predict(fit_duc, grid)
grid$y <- sin(pi * grid$x1) + grid$x2^2 - grid$x1 
grid$diff_mgvme_p <- (grid$mgcv - grid$my)/(grid$mgcv + 1e-16)
grid$diff_mgvme <- (grid$mgcv - grid$my)
grid$lm <- predict(fit_lm, as.data.frame(X_pred))
grid$mynopen <- as.vector(X_pred %*% fitmy_nopen$par )


plotdatg <- grid |>
  pivot_longer(cols = c(my, mgcv,lm, y),
               names_to = "model",
               values_to = "value")

ggplot() +
  geom_raster(plotdatg, mapping = aes(x = x1, y = x2, fill = value)) +
  facet_wrap(~ model) +
  scale_fill_viridis_c() +
  coord_equal() +
  theme_minimal()
ggplot() +
  geom_raster(grid, 
              mapping = aes(x = x1, y = x2, fill = diff_mgvme)) +
  scale_fill_viridis_c() +
  coord_equal() +
  theme_minimal()
ggplot() +
  geom_raster(grid, 
              mapping = aes(x = x1, y = x2, fill = diff_mgvme_p)) +
  scale_fill_viridis_c() +
  coord_equal() +
  theme_minimal()


##---- demo cross validation---------------------------------------------------
K = 5
folds <- sample(rep(1:K, length.out = n))
cv_score <- function(X, y, S, lambda, beta_init, folds) {
  
  total <- 0
  
  for (k in 1:K) {
    
    test_idx <- which(folds == k)
    train_idx <- setdiff(1:n, test_idx)
    
    X_train <- X[train_idx, , drop = FALSE]
    y_train <- y[train_idx]
    
    X_test <- X[test_idx, , drop = FALSE]
    y_test <- y[test_idx]
    
    beta_hat <- optim(beta_init,
                 negllk_norm,   # penalized
                 X = X_train,
                 y = y_train,
                 S = S,
                 lambda = lambda,
                 method = "BFGS")$par
    beta_init <- beta_hat
    
    #unpenalized likelihood for score
    ll <- negllk_norm(beta_hat, X_test, y_test, S, lambda, incl_pen = F)
    
    total <- total + ll  
  }
  
  return(total/length(y)) #mean
}

lambda_grid <- data.frame(lambda = exp(seq(-5, 5, length.out = 40)))

lambda_grid$cv_vals <- sapply(lambda_grid$lambda, function(lam) {
  cv_score(X, y, S, lam, beta_init, folds)
})
ggplot() +
  geom_line(data = lambda_grid, mapping = aes(x = log(lambda), y = cv_vals))

lambda_hat <- lambda_grid$lambda[which.min(lambda_grid$cv_vals)]
fitmy_lhat <- optim(beta_init, 
                     negllk_norm,
                     X = X, 
                     y = y,
                     S = (S),
                     lambda = lambda_hat,
                     method = "BFGS", control = list(maxit = 100))
grid$mylhat <- as.vector(X_pred %*% fitmy_lhat$par)
plotdatg <- grid |>
  pivot_longer(cols = c(mylhat, mgcv, y),
               names_to = "model",
               values_to = "value")

ggplot() +
  geom_raster(plotdatg, mapping = aes(x = x1, y = x2, fill = value)) +
  facet_wrap(~ model) +
  scale_fill_viridis_c() +
  coord_equal() +
  theme_minimal()


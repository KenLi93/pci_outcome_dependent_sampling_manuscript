naive_logit_reg <- function(data) {
  logit_reg <- glm(Y ~ A + X + Z + W, data = data, family = "binomial")
  
  coef_reg <- summary(logit_reg)$coefficients
  return(list(est = coef_reg["A", "Estimate"], 
              var = coef_reg["A", "Std. Error"] ^ 2, 
              se = coef_reg["A", "Std. Error"]))
}

oracle_logit_reg <- function(data) {
  logit_reg <- glm(Y ~ A + X + U, data = data, family = "binomial")
  
  coef_reg <- summary(logit_reg)$coefficients
  return(list(est = coef_reg["A", "Estimate"], 
              var = coef_reg["A", "Std. Error"] ^ 2, 
              se = coef_reg["A", "Std. Error"]))
}

qfunc <- function(tau, data){
  q_est <- with(data, c(1 + exp((1 - 2 * A) * (cbind(1, A, Z, X) %*% tau))))
  return(q_est)
}

hfunc <- function(eta, data){
  h_est <- with(data, c(exp(cbind(1, A, W, X) %*% eta)))
  return(h_est)
}



U_DR <- function(par, data) {
  beta <- par[1]
  tau <- par[2:5]
  eta <- par[6:9]
  
  q_est <- qfunc(tau, data)
  h_est <- hfunc(eta, data)
  h1 <- with(data, c(exp(cbind(1, 1, W, X) %*% eta)))
  h0 <- with(data, c(exp(cbind(1, 0, W, X) %*% eta)))
  U_q <- with(data,
              c((1 - Y) * q_est) * cbind(1, A, W, X, A * W, A * X, W * X, A * W * X) -
                (1 - Y) * cbind(2, 1, 2 * W, 2 * X, W, X, 2 * W * X, W * X))
  U_h <- with(data,
              c((1 - Y) * h_est - Y) * cbind(1, A, Z, X, A * Z, A * X, Z * X, A * Z * X))
  
  U_rr <- with(data,
               (1 - Y) * (h1 - h0 * exp(beta)) -
                 ((1-Y) * h_est - Y) * (2 * A - 1) * q_est * exp(beta * (1 - A)))
  U <- cbind(U_q, U_h, U_rr)
  return(U)
}

U_OR <- function(par, data){
  beta <- par[1]
  eta <- par[2:5]

  h_est <- hfunc(eta, data)
  h1 <- with(data, c(exp(cbind(1, 1, W, X) %*% eta)))
  h0 <- with(data, c(exp(cbind(1, 0, W, X) %*% eta)))

  U_h <- with(data,
              c((1 - Y) * h_est - Y) * cbind(1, A, Z, X, A * Z, A * X, Z * X, A * Z * X))
  
  U_rr <- with(data,
               (1 - Y) * (h1 - h0 * exp(beta)))
  U <- cbind(U_h, U_rr)
  return(U)
}





U_IPW <- function(par, data){
  beta <- par[1]
  tau <- par[2:5]
  
  q_est <- qfunc(tau, data)
  
  U_q <- with(data,
              c((1 - Y) * q_est) * cbind(1, A, W, X, A * W, A * X, W * X, A * W * X) -
                (1 - Y) * cbind(2, 1, 2 * W, 2 * X, W, X, 2 * W * X, W * X))
  
  U_rr <- with(data,
               (2 * A - 1) * Y * exp(- beta * A) * q_est)
  U <- cbind(U_q, U_rr)
  return(U)
}

GMMF <- function(mrf, param, data, ...){
  g0 <- mrf(par = param, data = data, ...)
  g <- apply(g0, 2, mean)
  gmmf <- sum(g ^ 2)
  
  return(gmmf)
}

G <- function(bfun, para, data){
  G <- numDeriv::jacobian(func = G1, bfun = bfun, x = para, data = data)
  return(G)
}

G1 <- function(bfun, param, data){
  G1 <- apply(bfun(param, data), 2, mean, na.rm=T)
  return(G1)
}

## sandwich variance
var.gmmf <- function(bfun, param, data){
  bG <- MASS::ginv(G(bfun, param, data))
  bg <- bfun(param, data)
  spsz <- dim(bg)[1]
  Omega <- t(bg) %*% bg / spsz
  Sigma <- bG %*% Omega %*% t(bG) / spsz
  return(Sigma)
}


nc_cont_dr <- function(data, beta_init = -1, qinit = rep(0, 4), hinit = rep(-0.5, 4)) {
  opt <- nlm(f = GMMF, p = c(beta_init, qinit, hinit), mrf = U_DR, data = data, iterlim = 200)
 
  param_all <- opt$estimate
  beta_est <- param_all[1]
  
  var_all <- var.gmmf(bfun = U_DR, para = param_all, data)
  var_beta <- var_all[1, 1]
  
  return(list(est = beta_est, var = var_beta, se = sqrt(var_beta)))
}

nc_cont_ipw <- function(data) {
  opt <- nlm(f = GMMF, p = rep(0, 5), mrf = U_IPW, data = data, iterlim = 200)
  
  param_all <- opt$estimate
  beta_est <- param_all[1]
  
  var_all <- var.gmmf(bfun = U_IPW, para = param_all, data)
  var_beta <- var_all[1, 1]
  
  return(list(est = beta_est, var = var_beta, se = sqrt(var_beta)))
}

nc_cont_or <- function(data, init = rep(-0.5, 5)) {
  opt <- nlm(f = GMMF, p = init, mrf = U_OR, data = data, iterlim = 200)
  
  param_all <- opt$estimate
  beta_est <- param_all[1]
  
  var_all <- var.gmmf(bfun = U_OR, para = param_all, data)
  var_beta <- var_all[1, 1]
  
  return(list(est = beta_est, var = var_beta, se = sqrt(var_beta)))
}


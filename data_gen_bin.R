
expit <- function(x) exp(x) / (1 + exp(x))
logit <- function(x) log(x / (1-x))

DGP_bin <- function(N, beta0) {
  U <- rbinom(N, 1, 0.5)  ## p_u = 0.5
  A <- rbinom(N, 1, 0.2 + 0.4 * U)  
  Z <- rbinom(N, 1, 0.2 + 0.1 * A + 0.4 * U + 0.2 * A * U)  ## p_0z = 0.4, p_uz = 0.2
  
  
  Y <- rbinom(N, 1, expit(logit(0.4) + beta0 * A - 0.7 * U))  
  
  W <- rbinom(N, 1, 0.2 + 0.4 * U)  ## p_0w = 0.4, p_uw = 0.2
  
  
  S <- rbinom(N, 1, exp(-1.7 + 0.2 * A + 0.4 * Y + 0.7 * U))
  
  
  full_data <- data.frame(A = A, Y = Y, Z = Z, W = W, S = S, U = U)
  
  study_data <- full_data[full_data$S==1, c("A", "Y", "Z", "W", "U")]
  
  return(study_data)
}
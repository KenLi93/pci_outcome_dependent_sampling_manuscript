
expit <- function(x) exp(x) / (1 + exp(x))
logit <- function(x) log(x / (1-x))

DGP_multinom <- function(N = 150000) {
  U <- runif(N, 0, 1)  ## p_u = 0.5
  X <- runif(N, 0, 1)
  
  eta0 <- logit(0.02)
  mu0A <- -1; muUA <- 1; muXA <- 1
  mu0Y <- log(0.02 / 0.98); muUY <- -2; muXY <- 2
  muXW <- 2; muUW <- 4
  muA1Z <- 0.25; muA2Z <- 0.15; muA3Z <- 0.2; muXZ <- 2; muUZ <- 4
  
  pA <- cbind(1,  
              expit(mu0A + muUA * U + muXA * X),  ## Pfizer 
              expit(mu0A + muUA * U + muXA * X),  ## Moderna
              expit(mu0A + muUA * U + muXA * X))  ## J & J
  
  A <- apply(pA, 1, function(pp) rmultinom(1, 1, pp)) ## multinomial treatment exposure
  A1 <- A[1,]
  A2 <- A[2,]
  A3 <- A[3,]
 
  Y <- rbinom(N, 1, expit(eta0 + log(0.85) * A1 + log(0.85) * A2 + log(0.6) * A3 + 
                            muUY * U + muXY * X))  
  
  W <- rnorm(N, muXW * X + muUW * U, 0.25)  ## p_0w = 0.4, p_uw = 0.2
  Z <- rnorm(N, muA1Z * A1 + muA2Z * A2 + muA3Z * A3 + muXZ * X + muUZ * U, 0.25)  ## p_0z = 0.4, p_uz = 0.2
  
  mu0S <- -5; muAS <- 0.4; muYS <- 4; muUS <- 0.3; muXS <- 0.2
  S <- rbinom(N, 1, exp(mu0S + muAS * (A1 + A2 + A3) + muYS * Y + muUS * U + muXS * X))
  
  
  full_data <- data.frame(A1 = A1, A2 = A2, A3 = A3, Y = Y, Z = Z, W = W, S = S, U = U, X = X)
  
  data <- full_data[full_data$S==1, c("A1", "A2", "A3", "Y", "Z", "W", "U", "X")]
  
  
  # psi_A <- log(0.2) 
  # psi_W <- muUY / muUW 
  # psi_X <- muXY - muXW * muUY / muUW
  # psi_0 <- mu0Y + muYS - 0.5 * psi_W ^ 2 * 0.25 ^ 2
  # hinit <- c(psi_0, psi_A, psi_W, psi_X)
  # 
  # print(hinit)
  # 
  # tau_Z <- muUZ / muUA; tau_X <- muXA - muXZ * muUZ / muUA
  # tau_A <- tau_Z ^ 2 * 0.25 ^ 2
  # tau_0 <- mu0A + muAS - 0.5 * tau_Z ^ 2 * 0.25 ^ 2
  # qinit <- c(tau_0, tau_A, tau_Z, tau_X)
  # 
  # print(qinit)
  
  return(data)
}

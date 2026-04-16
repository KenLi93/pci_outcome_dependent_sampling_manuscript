rm(list = ls())
library(numDeriv)
library(mvtnorm)
library(parallel)
source("data_gen_continuous.R")
source("nc_continuous_fun.R")



### One can choose different options below for replicating different cases of the simulation

NSIM <- 2000
n.cores <- ifelse(Sys.info()["sysname"] == "Windows", 1,
                  max(1, detectCores() - 1))

options(warn = 2)
param_grid <- expand.grid(N = c(100000, 200000),
                          beta0 = log(0.2),
                          eta0 = c(logit(0.02)))

ii <- as.numeric(Sys.getenv('LSB_JOBINDEX'))

N <- param_grid$N[ii]
beta0 <- param_grid$beta0[ii]
eta0 <- param_grid$eta0[ii]

sim_once <- function(N, beta0, eta0) {
  data <- DGP_continuous(N, eta0 = eta0, beta0 = beta0)
  
  ## candidates of hyperparameters
  naive_results <- try(naive_logit_reg(data), silent = T)
  oracle_results <- try(oracle_logit_reg(data), silent = T)
  ipw_results <- try(nc_cont_ipw(data), silent = T)
  dr_results <- try(nc_cont_dr(data), silent = T)
  or_results <- try(nc_cont_or(data), silent = T)
  
  if (class(naive_results) == "try-error" |
      class(oracle_results) == "try-error" |
      class(ipw_results) == "try-error" | 
      class(dr_results) == "try-error" | 
      class(or_results) == "try-error") {
    naive_est = NA; naive_se = NA; naive_var = NA
    oracle_est = NA; oracle_se = NA; oracle_var = NA
    IPW_est = NA; IPW_se = NA; IPW_var = NA
    DR_est = NA; DR_se = NA; DR_var = NA
    OR_est = NA; OR_se = NA; OR_var = NA
  } else {
    naive_est = naive_results$est; naive_se = naive_results$se; naive_var = naive_results$var
    oracle_est = oracle_results$est; oracle_se = oracle_results$se; oracle_var = oracle_results$var
    IPW_est = ipw_results$est; IPW_se = ipw_results$se; IPW_var = ipw_results$var
    DR_est = dr_results$est; DR_se = dr_results$se; DR_var = dr_results$var
    OR_est = or_results$est; OR_se = or_results$se; OR_var = or_results$var
  }
  
  c(nn = nrow(data),
    naive_est = naive_est, naive_se = naive_se, naive_var = naive_var,
    oracle_est = oracle_est, oracle_se = oracle_se, oracle_var = oracle_var,
    IPW_est = IPW_est, IPW_se = IPW_se, IPW_var = IPW_var,
    OR_est = OR_est, OR_se = OR_se, OR_var = OR_var,
    DR_est = DR_est, DR_se = DR_se, DR_var = DR_var)
  
}

result <- mclapply(1:NSIM, 
                   function(ii) {
                     set.seed(2026 + ii) 
                     sim_once(N = N, beta0 = beta0, eta0 = eta0)
                   }, 
                   mc.cores = n.cores)



saveRDS(result, 
        file = paste0("results/nc_cont_type_N_", N, 
                      "_beta0_", round(beta0 * 100), "_eta0_", round(eta0 * 100), ".rds"))






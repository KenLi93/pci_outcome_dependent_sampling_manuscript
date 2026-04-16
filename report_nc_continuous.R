library(dplyr)
library(xtable)
logit <- function(x) log(x / (1 - x))
param_grid <- expand.grid(N = c(100000, 200000),
                          beta0 = log(0.2),
                          eta0 = c(logit(0.02))) 

for (ii in 1:nrow(param_grid)) {
  N <- param_grid$N[ii]
  beta0 <- param_grid$beta0[ii]
  eta0 <- param_grid$eta0[ii]
  
  
  result <- readRDS(paste0("results/nc_cont_type_N_", N, 
                          "_beta0_", round(beta0 * 100), 
                          "_eta0_", round(eta0 * 100), ".rds")) %>% 
    bind_rows %>%
    as.data.frame
  
  param_grid$N_median[ii] <- median(result$nn, na.rm = T)
  param_grid$N_min[ii] <- min(result$nn, na.rm = T)
  param_grid$N_max[ii] <- max(result$nn, na.rm = T)
  
  param_grid$naive_bias[ii] <- mean(result$naive_est, na.rm = T) - beta0
  param_grid$naive_sd[ii] <- sd(result$naive_est, na.rm = T)
  param_grid$naive_se_mean[ii] <- mean(result$naive_se, na.rm = T)
  param_grid$naive_se_median[ii] <- median(result$naive_se, na.rm = T)
  param_grid$naive_se_q25[ii] <- quantile(result$naive_se, 0.25, na.rm = T)
  param_grid$naive_se_q75[ii] <- quantile(result$naive_se, 0.75, na.rm = T)
  param_grid$naive_se_min[ii] <- min(result$naive_se, na.rm = T)
  param_grid$naive_se_max[ii] <- max(result$naive_se, na.rm = T)
  param_grid$naive_coverage[ii] <- with(result, mean(naive_est - beta0 - qnorm(0.975) * naive_se < 0 &
                                                    naive_est - beta0 + qnorm(0.975) * naive_se > 0, na.rm = T))
  
  param_grid$oracle_bias[ii] <- mean(result$oracle_est, na.rm = T) - beta0
  param_grid$oracle_sd[ii] <- sd(result$oracle_est, na.rm = T)
  param_grid$oracle_se_mean[ii] <- mean(result$oracle_se, na.rm = T)
  param_grid$oracle_se_median[ii] <- median(result$oracle_se, na.rm = T)
  param_grid$oracle_se_q25[ii] <- quantile(result$oracle_se, 0.25, na.rm = T)
  param_grid$oracle_se_q75[ii] <- quantile(result$oracle_se, 0.75, na.rm = T)
  param_grid$oracle_se_min[ii] <- min(result$oracle_se, na.rm = T)
  param_grid$oracle_se_max[ii] <- max(result$oracle_se, na.rm = T)
  param_grid$oracle_coverage[ii] <- with(result, mean(oracle_est - beta0 - qnorm(0.975) * oracle_se < 0 &
                                                       oracle_est - beta0 + qnorm(0.975) * oracle_se > 0, na.rm = T))
  
  
  param_grid$OR_bias[ii] <- mean(result$OR_est, na.rm = T) - beta0
  param_grid$OR_sd[ii] <- sd(result$OR_est, na.rm = T)
  param_grid$OR_se_mean[ii] <- mean(result$OR_se, na.rm = T)
  param_grid$OR_se_median[ii] <- median(result$OR_se, na.rm = T)
  param_grid$OR_se_q25[ii] <- quantile(result$OR_se, 0.25, na.rm = T)
  param_grid$OR_se_q75[ii] <- quantile(result$OR_se, 0.75, na.rm = T)
  param_grid$OR_se_min[ii] <- min(result$OR_se, na.rm = T)
  param_grid$OR_se_max[ii] <- max(result$OR_se, na.rm = T)
  param_grid$OR_coverage[ii] <- with(result, mean(OR_est - beta0 - qnorm(0.975) * OR_se < 0 &
                                                    OR_est - beta0 + qnorm(0.975) * OR_se > 0, na.rm = T))
  
  param_grid$IPW_bias[ii] <- mean(result$IPW_est, na.rm = T) - beta0
  param_grid$IPW_sd[ii] <- sd(result$IPW_est, na.rm = T)
  param_grid$IPW_se_mean[ii] <- mean(result$IPW_se, na.rm = T)
  param_grid$IPW_se_median[ii] <- median(result$IPW_se, na.rm = T)
  param_grid$IPW_se_q25[ii] <- quantile(result$IPW_se, 0.25, na.rm = T)
  param_grid$IPW_se_q75[ii] <- quantile(result$IPW_se, 0.75, na.rm = T)
  param_grid$IPW_se_min[ii] <- min(result$IPW_se, na.rm = T)
  param_grid$IPW_se_max[ii] <- max(result$IPW_se, na.rm = T)
  param_grid$IPW_coverage[ii] <- with(result, mean(IPW_est - beta0 - qnorm(0.975) * IPW_se < 0&
                                                     IPW_est - beta0 + qnorm(0.975) * IPW_se > 0, na.rm = T))
  
  param_grid$DR_bias[ii] <- mean(result$DR_est, na.rm = T) - beta0
  param_grid$DR_sd[ii] <- sd(result$DR_est, na.rm = T)
  param_grid$DR_se_mean[ii] <- mean(result$DR_se, na.rm = T)
  param_grid$DR_se_median[ii] <- median(result$DR_se, na.rm = T)
  param_grid$DR_se_q25[ii] <- quantile(result$DR_se, 0.25, na.rm = T)
  param_grid$DR_se_q75[ii] <- quantile(result$DR_se, 0.75, na.rm = T)
  param_grid$DR_se_min[ii] <- min(result$DR_se, na.rm = T)
  param_grid$DR_se_max[ii] <- max(result$DR_se, na.rm = T)
  param_grid$DR_coverage[ii] <- with(result, mean(DR_est - beta0 - qnorm(0.975) * DR_se < 0 &
                                                    DR_est - beta0 + qnorm(0.975) * DR_se > 0, na.rm = T))
}


result_grid <- param_grid %>%
  transmute(N = N,
            n = paste0(N_median, "(", N_min, ",", N_max, ")"),
            OR_bias = format(round(OR_bias, 3), nsmall = 3),
            OR_sd = format(round(OR_sd, 3), nsmall = 3),
            OR_se = paste0(format(round(OR_se_median, 3), nsmall = 3),
                           "(", format(round(OR_se_q25, 3), nsmall = 3), 
                           ",", format(round(OR_se_q75, 3), nsmall = 3), ")"),
            OR_coverage = format(round(OR_coverage, 3), nsmall = 3),
            IPW_bias = format(round(IPW_bias, 3), nsmall = 3),
            IPW_sd = format(round(IPW_sd, 3), nsmall = 3),
            IPW_se = paste0(format(round(IPW_se_median, 3), nsmall = 3),
                            "(", format(round(IPW_se_q25, 3), nsmall = 3), 
                            ",", format(round(IPW_se_q75, 3), nsmall = 3), ")"),
            IPW_coverage = format(round(IPW_coverage, 3), nsmall = 3),
            DR_bias = format(round(DR_bias, 3), nsmall = 3),
            DR_sd = format(round(DR_sd, 3), nsmall = 3),
            DR_se = paste0(format(round(DR_se_median, 3), nsmall = 3),
                           "(", format(round(DR_se_q25, 3), nsmall = 3), 
                           ",", format(round(DR_se_q75, 3), nsmall = 3), ")"),
            DR_coverage = format(round(DR_coverage, 3), nsmall = 3))


result_grid_xtable <- xtable(result_grid,
                             digits = c(0, 0, 0, rep(3, 12)))

saveRDS(param_grid, file = "cont_v1_sim_results.rds")
print(result_grid_xtable, include.rownames = F,
      file = "cont_v1_xtable.txt")




result_grid2 <- param_grid %>%
  transmute(N = N,
            n = paste0(N_median, "(", N_min, ",", N_max, ")"),
            oracle_bias1 = format(round(oracle_bias, 3), nsmall = 3),
            oracle_sd1 = format(round(oracle_sd, 3), nsmall = 3),
            oracle_releff1 = format(1, nsmall = 0),
            oracle_se1 = paste0(format(round(oracle_se_median, 3), nsmall = 3),
                                "(", format(round(oracle_se_q25, 3), nsmall = 3), 
                                ",", format(round(oracle_se_q75, 3), nsmall = 3), ")"),
            oracle_coverage1 = format(round(oracle_coverage, 3), nsmall = 3),
            OR_bias1 = format(round(OR_bias, 3), nsmall = 3),
            OR_sd1 = format(round(OR_sd, 3), nsmall = 3),
            OR_releff1 = format(round(oracle_sd / OR_sd, 3), nsmall = 3),
            OR_se1 = paste0(format(round(OR_se_median, 3), nsmall = 3),
                            "(", format(round(OR_se_q25, 3), nsmall = 3), 
                            ",", format(round(OR_se_q75, 3), nsmall = 3), ")"),
            OR_coverage1 = format(round(OR_coverage, 3), nsmall = 3))

result_grid3 <- param_grid %>%
  transmute(N = N,
            n = paste0(N_median, "(", N_min, ",", N_max, ")"),
            IPW_bias1 = format(round(IPW_bias, 3), nsmall = 3),
            IPW_sd1 = format(round(IPW_sd, 3), nsmall = 3),
            IPW_releff1 = format(round(oracle_sd / IPW_sd, 3), nsmall = 3),
            IPW_se1 = paste0(format(round(IPW_se_median, 3), nsmall = 3),
                             "(", format(round(IPW_se_q25, 3), nsmall = 3), 
                             ",", format(round(IPW_se_q75, 3), nsmall = 3), ")"),
            IPW_coverage1 = format(round(IPW_coverage, 3), nsmall = 3),
            DR_bias1 = format(round(DR_bias, 3), nsmall = 3),
            DR_sd1 = format(round(DR_sd, 3), nsmall = 3),
            DR_releff1 = format(round(oracle_sd / DR_sd, 3), nsmall = 3),
            DR_se1 = paste0(format(round(DR_se_median, 3), nsmall = 3),
                            "(", format(round(DR_se_q25, 3), nsmall = 3), 
                            ",", format(round(DR_se_q75, 3), nsmall = 3), ")"),
            DR_coverage1 = format(round(DR_coverage, 3), nsmall = 3))

result_grid_xtable2 <- xtable(result_grid2, digits = c(0, 0, 0, rep(3, 10)))
result_grid_xtable3 <- xtable(result_grid3, digits = c(0, 0, 0, rep(3, 10)))

print(result_grid_xtable2, include.rownames = F, file = "complete_cont_v1_xtable_up.txt")
print(result_grid_xtable3, include.rownames = F, file = "complete_cont_v1_xtable_down.txt")
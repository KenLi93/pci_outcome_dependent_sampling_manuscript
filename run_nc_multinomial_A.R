rm(list = ls())
library(numDeriv)
library(mvtnorm)
library(dplyr)
library(parallel)
source("data_gen_multinomial_A.R")
source("nc_continuous_v1.R")

set.seed(2026)

### One can choose different options below for replicating different cases of the simulation


data <- DGP_multinom(N = 200000)

## estimate the VE of A1 relative to no treatment
## the proposed method in Web Appendix M for treatment of more than 2 categories is equivalent of
## the method in the main manuscript using a subset of the data with A = A1 or A = 0.
data_A1 <- data %>%
  as.data.frame %>%
  filter(A1 == 1 | (A1 == 0 & A2 == 0 & A3 == 0)) %>%
  rename(A = A1) %>%
  select(A, Y, Z, W, U, X)

naive_results_A1 <- try(naive_logit_reg(data_A1), silent = T)
ipw_results_A1 <- try(nc_cont_ipw(data_A1), silent = T)
or_results_A1 <- try(nc_cont_or(data_A1), silent = T)
dr_results_A1 <- try(nc_cont_dr(data_A1), silent = T)

## estimate the VE of A2 relative to no treatment

data_A2 <- data %>%
  as.data.frame %>%
  filter(A2 == 1 | (A1 == 0 & A2 == 0 & A3 == 0)) %>%
  rename(A = A2) %>%
  select(A, Y, Z, W, U, X)

naive_results_A2 <- try(naive_logit_reg(data_A2), silent = T)
ipw_results_A2 <- try(nc_cont_ipw(data_A2), silent = T)
or_results_A2 <- try(nc_cont_or(data_A2), silent = T)
dr_results_A2 <- try(nc_cont_dr(data_A2), silent = T)

## estimate the VE of A3 relative to no treatment

data_A3 <- data %>%
  as.data.frame %>%
  filter(A3 == 1 | (A1 == 0 & A2 == 0 & A3 == 0)) %>%
  rename(A = A3) %>%
  select(A, Y, Z, W, U, X)

naive_results_A3 <- try(naive_logit_reg(data_A3), silent = T)
ipw_results_A3 <- try(nc_cont_ipw(data_A3), silent = T)
or_results_A3 <- try(nc_cont_or(data_A3), silent = T)
dr_results_A3 <- try(nc_cont_dr(data_A3), silent = T)
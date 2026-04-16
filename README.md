# R code for “Doubly Robust Proximal Causal Inference under Confounded Outcome-Dependent Sampling”
This repository includes R Code for implementing the estimators and reproducing the simulation study in the manuscript "Doubly Robust Proximal Causal Inference under Confounded Outcome-Dependent Sampling".

Simulation with binary NCE, NCO, and unmeasured and measured confoundings (Scenario I in Section 7 of the manuscript):

- data_gen_bin.R: generating data.
- nc_bin_fun.R: implementing the oracle logistic regression estimator and the three proximal estimators with binary data.
- run_sim_nc_bin.R: main code to perform the simulation study.
- report_sim_bin.R: summarize and report the results of the simulation.

Simulation with condinuous NCE, NCO, and unmeasured and measured confoundings (Scenario II-V in Section 7 of the manuscript):

- data_gen_continuous.R: generating data.
- nc_continuous_fun.R: implementing the oracle logistic regression estimator and the three proximal estimators with continuous data (Scenario II).
- nc_continuous_fun_mis_trt_bridge.R: implementing the oracle logistic regression estimator and the three proximal estimators with continuous data, but the treatment confounding bridge function was misspecified (Scenario III).
- nc_continuous_fun_mis_outcome_bridge.R: implementing the oracle logistic regression estimator and the three proximal estimators with continuous data, but the outcome confounding bridge function was misspecified (Scenario IV).
- nc_continuous_fun_mis_both_bridge.R: implementing the oracle logistic regression estimator and the three proximal estimators with continuous data, but both the treatment and the outcome confounding bridge functions were misspecified (Scenario V).
- run_sim_nc_bin.R: main code to perform the simulation study.
- report_sim_bin.R: summarize and report the results of the simulation.

Preprint of the manuscript available at https://arxiv.org/abs/2208.01237.

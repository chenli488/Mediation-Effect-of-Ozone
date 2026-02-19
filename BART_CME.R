# This is the sample code for the manuscript: emperature and Respiratory Emergency Department Visits: 
# A Mediation Analysis with Ambient Ozone Exposure.

# Outcome: RESP; Mediator: q_O3; Exposure: w_tmax;

# About the data:
# 1. The analysis dataset is a simulated dataset from the real data, which will be in the Github respo as well.
# 2. Warning of perfect fitting will happen when running the program as the simulated data is generated 
# by the mediator model (linear) and outcome model (Quassi-Poisson) fitted on real data.

# Chen Li
# Feb 18th, 2026

library(cowplot)
library(gridExtra)
library('haven')
library(dplyr)
library(ggplot2)
library(scales)
library(naniar)
library('lubridate')
library(splines)
library(readxl)
library('rsample')
library('purrr')
library('numDeriv')
library('dvmisc')
library('Matrix')
library("predint")
library('gtools')
library('BART')

rm(list=ls(all.names=T))

## Read in simulated dataset -----
path <- '/Users/lichen/Library/CloudStorage/OneDrive-Emory/Documents/Research/Mediation analysis/JRSSA'
load(paste0(path, '/SimDS.RData'))

combine_la_warm <- SimDS

## Functions to Galculate CME ------
# 1. Simulate the coefficients
# Only simulate theta:
coef_sim <- function(data, outcome, temperature, pollutant, n){
  # data = combine_la_warm
  # outcome = 'RESPcnt'
  # temperature = 'w_tmax'
  # pollutant = 'q_O3'
  # n = 1000
  
  # The natural cubic spline for weighted tmax
  ncs <- ns(data[[temperature]], df = 6)
  
  ## The X + M ~ Y model -----
  # Quasi-Poisson regression
  qp_formula <- paste0(outcome, ' ~ ns(', temperature, ',df = 6) + ns(w_spfh, df = 6) + factor(HOLIDAY) + factor(WEEKDAY) + ns(day,df=6)*factor(year) + ', pollutant, ' + inter1 + inter2 + inter3')
  
  XM_Y <- glm(formula = qp_formula,
              family = 'quasipoisson',
              data = data)
  # summary(XM_Y)
  
  # Extract out the coefficients
  theta_coef <- summary(XM_Y)$coefficient[,c(1,2)]
  cov_XMY <- vcov(XM_Y)
  
  # Generate the simulation set.
  theta_simu <- mvrnorm(n, mu = theta_coef[,1], Sigma = cov_XMY)
  theta_simu <- as.data.frame(theta_simu)
  colnames(theta_simu) <- paste0("Theta_", rownames(theta_coef))
  
  simu_result <- list(theta = theta_simu)
  return(simu_result)
}

coef_simu_la <- coef_sim(data = combine_la_warm, outcome = 'RESPcnt', temperature = 'w_tmax', pollutant = 'q_O3', n = 20000)

# 2.Estimate the mediation formular
mf_simu <- function(data, outcome, temperature, pollutant, coeff, x_exp, x_star){
  # data = combine_la_warm
  # outcome = 'RESPcnt'
  # temperature = 'w_tmax'
  # pollutant = 'q_O3'
  # coeff = coef_simu_la
  # x_exp = 0.9
  # x_star = 0.5
  
  # The natural cubic spline for weighted tmax
  ncs <- ns(data[[temperature]], df = 6)
  
  # The basis function for exposed level
  fixed_exp_level <- quantile(data[[temperature]], x_exp)
  exp_matrix <- as.matrix(ns(fixed_exp_level, df = 6, knots = attr(ncs, "knots"), Boundary.knots = attr(ncs, "Boundary.knots")))
  exp_values <- as.numeric(exp_matrix)
  
  # The basis function for reference level
  fixed_star_level <- quantile(data[[temperature]], x_star)
  star_matrix <- as.matrix(ns(fixed_star_level, df = 6, knots = attr(ncs, "knots"), Boundary.knots = attr(ncs, "Boundary.knots")))
  star_values <- as.numeric(star_matrix)
  
  ## The X ~ M BART model
  # The formula
  xm_formula <- ~ w_tmax + w_spfh + factor(HOLIDAY) + factor(WEEKDAY) + day + factor(year)
  
  # covariates
  x <- model.matrix(xm_formula, data)
  x_train <- x[, -1] # Exclude intercept
  
  # outcome
  y <- data[,pollutant]
  y_train <- as.matrix(y)
  
  # BART model
  set.seed(666)
  nd <- 20000
  burn <- 5000
  post.xm <- wbart(x_train, y_train, nskip = burn, ndpost = nd, keepevery = 4)
  
  # Sigma
  post.sigma <- post.xm$sigma[-c(1:burn)]
  sigma <- post.sigma[seq(1, length(post.sigma), by = 4)]
  
  # Parameter for Dirichlet distribution 
  alpha <- rep(1, nrow(data))
  
  # The predicted posteriors for x_star
  x_test_star <- as.data.frame(x_train) %>%
    mutate(!!sym(temperature):= fixed_star_level)
  
  x_test_star <- as.matrix(x_test_star)
  M.hat_star <- predict(post.xm, x_test_star)
  
  # The predicted posteriors for x_exp
  x_test_exp <- as.data.frame(x_train) %>%
    mutate(!!sym(temperature):= fixed_exp_level)
  
  x_test_exp <- as.matrix(x_test_exp)
  M.hat_exp <- predict(post.xm, x_test_exp)
  
  ## The X + M ~ Y model 
  # Quasi-Poisson regression
  qp_formula <- paste0(outcome, ' ~ ns(', temperature, ',df = 6) + ns(w_spfh, df = 6) + factor(HOLIDAY) + factor(WEEKDAY) + ns(day,df=6)*factor(year) + ', pollutant, ' + inter1 + inter2 + inter3')
  
  simu_result <- data.frame()
  
  for (k in 1:nrow(coeff[['theta']])){
    # k = 1
    # extract the coefficients
    theta <- as.matrix(coeff[['theta']][k,])
    
    # extract the posterior prediction based on reference level exposure, exposed level exposure, and sigma
    pred.M_star <- as.matrix(M.hat_star[k,])
    pred.M_exp <- as.matrix(M.hat_exp[k,])
    sigma_k <- sigma[k]
    
    # Weight by a Dirichlet distribution
    weight <- rdirichlet(n = 1, alpha = alpha)
    
    ## Estimation of the expectation term
    # decide theta2 and theta3x
    theta2 <- theta[, grep(pollutant, colnames(theta))]
    
    theta3x_star <- case_when(x_star > 0 & x_star < 0.25 ~ 0,
                              x_star >= 0.25 & x_star < 0.50 ~ theta[, 'Theta_inter1'],
                              x_star >= 0.50 & x_star < 0.75 ~ theta[, 'Theta_inter2'],
                              x_star >= 0.75 & x_star <= 1.00 ~ theta[, 'Theta_inter3'],
                              TRUE ~ 0)
    
    theta3x_exp <- case_when(x_exp > 0 & x_exp < 0.25 ~ 0,
                             x_exp >= 0.25 & x_exp < 0.50 ~ theta[, 'Theta_inter1'],
                             x_exp >= 0.50 & x_exp < 0.75 ~ theta[, 'Theta_inter2'],
                             x_exp >= 0.75 & x_exp <= 1.00 ~ theta[, 'Theta_inter3'],
                             TRUE ~ 0)
    
    # get the expectation for every t for the kth iteration, for four different combinations.
    # 1. star & star
    expectation_star_star  <- exp((theta2 + theta3x_star)*pred.M_star + 0.5*sigma_k^2*(theta2 + theta3x_star)^2)
    
    # 2. star $ exp
    expectation_star_exp  <- exp((theta2 + theta3x_star)*pred.M_exp + 0.5*sigma_k^2*(theta2 + theta3x_star)^2)
      
    # 3. exp & star
    expectation_exp_star  <- exp((theta2 + theta3x_exp)*pred.M_star + 0.5*sigma_k^2*(theta2 + theta3x_exp)^2)
    
    # 4. exp & exp
    expectation_exp_exp  <- exp((theta2 + theta3x_exp)*pred.M_exp + 0.5*sigma_k^2*(theta2 + theta3x_exp)^2)
    
    ## Get the constant part of the mediation formular.
    # parameters
    theta0 <- as.matrix(theta[, grep('Intercept', colnames(theta))])
    theta1 <- as.matrix(theta[, grep(temperature, colnames(theta))])
    theta_indices <- unique(c(
      grep(temperature, colnames(theta)),
      grep(pollutant, colnames(theta)),
      grep('inter', colnames(theta)),
      grep("(Intercept)", colnames(theta))
    ))
    theta4 <- as.matrix(theta[ ,-theta_indices])
    
    # design matrix
    covariates <- model.matrix(~ ns(w_spfh, df = 6) + factor(HOLIDAY) + factor(WEEKDAY) + ns(day, df = 6) * factor(year), data = data)
    covariates <- covariates[,-1] # drop the intercept
    
    # calculation, for both expose and star level.
    # 1. star level
    cons_A_star <- matrix(rep(theta0 + t(theta1) %*% star_values, nrow(covariates)), nrow = nrow(covariates), byrow = TRUE)
    cons_star <- exp(cons_A_star + covariates %*% theta4)
    # 2. expose level
    cons_A_exp <- matrix(rep(theta0 + t(theta1) %*% exp_values, nrow(covariates)), nrow = nrow(covariates), byrow = TRUE)
    cons_exp <- exp(cons_A_exp + covariates %*% theta4)
    
    ## Get the results for every t and get the mean over time.
    results_star_star <- as.data.frame(cbind(cons_star, expectation_star_star)) %>%
      setNames(c('constant', 'expectation')) %>%
      mutate(formula_esti = constant * expectation)
    
    results_exp_star <- as.data.frame(cbind(cons_exp, expectation_exp_star)) %>%
      setNames(c('constant', 'expectation')) %>%
      mutate(formula_esti = constant * expectation)
    
    results_star_exp <- as.data.frame(cbind(cons_star, expectation_star_exp)) %>%
      setNames(c('constant', 'expectation')) %>%
      mutate(formula_esti = constant * expectation)
    
    results_exp_exp <- as.data.frame(cbind(cons_exp, expectation_exp_exp)) %>%
      setNames(c('constant', 'expectation')) %>%
      mutate(formula_esti = constant * expectation)
    
    results_esti_star_star <- weight %*% as.matrix(results_star_star$formula_esti)
    results_esti_exp_star  <- weight %*% as.matrix(results_exp_star$formula_esti)
    results_esti_star_exp  <- weight %*% as.matrix(results_star_exp$formula_esti)
    results_esti_exp_exp  <- weight %*% as.matrix(results_exp_exp$formula_esti)
    # This is the simulation estimation for the ith set of beta and theta.
    # We gonna have n(beta) of this simulation estimation.
    
    simu_result <- rbind(simu_result, c(results_esti_star_star, results_esti_exp_star, results_esti_star_exp, results_esti_exp_exp))
    print(k)
  }
  
  simu_result <- simu_result %>%
    setNames(c('star_star', 'exp_star', 'star_exp', 'exp_exp'))
  
  pnde <- data.frame(type = 'pnde', 
                     estimate = mean(simu_result$exp_star/simu_result$star_star),
                     upper = quantile(simu_result$exp_star/simu_result$star_star, 0.975),
                     lower = quantile(simu_result$exp_star/simu_result$star_star, 0.025))
  
  tnde <- data.frame(type = 'tnde', 
                     estimate = mean(simu_result$exp_exp/simu_result$star_exp),
                     upper = quantile(simu_result$exp_exp/simu_result$star_exp, 0.975),
                     lower = quantile(simu_result$exp_exp/simu_result$star_exp, 0.025))
  
  tnie <- data.frame(type = 'tnie', 
                     estimate = mean(simu_result$exp_exp/simu_result$exp_star),
                     upper = quantile(simu_result$exp_exp/simu_result$exp_star, 0.975),
                     lower = quantile(simu_result$exp_exp/simu_result$exp_star, 0.025))
  
  pnie <- data.frame(type = 'pnie',
                     estimate = mean(simu_result$star_exp/simu_result$star_star),
                     upper = quantile(simu_result$star_exp/simu_result$star_star, 0.975),
                     lower = quantile(simu_result$star_exp/simu_result$star_star, 0.025))
  
  te  <-  data.frame(type = 'te',
                     estimate = mean((simu_result$exp_star/simu_result$star_star)*(simu_result$exp_exp/simu_result$exp_star)),
                     upper = quantile((simu_result$exp_star/simu_result$star_star)*(simu_result$exp_exp/simu_result$exp_star), 0.975),
                     lower = quantile((simu_result$exp_star/simu_result$star_star)*(simu_result$exp_exp/simu_result$exp_star), 0.025))
  
  result <- rbind(pnde, tnde, tnie, pnie, te) %>%
    mutate(exp_level = x_exp)
  
  return(result)
}


## Main body of the program ------
coef_simu_la <- coef_sim(data = combine_la_warm, outcome = 'RESPcnt', temperature = 'w_tmax', pollutant = 'q_O3', n = 20000)

cme_95_50 <- mf_simu(combine_la_warm, 'RESPcnt', 'w_tmax', 'q_O3', coef_simu_la, 0.95, 0.5)
cme_90_50 <- mf_simu(combine_la_warm, 'RESPcnt', 'w_tmax', 'q_O3', coef_simu_la, 0.90, 0.5)
cme_85_50 <- mf_simu(combine_la_warm, 'RESPcnt', 'w_tmax', 'q_O3', coef_simu_la, 0.85, 0.5)
cme_80_50 <- mf_simu(combine_la_warm, 'RESPcnt', 'w_tmax', 'q_O3', coef_simu_la, 0.80, 0.5)
cme_75_50 <- mf_simu(combine_la_warm, 'RESPcnt', 'w_tmax', 'q_O3', coef_simu_la, 0.75, 0.5)
cme_70_50 <- mf_simu(combine_la_warm, 'RESPcnt', 'w_tmax', 'q_O3', coef_simu_la, 0.70, 0.5)
cme_65_50 <- mf_simu(combine_la_warm, 'RESPcnt', 'w_tmax', 'q_O3', coef_simu_la, 0.65, 0.5)
cme_60_50 <- mf_simu(combine_la_warm, 'RESPcnt', 'w_tmax', 'q_O3', coef_simu_la, 0.60, 0.5)
cme_55_50 <- mf_simu(combine_la_warm, 'RESPcnt', 'w_tmax', 'q_O3', coef_simu_la, 0.55, 0.5)

CME_results <- rbind(cme_95_50, cme_90_50, cme_85_50, cme_80_50, cme_75_50, cme_70_50, cme_65_50, cme_60_50, cme_55_50) %>%
  setNames(c('type', 'estimate', 'upper', 'lower', 'exp_level'))


source("R/worm_burden_distn_analysis.R")
source("R/Models/Deterministic/deterministic_age_stratified_model.R")
source("R/Models/model_helper_functions.R")
source("R/Models/model_parameters.R")
source("R/Models/Deterministic/deterministic_fit_functions.R")

require(deSolve)

# look at data from Senegal to fit to #########
comm_haem_sums %>% mutate(year = as.numeric(year)) %>% 
  ggplot(aes(x = year, y = w_est, col = school)) +
    geom_line() + 
    geom_point() + 
    geom_errorbar(aes(x = year, ymin = w_est - w_se_est, ymax = w_est + w_se_est), width = 0.05) +
    theme_bw()

# Let's try and fit to DT since worm burden there steadily declines like we would expect it to, so won't have to do anything hacky in the model (hopefully) to get it to fit
DT_dat <- comm_haem_sums %>% filter(school == "DT")  
  w_obs_SAC = DT_dat$w_est 
  w_se_obs_SAC = DT_dat$w_se_est
  
  w_obs_adult <- w_obs_SAC * 0.5   # Assume mean worm burden in adult pop is half of that in SAC pop
  w_se_obs_adult <- w_se_obs_SAC   # Assume same se since we don;t know anything about adult burden and want likelihood of fit to SAC to dominate
  
  n_H <- DT_dat$samp_size[1] * 3   #Assume total population is twice the size of SAC population

  area = 200
  H = n_H             #Total population
  prop_SAC = 0.3       #Percent of people that are school age children (SAC) that we treated/measured
  prop_adult = 1-prop_SAC #percent of people that are not school age children (assumed here to be adults)
  
  cvrg_SAC = 0.8  # Assume we treated 80% of SAC
  cvrg_adult = 0
  
  params = pars1
  
# Get somewhat arbitrary snail parameters
set.seed(212)  
  S_beta <- get_snail_beta_pars(40*area, .7, .1)
  E_beta <- get_snail_beta_pars(40*area, .2, .05)
  I_beta <- get_snail_beta_pars(40*area, .1, .025)
  
  S_shape1 = S_beta[1]
  S_shape2 = S_beta[2]
  E_shape1 = E_beta[1]
  E_shape2 = E_beta[2]
  I_shape1 = I_beta[1]
  I_shape2 = I_beta[2]
  
#Create array of values to feed to sim_fit function
  test_pars <- expand.grid(beta = seq(1e-9, 1e-5, length.out = 20),
                           lambda = seq(1e-9, 1e-5, length.out = 20)) %>% 
    mutate(negLL = map2_dbl(beta, lambda, sim_fit_SAC))
  
  beta_use <- test_pars$beta[which(test_pars$negLL == min(test_pars$negLL))]
  lambda_use <- test_pars$lambda[which(test_pars$negLL == min(test_pars$negLL))]
  
  #Run to equilibrium with best fit parameters   
    pars_DT_fit <- pars1
    pars_DT_fit["beta"] <- beta_use
    pars_DT_fit["lambda"] <- lambda_use
    
save(test_pars, pars_DT_fit, 
     area, w_obs_SAC, w_obs_adult, 
     cvrg_SAC, cvrg_adult, 
     prop_SAC, prop_adult, H,
     comm_haem_sums,
     file = "R/Models/Deterministic/deterministic_DT_fit.RData")

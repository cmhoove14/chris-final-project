source("R/worm_burden_distn_analysis.R")
source("R/Models/deterministic_age_stratified_model.R")
source("R/Models/model_helper_functions.R")
source("R/Models/model_parameters.R")

require(deSolve)

# Run to make sure model works ########
test <- sim_schisto_age_stratified(nstart = setNames(c(30*area, 10*area, 1*area,
                                               30, 30,
                                               15, 15),
                                             c('S', 'E', 'I', 
                                               'Wt_SAC', 'Wu_SAC',
                                               'Wt_adult', 'Wu_adult')),
                                   time = seq(0, 365*30, 30),
                                   parameters = pars1,
                                   events_df = NA)


# Functions, parameters, and data to help with fitting ##########
# Likelihood function for a single data point
  pointLL <- function(par1, par2, obs, dist){
    if(dist == "gaussian"){
      LL <- log(dnorm(obs, par1, par2))
    } else if (dist == "beta"){
      LL <- log(dbeta(obs, par1, par2))
    } else {
      LL <- NA
    }
    
    return(LL)
    
  }
  
#Set demographic parameters that vary by location 
  area = 1000
  H = 1000             #Total population
  prop_SAC = 0.5       #Percent of people that are school age children (SAC)
  prop_adult = 1-prop_SAC #percent of people that are not school age children (assumed here to be adults)
  
  cvrg_SAC = 0.6
  cvrg_adult = 0
  
# Time over which to run intervention simulation 
  sim_time <- c(0:(365*3))
  
#Time to pull worm burden estimates from model as observations to compare to distributions
  check_times <- c(0, 365, 365*2)
  
#Events data frame representing yearly mda of SAC  
  yrly_mdas <- data.frame(var = rep('Wt_SAC', 2),
                          time = c(1, 366),
                          value = rep(0.06, 2),     #94% efficacy
                          method = rep('mult', 2))
    
# Estimate somewhat arbitrary snail parameters to fit to
# Beta distns for snail infection classes  
  get_snail_beta_pars <- function(N, prop, sd){
    snail_dist <- fitdist(rnorm(1000, N*prop, N*sd)/N, "beta")[1]$estimate
      snail_shape1 <- snail_dist[[1]]
      snail_shape2 <- snail_dist[[2]]
    
    return(c(snail_shape1, snail_shape2))  
  }
  
# Function to run model with test parameter set and return negative log likelihood ########
# Function takes a test beta and lambda (the snail-to-man and man-to-snail transmission parameters) and updates the test parameter set (params) with them. It also takes a bunch of observed data including the observed estimated worm burden in SAC and adults and the observed distribution of snail infection densities relative to the total snail population  

  sim_fit <- function(beta_lambda){
    
    params["beta"] <- beta_lambda[1]
    params["lambda"] <- beta_lambda[2]
    
  #Run to equilibrium with parameters   
    run_to_eqbm <- sim_schisto_age_stratified(nstart = setNames(c(40*area*.75, 40*area*.22, 40*area*.02,
                                                                  w_obs_SAC[1], w_obs_SAC[1],
                                                                  w_obs_adult[1], w_obs_adult[1]),
                                                                c('S', 'E', 'I', 
                                                                  'Wt_SAC', 'Wu_SAC',
                                                                  'Wt_adult', 'Wu_adult')),
                                              time = seq(0, 365*30, 30),
                                              parameters = params,
                                              events_df = NA)
    
  #Simulate over period matching data collection  
    run_sim <- sim_schisto_age_stratified(nstart = setNames(as.numeric(run_to_eqbm[dim(run_to_eqbm)[1], c(2:8)]),
                                                            c('S', 'E', 'I', 
                                                              'Wt_SAC', 'Wu_SAC',
                                                              'Wt_adult', 'Wu_adult')), 
                                          time = sim_time,
                                          parameters = params,
                                          events_df = yrly_mdas)
    
    run_sim$N <- run_sim$S + run_sim$E + run_sim$I
    
  #Get log likelihood for individual data points  
    w_SAC_test <- run_sim$Wt_SAC[which(run_sim$time %in% check_times)]
      w_SAC_LL <- mapply(pointLL, par1 = w_obs_SAC, par2 = w_se_obs_SAC, obs = w_SAC_test, MoreArgs = list(dist = "gaussian"))
      
    w_adult_test <- run_sim$Wt_adult[which(run_sim$time %in% check_times)]
      w_adult_LL <- mapply(pointLL, par1 = w_obs_adult, par2 = w_se_obs_adult, obs = w_adult_test, MoreArgs = list(dist = "gaussian"))

    S_test <- run_sim$S[which(run_sim$time %in% check_times)] / run_sim$N[which(run_sim$time %in% check_times)]
      S_LL <- mapply(pointLL, par1 = S_shape1, par2 = S_shape2, obs = S_test, MoreArgs = list(dist = "beta"))
      
    E_test <- run_sim$E[which(run_sim$time %in% check_times)] / run_sim$N[which(run_sim$time %in% check_times)]
      E_LL <- mapply(pointLL, par1 = E_shape1, par2 = E_shape2, obs = E_test, MoreArgs = list(dist = "beta"))
      
    I_test <- run_sim$I[which(run_sim$time %in% check_times)] / run_sim$N[which(run_sim$time %in% check_times)]
      I_LL <- mapply(pointLL, par1 = I_shape1, par2 = I_shape2, obs = I_test, MoreArgs = list(dist = "beta"))
      
    all_negLL <- c(-w_SAC_LL, -w_adult_LL, -S_LL, -E_LL, -I_LL)
    
    all_negLL[!is.finite(all_negLL)] <- NA
      
    negLL <- sum(all_negLL, na.rm = T)
    LL_NAs <- sum(is.na(all_negLL))
    
    return(negLL)
  }

#Test fit function 
set.seed(212)  
  S_beta <- get_snail_beta_pars(40*area, .7, .1)
  E_beta <- get_snail_beta_pars(40*area, .2, .05)
  I_beta <- get_snail_beta_pars(40*area, .1, .025)
  
  params = pars1
  w_obs_SAC = 30 
  w_se_obs_SAC = 6
  w_obs_adult = 15
  w_se_obs_adult = 10
  S_shape1 = S_beta[1]
  S_shape2 = S_beta[2]
  E_shape1 = E_beta[1]
  E_shape2 = E_beta[2]
  I_shape1 = I_beta[1]
  I_shape2 = I_beta[2]
  
  sim_fit(beta_lambda = c(1e-5, 1e-7))
  
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
  n_H <- DT_dat$samp_size[1] * 2   #Assume total population is twice the size of SAC population

  area = 200
  H = n_H             #Total population
  prop_SAC = 0.5       #Percent of people that are school age children (SAC)
  prop_adult = 1-prop_SAC #percent of people that are not school age children (assumed here to be adults)
  
  cvrg_SAC = 0.6
  cvrg_adult = 0

  #opt_pars <- optim(pars = c(1e-5, 1e-7), sim_fit, lower = c(1e-9, 1e-9), upper = c(1e-3, 1e-3), trace = 1, method = "L-BFGS-B")
  
#Create array of values to feed to sim_fit function
  test_pars <- expand.grid(beta = seq(1e-9, 1e-6, length.out = 20),
                           lambda = seq(1e-9, 1e-5, length.out = 20))
  
  for(i in 1:nrow(test_pars)){
    test_pars[i,3] <- sim_fit(c(test_pars[i,1], test_pars[i,2]))
    print(c(i, test_pars[i,1], test_pars[i,2], test_pars[i,3]))
  }
  
  test_pars %>% rename(negLL = V3) %>% 
    ggplot(aes(x = beta, y = lambda, fill = negLL)) +
    geom_tile() + scale_fill_gradient(low = "white", high = "black")
  
  beta_use <- test_pars$beta[which(test_pars$V3 == min(test_pars$V3))]
  lambda_use <- test_pars$lambda[which(test_pars$V3 == min(test_pars$V3))]
  
  #Run to equilibrium with best fit parameters   
    pars1["beta"] <- beta_use
    pars1["lambda"] <- lambda_use

    best_fit_eqbm <- sim_schisto_age_stratified(nstart = setNames(c(40*area*.75, 40*area*.22, 40*area*.02,
                                                                  w_obs_SAC[1], w_obs_SAC[1],
                                                                  w_obs_adult[1], w_obs_adult[1]),
                                                                c('S', 'E', 'I', 
                                                                  'Wt_SAC', 'Wu_SAC',
                                                                  'Wt_adult', 'Wu_adult')),
                                              time = seq(0, 365*30, 30),
                                              parameters = pars1,
                                              events_df = NA)
    
  #Simulate over period matching data collection  
    best_fit_sim <- sim_schisto_age_stratified(nstart = setNames(as.numeric(best_fit_eqbm[dim(best_fit_eqbm)[1], c(2:8)]),
                                                            c('S', 'E', 'I', 
                                                              'Wt_SAC', 'Wu_SAC',
                                                              'Wt_adult', 'Wu_adult')), 
                                          time = sim_time,
                                          parameters = pars1,
                                          events_df = yrly_mdas)
    
    best_fit_sim <- best_fit_sim %>% 
      mutate(N = S+E+I,
             S_prop = S/N,
             E_prop = E/N,
             I_prop = I/N,
             time_year = (time/365)+2016)

comm_haem_sums %>% mutate(year = as.numeric(year)) %>% filter(school == "DT") %>% 
  ggplot(aes(x = year, y = w_est)) +
    geom_line() + 
    geom_point() + 
    geom_errorbar(aes(x = year, ymin = w_est - w_se_est, ymax = w_est + w_se_est), width = 0.05) +
    theme_bw() + ylim(c(0,40)) + xlim(c(2015.5, 2018.1)) +
    geom_line(data = best_fit_sim, aes(x = time_year, y = Wt_SAC))

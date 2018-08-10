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

#This function only fits to worm burden in SAC  
  sim_fit_SAC <- function(beta, lambda){
    
    params["beta"] <- beta
    params["lambda"] <- lambda
    
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
      
    w_SAC_LL[!is.finite(w_SAC_LL)] <- NA
      
    return(-sum(w_SAC_LL))
  }

#This function fits to all data (snail, SAC, adult)  
  sim_fit_all <- function(beta, lambda){
    
    params["beta"] <- beta
    params["lambda"] <- lambda
    
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

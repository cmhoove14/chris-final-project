#Load packages and other files ##########
require(adaptivetau)
require(deSolve)

#adaptivetau model ############
transitions = list(
  c(S = 1),             #New (susceptible) snail born
  c(S = -1),            #Susceptible snail dies
  c(S = -1, E = 1),     #Susceptible snail becomes exposed
  c(E = -1),            #Exposed snail dies
  c(E = -1, I = 1),     #Exposed snail becomes Infected
  c(I = -1),            #Infected snail dies
  c(Wt_SAC = 1),        #Adult worm in the treated SAC population
  c(Wt_SAC = -1),       #Adult worm in the treated SAC population dies
  c(Wu_SAC = 1),        #Adult worm in the untreated SAC population
  c(Wu_SAC = -1),        #Adult worm in the untreated SAC population dies
  c(Wt_adult = 1),      #Adult worm in the treated adult population
  c(Wt_adult = -1),      #Adult worm in the treated adult population dies
  c(Wu_adult = 1),      #Adult worm in the untreated adult population
  c(Wu_adult = -1))      #Adult worm in the untreated adult population dies

sfx <- function(x, p, t) {
  S = x['S']
  E = x['E']
  I = x['I']
  
  N = S + E + I
  
  Wt_SAC = x['Wt_SAC']
  Wu_SAC = x['Wu_SAC']
  Wt_adult = x['Wt_adult']
  Wu_adult = x['Wu_adult']
  
  W_SAC = (cvrg_SAC*Wt_SAC) + ((1-cvrg_SAC)*Wu_SAC)           #weighting treated and untreated populations among SAC
  W_adult = (cvrg_adult*Wt_adult) + ((1-cvrg_adult)*Wu_adult) #weighting treated and untreated populations among adults
    
  W_tot = W_SAC*prop_SAC + W_adult*prop_adult
  
    #Update clumping parameter, k from worm burden-dispersion relationship
      k1 = k_from_log_w(Wt_SAC)
      k2 = k_from_log_w(Wu_SAC)
      k3 = k_from_log_w(Wt_adult)
      k4 = k_from_log_w(Wu_adult)
    
    #Estimate mating probability within each strata 
      phi_Wt_SAC = phi_Wk(W = Wt_SAC, k = k1)  #Mating probability in treated SAC population
      phi_Wu_SAC = phi_Wk(W = Wu_SAC, k = k2)  #Mating probability in untreated SAC population
      phi_Wt_adult = phi_Wk(W = Wt_adult, k = k3)  #Mating probability in treated adult population
      phi_Wu_adult = phi_Wk(W = Wu_adult, k = k4)  #Mating probability in untreated adult population
      
    #Estimate number of mated female worms and their relative fecundity in each population
      fem_Wt_SAC = 0.5*(Wt_SAC*phi_Wt_SAC) * f_Wgk(Wt_SAC, gam, k1)
      
      fem_Wu_SAC = 0.5*(Wu_SAC*phi_Wu_SAC) * f_Wgk(Wu_SAC, gam, k2)
      
      fem_Wt_adult = 0.5*(Wt_adult*phi_Wt_adult) * f_Wgk(Wt_adult, gam, k3) 
      
      fem_Wu_adult = 0.5*(Wu_adult*phi_Wu_adult) * f_Wgk(Wu_adult, gam, k4)
      
      #multiplied by eggs per 10mL urine, mean urine produced by an individual in each population, number of inidivuals in the population, viability of eggs shed, and contamination factor (probability eggs shed end up in snail environment) below and summed across all populations to get miracidia (M) estimate used in estimating force of infection between humans and snails
 
  return(c(p$f_N * (1-N/p$C) * (S + E),   #Snail birth
           p$mu_N * S,                    #Susceptible snail death
           p$beta * (1-exp(-(sum(fem_Wt_SAC * p$m * ((p$H*p$prop_SAC)*p$cvrg_SAC) * p$v * p$u_SAC * p$rho_SAC,
                                 fem_Wu_SAC * p$m * ((p$H*p$prop_SAC)*(1-p$cvrg_SAC)) * p$v * p$u_SAC * p$rho_SAC,
                                 fem_Wt_adult * p$m * ((p$H*p$prop_adult)*p$cvrg_adult) * p$v * p$u_adult * p$rho_adult,
                                 fem_Wu_adult * p$m * ((p$H*p$prop_adult)*(1-p$cvrg_adult)) * p$v * p$u_adult * p$rho_adult))/N)) * S,  #Snail exposure
           p$mu_N * E,                    #Exposed snail dies
           p$sigma * E,                   #Exposed snail becomes infected
           (p$mu_N + p$mu_I) * I,         #Infected snail dies
           
           p$omega_SAC * p$lambda * I,    #infected snail produces adult worm in treated SAC population
           p$mu_W * Wt_SAC,               #Adult worm in treated SAC population dies
           
           p$omega_SAC * p$lambda * I,    #infected snail produces adult worm in untreated SAC population
           p$mu_W * Wu_SAC,               #Adult worm in untreated SAC population dies
           
           p$omega_adult * p$lambda * I,  #infected snail produces adult worm in treated adult population
           p$mu_W * Wt_adult,             #Adult worm in treated adult population dies
           
           p$omega_adult * p$lambda * I,  #infected snail produces adult worm in untreated adult population
           p$mu_W * Wu_adult))            #Adult worm in untreated adult population dies
           
}

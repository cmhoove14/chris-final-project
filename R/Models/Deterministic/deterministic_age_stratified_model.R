#Model with PDD ####################
schisto_age_stratified=function(t, n, parameters) { 
  with(as.list(parameters),{
    
    S=n[1]
    E=n[2]
    I=n[3]
    Wt_SAC=n[4]
    Wu_SAC=n[5]
    Wt_adult=n[6]
    Wu_adult=n[7]
    
    N=S+E+I
    
    W_SAC = (cvrg_SAC*Wt_SAC) + ((1-cvrg_SAC)*Wu_SAC)           #weighting treated and untreated populations among SAC
    W_adult = (cvrg_adult*Wt_adult) + ((1-cvrg_adult)*Wu_adult) #weighting treated and untreated populations among adults
    
    W_tot = W_SAC*prop_SAC + W_adult*prop_adult
    
    #Update clumping parameter, k from estimate of eggs burden per 10mL estimate
      k1 = k_from_log_w(Wt_SAC)
      k2 = k_from_log_w(Wu_SAC)
      k3 = k_from_log_w(Wt_adult)
      k4 = k_from_log_w(Wu_adult)
    
      #print(c(k1, k2, k3, k4))
    
    #Estimate mating probability within each strata 
      phi_Wt_SAC = phi_Wk(W = Wt_SAC, k = k1)  #Mating probability in treated SAC population
      phi_Wu_SAC = phi_Wk(W = Wu_SAC, k = k2)  #Mating probability in untreated SAC population
      phi_Wt_adult = phi_Wk(W = Wt_adult, k = k3)  #Mating probability in treated adult population
      phi_Wu_adult = phi_Wk(W = Wu_adult, k = k4)  #Mating probability in untreated adult population
      
      #print(c(phi_Wt_SAC, phi_Wu_SAC, phi_Wt_adult, phi_Wu_adult))
      
    #Estimate mean eggs produced per person in each strata as product of mated female worms, eggs produced per female worm per 10mL urine, and reduction in fecundity due to crowding
      eggs_Wt_SAC = 0.5*(Wt_SAC*phi_Wt_SAC) * m * f_Wgk(Wt_SAC, gam, k1)
      
      eggs_Wu_SAC = 0.5*(Wu_SAC*phi_Wu_SAC) * m * f_Wgk(Wu_SAC, gam, k2)
      
      eggs_Wt_adult = 0.5*(Wt_adult*phi_Wt_adult) * m * f_Wgk(Wt_adult, gam, k3) 
      
      eggs_Wu_adult = 0.5*(Wu_adult*phi_Wu_adult) * m * f_Wgk(Wu_adult, gam, k4)
 
      #print(c(eggs_Wt_SAC, eggs_Wu_SAC, eggs_Wt_adult, eggs_Wu_adult))
      
    #Estimate miracidia produced by each strata as product of mean eggs per 10 mL urine for individuals in each strata, mean mL urine produced by an average individual in each group/10, number of people in each strata, contamination coefficient for SAC/adults, and egg viability  
      M_Wt_SAC = eggs_Wt_SAC * ((H*prop_SAC)*cvrg_SAC) * v * u_SAC * rho_SAC
      
      M_Wu_SAC = eggs_Wu_SAC * ((H*prop_SAC)*(1-cvrg_SAC)) * v * u_SAC * rho_SAC
      
      M_Wt_adult = eggs_Wt_adult * ((H*prop_adult)*cvrg_adult) * v * u_adult * rho_adult
      
      M_Wu_adult =  eggs_Wu_adult * ((H*prop_adult)*(1-cvrg_adult)) * v * u_adult * rho_adult
       
    #Total miracidia snails are exposed to  
      M = M_Wt_SAC + M_Wu_SAC + M_Wt_adult + M_Wu_adult
      
    #Snail infection dynamics
      dSdt= f_N*(1-(N/C))*(S+E) - mu_N*S - beta*(1-exp(-M/N))*S #Susceptible snails
      
      dEdt= beta*M*S - (mu_N+sigma)*E #Exposed snails
      
      dIdt= sigma*E - (mu_N+mu_I)*I #Infected snails
    
    #worm burden in human populations
      dWt_SACdt= (omega_SAC*lambda*I) - (mu_W*Wt_SAC)
      dWu_SACdt= (omega_SAC*lambda*I) - (mu_W*Wu_SAC)
      dWt_adultdt= (omega_adult*lambda*I) - (mu_W*Wt_adult)
      dWu_adultdt= (omega_adult*lambda*I) - (mu_W*Wu_adult)
    
    return(list(c(dSdt,dEdt,dIdt,
                  dWt_SACdt,dWu_SACdt,
                  dWt_adultdt, dWu_adultdt)))
  }) 
} 

sim_schisto_age_stratified <- function(nstart, time, parameters, events_df = NA){
  
  if(is.na(events_df)){
    out <- ode(nstart, time, schisto_age_stratified, parameters)
  } else {
    out <- ode(nstart, time, schisto_age_stratified, parameters, events = list(data = events_df))
  }
  
  return(as.data.frame(out))
}

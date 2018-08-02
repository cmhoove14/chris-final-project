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
    
    W_SAC = (cvrg_SAC*Wt_SAC) + ((1-cvrg_SAC)*Wu_SAC)           #weighting treated and untreated populations
    W_adult = (cvrg_adult*Wt_adult) + ((1-cvrg_adult)*Wu_adult) #weighting treated and untreated populations
    
    W_tot = W_SAC*prop_SAC + W_adult*prop_adult
    
    #Estimate clumping parameter, k from worm burden data, this function needs to be developed
      #k_Wt_SAC = get_clump(Wt_SAC)
      #k_Wu_SAC = get_clump(Wu_SAC)
      #k_Wt_adult = get_clump(Wt_adult)
      #k_Wu_adult = get_clump(Wu_adult)
    
    #Estimate mating probability within each strata 
      phi_Wt_SAC = phi_Wk(W = Wt_SAC, k = k)  #Mating probability in treated SAC population
      phi_Wu_SAC = phi_Wk(W = Wu_SAC, k = k)  #Mating probability in untreated SAC population
      phi_Wt_adult = phi_Wk(W = Wt_adult, k = k)  #Mating probability in treated adult population
      phi_Wu_adult = phi_Wk(W = Wu_adult, k = k)  #Mating probability in untreated adult population

    #Estimate miracidia produced by each strata as function of mated female worms, eggs produced per female worm per 10mL urine, reduction in fecundity due to crowding, mL urine produced in each group/10, contamination coefficient for SAC/adults and egg viability  
      M_Wt_SAC = fem_worms(W = Wt_SAC,
                           phi = phi_Wt_SAC,
                           H = H,
                           prop = prop_SAC,
                           cvrg = cvrg_SAC) * m * f_Wgk(Wt_SAC, gam, k) * v * u_SAC * rho_SAC
      
      M_Wu_SAC = fem_worms(W = Wu_SAC,
                           phi = phi_Wu_SAC,
                           H = H,
                           prop = prop_SAC,
                           cvrg = (1-cvrg_SAC)) * m * f_Wgk(Wu_SAC, gam, k) * v * u_SAC * rho_SAC
      
      M_Wt_adult = fem_worms(W = Wt_adult,
                             phi = phi_Wt_adult,
                             H = H,
                             prop = prop_adult,
                             cvrg = cvrg_adult) * m * f_Wgk(Wt_adult, gam, k) * v * u_adult * rho_adult
      
      M_Wu_adult = fem_worms(W = Wu_adult,
                             phi = phi_Wu_adult,
                             H = H,
                             prop = prop_adult,
                             cvrg = (1-cvrg_adult)) * m * f_Wgk(Wu_adult, gam, k) * v * u_adult * rho_adult
    
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

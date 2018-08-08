#Small functions used in schisto modeling 

#mating function that implements PDD#######
  phi_Wk <- function(W, k) {
    if( k<= 0){
      return(1)
    } else {
      func <- function(x) {
      a <- ( W / (W + k) )
      b <- ((1-a)^(1+k))/(2*pi)
      return(( b*( 1-cos(x) ) / (( 1 + a*cos(x) )^(1+k)) ))
    }
    val = integrate(func, 0, 2*pi, subdivisions = 10000,
                    rel.tol = 1e-10, stop.on.error = FALSE)$value
    return(1-val)
    }
    
  }

#Parasite fecundity as a function of worm burden and clumping parameter #######
#represents reduced fecundity at high worm burdens due to crowding   
  f_Wgk <- function(W, gamma, k) {
    if(k <= 0){
      return(1)
    } else {
      return((1 + ((W*(1-(exp(-gamma))))/k))^(-k-1))
    }
    
  }


#Get number of female worms from mean worm burden in target population, mating probability, total number of people, proportion of people in target population, and coverage within target population (i.e. proportion in target population that are treated/not treated)
fem_worms <- function(W, phi, H, prop, cvrg){
  0.5*(W*phi)*((H*prop)*cvrg)
}

#Get clumping parameter from mean egg burden and prevalence
k_from_prev_M <- function(M, prev){
  uniroot(function(x) 1-prev-(1+M/x)^-x,
          interval = c(0,3))$`root`
}

#Get estimate of mean worm burden from mean egg burden taking into account assumption of mean eggs per mated female worm per 10mL urine, shape of mating probability function, and shape of female worm fecundity function  
w_from_m_k <- function(M, k){
  uniroot(function(w) 0.5*w*as.numeric(pars1["m"])*
            (1 + (w*(1-(exp(-as.numeric(pars1["gam"])))))/k)^(-k-1)*
            phi_Wk(w, k) - M,
          interval = c(0, 1000))$`root`
}

#Get clumping parameter from mean worm burden estimate
#this function estimated from our s. haematobium data in Senegal in script "worm_burden_distn_analysis.R"
#This is from log-linear relationship with mean worm burden
k_from_log_w <- function(w){
  0.03362+0.09033*log(w)
}


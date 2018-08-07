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

#Get clumping parameter from mean worm burden and prevalence
k_from_prev_W <- function(W, prev){
  uniroot(function(x) 1-prev-(1+W/x)^-x,
          interval = c(0,3))$`root`
}

#Get clumping parameter from egg output measured in egg/10mL urine
#this function estimated from our s. haematobium burden in Senegal in script "worm_burden_distn_analysis.R"
#This is from log-log relationship
k_from_log_w <- function(w){
  exp(-2.8079+0.39*log(50))
}

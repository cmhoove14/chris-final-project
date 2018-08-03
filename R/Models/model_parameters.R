#Set demographic parameters that vary by location ########
  area = 1000
  H = 1000             #Total population
  prop_SAC = 0.3       #Percent of people that are school age children (SAC)
  prop_adult = 1-prop_SAC #percent of people that are not school age children (assumed here to be adults)
  
  cvrg_SAC = 1
  cvrg_adult = 1

#Parameters used in elimination feasibility work #####################
pars1=c( 
  ##standard snail parameters 
    f_N=0.10, # recruitment rate (from sokolow et al)
    C=50*area, # carrying capacity corresponding to 50 snails per square meter
    mu_N=1/60, #Mean mortality rate of snails (assuming lifespan of 60days)
    sigma=1/40, #Transition rate from exposed to infected (assuming pre-patency period of 40 days)
    mu_I=1/10 - 1/60, # Increased mortality rate of infected snails
  
  #Adult Worm, Miracidia and Circariae Parameters
    mu_W = 1/(3.3*365), # death rate of adult worms
    m = 3.6,              # mean eggs shed per female worm per 10mL urine
    v = 0.08,           # mean egg viability of eggs shed into environment  
    gam = 0.001,        # parameter of fecundity reduction function 

  #Human parameters
    k_SAC=0.08, # clumping parameter of the negative binomial distribution in SAC
    k_adult=0.08,  # clumping parameter of the negative binomial distribution in adults
    u_SAC = 50, # mL urine per SAC/day/10mL assumed to be half of adult
    u_adult = 100, # mL urine per SAC/day/10mL (approximate, ranges from 80 - 200)
    rho_SAC = 0.4, # relative number of eggs shed by SAC that make it to snail habitat (~sanitation)
    rho_adult = 0.1, # relative number of eggs shed by adults that make it to snail habitat (~sanitation)
    omega_SAC = 1,     # relative infection risk of SAC (related to clean water access/education/water contact)
    omega_adult = 0.5,   # relative infection risk of adults (related to clean water access/education)
  
  #Transmission parameters
    lambda=1.2e-4, # snail-to-man transmission: p(infected snail sheds cercariae that infects human and matures into adult)
    beta=1.6e-6  # man-to-snail transmission: p(mated female worm produces a miracidia that infects a snail)
  
)
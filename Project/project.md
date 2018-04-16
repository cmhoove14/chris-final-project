Optimal control of neglected tropical diseases using a simple dynamic model
================
Chris Hoover
04-15-2018

------------------------------------------------------------------------

Neglected tropical diseases (NTDs) tend to be diseases of poverty that thrive in tropical, disadvantaged communities that often lack access to basic amenities such as clean water, sanitation, and healthcare. They are generally environmentally mediated meaning that some component of their transmission cycle depends on the environment (e.g. a vector, intermediate host, or free-living environmental stage). Control of NTDs generally relies on the administration of drugs that cure individuals from the disease and may confer short-term immunity. However, following drug administration, people are often reinfected by the environmental stage(s) of the pathogen that persist through mass drug administration (MDA) campaigns. Additional interventions that target these environmental reservoirs are often available, but underutilized due to their indirect effects on improving individuals' health.

Using a simple, generalizable model of NTD transmission presented in [Garchitorena et al](http://rstb.royalsocietypublishing.org/content/372/1722/20160128), two interventions will be considered: 1) drug administration, implemented as a pulse reduction in the state variable, *I*, that reduces the prevalence of the disease in the population and 2) environmental remediation (e.g. improvement in sanitation, vector control), implemented as a permanent alteration of a model parameter, that reduces the transmission of the disease.

Part 1: Reproduce Fig3c
-----------------------

The first goal is to code the model and reproduce Figure 3c (below) from [Garchitorena et al](http://rstb.royalsocietypublishing.org/content/372/1722/20160128) comparing the effects of interventions that only target infected people, such as mass drug administration (MDA), to the effects of combined interventions that target infected people and the environment, such as sanitation, mosquito spraying, mollusciciding, etc.

The model consists of two state variables, *I* and *W*, that correspond to the prevalence of infection in the human population (e.i. proportion infected at time=*t*) and the degree of contamination of the environment with the disease-causing agent, respectively. We make the simplifying assumption that individuals can only be susceptible, *S*, or infected, *I*, meaning *S* + *I* = 1 and eliminating the need for a recovered *R* compartment as is typical of SIR models but would complicate things here. The model equations then are:

$$\\frac{dI}{dt}=(\\beta\_EW+\\beta\_DI)(1-I)-\\gamma I $$

$$\\frac{dW}{dt}=\\Omega+V\\sigma\\lambda I-\\rho W $$

with parameter definitions and values in table 1 below

|     Parameter     |   Value  |                            Description                            |
|:-----------------:|:--------:|:-----------------------------------------------------------------:|
| *β*<sub>*E*</sub> | 9.26e-05 |       Transmission rate from environment to human population      |
| *β*<sub>*D*</sub> |  0.00555 |                  Human to human transmission rate                 |
|        *γ*        |   0.003  |         Rate of recovery from infected back to susceptible        |
|        *Ω*        |     0    |      Recruitment rate of infectious agents in the environment     |
|        *V*        |     1    |    Abundance of vectors/intermediate hosts/suitable environment   |
|        *λ*        |     1    |  Recruitment rate of infectious agents by infectious individuals  |
|        *σ*        |     1    | Fraction of infectious agents produced that reach the environment |
|  *ρ*<sup>\*</sup> |   0.017  |       Mortality rate of infectious agents in the environment      |

We make the simplifying assumptions that there is no exogenous production of infectious agents (*Ω* = 0) and that environment-human transmission is determined only by the transmission rate, *β*<sub>*E*</sub>, with no role of vectors or intermediate hosts (i.e. *V* = *σ* = *λ* = 1). This is approximately equivalent to a simple Cholera model with the form:

$$\\frac{dI}{dt}=(\\beta\_EW+\\beta\_DI)(1-I)-\\gamma I $$

$$\\frac{dW}{dt}=I-\\rho W $$

With this simplified model, we can get the equilibrium estimates of *W* and *I*:
$$\\frac{dW}{dt}=I-\\rho W$$
 At equilibrium:
$$W^\*=\\frac{I}{\\rho}$$
 Substituting we get the following quadratic with solutions at *I* = 0 and *I* = endemic equilibrium:
$$\\frac{dI}{dt}=\\Big(\\frac{\\beta\_EI}{\\rho}+\\beta\_DI\\Big)\\Big(1-I\\Big)-\\gamma I $$

Now we use these equilibrium estimates as starting values for the continuous time model and simulate ten years of interventions with the `deSolve` package to reproduce Fig 3c

``` r
#Parameter values from table 1 (some excluded) drawn from Fig S1 of the manuscript
pars <- c(beta_e = 9.26e-5, beta_d = 5.55e-3,
          gamma = 1/360, rho = 1/60)

#Function to get equilibrium prevalence (I) and equilibrium environmental contamination (W)
get_eq <- function(parameters){
  with(as.list(parameters),{
    init <- rootSolve::uniroot.all(function(I) ((beta_e*I / rho) + beta_d*I)*(1-I) - gamma*I,
                                   interval = c(0,1))
    eq_I <- init[which(init > 0)]
    eq_W <- eq_I/rho
    
    return(c(I = eq_I, W = eq_W))
  })  
}

#Function to run model `t` days, with starting conditions `n` and parameter set `p`
simp_mod = function(t, n, p){
  with(as.list(p),{
    I = n[1]
    W = n[2]
    
    dIdt = (beta_e*W + beta_d*I)*(1-I) - gamma*I
    dWdt = I - rho*W
    
    return(list(c(dIdt, dWdt)))
  })
} 
  
#Events representing interventions to implement in the model
n_years <- 10       #Simulate 10 years of annual intervention
drug_efficacy = 0.9 #90% reduction in prevalence at each drug treatment
run_time <- c(1:(365*n_years)) #Run model for 10 years

#data frame with event times for drug administration
drugs <- data.frame(var=rep('I', times = n_years),
                        time = c(0:(n_years-1))*365 + 1,
                        value = rep((1-drug_efficacy), times = n_years),
                        method = rep("mult", times = n_years))
    
#Run model with annual drug administration intervention
drug_only <- as.data.frame(ode(get_eq(pars), run_time, simp_mod, pars,
                               events = list(data = drugs)))

#Double mortality rate of environmental infectious agents for environmental intervention
pars_env <- pars
pars_env["rho"] <- pars["rho"] * 2 

#Run model with annual drug administration and environmental intervention 
drug_env <- as.data.frame(ode(get_eq(pars), run_time, simp_mod, pars_env,
                               events = list(data = drugs)))

#Plot results
rbind(drug_env, drug_only) %>% 
  mutate(Intervention = c(rep("Drug + Env", length(run_time)),
                          rep("Drug", length(run_time)))) %>% 
  ggplot(aes(x = time/365, y = I, lty = Intervention)) + 
    theme_bw() + theme(legend.position = c(0.85,0.85)) +
    geom_line(size = 0.75) + xlab("time (years)") + ylab("Prevalence") + ylim(c(0,0.75))
```

![](project_files/figure-markdown_github/mod_prep-1.png)

Looks a little bit different from Fig 3c in that there's a higher starting prevalence and more rebound between drug treatments in the `Drug + Env` group. This is likely due to transmission parameters that are too high, so let's check that out.

``` r
#Function to get R0 estimate given model parameters
get_r0 <- function(parameters){
  with(as.list(parameters),{
    r0_d <- beta_d / gamma
    r0_e <- beta_e / (gamma*rho)
    
    r0_d + r0_e
  })  
}

get_r0(pars)
```

    ## [1] 3.99816

So *R*<sub>0</sub> = 4 with the parameters pulled from the SI, but Fig3c says *R*<sub>0</sub> = 3. Let's see if we can bring the transmission parameters down a bit and get closer to Fig3c.

``` r
r0 <- 3
pars2 <- pars

#Use equations from Fig S1 to estimate transmission parameters given desired r0 of 3
pars2["beta_d"] <- pars["gamma"]*r0 / 2
pars2["beta_e"] <- pars["gamma"]*r0*pars["rho"] / 2

get_r0(pars2)
```

    ## [1] 3

Now that we have the same *R*<sub>0</sub>, let's try the figure again

``` r
# Rerun model with each intervention scenario and new parameter sets
drug_only2 <- as.data.frame(ode(get_eq(pars2), run_time, simp_mod, pars2,
                                events = list(data = drugs)))

pars_env2 <- pars2
pars_env2["rho"] <- pars2["rho"] * 2 

drug_env2 <- as.data.frame(ode(get_eq(pars2), run_time, simp_mod, pars_env2,
                               events = list(data = drugs)))

#Plot results 
rbind(drug_env2, drug_only2) %>% 
  mutate(intervention = c(rep("Drug + Env", length(run_time)),
                          rep("Drug", length(run_time)))) %>% 
  ggplot(aes(x = time/365, y = I, lty = intervention)) + 
    theme_bw() + theme(legend.position = c(0.85,0.85)) +
    geom_line(size = 0.75) + xlab("time (years)") + ylab("Prevalence") + ylim(c(0,0.75))
```

![](project_files/figure-markdown_github/unnamed-chunk-3-1.png)

Looks about right!

Part 2: Translate the model into a Markov Decision Process (MDP) solvable with stochastic dynamic programming (SDP)
-------------------------------------------------------------------------------------------------------------------

Next we want to translate our system of continuous time differential equations into a Markov Decision Process (MDP) consisting of **1)** a Markov chain in which the state of the system at time *t* + 1 is dependent only on the current state of the system (i.e. the state at *t*) and **2)** a decision or control action that is being made at each state transition (i.e. from *t* to *t* + 1). We therefore want a single, discrete time equation that captures about the same dynamics of our simple disease system.

Next, we work through the "Six steps of stochastic dynamic programming" as described in [Marescot et al](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12082)

1.  **Define the optimization objective of the problem.**
    In order to place infection and interventions on the same playing field (i.e. share the same units), we frame the problem as minimizing the costs, *C*, associated with 1) infection of individuals and 2) implementing interventions to reduce their infections:
    $$\\min\_{I,A}C=\\min\_{I,A}\\sum\_0^T{\\Pi(I\_t,A\_t)}$$
     where
    *Π*(*I*<sub>*t*</sub>, *A*<sub>*t*</sub>)=*d**I*<sub>*t*</sub> + *M**A*<sub>*t*</sub> + *M*(1 − *A*<sub>*t*</sub>)
     *d* is the cost associated with having prevalence of *I*<sub>*t*</sub>, *M* is the total capital available to allocate towards intervening, and *A*<sub>*t*</sub> is the proportion of that capital allocated towards drug administration (with the rest allocated towards environmental intervention).

*But *A*<sub>*t*</sub> isn't part of what we're trying to minimize, we're trying to minimize the total money spent over the time horizon (*T*). Could make two control variables, say *A*<sub>*t*</sub> and *B*<sub>*t*</sub> that represent capital allocated towards drug administration and capital allocated towards environmental intervention, respectively, but then we have to change the decision variable to be *A*<sub>*t*</sub> and *B*<sub>*t*</sub>, perhaps with a maximum of *M*. Then could include "leftover capital", *M* − (*A*<sub>*t*</sub> + *B*<sub>*t*</sub>), in the utility function somehow?*

1.  **Define the set of states that represent the configuration of the system at time *t***
    We define the state variable, *X*<sub>*t*</sub>, as the prevalence of infection in the human population, *I*<sub>*t*</sub> ∈ \[0, 1\]. Next, we assume that the dynamics of the infectious agents in the environment are faster than the dynamics of the prevalence in the human population, therefore they reach a steady state equilibrium:
    $$W^\*=\\frac{I}{\\rho}$$
     Again, substituting we get:
    $$\\frac{dI}{dt}=\\Big(\\frac{\\beta\_EI}{\\rho}+\\beta\_DI\\Big)\\Big(1-I\\Big)-\\gamma I $$
     which can be translated to a discrete time model as:
    $$I\_{t+1}=\\Big(\\frac{\\beta\_E(I\_t-MA\_t\\theta)}{\\rho+M(1-A\_t)\\mu}+\\beta\_D(I\_t-MA\_t\\theta)\\Big)\\Big(1-(I\_t-MA\_t\\theta)\\Big)-\\gamma(I\_t-MA\_t\\theta) $$
     Where *θ* is a constant that converts capital to units of prevalence and *μ* is a constant that converts capital to the same units as the mortality rate of infectious agents in the environment

2.  **Define the decision variable, *A*<sub>*t*</sub> that is the component of the system to be controlled to meet the objective**
    We define the decision variable *A*<sub>*t*</sub> as the proportion of capital, *M*, committed to drug administration

3.  **Build a transition model describing the system's dynamics as a function of the decision variable and the system state in the prior time step**
    See **Step 2**

4.  **Define the utility function *U*<sub>*t*</sub>(*X*<sub>*t*</sub>, *A*<sub>*t*</sub>) representing the desirability of acting in a given state**
    *TO DO: derive utility function, ideally incorporating terminal reward associated with permanent increase in mortality rate of infectious agents in the environment*

5.  **Determine the optimal solution of the optimization problem**

### Notes and random thoughts

-   If we want to incorporate permanent reductions in model parameters associated with non-MDA interventions, have to sum over the time horizon. i.e.:
    $$I\_{t+1}=\\Big(\\frac{\\beta\_E(I\_t-MA\_t\\theta)}{\\rho+\\sum\_0^TM(1-A\_t)\\mu}+\\beta\_D(I\_t-MA\_t\\theta)\\Big)\\Big(1-(I\_t-MA\_t\\theta)\\Big)-\\gamma(I\_t-MA\_t\\theta) $$
-   How to incorporate a discount rate (*δ*)? General idea is that you lose some value for waiting to spend capital, but we're not interested in maximizing revenue here, we're interested in minimizing cost
-   Should incorporate type of intervention as a function of disease prevalence: if prevalence is high, better to immediately reduce disease burden with drug administration; if prevalence is low, take advantage and allocate resources towards interventions that advance sustainable elimination

Optimal control of neglected tropical diseases using a simple dynamic model
================
Chris Hoover
04-04-2018

------------------------------------------------------------------------

Neglected tropical diseases (NTDs) tend to be diseases of poverty that thrive in tropical, disadvantaged communities that often lack access to basic amenities such as clean water, sanitation, and healthcare. They are generally environmentally mediated meaning that some component of their transmission cycle depends on the environment (e.g. a vector, intermediate host, or free-living environmental stage). Control of NTDs generally relies on the administration of drugs that cure individuals from the disease and may confer short-term immunity. However, following drug administration, people are often reinfected by the environmental stage(s) of the pathogen that persist through mass drug administration (MDA) campaigns. Additional interventions that target these environmental reservoirs are often available, but underutilized due to their indirect effects on improving individuals' health.

Using a simple, generalizable model of NTD transmission presented in [Garchitorena et al](http://rstb.royalsocietypublishing.org/content/372/1722/20160128), two interventions will be considered: 1) drug administration, implemented as a pulse reduction in the state variable, *I*, that reduces the prevalence of the disease in the population and 2) environmental remediation (e.g. improvement in sanitation, vector control), implemented as a permanent alteration of a model parameter, that reduces the transmission of the disease.

The first goal is to code the model and reproduce Figure 3c (below) from [Garchitorena et al](http://rstb.royalsocietypublishing.org/content/372/1722/20160128) comparing the effects of interventions that only target infected people, such as mass drug administration (MDA), to the effects of combined interventions that target infected people and the environment, such as sanitation, mosquito spraying, mollusciciding, etc.

![](SupplementaryFiles/Fig3c.png)

The model consists of two state variables, *I* and *W*, that correspond to the prevalence of infection in the human population (e.i. proportion infected at time=*t*) and the degree of contamination of the environment with the disease-causing agent, respectively. We make the simplifying assumption that individuals can only be susceptible, *S*, or infected, *I*, meaning *S* + *I* = 1 and eliminating the need for a recovered *R* compartment as is typical of SIR models but would complicate things here. The model equations then are:

$$\\frac{dI}{dt}=(\\beta\_EW+\\beta\_DI)(1-I)-\\gamma I $$

$$\\frac{dW}{dt}=\\Omega+V\\sigma\\lambda I-\\delta W $$

with parameter definitions and values in table 1 below

| Parameter         |     Value| Description                                                       |
|:------------------|---------:|:------------------------------------------------------------------|
| *β*<sub>*E*</sub> |  9.26e-05| Transmission rate from environment to human population            |
| *β*<sub>*D*</sub> |  5.55e-03| Human to human transmission rate                                  |
| *γ*               |  3.00e-03| Rate of recovery from infected back to susceptible                |
| *Ω*               |  0.00e+00| Recruitment rate of infectious agents in the environment          |
| *V*               |  1.00e+00| Abundance of vectors/intermediate hosts/suitable environment      |
| *λ*               |  1.00e+00| Recruitment rate of infectious agents by infectious individuals   |
| *σ*               |  1.00e+00| Fraction of infectious agents produced that reach the environment |
| *δ*               |  1.70e-02| Mortality rate of infectious agents in the environment            |

We make the simplifying assumptions that there is no exogenous production of infectious agents (*Ω* = 0) and that environment-human transmission is determined only by the transmission rate, *β*<sub>*E*</sub>, with no role of vectors or intermediate hosts (i.e. *V* = *σ* = *λ* = 1). This is approximately equivalent to a simple Cholera model.

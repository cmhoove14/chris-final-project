source("R/Models/deterministic_age_stratified_model.R")
source("R/Models/model_helper_functions.R")
source("R/Models/model_parameters.R")

require(deSolve)

pars1['k_SAC'] <- k_from_prev_W(115/3.6, 0.32)
pars1['k_adult'] <- k_from_prev_W(13/3.6, 0.12)

test <- sim_schisto_age_stratified(nstart = setNames(c(30*area, 10*area, 1*area,
                                               40, 40,
                                               20, 20),
                                             c('S', 'E', 'I', 
                                               'Wt_SAC', 'Wu_SAC',
                                               'Wt_adult', 'Wu_adult')),
                                   time = seq(0, 365*50, 50),
                                   parameters = pars1,
                                   events_df = NA)

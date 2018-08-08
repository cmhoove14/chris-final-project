source("~/RemaisWork/Schisto/R Codes/ag_schist/Human_parasitology/data_load_clean.R")
source("R/Models/model_helper_functions.R")
source("R/Models/model_parameters.R")

require(fitdistrplus)
require(tidyverse)

#Cleaning ##########
comm_haem_sum <- function(skewl, yr){
  egg_vec <- full_long %>% 
    filter(School == skewl & year == yr) %>% 
    pull(s_haem_mean_narm)
  
  n <- sum(!is.na(egg_vec))
  
  nb_sum <- fitdist(round(egg_vec[!is.na(egg_vec)]), "nbinom", method = "mme")[1]$estimate
  
  egg_mean = mean(full_long$s_haem_mean_narm[full_long$School == skewl & full_long$year == yr], na.rm = TRUE)
  haem_pos = length(which(full_long$s_haem_mean_narm[full_long$School == skewl & full_long$year == yr] > 0)) 
  haem_hvy = length(which(full_long$s_haem_mean_narm[full_long$School == skewl & full_long$year == yr] > 50)) 
  haem_prev = haem_pos / length(na.omit(full_long$s_haem_mean_narm[full_long$School == skewl & full_long$year == yr]))
  haem_hvy_prev = haem_hvy / length(na.omit(full_long$s_haem_mean_narm[full_long$School == skewl & full_long$year == yr]))
  
  return(c(school = as.character(skewl),
           year = yr,
           samp_size = as.numeric(n),
           nb_mu = as.numeric(nb_sum[[2]]),
           nb_k = as.numeric(nb_sum[[1]]),
           prev = as.numeric(haem_prev),
           hvy_prev = as.numeric(haem_hvy_prev)))
}

comm_yrs <- expand.grid(skewls = unique(full_long$School)[!is.na(unique(full_long$School))], 
                        yrs = unique(full_long$year))

comm_haem_sums <- data.frame(t(mapply(comm_haem_sum, skewl = comm_yrs[,1], yr = comm_yrs[,2])), stringsAsFactors = FALSE) %>% 
  mutate(samp_size = as.numeric(samp_size),
         nb_mu = as.numeric(nb_mu),
         nb_k = as.numeric(nb_k),
         prev = as.numeric(prev),
         hvy_prev = as.numeric(hvy_prev),
         prev_from_k_w = 1 - (1+nb_mu/nb_k)^-nb_k,
         k_from_prev_m = map2_dbl(nb_mu, prev, k_from_prev_M),
         w_est = map2_dbl(nb_mu, nb_k, w_from_m_k),
         w_sd_est = sqrt(w_est+(w_est^2)/nb_k),
         w_se_est = w_sd_est/samp_size)

rm(list = setdiff(ls(),"comm_haem_sums"))

#Exploratory ############
#Plot mean egg burden over time for each village
comm_haem_sums  %>% 
  ggplot(aes(x = as.numeric(as.character(year)), y = nb_mu, col = school)) +
    geom_line(size = 1.3) + theme_bw()

#Plot clumping parameter over time for each village
comm_haem_sums  %>% 
  ggplot(aes(x = as.numeric(as.character(year)), y = nb_k, col = school)) +
    geom_line(size = 1.3) + theme_bw()

#See how predicted prevalence from M and k matches observed prevalence
comm_haem_sums  %>% 
  ggplot(aes(x = prev_from_k_w, y = prev, size = samp_size)) +
    geom_point() + theme_bw() + stat_smooth(method = "lm")

#See how predicted k from W and prevalence matches observed k
comm_haem_sums  %>% 
  ggplot(aes(x = k_from_prev_m, y = nb_k, size = samp_size)) +
    geom_point() + theme_bw() + stat_smooth(method = "lm")

# PLot egg burden across prevalence
comm_haem_sums  %>% 
  ggplot(aes(x = prev, y = log(nb_mu), size = samp_size)) +
    geom_point() + theme_bw() + stat_smooth(method = "lm")

#Dispersion parameter models ##########
# PLot dispersion paraemter over mean worm burden
comm_haem_sums  %>% 
  ggplot(aes(x = w_est, y = nb_k, size = samp_size)) +
    geom_point() + theme_bw() +stat_smooth(method = "lm")

w_mod <- glm(nb_k ~ w_est, data = comm_haem_sums, weights = samp_size)
quad_w_mod <- nls(nb_k ~ a*w_est^2 + b*w_est + c, data = comm_haem_sums, start = list(a=1, b=1, c=1))
# PLot dispersion paraemter over log mean worm burden
comm_haem_sums  %>% 
  ggplot(aes(x = log(w_est), y = nb_k, size = samp_size)) +
    geom_point() + theme_bw() + stat_smooth(method = "lm")

log_w_mod <- glm(nb_k ~ log(w_est), data = comm_haem_sums, weights = samp_size)

# PLot log of dispersion parameter over mean worm burden
comm_haem_sums %>% 
  ggplot(aes(x = w_est, y = log(nb_k), size = samp_size)) +
    geom_point() + theme_bw() +stat_smooth(method = "lm")

w_mod_log <- glm(log(nb_k) ~ w_est, data = comm_haem_sums, weights = samp_size)

# PLot log-log dispersion parameter over mean worm burden
comm_haem_sums %>% 
  ggplot(aes(x = log(w_est), y = log(nb_k), size = samp_size)) +
    geom_point() + theme_bw() +stat_smooth(method = "lm")

log_w_mod_log <- glm(log(nb_k) ~ log(w_est), data = comm_haem_sums, weights = samp_size)

#Plot dispersion parameter across prevalence
comm_haem_sums %>% 
  ggplot(aes(x = prev, y = nb_k, size = samp_size)) +
    geom_point() + theme_bw() + stat_smooth(method = "lm")

prev_mod <- glm(nb_k ~ prev, data = comm_haem_sums, weights = samp_size)

#Plot log of dispersion parameter across prevalence
comm_haem_sums %>% 
  ggplot(aes(x = prev, y = log(nb_k), size = samp_size)) +
    geom_point() + theme_bw() + stat_smooth(method = "lm")

prev_mod_log <- glm(log(nb_k) ~ prev, data = comm_haem_sums, weights = samp_size)

# Plot dispersion parameter acrpss heavy worm burden prevalence
comm_haem_sums %>% 
  ggplot(aes(x = hvy_prev, y = nb_k, size = samp_size)) +
    geom_point() + theme_bw() + stat_smooth(method = "lm")

hvy_prev_mod <- glm(nb_k ~ hvy_prev, data = comm_haem_sums, weights = samp_size)

both_prev_mod <- glm(nb_k ~ prev + hvy_prev, data = comm_haem_sums, weights = samp_size)

mu_prev_mod <- glm(nb_k ~ prev + w_est, data = comm_haem_sums, weights = samp_size)

AIC(w_mod, quad_w_mod, log_w_mod, w_mod_log, log_w_mod_log, prev_mod, prev_mod_log, hvy_prev_mod, both_prev_mod, mu_prev_mod)

dev.off()
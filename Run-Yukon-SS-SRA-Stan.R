#==================================================================================================
# Project Name: YUKON RIVER CHINOOK IMEG WORKING GROUP - STAN version of age-structured state-space spawner-recruitment model with AR-1 process variation (Fleischman et al. CJFAS. 2013)
# Creator: Brendan Connors, Fisheries and Oceans Canada and Curry Cunningham, College of Fisheries and Ocean Sciences, UAF
# Date: 14.10.2020
#
# Purpose: Fit state-space spawner recruitment model to CDN-origin Yukon Chinook spawner, harvest and age composition data
#
#==================================================================================================
# NOTES:
#  1) starting simple with fixed maturity schedule and escapement and harvest from JTC report
#==================================================================================================

require(rstan)

# SA: added:
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

Esc <- read.csv('data/EscTab.csv') # escapement
Har <- read.csv('data/HarTab.csv') # harvest
S_cv <- rep(0.2,length(Esc$year)) # placeholder for CV on escapement observations
H_cv <- rep(0.05,length(Har$year)) # placeholder for CV on harvest observations
p_mat <- c(0.05,0.31,0.54,0.1) # placeholder fixed age at maturity
a_min <- 4
a_max <- 7
nyrs <- length(Esc$year)
A <- a_max - a_min + 1
nRyrs <- nyrs + A - 1

# Run Stan Model ====================================================

# STAN MODEL DATA
  stan.data <- list("nyrs" = nyrs,
                    "a_min" = a_min,
                    "a_max" = a_max,
                    "A" = A,
                    "nRyrs" = nyrs + A - 1,
                    "p" = p_mat,
                    "S_obs" = Esc$Esc/1000,
                    "H_obs" = Har$Harvest/1000,
                    "S_cv" = S_cv,
                    "H_cv" = H_cv)


# FIT STAN MODEL without reparameterization
  stan.fit.no_reparam_sigma_R0 <- stan(file = "SSSR_AR1.v1.stan",
                   model_name = "SSSR_AR1.v1",
                   data = stan.data,
                   chains = 4,
                   iter = 1000,
                   seed = 42,
                   control = list(adapt_delta = 0.99, max_treedepth = 20)
  )

pairs(stan.fit.no_reparam_sigma_R0, pars = c("lnalpha", "beta", "sigma_R", "sigma_R0", "phi", "mean_ln_R0","lnR[1]", "lnR[2]"))


# FIT STAN MODEL with reparameterization
stan.fit.reparam_sigma_R0 <- stan(file = "SSSR_AR1.sigmaR0-reparam.v1.stan",
                                     model_name = "SSSR_AR1.sigmaR0-reparam.v1",
                                     data = stan.data,
                                     chains = 4,
                                     iter = 1000,
                                     seed = 42,
                                     control = list(adapt_delta = 0.99, max_treedepth = 20)
)

pairs(stan.fit.reparam_sigma_R0, pars = c("lnalpha", "beta", "sigma_R", "sigma_R0", "phi", "mean_ln_R0","lnR[1]", "lnR[2]"))




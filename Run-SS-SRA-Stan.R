
require(rstan)
# require(plyr)
# require(tidyverse)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

Esc <- read.csv('data/EscTab.csv') # escapement
Har <- read.csv('data/HarTab.csv') # harvest
S_cv <- rep(0.2,length(Esc$year)) # placeholder for CV on escapement observations
H_cv <- rep(0.05,length(Har$year)) # placeholder for CV on harvest observations
A_obs <- read.csv('data/AgeComps.csv') # age composition
a_min <- 4
a_max <- 7
nyrs <- length(Esc$year)
A <- a_max - a_min + 1
nRyrs <- nyrs + A - 1

# Initialize Model Parameters =======================================

simdata <- data.frame(id = rep(1:nRyrs, each = 1),
                      x = runif(1,min=0.01, max=0.99),
                      y = runif(1,min=0.01, max=0.99),
                      z = runif(1,min=0.01, max=0.99),
                      zz = runif(1,min=0.01, max=0.99))

# init_fn <- function(chain_id=1) {
#   list(
#     "lnR"=abs(rnorm(nRyrs, mean=0, sd=5)),
#     "lnalpha"=abs(rnorm(1, mean=0, sd=1)),
#     "beta"=runif(1, min=0.01, max=9.99),
#     "sigma_R"=abs(rnorm(1, mean=0, sd=1)),
#     "sigma_R0"=abs(rnorm(1, mean=0, sd=1)),
#     "phi"=runif(1, min=-1, max=1),
#     "lnresid_0"=runif(1, min=-1, max=1),
#     "mean_ln_R0"=abs(rnorm(1, mean=0, sd=1)),
#     "U"=runif(nyrs, min=0.01, max=0.99),
#     "g"=plyr::daply(simdata %>% dplyr::mutate(id = as.integer(id)), "id",
#               function(df) df[1,c("x", "y", "z", "zz")]) %>% as.numeric %>% matrix(ncol=4)
#   )
# }
# Initial List of Lists for Multiple Chains
# init_ll <- lapply(1, function(id) init_fn(chain_id = id))

# SA: you can skip init_ll and just have init_fn() without chain_id (unless you want control over different inits in different chains)

# Run Stan Model ====================================================

# DATA
stan.data <- list("nyrs" = nyrs,
                  "a_min" = a_min,
                  "a_max" = a_max,
                  "A" = A,
                  "nRyrs" = nyrs + A - 1,
                  "A_obs" = as.matrix(A_obs[,2:5]),
                  "S_obs" = Esc$Esc/1000,
                  "H_obs" = Har$Harvest/1000,
                  "S_cv" = S_cv,
                  "H_cv" = H_cv)


# FIT
stan.fit <- stan(file = "SSSR_AR1-time-vary-mat.stan",
                 model_name = "SSSR_AR1-time-vary-mat",
                 data = stan.data,
                 chains = 4,
                 iter = 3000,
                 seed = 42,
                 init = init_ll,
                 thin = 2,
                 control = list(adapt_delta = 0.99, max_treedepth = 20))

# SUMMARIZE
shinystan::launch_shinystan(stan.fit) # there is still some autocorrelation in a few parameters, try increasing iterations and thinning more

pairs(stan.fit, pars = c("lnalpha", "beta", "sigma_R", "sigma_R0", "phi", "mean_ln_R0","lnR[1]", "lnR[2]"))

summary(stan.fit, pars= c("lnalpha", "beta", "sigma_R", "sigma_R0", "phi", "mean_ln_R0","S_max", "S_eq", "S_msy", "U_msy"))$summary


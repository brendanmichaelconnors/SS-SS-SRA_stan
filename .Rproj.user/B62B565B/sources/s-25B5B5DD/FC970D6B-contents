// Project Name: YUKON RIVER CHINOOK IMEG WORKING GROUP - STAN version of age-structured state-space spawner-recruitment model with AR-1 process variation (Fleischman et al. CJFAS. 2013)
// Creator: Brendan Connors, Fisheries and Oceans Canada and Curry Cunningham, College of Fisheries and Ocean Sciences, UAF
// Date: 14.10.2020
// Purpose: State-space spawner recruitment model
// Notes: (1) starting simple with fixed maturity schedule hence declaring p as fixed proportions

data{
  int nyrs;  // number of calender years
  int a_min;  // minimum age class
  int a_max;  // maximum age class
  int A; // number of age classes
  int nRyrs; // number of recruitment years
  vector[A] p; // age at maturity proportions, to start treating these as fixed
  vector[nyrs] S_obs; // observed spawners
  vector[nyrs] H_obs; // observed harvest
  vector[nyrs] S_cv; // spawner observation error CV
  vector[nyrs] H_cv; // harvest observation error CV
}

// SA: my guess is you don't want upper on beta and upper or lower on lnresid_0
// SA: it may not make a big difference
parameters{
  vector<lower=0>[nRyrs] lnR; // log recruitment states
  real<lower=0> lnalpha; // log Ricker a
  real<lower=0,upper=10> beta; // Ricker b
  real<lower=0> sigma_R; // process error
  real<lower=0> sigma_R0; // process error for first a.max years with no spawner link
  real<lower=-1,upper=1> phi; // lag-1 correlation in process error
  real<lower=-2,upper=2> lnresid_0;  // first residual for lag-1 correlation in process error
  real<lower=0> mean_ln_R0; // "true" mean log recruitment in first a.max years with no spawner link
  vector<lower=0.01,upper=0.99>[nyrs] U;  // harvest rate
}

transformed parameters{
  vector<lower=0>[nyrs] N;  // run size states
  vector<lower=0>[nyrs] S;  // spawner states
  vector[nyrs] lnS;  // log spawner states
  vector<lower=0>[nyrs] C;  // catch states
  vector[nyrs] lnC;  // log catch states
  vector<lower=0>[nRyrs] R;  // recruitment states
  vector[nRyrs] lnresid;  // log recruitment residuals
  vector[nRyrs] lnRm_1;  // log recruitment states in absence of lag-one correlation
  vector[nRyrs] lnRm_2;  // log recruitment states after accounting for lag-one correlation
  matrix<lower=0>[nyrs, A] N_ta; // returns by age matrix


  R = exp(lnR);

  for (t in 1:nyrs) {
    for(a in 1:A){
      N_ta[t,a] = R[t+A-a] * p[a];
    }
  }

  for(t in 1:nyrs) {
    N[t] = sum(N_ta[t,1:A]);
    S[t] = N[t] * (1 - U[t]);
    lnS[t] = log(S[t]);
    C[t] = N[t] * U[t];
    lnC[t] = log(C[t]);
  }

  // SA:
  for (i in 1:nRyrs) {
    lnresid[i] = 0.0;
    lnRm_1[i] = 0.0;
    lnRm_2[i] = 0.0;
  }
  for (y in (A+a_min):nRyrs) {
    lnRm_1[y] = lnS[y-a_max] + lnalpha - beta * S[y-a_max];
    lnresid[y] = lnR[y] - lnRm_1[y];
  }

  lnRm_2[A+a_min] =  lnRm_1[A+a_min] + phi * lnresid_0;

  for (y in (A+a_min+1):nRyrs) {
    lnRm_2[y] =  lnRm_1[y] + phi * lnresid[y-1];
  }
}

model{
  //  Define Priors
  lnalpha ~ normal(0,3);
  beta ~ normal(0,1);
  sigma_R ~ normal(0,2);
  phi ~ uniform(-1,1); // SA: not needed
  lnresid_0 ~ normal(0,20);

  // First `a.max` years of recruits, for which there is no spawner link
  // mean_ln_R0 ~ normal(0,10000); // SA: reduce scale?
  mean_ln_R0 ~ normal(0,20); // SA
  // sigma_R0 ~ inv_gamma(0.1,0.1); // SA: avoid if possible; see other .stan file
  sigma_R0 ~ student_t(7, 0, 3); // SA

  // SA: vectorize?
  for (i in 1:a_max){
    lnR[i] ~ normal(mean_ln_R0, sqrt(sigma_R0));
  }

  // State model
  lnR[(A+a_min):nRyrs] ~ normal(lnRm_2[(A+a_min):nRyrs], sigma_R);

  // SA: vectorize?
  // Observation model
  for(t in 1:nyrs){
    U[t] ~ uniform(0.01, 0.99); // SA: not needed; use Beta or normal?
    H_obs[t] ~ lognormal(lnC[t], sqrt(log((H_cv[t]^2) + 1)));
    S_obs[t] ~ lognormal(lnS[t], sqrt(log((S_cv[t]^2) + 1)));
  }
}



///
// This Stan program defines a model for estimating cholera positivity rates,
// i.e., the proportion of suspected cholera cases that are true Vc infections,
// using a meta-analysis across >100 studies, with covariates, random effect by study and
// adjusting for senstivity/specificity of PCR, culture, and RDT for identifying Vc infections.

// We have posterior distribution of estimates sensitivity/specificity from the latent class analysis,
// and data from >100 systematic review studies (from which we'd like to estimate positivity rate).
data {
  
  // positivity rates from meta-analysis by test
  int<lower=1> N_obs; // number of study-location-periods (observations) included in main analysis
  int<lower=1> J; // number of laboratory tests
  int<lower=0> num_test[N_obs, J]; // number of suspected cholera cases in each study that were tested by each test
  int<lower=0> num_pos[N_obs, J]; // number of participants from each study that tested positive by each test
  
  // random effect
  int<lower=1> N_re; // number of unique random effect categories
  int<lower=1, upper=N_re> re[N_obs]; // random effect corresponding to observation

  // covariates
  int<lower=1> p_vars; // number of variables to adjust for
  matrix[N_obs, p_vars] X; // covariate model matrix
  
  // sensitivity and specificity estimates
  int<lower=1> M; // number of posterior draws of sensitivity/specificity
  matrix[M, J] sens; // draws of sensitivity for each test
  matrix[M, J] spec; // draws of specificity for each test
}

parameters {
  // positivity overall
  real mu_logit_p;
  real<lower = 0> sigma_logit_p;

  vector[p_vars] beta; // fixed regression coefficients
  real<lower=0> sigma_re; // variability of random effect
  vector[N_re] eta_re; // standard normals for the random effect
}

transformed parameters {
  // probability of positive Vc result for each observation
  vector[N_obs] logit_p; 
  vector<lower=0, upper=1>[N_obs] p;
  logit_p = X * beta + sigma_re * eta_re[re];
  p = inv_logit(logit_p);
}

model {
  
  real lp[J, M];
  
  // proportion positive
  for (j in 1:J) {
    for (m in 1:M) {
      lp[j, m] = binomial_lpmf(num_pos[,j] | num_test[,j], p*sens[m,j]+(1-p)*(1-spec[m,j]));
    }
  }
  
  for (j in 1:J) {
    target += -log(M) + log_sum_exp(to_vector(lp[j,]));
  }

  // priors on coefficients and random effect
  target += normal_lpdf(beta | 0, 2); // increased variance
  target += normal_lpdf(eta_re | 0, 2); // increased variance

  // priors for pooling positivity rates
  target += normal_lpdf(logit_p | mu_logit_p, sigma_logit_p);
  target += normal_lpdf(sigma_logit_p | 0, 2); // increased variance
  target += normal_lpdf(mu_logit_p | 0, 2); // increased variance
}

generated quantities {

  // posterior predictive distribution for underlying true positive
  real logit_p_pred = normal_rng(mu_logit_p, sigma_logit_p);
  real p_pred = inv_logit(logit_p_pred);
}


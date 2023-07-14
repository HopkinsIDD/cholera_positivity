///
// This Stan program defines a model for estimating cholera positivity rates,
// i.e., the proportion of suspected cholera cases that are true Vc infections,
// using a meta-analysis across >100 studies, with covariates, random effect by study and
// adjusting for senstivity/specificity of PCR, culture, and RDT for identifying Vc infections.

// We have posterior distribution of estimates sensitivity/specificity from the latent class analysis,
// and data from >100 systematic review studies (from which we'd like to estimate positivity rate).
data {
  
  // positivity rates from meta-analysis by test
  int<lower=1> N_obs; // number of observations included in main analysis
  int<lower=1> J; // number of laboratory tests
  int<lower=1> N_obs_tests; // number of observations and tests included in the analysis
  int<lower=1> num_test[N_obs_tests]; // number of suspected cholera cases in each study that were tested
  int<lower=0> num_pos[N_obs_tests]; // number of participants from each study that tested positive
  
  // indices for test type used
  int<lower=1> N_culture; // number of observations that used culture
  int<lower=1> N_pcr; // number of observations that used pcr
  int<lower=1> N_rdt; // number of observations that used rdt
  int<lower=1> culture_id_long[N_culture]; // row ids for where culture was used for long dataset
  int<lower=1> pcr_id_long[N_pcr]; // row ids for where pcr was used for long dataset
  int<lower=1> rdt_id_long[N_rdt]; // row ids for where rdt was used for long dataset
  int<lower=1> culture_id_wide[N_culture]; // row ids for where culture was used for wide dataset
  int<lower=1> pcr_id_wide[N_pcr]; // row ids for where pcr was used for wide dataset
  int<lower=1> rdt_id_wide[N_rdt]; // row ids for where rdt was used for wide dataset
  
  // random effect
  int<lower=1> N_re; // number of unique random effect categories
  int<lower=1, upper=N_re> re[N_obs]; // random effect corresponding to observation for wide dataset

  // covariates
  int<lower=1> p_vars; // number of variables to adjust for
  matrix[N_obs, p_vars] X; // covariate model matrix for wide dataset
  
  // sensitivity and specificity estimates
  int<lower=1> M; // number of posterior draws of sensitivity/specificity
  matrix[M, J] sens; // draws of sensitivity for each test
  matrix[M, J] spec; // draws of specificity for each test
}

transformed data {
  // sample one of the sensitivity/specificity draws;
  int<lower = 1, upper = M> m;
  simplex[M] uniform = rep_vector(1.0 / M, M);
  m = categorical_rng(uniform); 
}

parameters {
  real alpha; // intercept
  vector[p_vars] beta; // fixed regression coefficients
  real<lower=0> sigma_re; // variability of random effect
  vector[N_re] eta_re; // standard normals for the random effect
}

transformed parameters {
  // probability of positive Vc result for each observation
  vector[N_obs] logit_p; 
  vector<lower=0, upper=1>[N_obs] p;
  logit_p = alpha + X * beta + sigma_re * eta_re[re];
  p = inv_logit(logit_p);
}

model {
  
  // proportion positive
  // culture
  target += binomial_lpmf(num_pos[culture_id_long] | num_test[culture_id_long], p[culture_id_wide]*sens[m,1]+(1-p[culture_id_wide])*(1-spec[m,1]));
  // pcr
  target += binomial_lpmf(num_pos[pcr_id_long] | num_test[pcr_id_long], p[pcr_id_wide]*sens[m,2]+(1-p[pcr_id_wide])*(1-spec[m,2]));
  // rdt
  target += binomial_lpmf(num_pos[rdt_id_long] | num_test[rdt_id_long], p[rdt_id_wide]*sens[m,3]+(1-p[rdt_id_wide])*(1-spec[m,3]));

  // priors on intercept, coefficients, random effect, and sigma_re
  target += normal_lpdf(beta | 0, 2);
  target += normal_lpdf(eta_re | 0, 1);
  target += normal_lpdf(sigma_re | 0, 1);
  target += normal_lpdf(alpha | 0.9, 2); // shift
}

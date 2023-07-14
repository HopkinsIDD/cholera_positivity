### Model settings ------------------------

# JAGS settings
jags_redo <- FALSE # set this to TRUE to re-run the jags_run_id, FALSE to load previous model fits
jags_run_id <- 'jags-cd'

# Stan settings
stan_redo <- TRUE # set this to TRUE to re-run the stan_run_id corresponding to the stan_model and draws
stan_model <- 'mod1' # which Stan model to run
adjust_tests <- TRUE # whether or not to adjust for sensitivity/specificity of the tests
sa <- TRUE # whether or not to run sensitivity analysis where we shift prior on alpha
test <- FALSE # add a 'test' flag to run_id


### Setup ------------------------

# packages
pacman::p_load('tidyverse', 'data.table', 'readxl',
               'stringr', 'lubridate', 'reshape2',
               'cmdstanr', 'R2jags', 'boot',
               'sf', 'raster', 'countrycode')

# custom functions
source(here::here('code', 'utils.R'))

# Coefficient equations corresponding to each model
stan_mod_eqns <- list(
  'mod1' = '1',
  'mod2' = 'sampling_quality + cases_min_age_bin',
  'mod3' = 'sampling_quality + cases_min_age_bin + surveillance_type',
  'mod4' = 'sampling_quality + surveillance_type'
)
stan_coef_eqn <- stan_mod_eqns[[stan_model]]
stan_run_id <- paste0(stan_model, 
                      ifelse(adjust_tests, '', '-unadj'),
                      ifelse(sa, '-sa', ''),
                      ifelse(test, '-test', ''))

# use all cores
options(mc.cores = parallel::detectCores())

# set random seed
set.seed(7341)


### Load and clean data ------------------------

#  load data
df <- read.csv(here::here('data', 'extracted_data', 'cc_extractions.csv'))
vars <- read.csv(here::here('data', 'extracted_data', 'extraction_variables.csv'))

# set column variable names
names(df) <- gsub('\\.', ' ', names(df))
df <- relocate(df, vars$Variable.label)
names(df) <- vars$Variable
df <- df[!is.na(names(df))]

# format columns
df <- df %>%
  # primary dataset
  filter(primary==1) %>%
  # dates
  mutate(outbreak_L = as.Date(outbreak_L, format = c('%d-%B-%Y')),
         outbreak_R = as.Date(outbreak_R, format = c('%d-%B-%Y')),
         TL = as.Date(TL, format = c('%d-%B-%Y')),
         TR = as.Date(TR, format = c('%d-%B-%Y'))) %>%
  # remove variable name text in parentheses and set factor levels
  mutate(study_design = gsub('\\s*\\([^\\)]+\\)', '', study_design),
         sampling_approach = gsub('\\s*\\([^\\)]+\\)', '', sampling_approach),
         study_design = factor(study_design, 
                               levels = c('Diagnostic test accuracy',
                                          'Vaccine effectiveness',
                                          'Surveillance', 'Other')),
         sampling_approach = factor(sampling_approach, 
                                    levels = c('All suspected cases', 'Systematic',
                                               'Simple random', 'Other', 
                                               'Convenience', 'Unknown')),
         sampling_quality = ifelse(sampling_approach %in% c('All suspected cases', 
                                                            'Systematic',
                                                            'Simple random'),
                                      'High', 'Low'),
         cases_dehydrated = ifelse(cases_dehydrated == 1, 'Yes', 'No')) %>%
  # region as defined by World Bank Development Indicators
  mutate(region = countrycode(country_iso3, 
                              origin = 'iso3c', 
                              destination = 'region')) %>%
  mutate(region = ifelse(region != 'Sub-Saharan Africa',
                         region, 
                         countrycode(country_iso3, 
                                     origin = 'iso3c', 
                                     destination = 'region23'))) %>%
  # surveillance is outbreak or non-outbreak
  mutate(surveillance_type = ifelse(surveillance_type == 'Outbreak surveillance',
                                    'Outbreak', 'Non-outbreak')) %>%
  # binary version of age variable
  mutate(cases_min_age_bin = ifelse(cases_min_age == 0, 0, 1)) %>%
  # remove incomplete extractions
  filter(!is.na(sampling_approach))


### Prepare validation data for JAGS model ------------------

# Re-create dataset using validation data from Fig 2, Table 4, and text of Debes et al. 2021
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8176501/
df_ken <- rbind(
  data.frame(id = 1:8, culture=1, pcr=0, rdt=1),
  data.frame(id = 1:69, culture=1, pcr=1, rdt=1),
  data.frame(id = 1:2, culture=0, pcr=1, rdt=1),
  data.frame(id = 1:2, culture=0, pcr=1, rdt=0),
  data.frame(id = 1:149, culture=0, pcr=0, rdt=0)
)
df_ken$id <- NULL

# Re-create dataset using validation data from Table 2 and Fig 3 of Sayeed et al. 2018
# The results for Cholkit and Crystal VC are the same comapred to the other tests
# https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0006286
df_bgd <- rbind(
  data.frame(id = 1:15, culture=1, pcr=1, rdt=1),
  data.frame(id = 1:6, culture=0, pcr=1, rdt=1),
  data.frame(id = 1:4, culture=1, pcr=0, rdt=1),
  data.frame(id = 1:5, culture=0, pcr=0, rdt=1),
  data.frame(id = 1:47, culture=0, pcr=0, rdt=0)
)
df_bgd$id <- NULL

# Validation data from Table 1 of Ontweka et al. 2016 PLOS ONE
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0168257
df_ssd <- read_xlsx(here::here('data', 'raw_data', 'pone.0168257.s002.xlsx')) %>%
  dplyr::select(`...12`, `...14`, `...9`) 

names(df_ssd) <- c('culture', 'pcr', 'rdt')
df_ssd <- df_ssd[-c(1,2),]

df_ssd <- df_ssd %>%
  mutate(culture = ifelse(culture=='Negative', 0, 1),
         rdt = as.numeric(rdt),
         rdt = ifelse(rdt==2, 1, rdt),
         pcr = as.numeric(pcr)) %>%
  na.omit()

# Validation data from Mwaba et al. 2018
# https://onlinelibrary.wiley.com/doi/full/10.1111/tmi.13084
df_zmb_pcr <- read_xlsx(here::here('data', 'raw_data', 'Lusaka_2016_PCR.xlsx')) %>%
 dplyr::select(Sample, `Final result`) %>%
 na.omit() %>%
 mutate(`Final result` = ifelse(`Final result` == '+', 1, 0))

df_zmb_rdt <- read_xlsx(here::here('data', 'raw_data', 'Lusaka_2016_RDT.xlsx'),
                       range = 'C4:L237') %>%
 # Included in the paper
 dplyr::filter(Include==1) %>%
 mutate(Sample = as.numeric(gsub('o', '', ID))) %>%
 mutate(`Direct RDT` = ifelse(`APW RDT` == '+' | `Direct RDT` == 'Pos', 1, 0),
        Culture = ifelse(Culture == '+' | Culture == 'Pos', 1, 0))

df_zmb <- merge(df_zmb_rdt%>%dplyr::select(Sample, Culture, `Direct RDT`),
               df_zmb_pcr) %>%
 dplyr::select(Culture, `Final result`, `Direct RDT`)

names(df_zmb) <- c('culture', 'pcr', 'rdt')
rm(df_zmb_pcr, df_zmb_rdt)

# create parent dataset and vector with study label
df_val <- rbind(df_ken, df_bgd, 
                df_ssd, df_zmb
                )
val_id <- c(rep(1, nrow(df_ken)),
            rep(2, nrow(df_bgd)),
            rep(3, nrow(df_ssd)),
            rep(4, nrow(df_zmb))
            )


### Run JAGS model to get sensitivity and specificity --------------

# set run id
message('Running analysis: ', jags_run_id)

# file to save to
jags_est_file <- here::here('data', 'generated_data', 'sens_spec_estimates', 
                            paste0('sens-spec-', jags_run_id, '.rds'))

# run model
if (!file.exists(jags_est_file) | jags_redo) {
  
  jags_fit <- jags(data = list('n'=nrow(df_val), # number of observations
                               'K'=ncol(df_val), # number of tests
                               'S'=length(unique(val_id)), # number of studies
                               'sid'=val_id, # vector with study ids
                               'x'=df_val, # dataset
                               'z'=rep(0, nrow(df_val))),
                   parameters.to.save = c('pi', 'mu_se', 'mu_sp', 'se', 'sp'), 
                   n.chains = 4, 
                   n.iter = 200000,
                   model.file = here::here('code', 'latent_class_cd_meta.jags'))
  
  saveRDS(jags_fit, jags_est_file)
  
} else {
  message('Loading pre-computed estimates')
  jags_fit <- readRDS(jags_est_file)
}

print(jags_fit)

# draws of sensitivity and specificity
se_draws <- jags_fit$BUGSoutput$sims.list$mu_se
sp_draws <- jags_fit$BUGSoutput$sims.list$mu_sp
colnames(se_draws) <- colnames(sp_draws) <- c('Culture', 'PCR', 'RDT')


### Prepare epi data for Stan model ------------------------

# set covariates
covs <- c('sampling_quality', 'cases_min_age_bin', 'cases_dehydrated', # adjust for these
          'surveillance_type') # examine differences by these (and age categories in separate analysis)

# positivity rates from meta-analysis primary dataset
sr_dat <- df %>%
  group_by_at(c('study_id', 'primary', 'country_iso3', 'region', covs,
                'diagnostic_test_type', 'diagnostic_test_name')) %>%
  summarize(n_cases_tested = sum(n_cases_tested),
            n_cases_positive = sum(n_cases_positive)) %>%
  ungroup() %>%
  # subset to primary dataset
  filter(primary==1)

# cast number tested by test type
test_dat <- sr_dat %>%
  reshape2::dcast(as.formula(paste0(paste0(c('study_id', 'country_iso3', 'region', covs),
                                           collapse = '+'),
                                    '~ diagnostic_test_type')),
                  value.var = 'n_cases_tested', fill = NA)

# cast number positive by test type
pos_dat <- sr_dat %>%
  reshape2::dcast(as.formula(paste0(paste0(c('study_id', 'country_iso3', 'region', covs),
                                           collapse = '+'),
                                    '~ diagnostic_test_type')),
                  value.var = 'n_cases_positive', fill = NA)

# add random effect label by observation
test_dat <- test_dat %>% mutate(study_lab = 1:nrow(test_dat))
pos_dat <- pos_dat %>% mutate(study_lab = 1:nrow(test_dat))

# save covariate matrix
dat_covs <- pos_dat

# reshape long by test and remove NA
test_dat <- reshape2::melt(test_dat, id.vars = c('study_id', 'study_lab', 'country_iso3', 'region', covs),
                           variable.name = 'diagnostic_test_type', value.name = 'n_cases_tested') %>%
  filter(!is.na(n_cases_tested)) %>%
  mutate(test_id = ifelse(diagnostic_test_type == 'Culture', 1,
                          ifelse(diagnostic_test_type == 'PCR', 2, 3)))

pos_dat <- reshape2::melt(pos_dat, id.vars = c('study_id', 'study_lab', 'country_iso3', 'region', covs),
                          variable.name = 'diagnostic_test_type', value.name = 'n_cases_positive') %>%
  filter(!is.na(n_cases_positive)) %>%
  mutate(test_id = ifelse(diagnostic_test_type == 'Culture', 1,
                          ifelse(diagnostic_test_type == 'PCR', 2, 3)))


### Model with all data and covariates ------------------------

# set run id
message('Running analysis: ', stan_run_id)

# file to save to
out_est_file <- here::here('data', 'generated_data', 'adjusted_positivity_rates',
                           paste0('posrate-', stan_run_id, '.rds'))

# model script
model_script <- ifelse(adjust_tests, 
                       here::here('code', 'meta_analysis_re.stan'),
                       here::here('code', 'meta_analysis_re_unadj.stan'))

model_script <- ifelse(sa, 
                       here::here('code', 'meta_analysis_re_sa.stan'),
                       model_script)

# run model
if (!file.exists(out_est_file) | stan_redo) {
  posrate <- run_analysis_stan(
    model_script = model_script,
    dat_suspected = test_dat,
    dat_confirmed = pos_dat,
    dat_covars = dat_covs,
    re_ids = dat_covs$study_lab,
    test_ids = test_dat$test_id,
    run_id = stan_run_id,
    coef_eqn = stan_coef_eqn,
    sens_draws = se_draws,
    spec_draws = sp_draws,
    redo = stan_redo,
    adjust_tests = adjust_tests
  )

  saveRDS(posrate, out_est_file)

} else {
  posrate <- readRDS(out_est_file)
}

### Check model convergence and other diagnostics ------------------------------

stan_out_file <- here::here("data", "generated_data", "model_fits", 
                            paste0("stan_fit_", stan_run_id,".rds"))

stan_est <- readRDS(stan_out_file)

stan_est$cmdstan_diagnose()

#shinystan::launch_shinystan(stan_est)

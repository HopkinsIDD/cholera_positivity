### Model settings ------------------------

# JAGS settings
jags_redo <- FALSE # set this to TRUE to re-run the jags_run_id, FALSE to load previous model fits
jags_run_id <- 'jags-cd'

# Stan settings
stan_redo <- FALSE # set this to TRUE to re-run the stan_run_id corresponding to the stan_model and draws
stan_model <- 'mod1'# which Stan model to run
se_sp_draws <- 1000 # number of draws of sens/spec to use from the JAGS model output
data_sub <- 'full' # subset of the data to run; options = c('full', 'outbreak', 'non-outbreak', 'age0', 'age1', 'age2', 'age5', 'high', 'low')
sa <- FALSE # whether or not we're running a sensitivity analysis on the priors (increase variance)
test <- FALSE # add a 'test' flag to run_id


### Setup ------------------------

# packages
pacman::p_load('tidyverse', 'data.table', 'readxl',
               'stringr', 'lubridate', 'reshape2',
               'cmdstanr', 'R2jags', 
               'sf', 'raster', 'countrycode')

# custom functions
source(here::here('code', 'utils.R'))

# Coefficient equations corresponding to each model
stan_mod_eqns <- list(
  'mod1' = '1',
  'mod2' = 'sampling_quality + cases_min_age',
  'mod2c' = 'sampling_quality + cases_min_age_categ',
  'mod3' = 'sampling_quality + cases_min_age + surveillance_type',
  'mod3c' = 'sampling_quality + cases_min_age_categ + surveillance_type'
)
stan_coef_eqn <- stan_mod_eqns[[stan_model]]
stan_run_id <- paste0(stan_model, ifelse(!sa, '', '-sa'), '-d', se_sp_draws, '-', data_sub, ifelse(test, '-test', ''))

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
  # categorical version of age variable
  mutate(cases_min_age_categ = paste0('age', cases_min_age)) %>%
  # remove incomplete extractions
  filter(!is.na(sampling_approach))


### Sub-group analysis ----------------------------------------

# subset, if applicable
if (data_sub == 'outbreak') {
  df <- df %>% filter(surveillance_type == 'Outbreak')
}
if (data_sub == 'non-outbreak') {
  df <- df %>% filter(surveillance_type == 'Non-outbreak')
}
if (data_sub == 'high') {
  df <- df %>% filter(sampling_quality == 'High')
}
if (data_sub == 'low') {
  df <- df %>% filter(sampling_quality == 'Low')
}
if (data_sub == 'age0') {
  df <- df %>% filter(cases_min_age == 0)
}
if (data_sub == 'age1') {
  df <- df %>% filter(cases_min_age == 1)
}
if (data_sub == 'age2') {
  df <- df %>% filter(cases_min_age == 2)
}
if (data_sub == 'age5') {
  df <- df %>% filter(cases_min_age == 5)
}
if (! data_sub %in% c('full', 'outbreak', 'non-outbreak', 'age0', 'age1', 'age2', 'age5', 'high', 'low')) {
  stop(data_sub, ' is not an option for sub-group analysis')
}


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

# randomly select draws to run in stan
se_draws <- jags_fit$BUGSoutput$sims.list$mu_se[round(runif(se_sp_draws, 1, 4000)), ]
sp_draws <- jags_fit$BUGSoutput$sims.list$mu_sp[round(runif(se_sp_draws, 1, 4000)), ]
colnames(se_draws) <- colnames(sp_draws) <- c('Culture', 'PCR', 'RDT')


### Prepare epi data for Stan model ------------------------

# set covariates
covs <- c('sampling_quality', 'cases_min_age', 'cases_min_age_categ', 'cases_dehydrated', # adjust for these
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
                  value.var = 'n_cases_tested', fill = 0)

# cast number positive by test type
pos_dat <- sr_dat %>%
  reshape2::dcast(as.formula(paste0(paste0(c('study_id', 'country_iso3', 'region', covs),
                                           collapse = '+'),
                                    '~ diagnostic_test_type')),
                  value.var = 'n_cases_positive', fill = 0)

# add random effect label by observation
test_dat <- test_dat %>% mutate(study_lab = 1:nrow(test_dat))
pos_dat <- pos_dat %>% mutate(study_lab = 1:nrow(test_dat))

# add PCR column as 0's if we don't have any of this data by the sub-set
if (!'PCR' %in% names(test_dat)) {
  test_dat$PCR <- 0
  pos_dat$PCR <- 0
}


### Model with all data and covariates ------------------------

# set run id
message('Running analysis: ', stan_run_id)

# file to save to
out_est_file <- here::here('data', 'generated_data', 'adjusted_positivity_rates',
                           paste0('posrate-', stan_run_id, '.rds'))

# model script
model_script <- ifelse(!sa, here::here('code', 'meta_analysis_re.stan'),
                       here::here('code', 'meta_analysis_re_sa.stan'))

# run model
if (!file.exists(out_est_file) | stan_redo) {
  posrate <- run_analysis_stan(
    model_script = model_script,
    dat_suspected = test_dat,
    dat_confirmed = pos_dat,
    re_ids = test_dat$study_lab,
    run_id = stan_run_id,
    coef_eqn = stan_coef_eqn,
    sens_draws = se_draws,
    spec_draws = sp_draws,
    redo = stan_redo
  )

  saveRDS(posrate, out_est_file)

} else {
  posrate <- readRDS(out_est_file)
}

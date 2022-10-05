## Source code and data for _Towards estimating true cholera burden: a systematic review and meta-analysis of Vibrio cholerae positivity_
#### Kirsten E. Wiens, Hanmeng Xu, Kaiyue Zou, John Mwaba, Justin Lessler, Espoir B. Malembaka, Maya N. Demby, Godfrey Bwire, Firdausi Qadri, Elizabeth C. Lee, Andrew S. Azman



### Instructions

1. Find data extracted from studies in the systematic review in data/extracted_data.
2. Download shapefile from the links below and save in folder data/shapefiles.

    + Global admin 0 shapefile with iso3 codes: https://www.naturalearthdata.com/downloads/50m-cultural-vectors/50m-admin-0-countries-2/

3. If running latent class analysis to generate new sensitivity/specificity estimates by test type with the JAGS program, use "launch_rstan_models.R". Instructions and JAGS settings are at the top of the script.
4. If running hierarchical meta-analysis to generate new cholera positivity rate estimates with the Stan program, also use "launch_rstan_models.R". Settings for Stan models are set at the top of that script. If running the analysis with >50 draws from the posterior estimates of sensitivity/specificity, suggest running this on a high performance computing cluster.
5. Use "main.Rmd" to run analyses that correspond to figures and tables in manuscript.


### Scripts in code folder

| File                       | Description                                                                                                                |
| :------------------------- |:---------------------------------------------------------------------------------------------------------------------------|
| latent_class_cd_meta.jags  | Latent class analysis program with conditional dependence adapted from code provided in Wang et al. (2020)                 |
| launch_rstan_models.R      | Script to prep data and run JAGS and/or Stan models with different settings indicated at top of script                     |
| main.Rmd                   | Notebook that runs manuscript analyses and creates tables/figures using previously-generated positivity rate estimates          |
| meta_analysis_re.stan      | Stan program for estimating postivity rate using hierarchical meta-analysis accounting for test performance and covariates |    
| meta_analysis_re_sa.stan   | Stan program identical to above but with increased variance on all priors, used for sensitivity analysis                 |  
| utils.R                    | Helper functions used in R scripts and Rmd notebooks                                                                      |

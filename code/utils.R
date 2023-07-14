# summarize variable categories for table
summarize_cats <- function(indat) {
  
  cats <- indat %>%
    # melt long
    reshape2::melt(measure.vars = names(indat)) %>%
    # calculate number and percent
    group_by(variable, value) %>%
    mutate(Total = length(value)) %>%
    ungroup() %>%
    group_by(variable) %>%
    mutate(Percent = round(Total/length(variable)*100,1)) %>%
    ungroup() %>%
    unique() %>%
    # prettier variable names
    mutate(variable = gsub('_', ' ', variable),
           variable = str_to_sentence(variable),
           value = ifelse(is.na(value), 'Not reported', value))
  
  # order and remove row names
  setorderv(cats, c('variable', 'value'), order = c(-1, 1))
  rownames(cats) <- c()
  
  # return
  return(cats)
}


# make simple kable table with merged first column
pretty_table <- function(intab) {
  
  tab <- intab %>%
    # make kable
    kable(escape = F, font = 12, row.names = FALSE) %>%
    # simple striped, bordered table
    kable_styling(bootstrap_options = c('bordered', 'striped', 'condensed'), 
                  position = 'left', full_width = FALSE) %>%
    # rows are black
    row_spec(1:nrow(intab), color='black') %>%
    # collapse rows
    collapse_rows(columns = 1, valign = 'top')
  
  # return
  return(tab)
}


# make a custom forest plot
forest_plot <- function(plot_dt) {
  
  # set colors
  pal <- brewer.pal(n=5, name='Dark2')
  cols <- c('Culture' = pal[1], 'PCR' = pal[2], 'RDT' = pal[3], 'Adjusted' = 'black')
  
  # make forest plot
  ggplot(plot_dt) +
    # study-level proportion by culture
    geom_point(aes(y = id, x = Culture_p, color = 'Culture'), shape = 1, size = 0.85) +
    geom_errorbarh(aes(y = id, xmin = Culture_lower, xmax = Culture_upper, color = 'Culture'), height = 0.1) +
    # study-level proportion by PCR
    geom_point(aes(y = id, x = PCR_p, color = 'PCR'), shape = 1, size = 0.85) +
    geom_errorbarh(aes(y = id, xmin = PCR_lower, xmax = PCR_upper, color = 'PCR'), height = 0.1) +
    # study-level proportion by RDT
    geom_point(aes(y = id, x = RDT_p, color = 'RDT'), shape = 1, size = 0.85) +
    geom_errorbarh(aes(y = id, xmin = RDT_lower, xmax = RDT_upper, color = 'RDT'), height = 0.1) +
    # study-level underlying positivity rates
    geom_point(aes(y = id, x = mean, color = 'Adjusted'), size = 0.95, alpha = 0.7) +
    geom_errorbarh(aes(y = id, xmin = lower, xmax = upper), height = 0.1, alpha = 0.7) +
    # scale axes
    scale_x_continuous(limits = c(0,1), breaks = c(0,0.25,0.5,0.75,1), name = 'Proportion positive') +
    scale_y_continuous(name = '', breaks = 1:(length(plot_dt$study_lab)),
                       labels = plot_dt$study_lab, 
                       trans = 'reverse') +
    # legend
    scale_color_manual(name = NULL, values = cols) +
    # rotate labels
    theme_minimal() +
    theme(panel.spacing = unit(1, 'lines'),
          legend.position = c(0.86, 0.9)) +
    theme(strip.text.y = element_text(angle = 0),
          axis.text.y = element_text(size = 8))
}


# proportion positive by grouping variable of interest
plot_posrate <- function(indat, x_col, 
                         y_col = 'prop_pos',
                         rotate_x_labs = FALSE,
                         color_col = x_col,
                         palette = 'Dark2',
                         num_colors = length(unique(indat[,color_col]))) {
  # plot
  ggplot(indat, aes(x = as.factor(get(x_col)), y = get(y_col))) +
    # boxplot without outliers
    geom_boxplot(outlier.shape = NA) +
    # points (will include the outliers)
    geom_jitter(aes(color = as.factor(get(color_col))),
                alpha = 0.75) +
    theme_bw() +
    # wrap x axis labels
    scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
    # color by grouping variable
    scale_color_manual(values=brewer.pal(n=num_colors, name=palette)) +
    # format axes
    theme(axis.text.x = element_text(angle = ifelse(rotate_x_labs, 45, 0), 
                                     hjust = ifelse(rotate_x_labs, 1, 0.5)),
          axis.text = element_text(size = 10),
          legend.position = 'none') +
    ylab('Proportion positive') + xlab('') +
    # title based on grouping variable
    ggtitle(str_to_sentence(gsub('_', ' ', x_col, )))
}


# country-level data coverage plot
data_map <- function(dt_plot, shp_plot, ad, 
                     num_max = 10,
                     plotting_prop = FALSE) {
  
  dt_plot <- dt_plot %>%
    # subset to gegraphic level
    filter(variable==ad) %>%
    # set an numbers above max num to max
    mutate(plot_var = ifelse(observations > num_max, num_max, observations))
  
  # add to shapefile
  dt_shp <- merge(shp_plot, dt_plot, by.x = 'ISO_A3', by.y = 'country_iso3', all = T)
  
  # plot
  myplot <- 
    ggplot(data = dt_shp) +
    geom_sf(colour = 'grey20', aes(fill = plot_var), show.legend = T) +
    coord_sf(xlim = c(-112, 135), ylim = c(-35, 62)) + 
    ggtitle(str_to_title(ad)) +
    theme(text = element_text(size = 14),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          rect = element_blank(),
          panel.border = element_rect(colour = 'black', fill=NA, size = 0.5)) +
    if (plotting_prop == TRUE) {
      scale_fill_viridis_c('Proportion', na.value = 'grey80', limits = c(0, 1))
    } else {
      scale_fill_viridis_c('Number', na.value = 'grey80', 
                           option = 'plasma', limits = c(1, num_max))
    } 
  return(myplot)
}


# function to make small changes to admin names, to make the name match with names in shpfiles
edit_admin_names <- function(df){
  
  df <- df %>%
    # names from VNM
    mutate(admin1 = str_replace_all(admin1, "Ha Noi", "Hà Nội"),
           admin1 = str_replace_all(admin1, "Thua Thien Hue", "Thừa Thiên Huế"),
           admin2 = str_replace_all(admin2, "Hue", "Huế"),
           admin2 = str_replace_all(admin2, "Hoang Mai", "Hoàng Mai"),
           admin2 = str_replace_all(admin2, "Thanh Xuan", "Thanh Xuân"),
           admin2 = str_replace_all(admin2, "Dong Da", "Đống Đa"),
           admin2 = str_replace_all(admin2, "Cau Giay", "Cầu Giấy")) %>%
    # admin names that were slightly different from names in shpfiles 
    mutate(admin1 = str_replace_all(admin1, "Sekong", "Xékong"),
           admin1 = str_replace_all(admin1, "Karbala", "Karbala'"),
           admin1 = str_replace_all(admin1, "al-Ta'mim", "At-Ta'mim"),
           admin1 = str_replace_all(admin1, "Kaijado", "Kajiado"),
           admin1 = str_replace_all(admin1, "Extreme-North", "Extrême-Nord"),
           admin2 = str_replace_all(admin2, "Limete", "Kinshasa"), 
           admin2 = str_replace_all(admin2, "Anse-à-Veau", "l'Anse-à-Veau"),
           admin2 = str_replace_all(admin2, "Les Gonaïves", "les Gonaïves"),
           admin2 = str_replace_all(admin2, "Grand Rivière du Nord", "Grande-Rivière du Nord"),
           admin2 = str_replace_all(admin2, "Trou du Nord", "le Trou-du-Nord"),
           admin2 = str_replace_all(admin2, "Golfe", "Golfe (incl Lomé)"),
           admin2 = str_replace_all(admin2, "Lome", "Golfe (incl Lomé)"),
           admin2 = str_replace_all(admin2, "Forecariah", "Forécariah"),
           admin2 = str_replace_all(admin2, "Dubreka,", "Dubréka"),
           admin2 = str_replace_all(admin2, "Calcutta", "Kolkata"))
  
  return(df)
}


# function to get smallest admin level for each entry
get_smallest_admin <- function(df){
  
  df_admin <- df %>%
    # which admin level is the smallest
    mutate(admin_min_level = case_when(admin1 == "" ~ 0,
                                   admin1 != "" & admin2 == "" ~ 1,
                                   admin2 != "" ~ 2)) %>%
    # get the name of the smallest admin
    mutate(admin_min = case_when(admin1 == "" ~ country_iso3,
                                 admin1 != "" & admin2 == "" ~ admin1,
                                 admin2 != "" ~ admin2))
  return(df_admin)
  
}
  

# function: read in df and return a list of shapefile which contain:
# three shapefiles for all three admin levels
match_polygon <- function(df){
  
  # country codes:
  cc <- unique(df$country_iso3)[nchar(unique(df$country_iso3)) == 3]

  # download shapefiles
  shp0 <- gadm_sf_loadCountries(cc, level=0)$sf
  shp1 <- gadm_sf_loadCountries(cc, level=1)$sf
  shp2 <- gadm_sf_loadCountries(cc, level=2)$sf
  
  # get the smallest admin level 
  df_admin <- get_smallest_admin(df)
  
  # admin level 0
  # if the entry contains more than 1 locations, separate into different rows (each row will then contains only one country)
  # when the min admin level is 0
  shp_admin0 <- right_join(shp0,
                    df_admin %>% 
                      filter(admin_min_level == 0) %>% 
                      mutate(admin_min = strsplit(as.character(admin_min), ", ")) %>%
                      unnest(admin_min) %>% as.data.frame(),
                    by = c("ISO" = "admin_min"))
  
  
  
  # when min admin level is 1
  shp_admin1 <- df_admin %>% 
    filter(admin_min_level == 1) %>%
    mutate(admin_min = strsplit(as.character(admin_min), ", ")) %>%
    unnest(admin_min)

  # find those that are not in the admin_list
  # unique(shp_admin1$admin_min)[(unique(shp_admin1$admin_min) %in% unique(shp1$NAME_1) == FALSE)]
  
  # merge at admin1
  shp_admin1 <- right_join(shp1, shp_admin1, by = c("NAME_1" = "admin_min","ISO"= "country_iso3")) %>%
    mutate(admin_min = NAME_1)
  
  # when min admin level is 2
  shp_admin2 <- df_admin %>% 
    filter(admin_min_level == 2) %>%
    mutate(admin_min = strsplit(as.character(admin_min), ", ")) %>%
    unnest(admin_min)
  
  # find those that are not in the shapefile admin2 list
  # unique(shp_admin2$admin_min)[(unique(shp_admin2$admin_min) %in% unique(shp2$NAME_2) == FALSE)]
  # admin2 that couldn't be matched: 
  # 1) Limete (in Kinshasa, changed admin2 to Kinshasa),  
  # 2) Muambe (couldn't find the matched polygon, just left this)
  # 3）Dubreka
  
  # merge at admin2
  shp_admin2 <- right_join(shp2, shp_admin2, by = c("NAME_2" = "admin_min", "ISO" = "country_iso3")) %>%
    mutate(admin_min = NAME_2)
  
  
  shp_list <- list("shp_admin0" = shp_admin0,
                   "shp_admin1" = shp_admin1,
                   "shp_admin2" = shp_admin2)
  
  # save the shp_list (need to re-run and update this file if new data come)
  saveRDS(shp_list, "data/shapefiles/shp_list.RData")

}


# this function adds a column of incidence to df (calculated from get_incidence() function above)
load_incidence <- function(df){

  df <- get_smallest_admin(df)
  df_admin0 <- df %>% filter(admin_min_level == 0) %>% mutate(inc_calc = NA)
  df_admin1 <- df %>% filter(admin_min_level == 1) %>% mutate(inc_calc = NA)
  df_admin2 <- df %>% filter(admin_min_level == 2) %>% mutate(inc_calc = NA)
  
  # load file with incidence
  inc_list <- readRDS(here::here("data", "shapefiles", "inc_list.RData"))
  inc_admin0 <- inc_list[["inc_admin0"]] 
  inc_admin1 <- inc_list[["inc_admin1"]]
  inc_admin2 <- inc_list[["inc_admin2"]]
  
  # fill in inc_calc column in df_admin0
  for(i in 1:nrow(df_admin0)){
    message(i)
    
    # if this entry has only one country
    if(str_detect(df_admin0$admin_min[i], ",") == FALSE){
      
      df_admin0$inc_calc[i] <- inc_admin0[which(inc_admin0$ISO == df_admin0$admin_min[i]),]$inc_admin0
      
    }else{ # when this entry has more than one country, get the mean inc of all places
      
      places <- strsplit(df_admin0$admin_min[i], ", ")[[1]]
      temp <- inc_admin0 %>% filter(ISO %in% places) %>% filter(inc_admin0 != 0) # remove those places that were not covered by incidence raster
      if(nrow(temp) == 0){inc = 0}
      if(nrow(temp) != 0){inc <- mean(temp$inc_admin0, na.rm = T)}
      df_admin0$inc_calc[i] <- inc
      
    }
    i = i + 1
  }
  
  # fill in inc_calc column in df_admin1
  for(i in 1:nrow(df_admin1)){
    
    # if this entry has only one admin1 area
    if(str_detect(df_admin1$admin_min[i], ",") == FALSE){
      
      temp <- inc_admin1[which(inc_admin1$NAME_1 == df_admin1$admin_min[i] & inc_admin1$ISO == df_admin1$country_iso3[i]),]
      if(nrow(temp) == 0){inc = 0}
      if(nrow(temp) != 0){inc <- mean(temp$inc_admin1, na.rm = T)}
      df_admin1$inc_calc[i] <- inc
      
    }else{ # when this entry has more than one admin1, get the mean inc of all places
      
      places <- strsplit(df_admin1$admin_min[i], ", ")[[1]]
      temp <- inc_admin1 %>% filter(NAME_1 %in% places) %>% filter(inc_admin1 != 0)
      if(nrow(temp) == 0){inc = 0}
      if(nrow(temp) != 0){inc <- mean(temp$inc_admin1, na.rm = T)}
      df_admin1$inc_calc[i] <- inc
      
    }
    i = i + 1
  }
  
  # fill in inc_calc column in df_admin2
  for(i in 1:nrow(df_admin2)){
    
    # if this entry has only one admin2 area
    if(str_detect(df_admin2$admin_min[i], ",") == FALSE){
      
      temp <- inc_admin2[which(inc_admin2$NAME_2 == df_admin2$admin_min[i] & inc_admin2$ISO == df_admin2$country_iso3[i]),]
      if(nrow(temp) == 0){inc = 0}
      if(nrow(temp) != 0){inc <- mean(temp$inc_admin2, na.rm = T)}
      df_admin2$inc_calc[i] <- inc
      
    }else{ # when this entry has more than one admin2, get the mean inc of all places
      
      places <- strsplit(df_admin2$admin_min[i], ", ")[[1]]
      temp <- inc_admin2 %>% filter(NAME_2 %in% places) %>% filter(inc_admin2 != 0)
      if(nrow(temp) == 0){inc = 0}
      if(nrow(temp) != 0){inc <- mean(temp$inc_admin2, na.rm = T)}
      df_admin2$inc_calc[i] <- inc
      
    }
    i = i + 1
  }
  
  # bind all admin levels
  df <- rbind(df_admin0, df_admin1, df_admin2) %>% arrange(entry_id)
  
  # replace all those 0 (not covered by incidence raster) with NA
  df <- df %>% mutate(inc_calc = ifelse(inc_calc == 0, NA, inc_calc))
  
  return(df)
}


#' @title Run Stan analysis
#' @description Runs the Stan meta-analysis model
#'
#' @param model_script Stan model script
#' @param dat_suspected number of suspected cases tested by test
#' @param dat_confirmed number of suspected cases confirmed by test
#' @param dat_covars covariates by study, wide format
#' @param re_ids vector of ids corresponding to the variable to use as a random effect
#' @param test_ids vector of ids corresponding to which test was used
#' @param run_id unique id for run, will be in saved file name
#' @param coef_eqn formula in character format expressing the probability of positivity on the logit scale
#' @param sens_draws vector of posterior draws of sensitivity from latent class analysis
#' @param spec_draws vector of posterior draws of specificity from latent class analysis
#' @param redo redo fit or load pre-computed posteriors if available
#' @param chains number of chains to run in sampling
#' 
#' @return a list with parameter posteriors and results
run_analysis_stan <- function(model_script,
                              dat_suspected,
                              dat_confirmed,
                              dat_covars,
                              re_ids,
                              test_ids,
                              run_id,
                              coef_eqn,
                              sens_draws,
                              spec_draws,
                              redo,
                              chains,
                              adjust_tests,
                              ...) {
  
  # Set analysis data
  dat_suspected <- as_tibble(dat_suspected)
  dat_confirmed <- as_tibble(dat_confirmed)
  dat_covars <- as_tibble(dat_covars)

  # Get test names and number
  lab_tests <- colnames(sens_draws)
  n_tests <- length(lab_tests)
  
  # Set model matrix
  X <- model.matrix(as.formula(paste("~", coef_eqn)), data = dat_covars)
  
  # Unique study ids from 1 to N
  u_re_ids <- unique(re_ids)
  re_ids <- map_dbl(re_ids, ~which(u_re_ids == .))
  
  # name of the stan output file 
  # just in case results directory isn't created already
  if(!dir.exists(here::here("data", "generated_data", "model_fits"))){
    dir.create(here::here("data", "generated_data", "model_fits"))
  }
  
  stan_out_file <- here::here("data", "generated_data", "model_fits", 
                              paste0("stan_fit_", run_id,".rds"))
  
  if (!file.exists(stan_out_file) | redo) {
    
    mod <- cmdstanr::cmdstan_model(model_script)
    
    # Run stan model with adjustments for test accuracy
    if (adjust_tests) {
      
      stan_est <- mod$sample(data = list(
        N_obs = length(u_re_ids), # random effect on each observation
        N_obs_tests = nrow(dat_suspected),
        J = n_tests,
        N_re = length(u_re_ids), # random effect on each observation
        re = re_ids,
        p_vars = ncol(X),
        X = X,
        num_test = dat_suspected$n_cases_tested,
        num_pos = dat_confirmed$n_cases_positive,
        # number of studies that used each type of test
        N_culture = length(which(test_ids == 1)), 
        N_pcr = length(which(test_ids == 2)),
        N_rdt = length(which(test_ids == 3)),
        # row IDs for the dataset that is long by test
        culture_id_long = which(test_ids == 1),
        pcr_id_long = which(test_ids == 2),
        rdt_id_long = which(test_ids == 3),
        # row IDs for the dataset that is wide by test
        culture_id_wide = which(!is.na(dat_covs$Culture)),
        pcr_id_wide = which(!is.na(dat_covs$PCR)),
        rdt_id_wide = which(!is.na(dat_covs$RDT)),
        # draws of sensitivity and specificity
        M = nrow(sens_draws),
        sens = sens_draws,
        spec = spec_draws
      ),
      chains = 4,
      parallel_chains = 4,
      iter_sampling = 2000,
      iter_warmup = 1000,
      step_size = 0.01,
      adapt_delta = 0.99,
      max_treedepth = 15,
      refresh = 100
      )
    }
    
    # Run stan model without adjustments for test accuracy
    if (!adjust_tests) {
      
      stan_est <- mod$sample(data = list(
        N_obs = length(u_re_ids), # random effect on each observation
        N_obs_tests = nrow(dat_suspected),
        N_re = length(u_re_ids), # random effect on each observation
        re = re_ids,
        p_vars = ncol(X),
        X = X,
        num_test = dat_suspected$n_cases_tested,
        num_pos = dat_confirmed$n_cases_positive,
        # number of studies that used each type of test
        N_culture = length(which(test_ids == 1)), 
        N_pcr = length(which(test_ids == 2)),
        N_rdt = length(which(test_ids == 3)),
        # row IDs for the dataset that is long by test
        culture_id_long = which(test_ids == 1),
        pcr_id_long = which(test_ids == 2),
        rdt_id_long = which(test_ids == 3),
        # row IDs for the dataset that is wide by test
        culture_id_wide = which(!is.na(dat_covs$Culture)),
        pcr_id_wide = which(!is.na(dat_covs$PCR)),
        rdt_id_wide = which(!is.na(dat_covs$RDT))
      ),
      chains = 4,
      parallel_chains = 4,
      iter_sampling = 2000,
      iter_warmup = 1000,
      step_size = 0.01,
      adapt_delta = 0.99,
      max_treedepth = 15,
      refresh = 100
      )
    }
    
    # print summary
    stan_est$cmdstan_summary()
    # save model objects
    stan_est$draws()
    stan_est$sampler_diagnostics()
    saveRDS(stan_est, stan_out_file)
  } else {
    # load pre-computed posteriors
    cat("Model not re-run. Loading pre-computed posteriors from ", stan_out_file, "\n")
    stan_est <- readRDS(stan_out_file)
  }

  # extract draws for all parameters
  draws <- stan_est$draws(format = "df")
  
  # covariate table
  cov_tbl <- draws[grep('beta', names(draws))]
  names(cov_tbl) <- colnames(X)
  cov_tbl <- exp(cov_tbl)
  
  # pos table
  p_tbl <- draws[grep('^p\\[', names(draws))]
  
  # re table
  e_tbl <- draws[grep('^eta_re\\[', names(draws))]

  # store results for all parameters
  res <- list(
    beta = cov_tbl,
    model_mtx = X,
    num_test = dat_suspected$n_cases_tested,
    num_pos = dat_confirmed$n_cases_positive,
    sens = sens_draws,
    spec = spec_draws,
    p = p_tbl,
    e = e_tbl,
    sigma_re = draws['sigma_re'],
    alpha = draws['alpha']
  )
  
  return(res)
}


# plot prior and posterior for sensitivity/specificity draws from the JAGS model
plot_prior_post <- function(post_draws, 
                            x_lab,
                            arg1 = 1,
                            arg2 = 1,
                            dprior = 'dbeta') {
  # purples and greens
  pal <- brewer.pal(11, 'PRGn')
  
  # add plot
  ggplot() + 
    stat_function(fun = get(dprior), args = list(arg1, arg2), geom = "area", 
                  fill = pal[9], alpha = 0.2) + 
    stat_function(fun = get(dprior), args = list(arg1, arg2), 
                  aes(color = pal[9])) +
    geom_density(aes(post_draws, color = pal[2])) + 
    scale_color_manual(values = c(pal[9], pal[2]), 
                       name = NULL, labels = c('Prior', 'Posterior')) +
    theme_classic() +
    xlab(x_lab) + ylab('Density') + 
    theme(legend.position = c(0.2, 0.7),
          legend.background=element_blank())
}


# calculate 95% confidence interval for a proportion using binomial probability
# formula here: https://www.statology.org/binomial-confidence-interval-r/
binomial_ci <- function(p, n, a = 0.05) {
  ci <- p + c(-qnorm(1-a/2), qnorm(1-a/2))*sqrt((1/n)*p*(1-p))
  return(ci)
}

# calculate variance of a proportion
prop_var <- function(p, n) {
  p*(1-p)/n
}

# create a data frame with the number of observations by study methods
# for use in post-stratification later
obs_by_cat <- function(dat, covs) {
  
  # get study observations
  dat <- dat %>%
    dplyr::select(c('study_id', 'country_iso3', 'cases_min_age_bin',
                    'sampling_quality', 'surveillance_type')) %>%
    unique()

  # count number by covariates
  props <- dat %>%
    dplyr::count(across(all_of(covs))) %>%
    mutate(prop = n/nrow(dat)) %>% 
    replace(is.na(.), 0)
  
  # arrange alphabetically to make sure in same order as covs
  props <- props %>%
    dplyr::select(c(sort(covs), 'n', 'prop'))
  
  # convert to 0/1 values
  if ('sampling_quality' %in% covs) {
  props <- props %>%
    mutate(sampling_quality = ifelse(sampling_quality == 'High', 0, 1))
  }
  if ('surveillance_type' %in% covs) {
    props <- props %>%
      mutate(surveillance_type = ifelse(surveillance_type == 'Outbreak', 1, 0))
  }
  
  return(props)
}

# stratify estimates of V. cholerae positivity
# using alpha, beta, and integrating over random effects,
# nb: also post-stratifies
# assuming that the proportion of all potential studies that use different methods
# and case definitions match the proportions we found in the systematic review
# but we do not include the post-stratified results because they are not meaningful
stratify <- function(mod, cov_cats) {
  
  # load model draws
  draws <- readRDS(here::here('data', 'generated_data', 'adjusted_positivity_rates',
                              paste0('posrate-', mod, '.rds')))
  
  # extract parameters
  beta <- cbind(draws[grep('alpha', names(draws))],
                draws[grep('beta', names(draws))][[1]])
  beta$`(Intercept)` <- NULL # using alpha intercept
  sigma <- draws[grep('sigma_re', names(draws))][[1]][[1]]
  
  # covariate matrix
  mod_eqns <- list(
    'mod1' = '1',
    'mod1-unadj' = '1',
    'mod1-sa' = '1',
    'mod2' = 'sampling_quality + cases_min_age_bin',
    'mod3' = 'sampling_quality + cases_min_age_bin + surveillance_type',
    'mod4' = 'sampling_quality + surveillance_type'
  )
  coef_eqn <- mod_eqns[[mod]]
  cat_mat <- cov_cats %>%
    model.matrix(as.formula(paste("~", coef_eqn)), data = .)
  
  # compute positivity estimates by age and sampling categories
  cl <- parallel::makeCluster(8)
  doParallel::registerDoParallel(cl)
  
  cat_p <- foreach(i = 1:nrow(cat_mat),
                   .combine = rbind,
                   .inorder = F,
                   .packages = c("tidyverse", "foreach")) %dopar%
    {
      foreach(j = 1:nrow(beta),
              .combine = rbind,
              .inorder = T) %do%
        {
          
          # compute positivity integrating across in observation random effects
          pos <- integrate(function(x) {
            plogis(qnorm(
              x, as.matrix(beta[j, , drop = F]) %*% t(cat_mat[i, , drop = F]),
              sigma[j]
            ))
          }, 0, 1)[[1]]
          
          # save results
          tibble(
            strat = i,
            prop = cov_cats$prop[i],
            pos = pos
          ) %>%
            mutate(sim = j)
          
      }
  }
  
  parallel::stopCluster(cl)
  
  # overall positivity by draw, weighted by proportion of studies in each strata
  p_est <- cat_p %>%
    group_by(sim) %>%
    summarize(pos = weighted.mean(pos, w = prop)) %>%
    ungroup()
  
  # mean positivity by strata
  strata_summary <- cat_p %>%
    group_by(strat) %>%
    summarize(mean = mean(pos),
              lower = quantile(pos, 0.025),
              upper = quantile(pos, 0.975)) %>%
    ungroup() %>%
    cbind(cat_mat)
  
  # end function
  return(list(p_est, strata_summary))
}

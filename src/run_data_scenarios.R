set.seed(453453)


path.export <- 'figures/'

if (!dir.exists(path.export)) {
  dir.create(path.export)
}

# options(mc.cores = 4)
options(clustermq.scheduler="multiprocess")
sampling_options <- list(
  cores = 4,
  chains = 4,
  iter = 10000,
  warmup = 4000
)



# Define trial data and prior -----------------------------------------------------------------
source("src/example_model.R")
source("src/combo2-thall2003.R")


# Define blrm_trial objects -------------------------------------------------------------------
trial_independent_template <- blrm_trial(
  data = all_hist_trial_data,
  dose_info = trial_dose_info,
  drug_info = drug_info,
  formula_generator = function(ref_doses) {
    # Only include 1-st order terms (no interactions)
    blrm_formula_linear(ref_doses, max_interaction_level = 1) 
  }
)

trial_linear_template <- blrm_trial(
  data = all_hist_trial_data,
  dose_info = trial_dose_info,
  drug_info = drug_info,
  
  # Default (second-order) linear interaction model
  formula_generator = blrm_formula_linear
)

trial_saturating_template <- blrm_trial(
  data = all_hist_trial_data,
  dose_info = trial_dose_info,
  drug_info = drug_info,
  
  # Proposed (second-order) saturating interaction model
  formula_generator = blrm_formula_saturating
)


# Define starting doses -----------------------------------------------------------------------
# Dosings with id 0 are not started yet (negative ones have finished)
starting_doses <- filter(
  summary(trial_independent_template, "dose_info"), 
  starting_dose == TRUE
)

# Use 4 chains and 4 cores




# Define combo doses --------------------------------------------------------


q_fun <- function(dose, mu_sd_inter, drug_info, trial_independent_template, trial_linear_template, trial_saturating_template, sampling_options, path.export) {
  source("src/install_packages.R")
  source("src/example_model.R")
  source("src/combo2-thall2003.R")
  source("src/plotting_functions.R")
  
  set.seed(453453)
  dose_A <- dose
  dose_B <- dose
  
  sd_inter_slope_thall <- ifelse(mu_sd_inter == 0.5, 0.5, 1.0)
  trial_thall2003_prior <- set_thall2003_prior(
    sd_inter_slope = sd_inter_slope_thall, 
    sd_inter_inter = sqrt(2*2^2)*sd_inter_slope_thall,
    drug_info = drug_info,
    sampling_options = sampling_options
  ) 
  
  # Define prior ------------------------------------------------------------
  trial_independent <- set_prior(trial_independent_template, mu_sd_inter, sampling_options) # , chains = chains
  trial_linear <- set_prior(trial_linear_template, mu_sd_inter, sampling_options) # , chains = chains
  trial_saturating <- set_prior(trial_saturating_template, mu_sd_inter, sampling_options) # , chains = chains
  trial_thall2003 <- update_thall2003(trial_thall2003_prior, add_data = all_hist_trial_data, drug_info = drug_info)
  
  # Trial objects that simulate from prior ----------------------------------
  trial_independent_prior <- update(trial_independent, prior_PD = TRUE)
  trial_linear_prior <- update(trial_linear, prior_PD = TRUE)
  trial_saturating_prior <- update(trial_saturating, prior_PD = TRUE)
  
  # Define a toxic doublet data scenario ------------------------------------
  new_toxic_data <- filter(summary(trial_independent, "dose_info"), A == dose_A, B == dose_B)
  new_toxic_data$num_patients <- 5
  new_toxic_data$num_toxicities <- 5
  
  trial_independent_toxic <- update(trial_independent, add_data = new_toxic_data)
  trial_linear_toxic <- update(trial_linear, add_data = new_toxic_data)
  trial_saturating_toxic <- update(trial_saturating, add_data = new_toxic_data)
  trial_thall2003_toxic <- update_thall2003(trial_thall2003, add_data = new_toxic_data, drug_info = drug_info)
  
  # Define a non-toxic doublet data scenario --------------------------------
  new_non_toxic_data <- new_toxic_data
  new_non_toxic_data$num_toxicities <- 0
  
  trial_independent_non_toxic <- update(trial_independent, add_data = new_non_toxic_data)
  trial_linear_non_toxic <- update(trial_linear, add_data = new_non_toxic_data)
  trial_saturating_non_toxic <- update(trial_saturating, add_data = new_non_toxic_data)
  trial_thall2003_non_toxic <- update_thall2003(trial_thall2003, add_data = new_non_toxic_data, drug_info = drug_info)
  
  
  # Summarize data ----------------------------------------------------------
  data_prior <- 
    add_column(summary(trial_independent_prior, "dose_prediction", prob = c(0.5, 0.95)), type="independent") %>%
    bind_rows(add_column(summarize_thall2003(trial_thall2003_prior, candidate_doses = trial_dose_info, drug_info = drug_info), type="Thall (2003)")) %>% 
    bind_rows(add_column(summary(trial_linear_prior, "dose_prediction", prob = c(0.5, 0.95)), type="linear")) %>%
    bind_rows(add_column(summary(trial_saturating_prior, "dose_prediction", prob = c(0.5, 0.95)), type="saturating"))
  
  
  data_trial <- 
    add_column(summary(trial_independent, "dose_prediction", prob = c(0.5, 0.95)), type="independent") %>%
    bind_rows(add_column(summarize_thall2003(trial_thall2003, candidate_doses = trial_dose_info, drug_info = drug_info), type="Thall (2003)")) %>%
    bind_rows(add_column(summary(trial_linear, "dose_prediction", prob = c(0.5, 0.95)), type="linear")) %>%
    bind_rows(add_column(summary(trial_saturating, "dose_prediction", prob = c(0.5, 0.95)), type="saturating")) 
  
  
  data_toxic <- 
    add_column(summary(trial_independent_toxic, "dose_prediction", prob = c(0.5, 0.95)), type="independent") %>%
    bind_rows(add_column(summarize_thall2003(trial_thall2003_toxic, candidate_doses = trial_dose_info, drug_info = drug_info), type="Thall (2003)")) %>%
    bind_rows(add_column(summary(trial_linear_toxic, "dose_prediction", prob = c(0.5, 0.95)), type="linear")) %>%
    bind_rows(add_column(summary(trial_saturating_toxic, "dose_prediction", prob = c(0.5, 0.95)),  type="saturating")) 
  
  
  data_non_toxic <- 
    add_column(summary(trial_independent_non_toxic, "dose_prediction", prob = c(0.5, 0.95)), type="independent") %>%
    bind_rows(add_column(summarize_thall2003(trial_thall2003_non_toxic, candidate_doses = trial_dose_info, drug_info = drug_info), type="Thall (2003)")) %>%
    bind_rows(add_column(summary(trial_linear_non_toxic, "dose_prediction", prob = c(0.5, 0.95)), type="linear")) %>%
    bind_rows(add_column(summary(trial_saturating_non_toxic, "dose_prediction", prob = c(0.5, 0.95)),  type="saturating")) 
  
  
  # Plot level plots after updating -------------------------------------------------------------
  make_levelplots_ggplot(data_prior, "Prior", mu_sd_inter, path.export)
  make_levelplots_ggplot(data_trial, "Historical data", mu_sd_inter, path.export)
  make_levelplots_ggplot(data_toxic, paste0("Historical data + 5/5 DLTs at A=", dose_A, ", B=", dose_B), mu_sd_inter, path.export)
  make_levelplots_ggplot(data_non_toxic, paste0("Historical data + 0/5 DLTs at A=", dose_A, ", B=", dose_B), mu_sd_inter, path.export)
  
  # Inspect marginals after updating ------------------------------------------------------------
  plot_marginals(data_prior, "Prior", mu_sd_inter, path.export)
  plot_marginals(data_trial, "Historical data", mu_sd_inter, path.export)
  plot_marginals(data_toxic, paste0("Historical data + 5/5 DLTs at A=", dose_A, ", B=", dose_B), mu_sd_inter, path.export)
  plot_marginals(data_non_toxic, paste0("Historical data + 0/5 DLTs at A=", dose_A, ", B=", dose_B), mu_sd_inter, path.export)
}

scenarios <- expand_grid(dose = c(100, 200), mu_sd_inter = c(0.5, 1.5) )
clustermq::Q_rows(
  df=scenarios, fun = q_fun, n_jobs = as.numeric(ceiling(parallelly::availableCores() / 4)), 
  const = list(
    trial_independent_template = trial_independent_template, 
    trial_linear_template = trial_linear_template,
    trial_saturating_template = trial_saturating_template,
    drug_info = drug_info, sampling_options = sampling_options,
    path.export = path.export
  )
)

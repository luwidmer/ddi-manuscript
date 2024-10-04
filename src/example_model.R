sort_drug_columns <- function(df) {
  initial_cols <- c("group_id")
  if (any(colnames(df) == "stratum_id")) {
    initial_cols <- c(initial_cols, "stratum_id")
  }
  
  bind_cols(
    select(df, one_of(initial_cols)), 
    select(df, -one_of(c(initial_cols, "num_patients", "num_toxicities"))), 
    select(df, "num_patients", "num_toxicities")
  )
}



drug_info <- tibble(
  drug_name = c("A", "B"),
  dose_ref  = c(200, 200),
  dose_unit = "mg/day"
)

# Define historical data
hist_data <- function()
{
  # Define historical trials ------------------------------------------------
  # Hypothetical single drug A study
  hist_trial_A <- tibble(
    group_id       = "single_trial",
    A              = c(50, 100, 200, 300, 400, 600),
    num_patients   = 10,
    num_toxicities = c( 0,   1,   1,   2,   3,   6)
  )
  hist_trial_A
  
  # Hypothetical single drug B study
  hist_trial_B <- tibble(
    group_id       = hist_trial_A$group_id,
    B              = hist_trial_A$A,
    num_patients   = 5,
    num_toxicities = c( 0,   0,   1,   1,   1,   3)
  )
  hist_trial_B
  
  # Build data frame for historical data, and 
  # re-order so all doses come before num_patients / num_toxicities
  sort_drug_columns(bind_rows_0(hist_trial_A, hist_trial_B))
}
all_hist_trial_data <- hist_data()



# Define new trial --------------------------------------------------------
dose_info <- function()
{
  trial_test_AB <- expand_grid(
    group_id = "single_trial",
    A = c(0, 50, 100, 200, 300, 400, 500, 600),
    B = c(0, 50, 100, 200, 300, 400, 500, 600)
  )
  
  # Create dosing tibble --------------------------------------------------
  trial_dose_info <- trial_test_AB
  
  # Define starting doses for the trial -----------------------------------
  trial_dose_info <- mutate(
    trial_dose_info, 
    starting_dose = (A == 100 & B == 100)
  )
  
  trial_dose_info
}
trial_dose_info <- dose_info()


set_prior <- function(trial, mu_sd_inter, sampling_options) 
{
  dims <- summary(trial, "dimensionality")
  
  do.call(
    update, 
    args = c(list(
      trial,
      # Prior mean and sd on mu_{alpha_i}, mu_{beta_i}
      prior_EX_mu_mean_comp  = matrix(c(logit(0.10), 0), nrow = dims$num_components, ncol = 2, TRUE),
      prior_EX_mu_sd_comp    = matrix(c(2, 1), nrow = dims$num_components, ncol = 2, TRUE),
      
      # Prior mean and sd on tau_{alpha_{s,i}}, tau{beta_{s,i}} 
      prior_EX_tau_mean_comp = abind(
        matrix(c(0, 0), nrow = dims$num_components, ncol = 2, TRUE), # add more strata here if needed
        along = 0
      ),
      prior_EX_tau_sd_comp = abind(
        matrix(0, nrow = dims$num_components, ncol = 2, TRUE), # add more strata here if needed
        along = 0
      ),
      
      # Prior mean and sd on mu_{eta}
      prior_EX_mu_mean_inter  = rep(0,   dims$num_interaction_terms),
      prior_EX_mu_sd_inter    = rep(mu_sd_inter, dims$num_interaction_terms), 
      prior_EX_tau_mean_inter = matrix(0, nrow = dims$num_strata, ncol = dims$num_interaction_terms),
      prior_EX_tau_sd_inter   = matrix(0, nrow = dims$num_strata, ncol = dims$num_interaction_terms),
      
      prior_is_EXNEX_comp = rep(FALSE, dims$num_components),
      prior_EX_prob_comp = matrix(1, nrow = dims$num_groups, ncol = dims$num_components),
      
      prior_is_EXNEX_inter = rep(FALSE, dims$num_interaction_terms),
      prior_EX_prob_inter = matrix(1, nrow = dims$num_groups, ncol = dims$num_interaction_terms),
      
      prior_tau_dist = 0), sampling_options)
    )
}
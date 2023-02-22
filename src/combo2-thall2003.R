
prepare_dummy_data <- function(drug_info) 
{
  dref_A <- filter(drug_info, drug_name == "A")$dose_ref
  dref_B <- filter(drug_info, drug_name == "B")$dose_ref

  drug_A <- 0.01
  drug_B <- 0.01
  
  dummy_data <- tibble(
    group_id = factor("trial_nogroup"),
    drug_A = c(drug_A),
    drug_B = c(drug_B),
    std_drug_A = drug_A/dref_A,
    std_drug_B = drug_B/dref_B,
    num_patients = c(0), 
    num_toxicities = c(0),
    cohort_time = 0
  )

  return(dummy_data)
}

thall2003_model <-
  bf(
    num_toxicities | trials(num_patients) ~ 
    log(
      exp(interA) * (std_drug_A ^ slopeA) + 
      exp(interB) * (std_drug_B ^ slopeB) + 
      #exp(interA + interB) * ((std_drug_A ^ slopeA) * (std_drug_B ^ slopeB)) + 0 * interAB * slopeAB # no interaction case
      exp(interAB) * ((std_drug_A ^ slopeA) * (std_drug_B ^ slopeB)) ^ slopeAB
    )
    , nl = TRUE
  ) +
  lf(interA ~ 1) +
  lf(slopeA ~ 1) +
  lf(interB ~ 1) +
  lf(slopeB ~ 1) +
  lf(interAB ~ 1) +
  lf(slopeAB ~ 1) 

thall2003_prior <- function(sd_inter_inter = sqrt(2*2^2), sd_inter_slope = 0.5) {
  prior_string(paste0("normal(", logit(0.10), ", 2)"), nlpar="interA") +
  prior(lognormal(0, 1), nlpar="slopeA", lb = 0) +
  prior_string(paste0("normal(", logit(0.10), ", 2)"), nlpar="interB") +
  prior(lognormal(0, 1), nlpar="slopeB", lb = 0) +
  prior_string(paste0("normal(", 2*logit(0.10), ", ", sd_inter_inter, ")"), nlpar="interAB") +
  prior_string(paste0("lognormal(0, ", sd_inter_slope, ")"), nlpar="slopeAB", lb = 0)
}

# combo2_fit <- brm(
#   thall2003_model,
#   data = codata_combo2_nogroup_std,
#   family = binomial,
#   prior = thall2003_prior(),
#   chains=4,
#   file="brms_thal2003",
#   file_refit="on_change",
#   iter = 10000
# )

set_thall2003_prior <- function(drug_info, iter, warmup, chains, ...) {
  brm(
    thall2003_model,
    data = prepare_dummy_data(drug_info),
    family = binomial,
    prior = thall2003_prior(...),
    file="brms_thall2003",
    file_refit="always",
    chains = chains,
    iter = iter,
    warmup = warmup,
    control = list(adapt_delta = 0.9)
  )
}

# thall2003_trial <- function(data, dose_info, drug_info) {
#   dref_A <- filter(drug_info, drug_name == "A")$dose_ref
#   dref_B <- filter(drug_info, drug_name == "B")$dose_ref
#   
#   data_brms <- mutate(
#     data, 
#     std_drug_A = A/dref_A,
#     std_drug_B = B/dref_B
#   )
#   
#   brm(
#     thall2003_model,
#     data = data,
#     family = binomial,
#     prior = thall2003_prior(...),
#     chains = chains,
#     file="brms_thall2003",
#     file_refit="on_change",
#     iter = iter,
#     warmup = warmup
#   )
# }

update_thall2003 <- function(brmfit_object, add_data, drug_info)
{
  dref_A <- filter(drug_info, drug_name == "A")$dose_ref
  dref_B <- filter(drug_info, drug_name == "B")$dose_ref
  
  add_data_brms <- mutate(
    add_data, 
    std_drug_A = A/dref_A,
    std_drug_B = B/dref_B
  )
  
  update(brmfit_object, newdata = bind_rows(brmfit_object$data, add_data_brms))
}

summarize_thall2003 <- function(brmfit_object, candidate_doses, drug_info, prob_intervals = c(0, 0.16, 0.33, 1), probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) 
{
  dref_A <- filter(drug_info, drug_name == "A")$dose_ref
  dref_B <- filter(drug_info, drug_name == "B")$dose_ref
  
  candidate_doses_no_0 <- candidate_doses %>% mutate(
    A = pmax(A, 1e-6),
    B = pmax(B, 1e-6)
  )
  
  candidate_doses_brms <- candidate_doses_no_0 %>% mutate(
    std_drug_A = A/dref_A,
    std_drug_B = B/dref_B,
    num_toxicities = 0,
    num_patients = 0
  )
  
  raw_draws <- posterior_linpred(brmfit_object, newdata = candidate_doses_brms, transform = TRUE)
  
  output <- candidate_doses %>% bind_cols(
    as_tibble(posterior_summary(raw_draws, probs = probs))
  )
  
  draws_matrix <- posterior::as_draws_matrix(raw_draws)
  
  interval_probs <- posterior::summarize_draws(
    draws_matrix,
    list(
      function(x){
        prop.table(table(cut(x, breaks = prob_intervals)))
      }
    )
  )
  
  interval_probs <- interval_probs %>% rename(
    "prob_underdose" = `(0,0.16]`,
    "prob_target" = `(0.16,0.33]`,
    "prob_overdose" = `(0.33,1]`
  ) %>% select(-variable)
  
  
  return(
    output %>%
    bind_cols(interval_probs) %>%
    mutate(
      ewoc_ok = prob_overdose <= 0.25
    ) %>% 
    rename(
      "mean" = "Estimate",
      "sd" = "Est.Error",
      "2.5%" = "Q2.5",
      "25%" = "Q25",
      "50%" = "Q50",
      "75%" = "Q75",
      "97.5%" = "Q97.5"
    )
  )
  
}

# foo <- summarize_thall2003(combo2_fit, codata_combo2_nogroup_std)

# combo2_fit
# 
# posterior_summary(combo2_fit)
# 
# mcplot <- brms::mcmc_plot(combo2_fit)
# print(mcplot)
# 
# print(
#   codata_combo2_nogroup_std %>%
#     bind_cols(
#       as_tibble(posterior_summary(posterior_linpred(combo2_fit, transform = T), ))
#     ) %>%
#     mutate(ewoc_ok = Q75 < 0.33) %>%
#     kable(digits = 3)
# )

library(OncoBayes2)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(posterior)
library(bayesplot)
library(tidybayes)
library(egg)
dose_info <- tibble(
  group_id = "single_group",
  stratum_id = "all",
  drug_A = c(10, 25, 50, 100, 200, 400, 800)
)

drug_info <- tibble(
  drug_name="drug_A",
  dose_ref = 120,
  dose_unit = "AU"
)

demo_blrm_trial_init <- blrm_trial(
  data = NULL, dose_info = dose_info, drug_info = drug_info
)

demo_blrm_trial <- demo_blrm_trial_init |> update(
  prior_EX_mu_mean_comp = matrix(
    c(logit(0.2), # mean of intercept on logit scale
      log(1)),    # mean of log-slope on logit scale
    nrow = 1,
    ncol = 2
  ),
  prior_EX_mu_sd_comp = matrix(
    c(2,  # sd of intercept
      log(4)/1.96), # sd of log-slope
    nrow = 1,
    ncol = 2
  ),
  ## Here we take tau as known and as zero.
  ## This disables the hierarchical prior which is
  ## not required in this example as we analyze a
  ## single trial.
  prior_EX_tau_mean_comp = matrix(
    c(0, 0),
    nrow = 1,
    ncol = 2
  ),
  prior_EX_tau_sd_comp = matrix(
    c(1, 1),
    nrow = 1,
    ncol = 2
  ),
  prior_EX_prob_comp = matrix(1, nrow = 1, ncol = 1),
  prior_tau_dist = 0,
  prior_PD = FALSE
)

demo_blrm_trial_wide <- demo_blrm_trial_init |> update(
  prior_EX_mu_mean_comp = matrix(
    c(logit(0.2), # mean of intercept on logit scale
      log(1)),    # mean of log-slope on logit scale
    nrow = 1,
    ncol = 2
  ),
  prior_EX_mu_sd_comp = matrix(
    c(2,  # sd of intercept
      2), # sd of log-slope
    nrow = 1,
    ncol = 2
  ),
  ## Here we take tau as known and as zero.
  ## This disables the hierarchical prior which is
  ## not required in this example as we analyze a
  ## single trial.
  prior_EX_tau_mean_comp = matrix(
    c(0, 0),
    nrow = 1,
    ncol = 2
  ),
  prior_EX_tau_sd_comp = matrix(
    c(1, 1),
    nrow = 1,
    ncol = 2
  ),
  prior_EX_prob_comp = matrix(1, nrow = 1, ncol = 1),
  prior_tau_dist = 0,
  prior_PD = FALSE
)

add_risk_rvar <- function(newdata, model) {
  newdata |> mutate(
    risk = rvar(posterior_linpred(model$blrmfit,
                                  newdata=newdata,
                                  allow_new_levels=FALSE,
                                  transform = T))
  )
}

label_percent <- function (accuracy = NULL, scale = 100, prefix = "", suffix = "%",
                           big.mark = " ", decimal.mark = ".", trim = TRUE,
                           ...) {
  scales::number_format(accuracy = accuracy, scale = scale, prefix = prefix,
                        suffix = suffix, big.mark = big.mark, decimal.mark = decimal.mark,
                        trim = trim, ...)
  
}


plot_curve_dist <- function(newdata, dose_info, model) {
  doses <- dose_info$drug_A
  newdata |> 
    add_risk_rvar(model) |>
    unnest_rvars() |>
    ggplot(aes(x = drug_A, y = risk, group = .draw)) +
    geom_line( alpha=0.01) +
    scale_x_continuous(breaks = doses, trans = "log10") +
    hline_at(c(0.16, 0.33), linetype = I(2)) +
    xlab("Dose [mg]") +
    ylab("P(DLT)") +
    bayesplot::bayesplot_theme_get() +
    scale_y_continuous(breaks = sort(c((0:10) / 10, 0.16, 0.33)),
                         minor_breaks = NULL,
                         limits = c(0, 1),
                         labels = label_percent(accuracy = 1),
                         oob = scales::squish)
  
}

plot_toxicity_curve_nolabel <- function(..., dose_info, additional_text = "") {
  plot_toxicity_curve(...) + theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) + 
    scale_x_continuous(breaks = dose_info$drug_A, trans = "log10") +
    xlab("Dose [mg]") +
    ggtitle(
      label = "DLT Risk",
      subtitle = additional_text
    ) + legend_none()
}


newdata <- tibble(drug_A = exp(seq(from = log(10), to = log(800), length.out = 101)), group_id = "single_group", stratum_id="all")


pdist_demo <- plot_curve_dist(newdata, dose_info, demo_blrm_trial) 
pdist_demo_wide <- plot_curve_dist(newdata, dose_info, demo_blrm_trial_wide)

p_demo <- plot_toxicity_curve_nolabel(demo_blrm_trial, dose_info = dose_info, additional_text = "log(β) ~ N(1, 0.7)\n ")
p_demo_wide <- plot_toxicity_curve_nolabel(demo_blrm_trial_wide, dose_info = dose_info, additional_text = "log(β) ~ N(1, 2)\n ")

new_data <- tibble(group_id = "single_group", drug_A = c(10, 25, 50, 100), num_patients = 3, num_toxicities=0)

demo_blrm_trial_wide_upd <- demo_blrm_trial_wide |> update(add_data = new_data)
demo_blrm_trial_upd <- demo_blrm_trial |> update(add_data = new_data)

p_demo_upd <- demo_blrm_trial_upd |> plot_toxicity_curve_nolabel(additional_text = "log(β) ~ N(1, 0.7),\n0/3 DLTs at 10, 25, 50, 100 mg", dose_info = dose_info)  
pdist_demo_upd <- plot_curve_dist(newdata, dose_info, demo_blrm_trial_upd)

p_demo_wide_upd <- plot_toxicity_curve_nolabel(demo_blrm_trial_wide_upd, additional_text = "log(β) ~ N(1, 2),\n0/3 DLTs at 10, 25, 50, 100 mg", dose_info = dose_info)  
pdist_demo_wide_upd <- plot_curve_dist(newdata, dose_info, demo_blrm_trial_wide_upd)
summary(demo_blrm_trial_wide_upd, summarize="dose_prediction")

result <- grid.arrange(p_demo, p_demo_upd, p_demo_wide, p_demo_wide_upd, pdist_demo, pdist_demo_upd, pdist_demo_wide, pdist_demo_wide_upd, nrow = 2, ncol = 4)
ggsave("figures/Zhang2022.png", result, scale = 1)

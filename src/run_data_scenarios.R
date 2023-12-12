set.seed(453453)


path.export <- 'figures/'

if (!dir.exists(path.export)) {
  dir.create(path.export)
}

options(mc.cores = 4)
chains <- 4
iter <- 10000
warmup <- 4000


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



make_levelplots_ggplot <- function (data, title)
{
  title_text <- title
  # * `prob_overdose` * `ewoc_ok`
  data$A <- as.factor(data$A)
  data$B <- as.factor(data$B)
  data$ewoc_fail <- !data$ewoc_ok
  data$type <- factor(data$type, levels=sort(unique(data$type))[c(1,4,2,3)])
  
  vars <- c("ewoc_fail", "prob_overdose", "prob_target", "prob_underdose")
  var_labels <- c("EWOC failure", "P(overdose)", "P(target)", "P(underdose)")
  data <- select(data, "A", "B", "type", vars)

  data_long <- data %>% 
    pivot_longer(cols=vars, names_to="variable", values_to="Probability")
  
  data_long$variable <- case_when(
    data_long$variable == "ewoc_fail" ~ "EWOC failure",
    data_long$variable == "prob_overdose" ~ "P(overdose)",
    data_long$variable == "prob_target" ~ "P(target)",
    data_long$variable == "prob_underdose" ~ "P(underdose)",
    TRUE ~ "ERROR"
  )
  
  data_long$variable <- factor(data_long$variable, levels=var_labels)
  data_long <- arrange(data_long, variable)    
  
  foo <- ggplot(data_long, aes(A, B, fill = `Probability`)) + 
    theme_minimal() +
    facet_grid(variable ~ type) +
    geom_raster() + 
    coord_fixed() +
    # scale_fill_viridis(option="B") +
    scale_fill_gradientn(colours=c("#2166acFF","#f7f7f7FF","#b2182bFF")) +
    #geom_contour(colour = "black", breaks=c(0,0.25,1), aes(z = `Probability`)) +
    theme(axis.text.x = element_text(angle = 45)) +
    xlab("Drug 1 (dose)") + ylab("Drug 2 (dose)") + labs(title=title_text)
  
  filename <- paste0('ggplot_mu_sd_inter_', mu_sd_inter, "_", title)
  
  ggsave(
    filename = paste0(path.export, str_replace_all(filename,"[[:punct:]\\s]+","_"), ".pdf"),
    plot = foo
  )
  graphics.off()
  
  return(NULL)
}


plot_marginals <- function(data, title)
{
  data$type <- factor(data$type, levels=sort(unique(data$type))[c(1,4,2,3)])
  
  p1 <- ggplot(filter(data, A == 0), aes(x=factor(B), colour=ewoc_ok)) +
    facet_wrap(~type, labeller=label_both) +
    scale_y_continuous(breaks=c(0, 0.16, 0.33, 0.4, 0.6, 0.8, 1.0)) +
    coord_cartesian(ylim=c(0,0.8)) +
    geom_hline(yintercept = c(0.16, 0.33), linetype = "dotted") +
    geom_pointrange(aes(y=`50%`, ymin=`2.5%`, ymax=`97.5%`)) +
    geom_linerange(aes(ymin=`25%`, ymax=`75%`), linewidth=1.5) +
    #ggtitle(paste0("DLT Probability, ", title), "Shown is the median (dot), 50% CrI (thick line) and 95% CrI (thin line)") +
    ylab(NULL) + 
    xlab("Dose Drug 2 [mg]")
  
  p2 <- ggplot(filter(data, B == 0), aes(x=factor(A), colour=ewoc_ok)) +
    facet_wrap(~type, labeller=label_both) +
    scale_y_continuous(breaks=c(0, 0.16, 0.33, 0.4, 0.6, 0.8, 1.0)) +
    coord_cartesian(ylim=c(0,0.8)) +
    geom_hline(yintercept = c(0.16, 0.33), linetype = "dotted") +
    geom_pointrange(aes(y=`50%`, ymin=`2.5%`, ymax=`97.5%`)) +
    geom_linerange(aes(ymin=`25%`, ymax=`75%`), linewidth=1.5) +
    #ggtitle(paste0("DLT Probability, ", title), "Shown is the median (dot), 50% CrI (thick line) and 95% CrI (thin line)") +
    ylab(NULL) + 
    xlab("Dose Drug 1 [mg]")
  
  filename <- paste0('marginal_mu_sd_inter_', mu_sd_inter, "_", title)
  #png(paste0(path.export, str_replace_all(filename,"[[:punct:]\\s]+","_"), ".png"), width=1000, height=1000, res=100)

  plot <- arrangeGrob(p1, p2, nrow = 1)
  ggsave(
    filename = paste0(path.export, str_replace_all(filename,"[[:punct:]\\s]+","_"), ".pdf"),
    plot = plot,#    scale = 2
    width = 16,
    height = 5
  )
  graphics.off()
  
  # grid.draw(plot)
  return(NULL)
}


# Define combo doses --------------------------------------------------------


for (dose in c(100, 200)) {
  dose_A <- dose
  dose_B <- dose
  for (mu_sd_inter in c(0.5, 1.5)) 
  {    
    sd_inter_slope_thall <- ifelse(mu_sd_inter == 0.5, 0.5, 1.0)
    trial_thall2003_prior <- set_thall2003_prior(
      sd_inter_slope = sd_inter_slope_thall, 
      sd_inter_inter = sqrt(2*2^2)*sd_inter_slope_thall,
      drug_info = drug_info,
      iter = iter,
      warmup = warmup,
      chains = chains
    ) 
    
    
    
    
    # Define prior ------------------------------------------------------------
    trial_independent <- set_prior(trial_independent_template, mu_sd_inter, iter = iter, warmup = warmup, chains = chains)
    trial_linear <- set_prior(trial_linear_template, mu_sd_inter, iter = iter, warmup = warmup, chains = chains)
    trial_saturating <- set_prior(trial_saturating_template, mu_sd_inter, iter = iter, warmup = warmup, chains = chains)
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
    make_levelplots_ggplot(data_prior, "Prior")
    make_levelplots_ggplot(data_trial, "Historical data")
    make_levelplots_ggplot(data_toxic, paste0("Historical data + 5/5 DLTs at A=", dose_A, ", B=", dose_B))
    make_levelplots_ggplot(data_non_toxic, paste0("Historical data + 0/5 DLTs at A=", dose_A, ", B=", dose_B))
    
    # Inspect marginals after updating ------------------------------------------------------------
    plot_marginals(data_prior, "Prior")
    plot_marginals(data_trial, "Historical data")
    plot_marginals(data_toxic, paste0("Historical data + 5/5 DLTs at A=", dose_A, ", B=", dose_B))
    plot_marginals(data_non_toxic, paste0("Historical data + 0/5 DLTs at A=", dose_A, ", B=", dose_B))
  }
}
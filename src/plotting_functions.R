

make_levelplots_ggplot <- function (data, title, mu_sd_inter, path.export)
{
  title_text <- title
  # * `prob_overdose` * `ewoc_ok`
  data$A <- as.factor(data$A)
  data$B <- as.factor(data$B)
  data$ewoc_fail <- !data$ewoc_ok
  data$type <- factor(data$type, levels=sort(unique(data$type))[c(1,4,2,3)])
  
  vars <- c("ewoc_fail", "prob_overdose", "prob_target", "prob_underdose")
  var_labels <- c("EWOC failure", "P(overdose)", "P(target)", "P(underdose)")
  data <- select(data, all_of(c("A", "B", "type", vars)))
  
  data_long <- data %>% 
    pivot_longer(cols=all_of(vars), names_to="variable", values_to="Probability")
  
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


plot_marginals <- function(data, title, mu_sd_inter, path.export)
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
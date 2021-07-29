library(tidyverse)
library(lme4)

dat1 <- read.table("~/Dropbox/oil_envr/methods/germination/markdown/ps_table1-paste.txt",T,'\t')

# Create germination percent
dat <- dat1 %>%
  slice(-c(10, 11)) %>%
  select(-code, -notes, -TP_11, -TP_12, -TP_13, -TP_14, 
         -TP_15, -TP_16, -TP_17, -TP_18, -TP_19, -TP_20) %>%
  pivot_longer(
    cols = paste("TP", 01:10, sep="_"),
    names_to = "time_point",
    values_to = "germ_count"
  ) %>%
  mutate(germ_perc = germ_count / N_seeds,
         sand1 = ifelse(sand == 0, 'no', 'yes'),
         parafilm1 = ifelse(parafilm == 0, 'no', 'yes'))

# Fit the model
mod0 <- lm(germ_perc ~ time_point + sand + parafilm + N_seed_weight_g,
             data = dat )

mod1 <- lm(germ_perc ~ sand + parafilm + N_seed_weight_g,
             data = dat )

# Fit with brms
library(brms)
mod00 <- brm(germ_perc ~ time_point + sand1 + parafilm1 + N_seed_weight_g,
             data = dat,
             family = "gaussian",
             seed = 2398, chains = 2L, iter = 3000) %>%
  add_criterion("loo", reloo = TRUE)

mod11 <- brm(germ_perc ~ sand1 + parafilm1 + N_seed_weight_g,
             data = dat,
             family = "gaussian",
             seed = 2398, chains = 2L, iter = 3000) %>%
  add_criterion("loo", reloo = TRUE)

             
# View the effects of sand and parafilm on germination percent
conditional_effects(mod00, effects = "parafilm1:sand1")

# View the effects of parafilm and time_point on germination
conditional_effects(mod00, effects = "time_point:parafilm1")

# Make the plots a bit nicer
# Treat disease as continuous
ce_mod00_ps <- plot(conditional_effects(mod00, effects = "parafilm1:sand1"),
     plot = FALSE,
     cat_args = list(size=3))[[1]] + # car_args adjusts the size of the points
  labs(x="Parafilm Added", y = "Scaled percent germination") + 
  scale_fill_discrete(guide = FALSE) +
  guides(color=guide_legend(title="Sand Added"))

# Plot for the time_points
# First make an object to order the timepoints correctly
tp_order = c("TP_1", "TP_2", "TP_3", "TP_4", "TP_5", "TP_6", "TP_7", "TP_8", "TP_9", "TP_10")
tp_labels = seq(from = 12, to = 120, by = 12)

# Plot it
# Just the time points
ce_mod00_t <- plot(conditional_effects(mod00, effects = "time_point"),
     plot = FALSE,
     cat_args = list(size=3))[[1]] + # car_args adjusts the size of the points
  labs(x="", y = "Percent germination") + 
  scale_x_discrete(limits = tp_order, labels = tp_labels)

# Time points with parafilm  
ce_mod00_tp <- plot(conditional_effects(mod00, effects = "time_point:parafilm1"),
     plot = FALSE,
     cat_args = list(size=3))[[1]] + # car_args adjusts the size of the points
  labs(x="Hours", y = "Scaled percent germination") + 
  scale_fill_discrete(guide = FALSE) +
  guides(color = guide_legend(title="Parafilm added")) +
  scale_x_discrete(limits = tp_order, labels = tp_labels) +
  theme(legend.position = c(0.8, 0.25))
  
  
# view conditional effects plots together
#library(cowplot)
pilot_26C_plots <- plot_grid(
  plot_grid(ce_mod00_t, ce_mod00_tp, ncol = 1, labels = c("A", "B")),
  ce_mod00_ps,
  ncol = 2, labels = c("", "C"))

# Save the plot as a pdf
ggsave("pilot_study_results_26C.pdf",
       plot = pilot_26C_plots,
       width = 8, 
       height = 5.5,
       units = "in")

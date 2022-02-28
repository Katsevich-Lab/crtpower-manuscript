### load libraries

library(glmnet)     # ridge regression
library(qqplotr)    # QQ plots
library(gtable)     # plot layout
library(grid)       # plot layout
library(tidyverse)  # data manipulation and plotting

### load simulation functions

source("simulation_functions.R")

### specify fixed parameters

n_train = 100     # number of samples for training g-hat
pi_init = 0.1     # probability Markov chain will start at 1
pi_flip = 0.1     # probability Markov chain will flip
sd_eps  = 1       # noise level
reps    = 1000    # number of random samples

### specify variable parameter values

n_test_vals = c(10, 25, 100)  # number of samples for CRT
SNR_vals = c(0,1,5)           # signal-to-noise ratio
p_vals = c(20,100,500)        # dimension of Z

### specify default variable parameter values

n_test_default = 100
SNR_default = 1
p_default = 500

### assemble parameter settings

variable_parameters = bind_rows(
  tibble(n_test = n_test_vals, SNR = SNR_default, p = p_default, 
         variable_setting = paste0("n = ", n_test_vals),
         fixed_setting = sprintf("SNR = %d, p = %d", SNR_default, p_default)),
  tibble(n_test = n_test_default, SNR = SNR_vals, p = p_default,
         variable_setting = paste0("SNR = ", SNR_vals),
         fixed_setting = sprintf("n = %d, p = %d", n_test_default, p_default)),
  tibble(n_test = n_test_default, SNR = SNR_default, p = p_vals,
         variable_setting = paste0("p = ", p_vals),
         fixed_setting = sprintf("n = %d, SNR = %d", n_test_default, SNR_default))) %>%
  mutate(setting = row_number(), 
         variable_setting = factor(variable_setting),
         fixed_setting = factor(fixed_setting))

### run the simulation

num_settings = nrow(variable_parameters)
unconditional_zvalues = matrix(NA, num_settings, reps, 
                               dimnames = list(NULL, paste0("rep_", 1:reps)))
conditional_zvalues = matrix(NA, num_settings, reps, 
                             dimnames = list(NULL, paste0("rep_", 1:reps)))

for(setting in 1:num_settings){
  cat(sprintf("Working on setting %d out of %d...\n", setting, num_settings))
  
  # extract variable parameters
  n_test = variable_parameters %>% filter(setting == !!setting) %>% pull(n_test)
  SNR = variable_parameters %>% filter(setting == !!setting) %>% pull(SNR)
  p = variable_parameters %>% filter(setting == !!setting) %>% pull(p)
  
  # training
  set.seed(1) # set seed for reproducibility
  Z_train = sample_Z(n_train, p, pi_init, pi_flip)
  Y_train = sample_Y_given_Z(Z_train, SNR, sd_eps)
  g_hat = cv.glmnet(x = Z_train, y = Y_train, alpha = 0)
  
  # CRT unconditional
  set.seed(1) # set seed for reproducibility
  for(rep in 1:reps){
    X_test = sample_X(n_test, pi_init)
    Z_test = sample_Z_given_X(X_test, p, pi_flip)
    Y_test = sample_Y_given_Z(Z_test, SNR, sd_eps)
    unconditional_zvalues[setting,rep] = 
      compute_U(X_test, Y_test, Z_test, g_hat, pi_init, pi_flip)
  }
  
  # CRT conditional on Y,Z
  set.seed(1) # set seed for reproducibility
  Z_test = sample_Z(n_test, p, pi_init, pi_flip)
  Y_test = sample_Y_given_Z(Z_test, SNR, sd_eps)
  for(rep in 1:reps){
    X_test = sample_X_given_Z(Z_test, pi_init, pi_flip)
    conditional_zvalues[setting,rep] = 
      compute_U(X_test, Y_test, Z_test, g_hat, pi_init, pi_flip)
  }
}

unconditional_results = variable_parameters %>%
  bind_cols(unconditional_zvalues %>% as_tibble()) %>%
  pivot_longer(cols = starts_with("rep_"), 
               names_to = "rep",
               names_prefix = "rep_",
               names_transform = list(rep = as.integer),
               values_to = "unconditional_zvalue")

# plot unconditional results
p1 = unconditional_results %>%
  filter(setting %in% c(1,2,3)) %>%
  mutate(variable_setting = 
           factor(variable_setting, levels = paste0("n = ", n_test_vals))) %>%
  ggplot(aes(sample = unconditional_zvalue)) + 
  geom_qq_band() + 
  stat_qq_line() + 
  stat_qq() +
  scale_x_continuous(expand = c(0.01,0.01)) +
  facet_grid(fixed_setting ~ variable_setting) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())

p2 = unconditional_results %>%
  filter(setting %in% c(4,5,6)) %>%
  mutate(variable_setting = 
           factor(variable_setting, levels = paste0("SNR = ", SNR_vals))) %>%
  ggplot(aes(sample = unconditional_zvalue)) + 
  geom_qq_band() + 
  stat_qq_line() + 
  stat_qq() +
  scale_x_continuous(expand = c(0.01,0.01)) +
  facet_grid(fixed_setting ~ variable_setting) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  ylab("Observed Quantile")

p3 = unconditional_results %>%
  filter(setting %in% c(7,8,9)) %>%
  mutate(variable_setting = 
           factor(variable_setting, levels = paste0("p = ", p_vals))) %>%
  ggplot(aes(sample = unconditional_zvalue)) + 
  geom_qq_band() + 
  stat_qq_line() + 
  stat_qq() +
  scale_x_continuous(expand = c(0.01,0.01)) +
  facet_grid(fixed_setting ~ variable_setting) +
  theme_bw() +
  theme(axis.title.y = element_blank()) + 
  xlab("Expected Quantile")

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g3 <- ggplotGrob(p3)
g <- rbind(g1, g2, g3, size = "first")
g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths)
grid.newpage()
grid.draw(g)
ggsave(plot = g,
       path = "../../writeups/simulation",
       filename = "unconditional.png", 
       device = "png",
       width = 6, 
       height = 6)



conditional_results = variable_parameters %>%
  bind_cols(conditional_zvalues %>% as_tibble()) %>%
  pivot_longer(cols = starts_with("rep_"), 
               names_to = "rep",
               names_prefix = "rep_",
               names_transform = list(rep = as.integer),
               values_to = "conditional_zvalue")

# plot conditional results
p1 = conditional_results %>%
  filter(setting %in% c(1,2,3)) %>%
  mutate(variable_setting = 
           factor(variable_setting, levels = paste0("n = ", n_test_vals))) %>%
  ggplot(aes(sample = conditional_zvalue)) + 
  geom_qq_band() + 
  stat_qq_line() + 
  stat_qq() +
  scale_x_continuous(expand = c(0.01,0.01)) +
  facet_grid(fixed_setting ~ variable_setting) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())

p2 = conditional_results %>%
  filter(setting %in% c(4,5,6)) %>%
  mutate(variable_setting = 
           factor(variable_setting, levels = paste0("SNR = ", SNR_vals))) %>%
  ggplot(aes(sample = conditional_zvalue)) + 
  geom_qq_band() + 
  stat_qq_line() + 
  stat_qq() +
  scale_x_continuous(expand = c(0.01,0.01)) +
  facet_grid(fixed_setting ~ variable_setting) +
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) + 
  ylab("Observed Quantile")

p3 = conditional_results %>%
  filter(setting %in% c(7,8,9)) %>%
  mutate(variable_setting = 
           factor(variable_setting, levels = paste0("p = ", p_vals))) %>%
  ggplot(aes(sample = conditional_zvalue)) + 
  geom_qq_band() + 
  stat_qq_line() + 
  stat_qq() +
  scale_x_continuous(expand = c(0.01,0.01)) +
  facet_grid(fixed_setting ~ variable_setting) +
  theme_bw() +
  theme(axis.title.y = element_blank()) + 
  xlab("Expected Quantile")

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g3 <- ggplotGrob(p3)
g <- rbind(g1, g2, g3, size = "first")
g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths)
grid.newpage()
grid.draw(g)
ggsave(plot = g,
       path = "../../writeups/simulation",
       filename = "conditional.png", 
       device = "png",
       width = 6, 
       height = 6)
library(tidyverse)
library(opdesmixr)
library(here)


beta_values_flies_01 = tibble(
  var = c("R", "O", "Y", "G", "B", "P", "R*O", "R*Y", "R*G", "R*B", "R*P", "R*UV", "O*Y", "O*G", "O*B", "O*P", "O*UV", "Y*G", "Y*B", "Y*P", "Y*UV", "G*B", "G*P", "G*UV", "B*P", "B*UV", "P*UV", "R*intensity", "O*intensity", "Y*intensity", "G*intensity", "B*intensity", "P*intensity", "UV*intensity", "intensity^2", "no_choice"),
  mean = c(-0.565, -0.794, -1.279, -1.719, -1.567, 0.194, -0.803, -2.255, -1.168, -0.845, -1.220, 3.177, 1.357, -0.655, -1.379, -3.561, 3.626, -0.030, 1.746, 1.690, 2.852, 1.283, -2.660, 1.169, -3.259, 2.734, -2.737, -0.404, -0.202, 0.091, 0.165, -0.112, 0.253, 0.171, 0.322, 2.242),
  sd = c(0.1743, 0.1744, 0.1465, 0.1689, 0.1757, 0.1445, 0.8643, 1.0034, 1.0645, 0.9793, 0.6961, 0.5002, 0.8808, 1.1355, 0.8769, 0.7177, 0.4771, 1.1870, 0.9171, 1.0182, 0.5404, 1.2326, 0.9853, 0.5469, 1.1914, 0.5641, 0.4595, 0.0986, 0.1188, 0.1060, 0.1249, 0.1144, 0.1058, 0.0754, 0.0240, 0.0447)
)

length_beta_flies_01 = nrow(beta_values_flies_01)

q_flies_01 = 7
J_flies_01 = 2
S_flies_01 = 60
order_flies_01 = 2
n_pv_flies_01 = 1
no_choice_flies_01 = T

set.seed(2023)

n_draws_flies_01 = 128

beta_flies_01_draws = get_halton_draws(
  beta = beta_values_flies_01$mean, 
  sd = beta_values_flies_01$sd, 
  n_draws = n_draws_flies_01)

# beta_flies_01_draws = get_correlated_halton_draws(
#   beta = beta_values_flies_01$mean, 
#   sigma = diag(beta_values_flies_01$sd^2), 
#   n_draws = 120)
# 
# beta_flies_01_draws = matrix(0.0, nrow = n_draws_flies_01, ncol = length_beta_flies_01)
# for(i in 1:length_beta_flies_01){
#   beta_flies_01_draws[, i] = as.numeric(
#     get_halton_draws(
#       beta = beta_values_flies_01$mean[i], 
#       sd = beta_values_flies_01$sd[i],
#       ndraws = n_draws_flies_01)
#   )
# }


# To see how different the means and SDs from the draws are from the theoretical values
data.frame(
  variable = beta_values_flies_01$var,
  real_mean = beta_values_flies_01$mean,
  mean_draws = apply(beta_flies_01_draws, 2, mean),
  real_sd = beta_values_flies_01$sd,
  sd_draws = apply(beta_flies_01_draws, 2, sd)
  ) %>% 
  mutate(
    mean_ratio = mean_draws/real_mean,
    sd_ratio = sd_draws/real_sd
    )


# Plot the real betas versus the draws
data.frame(
  variable = beta_values_flies_01$var,
  mean = beta_values_flies_01$mean,
  sd = beta_values_flies_01$sd,
  type = "real"
) %>% 
  bind_rows(
    data.frame(
      variable = beta_values_flies_01$var,
      mean = apply(beta_flies_01_draws, 2, mean),
      sd = apply(beta_flies_01_draws, 2, sd),
      type = "draws"
    )
  ) %>% 
  ggplot(aes(x = variable, color = type)) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd), size = 0.6, position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd), size = 1.2, position = position_dodge(width = 0.5)) +
  coord_flip() +
  theme_bw() +
  xlab("Beta") +
  ylab("Value")




i_opt_design_filename_flies_01 = paste0(here("out/second_i_optimal_design.rds"))

# Around 22 minutes per iteration
# Hence, 15 iterations would take around 5.5 hours
if(file.exists(i_opt_design_filename_flies_01)){
  cat("I_B optimal design already exists.\n")
  i_opt_design_flies_01 = readRDS(i_opt_design_filename_flies_01)
} else{
  cat("Doing I_B optimal design for insect experiment.\n")
  (t1I = Sys.time())
  i_opt_design_flies_01 =  mnl_mixture_coord_exch(
    n_random_starts = 4,
    q = q_flies_01, 
    J = J_flies_01, 
    S = S_flies_01, 
    n_pv = n_pv_flies_01,
    order = order_flies_01, 
    beta = beta_flies_01_draws,
    transform_beta = F,
    opt_method = "B",
    opt_crit = "I",
    max_it = 15,
    verbose = 1,
    plot_designs = F,
    seed = 2023,
    n_cores = 4,
    save_all_designs = T,
    no_choice = T
  )
  (t2I = Sys.time())
  t2I - t1I
  
  saveRDS(i_opt_design_flies_01, i_opt_design_filename_flies_01)
}



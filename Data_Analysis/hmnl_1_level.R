library(tidyverse)
library(rstan)

init_fun <- function() {
  out_list = list(
    # tau_level_0 = rep(1, n_mixture_cols+1), 
    # Omega_level_0 = diag(1, ncol = n_mixture_cols+1, nrow = n_mixture_cols+1)
    # tau_level_0 = rnorm(n_mixture_cols+1, mean = 1, sd = 0.1), 
    tau_level_0 = rgamma(n_mixture_cols+1, shape = 4, scale = 0.25),
    Omega_level_0 = rlkjcorr(K = n_mixture_cols+1, eta = 30),
    L_Omega_level_0 = chol(rlkjcorr(K = n_mixture_cols+1, eta = 30))
  )
  
  return(out_list)
}

# To stop Stan files from crashing RStudio
rstan_options(javascript=FALSE)

source("/Users/apple/Desktop/硕士论文的code/6.9/stan files flies/hmnl_reparameterized/japanese_flies_hmnl_1_level_utils.R")

counts_japanese_sample = read_csv("/Users/apple/Desktop/硕士论文的code/6.9/maxcount.csv")
counts_japanese_sample$image=1
n_experiments = 7
# Stan data ------------
unique_experiment_choice_set = counts %>%
  select(experiment, folder) %>%
  distinct()


index_dataframe = lapply(1:nrow(unique_experiment_choice_set), function(i){
  exp_i = unique_experiment_choice_set$experiment[i]
  cs_i =  unique_experiment_choice_set$folder[i]
  
  indices_i = which(counts_japanese_sample$folder == cs_i & counts_japanese_sample$experiment == exp_i)
  out = tibble(
    choice_set = cs_i,
    start = indices_i[1],
    end = indices_i[length(indices_i)],
    exp = exp_i
  )
}) %>%
  bind_rows()


unique_experiment_choice_set_image = counts_japanese_sample %>%
  select(experiment, folder, image) %>%
  distinct()

index_dataframe_image = lapply(1:nrow(unique_experiment_choice_set_image), function(i){
  exp_i = unique_experiment_choice_set_image$experiment[i]
  cs_i =  unique_experiment_choice_set_image$folder[i]
  image_i =  unique_experiment_choice_set_image$image[i]
  
  indices_i = which(counts_japanese_sample$folder == cs_i & counts_japanese_sample$experiment == exp_i & counts_japanese_sample$image == image_i)
  out = tibble(
    exp = exp_i,
    choice_set = cs_i,
    image = image_i,
    start = indices_i[1],
    end = indices_i[length(indices_i)]
  )
}) %>%
  bind_rows() %>% 
  mutate(choice_set_2 = as.integer(fct(paste(exp, choice_set))))



df_to_fit <- counts_japanese_sample %>% 
  select(all_of(c('R','O','Y','G','B','P','UV','intensity', 'no_choice'))) %>%
  as.data.frame()


X_stan_list = create_model_matrix_second_order_scheffe(df_to_fit)




n_mixture_cols = dim(X_stan_list$X)[2] - 1



stan_data <- list(
  n_alt = 3, 
  n_exp = n_experiments,
  n_var = ncol(X_stan_list$X), 
  n_obs = nrow(X_stan_list$X),
  Ny = counts_japanese_sample$max_count,
  X = X_stan_list$X,
  start = index_dataframe_image$start, # the starting observation for each image
  end = index_dataframe_image$end,  # the ending observation for each image
  n_images = nrow(index_dataframe_image),
  exp_index = index_dataframe_image$exp,
  
  # prior
  tau_level_0_mean = 1,
  tau_level_0_sd = 0.5,
  L_Omega_level_0_param = 5,
   beta_level_0_mean = rep(0, ncol(X_stan_list$X)),
   beta_level_0_sd = rep(2, ncol(X_stan_list$X))
  #beta_level_0_mean = prior_tibble$mean,
  
  # 3 times the estimated SD as a prior
  #beta_level_0_sd = 3*prior_tibble$sd
)

model_txt="data {
  int<lower=2> n_alt; // number of alternatives/outcomes
  int<lower=1> n_var; // of covariates
  int<lower=1> n_obs; // of observations
  int<lower=1> n_exp; // Number of experiments
  vector[n_obs] Ny; // counts of the response variable
  matrix[n_obs, n_var] X; // attribute matrix
  int n_images; // Number of images
  int start[n_images]; // the starting observation for each image
  int end[n_images]; // the ending observation for each image
  int exp_index[n_images];
  
  // priors
  real tau_level_0_mean;
  real<lower=0> tau_level_0_sd;
  real L_Omega_level_0_param;
  vector[n_var] beta_level_0_mean;
  vector[n_var] beta_level_0_sd;
}

parameters {
  vector[n_var] beta_level_0;
  matrix[n_exp, n_var] z;
  // corr_matrix[n_var] Omega_level_0;
  cholesky_factor_corr[n_var] L_Omega_level_0;
  vector<lower=0>[n_var] tau_level_0;
}

transformed parameters {
  // reparametrization of beta
  matrix[n_exp, n_var] beta_level_1 = rep_matrix(beta_level_0', n_exp) + z*diag_pre_multiply(tau_level_0, L_Omega_level_0);
}

model {
  vector[n_alt] utilities;

  to_vector(z) ~ normal(0, 1);
  // In LKJ:
  // if eta = 1, then the density is uniform over correlation matrices of a given order K (the number of row/column).
  // if eta > 1, the correlation values in correlation matrices are going to centered around 0. higher eta indicate no correlations (converge to identity correlation matrix).
  // https://yingqijing.medium.com/lkj-correlation-distribution-in-stan-29927b69e9be
  
  tau_level_0 ~ normal(tau_level_0_mean, tau_level_0_sd); 
  L_Omega_level_0 ~ lkj_corr_cholesky(L_Omega_level_0_param);
  
  beta_level_0 ~ normal(beta_level_0_mean, beta_level_0_sd);
  

  for(i in 1:n_images){
    utilities = X[start[i]:end[i], 1:n_var]*beta_level_1[exp_index[i]]';
    target += sum(log_softmax(utilities) .* Ny[start[i]:end[i]]);
  }
}


generated quantities{
  matrix[n_var, n_var] Sigma_level_0 = diag_pre_multiply(tau_level_0, L_Omega_level_0) * diag_pre_multiply(tau_level_0, L_Omega_level_0)';
  
 
   vector[n_images] log_lik;
 
  for(i in 1:n_images){
      log_lik[i] = sum(log_softmax(X[start[i]:end[i], 1:n_var]*beta_level_1[exp_index[i]]') .* Ny[start[i]:end[i]]);
    }
  
}"


###
model_hhml <- stan(
  model_code = model_txt, 
  data = stan_data, 
  warmup = 2000,
  iter = 5000, 
  chains = 4, 
  seed = 2023,
  cores = 4)
model_hhml

saveRDS(model_hhml,'/Users/apple/Desktop/result_hierarchical_multinomial_logit_model_with_max_count.rds')

model_hmnl= readRDS('/Users/apple/Desktop/result_hierarchical_multinomial_logit_model_with_max_count.rds')




###

betas_level_0_summary_01 = summary(model_hmnl, pars = c("beta_level_0"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "sd", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = colnames(X_stan_list$X),
         ix = 1:n()) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = F)) 




betas_level_0_summary_01 %>% 
  mutate(variable = fct_reorder(variable, ix, .desc = F)) %>% 
  ggplot(aes(x = variable)) +
  geom_hline(yintercept = 0) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle(paste0("Overall-level parameters, experiments 1 to ", n_experiments))




names_betas_level_1 = X_stan_list$names_betas_level_1


betas_level_1_summary_01 = summary(model_hhml, pars = c("beta_level_1"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean","sd", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(
    exp = paste0(
      "Exp ",
      rep(1:n_experiments, each = length(names_betas_level_1))
    ),
    variable = rep(names_betas_level_1, n_experiments)
  ) %>% 
  group_by(exp) %>% 
  mutate(ix = 1:n()) %>%
  ungroup() %>% 
  # mutate(variable = paste0(exp, ", ", var)) %>% 
  mutate(variable = fct_reorder(variable, ix, .desc = F)) 



betas_level_1_summary_01 %>% 
  ggplot(aes(x = variable)) +
  geom_hline(yintercept = 0) +
  geom_point(aes(y = mean), color = "black") +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), color = "dark grey") +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), color = "black", size = 1) +
  facet_wrap(~exp) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  ggtitle(paste0("Batch-level parameters of experiments 1 to ", n_experiments)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



betas_level_1_summary_01 %>% 
  ggplot(aes(x = variable, color = exp)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.6), size = 0.9) +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), size = 0.8, position = position_dodge(width = 0.6)) +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), size = 1.2, position = position_dodge(width = 0.6)) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  geom_vline(xintercept = 1:36 + 0.5, size = 0.3, color = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ggtitle(paste0("Batch-level parameters of experiments 1 to ", n_experiments))


betas_level_1_summary_01 %>% 
  ggplot(aes(x = variable, color = exp)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.6), size = 0.9) +
  geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd), size = 0.8, position = position_dodge(width = 0.6)) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd), size = 1.2, position = position_dodge(width = 0.6)) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  geom_vline(xintercept = 1:36 + 0.5, size = 0.3, color = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ggtitle(paste0("Batch-level parameters of experiments 1 to ", n_experiments))



betas_level_0_summary_01 %>% 
  mutate(exp = "All experiments", line_thickness_1 = 2, line_thickness_2 = 1) %>% 
  bind_rows(betas_level_1_summary_01 %>% 
              mutate(line_thickness_1 = 0.6, line_thickness_2 = 0.4)) %>% 
  ggplot(aes(x = variable, group = exp)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.7), size = 0.9) +
  geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd, size = line_thickness_2), 
                 position = position_dodge(width = 0.7)) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd, size = line_thickness_1),
                 position = position_dodge(width = 0.7)) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  geom_vline(xintercept = 1:36 + 0.5, size = 0.3, color = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ggtitle(paste0("Batch-level parameters of experiments 1 to ", n_experiments)) +
  scale_size_continuous(range = c(0.5, 2.5))



betas_level_0_summary_01 %>% 
  mutate(exp = "All experiments", line_thickness_1 = 2, line_thickness_2 = 1) %>% 
  bind_rows(betas_level_1_summary_01 %>% 
              mutate(line_thickness_1 = 0.8, line_thickness_2 = 0.5)) %>% 
  ggplot(aes(x = variable, color = exp)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.7), size = 0.9) +
  geom_linerange(aes(ymin = mean - 2*sd, ymax = mean + 2*sd, size = line_thickness_2), 
                 position = position_dodge(width = 0.7), show.legend = FALSE) +
  geom_linerange(aes(ymin = mean - sd, ymax = mean + sd, size = line_thickness_1),
                 position = position_dodge(width = 0.7), show.legend = FALSE) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  geom_vline(xintercept = 1:36 + 0.5, size = 0.3, color = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ggtitle(paste0("Batch-level parameters of experiments 1 to ", n_experiments)) +
  scale_size_continuous(range = c(0.5, 2)) +
  scale_color_manual(values = c("red", "#381532", "#4b1b42", "#5d2252", "#702963",
                                "#833074", "#953784", "#a83e95"))



Sigma_level_0_posterior_median_01 = matrix(
  as.data.frame(summary(model_hhml, pars = c("Sigma_level_0"), probs = c(0.5))$summary)$`50%`, 
  ncol = stan_data$n_var)



diag(Sigma_level_0_posterior_median_01)

























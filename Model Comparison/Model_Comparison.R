library(tidyverse)
library(rstan)

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
  exp_index = index_dataframe_image$exp
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
}

parameters {
    vector[n_var] beta;  // attribute effects 
}


model {
  vector[n_alt] utilities;


  for(i in 1:n_images){
    utilities = X[start[i]:end[i], 1:n_var]*beta;
    target += sum(log_softmax(utilities) .* Ny[start[i]:end[i]]);
  }
}


generated quantities{
  
  
 
   vector[n_images] log_lik;
 
  for(i in 1:n_images){
      log_lik[i] = sum(log_softmax(X[start[i]:end[i], 1:n_var]*beta) .* Ny[start[i]:end[i]]);
    }
  
}"


###
model_mnl <- stan(
  model_code = model_txt, 
  data = stan_data, 
  warmup = 2000,
  iter = 5000, 
  chains = 4, 
  seed = 2023,
  cores = 4)
model_mnl

saveRDS(model_mnl,'/Users/apple/Desktop/result_simple_multinomial_logit_model_with_max_count.rds')

##
betas_model_mnl_summary= summary(model_mnl, pars = c("beta"), probs = c(0.025, 0.1, 0.5, 0.9, 0.975))$summary %>% 
  as.data.frame() %>% 
  select("mean", "sd", "2.5%", "10%", "90%", "97.5%") %>% 
  mutate(variable = colnames(X_stan_list$X),
         ix = 1:n()) %>%
  mutate(variable = fct_reorder(variable, ix, .desc = F)) 
betas_model_mnl_summary$model='MNLM'
betas_level_0_summary_01$model='MLM'
beta_compare=rbind(betas_model_mnl_summary,betas_level_0_summary_01)



beta_compare %>% 
  ggplot(aes(x = variable, color = model)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 0, size = 0.3) +
  geom_point(aes(y = mean), position = position_dodge(width = 0.6), size = 2) +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), size = 0.8, position = position_dodge(width = 0.6)) +
  geom_linerange(aes(ymin = `10%`, ymax = `90%`), size = 1.5, position = position_dodge(width = 0.6)) +
  theme_bw() +
  xlab("Parameter") +
  ylab("Value") +
  geom_vline(xintercept = 1:36 + 0.5, size = 0.3, color = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  ggtitle("MLM VS MNLM ")

##loo

model_hmnl=readRDS('/Users/apple/Desktop/result_hierarchical_multinomial_logit_model_with_max_count.rds')
library(loo)
loo_fit_hmnl <- loo(model_hmnl)
print(loo_fit_hmnl)
loo_fit_mnl <- loo(model_mnl)
print(loo_fit_mnl)
compare_loo <- loo_compare(loo_fit_mnl, loo_fit_hmnl)
compare_loo

loo_fit_hmnl_420=loo_fit_hmnl$pointwise[,1]
loo_fit_mnl_420=loo_fit_mnl$pointwise[,1]
loo_420=as.data.frame(cbind(loo_fit_hmnl_420,loo_fit_mnl_420))
#t_loo <- t.test(loo_fit_hmnl_420, loo_fit_mnl_420, paired = TRUE,alternative = "greater")

x <- seq(1, nrow(loo_420))
plot <- ggplot(loo_420, aes(x = x))
plot <- plot + geom_point(aes(y = loo_fit_hmnl_420, color ='MLM' ))
plot <- plot + geom_point(aes(y = loo_fit_mnl_420, color = "MNLM"))
plot <- plot + labs(title = "Individual Contriburion for elpd_loo", x = "Index", y = "elpd_loo")
plot <- plot +geom_vline(xintercept = c(60, 120, 180, 240, 300, 360),
             linetype = "dashed", color = "black")

##waic
log_lik_hmnl= extract_log_lik(model_hmnl)
waic_fit_hmnl=waic(log_lik_hmnl)
print(waic_fit_hmnl)

log_lik_mnl= extract_log_lik(model_mnl)
waic_fit_mnl=waic(log_lik_mnl)
print(waic_fit_mnl)
compare_waic <- loo_compare(waic_fit_mnl, waic_fit_hmnl)

waic_fit_hmnl_420=waic_fit_hmnl$pointwise[,1]
waic_fit_mnl_420=waic_fit_mnl$pointwise[,1]
waic_420=as.data.frame(cbind(waic_fit_hmnl_420,waic_fit_mnl_420))
#t_waic <- t.test(waic_fit_hmnl_420, waic_fit_mnl_420, paired = TRUE,alternative = "greater")

x <- seq(1, nrow(waic_420))
plot <- ggplot(waic_420, aes(x = x))
plot <- plot + geom_point(aes(y = waic_fit_hmnl_420, color ='MLM' ))
plot <- plot + geom_point(aes(y = waic_fit_mnl_420, color = "MNLM"))
plot <- plot + labs(title = "Individual Contriburion for elpd_waic", x = "Index", y = "elpd_waic")
plot <- plot +geom_vline(xintercept = c(60, 120, 180, 240, 300, 360),
                         linetype = "dashed", color = "black")

diff_loo=loo_fit_hmnl_420-loo_fit_mnl_420
diff_waic=waic_fit_hmnl_420-waic_fit_mnl_420
diff=as.data.frame(cbind(diff_loo,diff_waic))

x <- seq(1, nrow(diff))
plot <- ggplot(diff, aes(x = x))
plot <- plot + geom_point(aes(y = diff_loo))
plot <- plot + labs(title = "Difference of elpd_loo between MLM and MNLM", x = "Index", y = "diff_elpd_loo")
plot <- plot +geom_vline(xintercept = c(60, 120, 180, 240, 300, 360),
                         linetype = "dashed", color = "black")
plot <- plot +geom_hline(yintercept = 0,
                        color = "red")

x <- seq(1, nrow(diff))
plot <- ggplot(diff, aes(x = x))
plot <- plot + geom_point(aes(y = diff_waic))
plot <- plot + labs(title = "Difference of elpd_waic between MLM and MNLM", x = "Index", y = "diff_elpd_waic")
plot <- plot +geom_vline(xintercept = c(60, 120, 180, 240, 300, 360),
                         linetype = "dashed", color = "black")
plot <- plot +geom_hline(yintercept = 0,
                         color = "red")
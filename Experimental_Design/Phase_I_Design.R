library(opdesmixr)
library(tidyverse)
library(here)


designs_folder = here("out/insect_pv_designs/")
dir.create(designs_folder, showWarnings = F)

n_cores = parallel::detectCores()




q = 7
J = 2
S = 60 #total number of choice sets
n_pv = 1

beta_vec = c(rep(0,35))
names(beta_vec) = c(
  "R","O","Y","G","B","P",
  
  "R:O","R:B","R:Y","R:P","R:G","R:UV", 
  "O:B","O:Y","O:P","O:G","O:UV",
  "B:Y","B:P","B:G","B:UV",
  "Y:P","Y:G","Y:UV",
  "P:G","P:UV",
  "G:UV",
  
  "R:z","O:z","Y:z","G:z","B:z","P:z","UV:z",
  "I(z^2)"
)

SDs = c(rep(100,35))
names(SDs) = names(beta_vec)
var_cov_mat = diag(SDs^2)


n_draws = 128
beta_correlated_draws_insect = get_correlated_halton_draws(beta_vec, var_cov_mat, n_draws)


# n_random_initial_starts = 80
n_random_initial_starts = 48
max_it_insect = 15
# max_it_insect = 3
seed = 2022


insect_pv_d_opt_filename = paste0(designs_folder, "insect_pv_d_optimal_", max_it_insect, "iter.rds")
insect_pv_i_opt_filename = paste0(designs_folder, "insect_pv_i_optimal_", max_it_insect, "iter.rds")


if(file.exists(insect_pv_d_opt_filename)){
  cat("D_B optimal design already exists.\n")
} else{
  # 40 seconds in 4 cores with 4 random starts and 4 iterations and point estimate
  # 4 minutes in 4 cores with 4 random starts and 3 iterations and 12 draws
  # 42 minutes in 4 cores with 4 random starts and 3 iterations and 128 draws
  # 8643.5 seconds (2.5 hours) in 4 cores with 4 random starts and 10 iterations and 128 draws
  # 10174 seconds in 4 cores with 4 random starts and 15 iterations and 128 draws
  # 45,000 seconds in 12 cores with 48 random starts and 15 iterations and 128 draws
  cat("Doing D_B optimal design for insect experiment.\n")
  (t1D = Sys.time())
  insect_pv_D_opt = mnl_mixture_coord_exch(
    n_random_starts = n_random_initial_starts,
    q = q,
    J = J,
    S = S,
    n_pv = n_pv,
    order = 2,
    beta = beta_correlated_draws_insect,
    transform_beta = F,
    opt_method = "D",
    n_cox_points = 10,
    opt_crit = "D",
    max_it = max_it_insect,
    verbose = 1,
    plot_designs = F,
    seed = seed,
    n_cores = n_cores,
    save_all_designs = T
  )
  
  (t2D = Sys.time())
  t2D - t1D
  
  saveRDS(insect_pv_D_opt, insect_pv_d_opt_filename)
}



if(file.exists(insect_pv_i_opt_filename)){
  cat("I_B optimal design already exists.\n")
} else{
  # 9804.4 seconds in 4 cores with 4 random starts and 15 iterations and 128 draws
  # 45,000 seconds in 12 cores with 48 random starts and 15 iterations and 128 draws
  cat("Doing I_B optimal design for insect experiment.\n")
  (t1I = Sys.time())
  insect_pv_I_opt =  mnl_mixture_coord_exch(
    n_random_starts = n_random_initial_starts,
    q = q,
    J = J,
    S = S,
    n_pv = n_pv,
    order = 2,
    beta = beta_correlated_draws_insect,
    transform_beta = F,
    opt_method = "D",
    n_cox_points = 10,
    opt_crit = "I",
    max_it = max_it_insect,
    verbose = 1,
    plot_designs = F,
    seed = seed,
    n_cores = n_cores,
    save_all_designs = T
  )
  (t2I = Sys.time())
  t2I - t1I
  
  saveRDS(insect_pv_I_opt, insect_pv_i_opt_filename)
}







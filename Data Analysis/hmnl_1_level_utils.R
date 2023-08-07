rlkjcorr = function(K, eta = 1)  {
  stopifnot(is.numeric(K), K >= 2, K == as.integer(K))
  stopifnot(eta > 0)
  
  alpha <- eta + (K - 2)/2
  r12 <- 2 * rbeta(1, alpha, alpha) - 1
  R <- matrix(0, K, K)
  R[1, 1] <- 1
  R[1, 2] <- r12
  R[2, 2] <- sqrt(1 - r12^2)
  if (K > 2) 
    for (m in 2:(K - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m/2, alpha)
      z <- rnorm(m, 0, 1)
      z <- z/sqrt(crossprod(z)[1])
      R[1:m, m + 1] <- sqrt(y) * z
      R[m + 1, m + 1] <- sqrt(1 - y)
    }
  return(crossprod(R))
}





create_model_matrix_second_order_scheffe = function(df_to_convert){
  
  correct_var_names = c("R", "O", "Y", "G", "B", "P", "UV", "intensity", "no_choice")
  # correct_var_names = c("R", "O", "Y", "G", "B", "P", "UV", "intensity", "is_right", "is_daytime", "no_choice")
  
  if(is.character(all.equal(names(df_to_convert), correct_var_names))) {
    stop("Names of df_to_convert must be c(", paste(correct_var_names, collapse = ", "), ")")
  }
  
  
  
  mixture_variable_names_no_UV = c('R','O','Y','G','B','P')
  mixture_variable_names = c(mixture_variable_names_no_UV, 'UV')
  
  mixture_pairwise_interaction_names = c(
    'R*O','R*Y','R*G','R*B','R*P','R*UV',
    'O*Y','O*G','O*B','O*P','O*UV',
    'Y*G','Y*B','Y*P','Y*UV',
    'G*B','G*P','G*UV',
    'B*P','B*UV',
    'P*UV')
  
  mixture_intensity_interaction_names = c('R*intensity','O*intensity','Y*intensity','G*intensity',
                                          'B*intensity','P*intensity','UV*intensity')
  
  
  X_other_vars_nc <- df_to_convert %>% 
    # dplyr::select(all_of(c('no_choice','is_right', 'is_daytime'))) %>%
    dplyr::select(all_of(c('no_choice'))) %>%
    as.data.frame()
  
  X_main <- df_to_convert %>% 
    dplyr::select(all_of(c(mixture_variable_names,'intensity'))) %>%
    as.data.frame()
  
  
  # X_color_pairwise <- combn(mixture_variable_names, 2, function(x) X_main[x[1]]*X_main[x[2]], simplify = FALSE) %>%
  #   bind_cols() %>%
  #   set_names(mixture_pairwise_interaction_names)
  
  X_color_pairwise = combn(seq_along(mixture_variable_names), 2, function(x) X_main[, x[1]]*X_main[, x[2]], simplify = T) %>% 
    matrix(ncol = length(mixture_pairwise_interaction_names)) %>% 
    as.data.frame() %>% 
    set_names(mixture_pairwise_interaction_names)
  
  
  X_color_intensity <- X_main %>%
    dplyr::select(all_of(c(mixture_variable_names,'intensity'))) %>%
    mutate(across(c(mixture_variable_names,'intensity'), ~.*intensity)) %>%
    set_names(c(mixture_intensity_interaction_names,'intensity^2'))
  
  
  X <- X_main %>% 
    dplyr::select(-UV, -intensity) %>% 
    bind_cols(X_color_pairwise) %>% 
    bind_cols(X_color_intensity) %>%
    bind_cols(X_other_vars_nc) %>% 
    as.matrix()
  
  names_betas_level_1 = c(mixture_variable_names_no_UV, mixture_pairwise_interaction_names, mixture_intensity_interaction_names,'intensity^2', "no_choice")
  
  return(list(
    X = X,
    mixture_variable_names_no_UV = mixture_variable_names_no_UV,
    mixture_variable_names = mixture_variable_names,
    mixture_pairwise_interaction_names = mixture_pairwise_interaction_names,
    mixture_intensity_interaction_names = mixture_intensity_interaction_names,
    names_betas_level_1 = names_betas_level_1
  ))
}




create_lattice_design = function(n_var = 3, n_levels = 5, limits = NULL){
  # Examples: 
  # create_lattice_design(3, 10)
  #
  # limits = cbind(c(0.4, 0.55), c(0, 0.5), c(0.1, 0.45))
  # create_lattice_design(3, 10, limits)
  
  # create_lattice_design(3, 20) %>%
  #   ggtern::ggtern(ggtern::aes(x1, x2, x3)) +
  #   geom_point(shape = "x", size = 4) +
  #   theme_minimal() +
  #   ggtern::theme_nomask()
  
  if(!is.null(limits)){
    if(ncol(limits) != n_var) stop("limits should have ", n_var, " columns because n_var = ", n_var)
  } else{
    limits = matrix(replicate(n_var, c(0, 1)), ncol = n_var)
  }
  
  lattice_df = lapply(1:n_var, function(i){
    if(limits[1, i] != limits[2, i]){
      out = tibble(x = seq(limits[1, i], limits[2, i], length.out = n_levels)) %>% 
        set_names(paste0("x", i))
    } else{
      out = tibble(x = limits[1, i]) %>% 
        set_names(paste0("x", i))
    }
    return(out)
  }) %>% 
    bind_cols() %>% 
    expand.grid() %>% 
    as_tibble() %>% 
    distinct() %>% 
    slice(which(abs(apply(., 1, sum) - 1) < 1e-15))
  
  return(lattice_df)
}


# create_lattice_design = function(n_var = 3, n_levels = 5){
#   # Example: create_lattice_design(3, 10)
#   
#   # create_lattice_design(q, 20) %>% 
#   #   ggtern::ggtern(ggtern::aes(x1, x2, x3)) +
#   #   geom_point(shape = "x", size = 4) +
#   #   theme_minimal() +
#   #   ggtern::theme_nomask()
#   
#   lattice_df = lapply(1:n_var, function(i){
#     tibble(x = seq(0, 1, by = (1/n_levels))) %>% 
#       set_names(paste0("x", i))
#   }) %>% 
#     bind_cols() %>% 
#     expand.grid() %>% 
#     as_tibble() %>% 
#     slice(which(abs(apply(., 1, sum) - 1) < 1e-12))
#   
#   return(lattice_df)
# }
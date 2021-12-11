#'
#'
update_nu_gamma <- function(nu_int, gamma_int, beta_int, Lamb_int,
                            lambda2 = 0, lambda3 = 0, kappa = 1,
                            p_x, K, coeff_group_numbers, diff_tol = 1e-5, iter_max = 500){

  beta_int_mat <- matrix(beta_int, ncol = p_x, nrow = K, byrow = T)
  nu_int_mat <- matrix(nu_int, ncol = p_x, nrow = K, byrow = T)
  Lamb_int_mat <- matrix(Lamb_int, ncol = p_x, nrow = K, byrow = T)

  nu_est <- c()
  gamma_est <- list()

  for(j in 1:p_x){
    nu_int_p <- nu_int_mat[, j]
    beta_int_p <- beta_int_mat[, j]
    Lamb_int_p <- Lamb_int_mat[, j]
    gamma_int_p <- gamma_int[[j]]

    f_val <- pen_obj_by_p(nu_int_p, gamma_int_p, beta_int_p, Lamb_int_p, kappa, lambda2, lambda3)

    i <- diff <- 1
    while(i < iter_max & diff > diff_tol){
      # Step1: given gamma, update nu
      beta_tilde <- beta_int_p + Lamb_int_p/kappa

      nu_int_p_new <- sapply(1:length(beta_int_p), function(s){
          beta_i <- beta_tilde[s]
          centers <- gamma_int_p
          opt_center <- centers[which.min((beta_i - centers)^2)]
          (2 * opt_center * lambda2 + kappa * beta_i)/(2 * lambda2 + kappa)
        })

      # Step2: given nu, update gamma
        if(lambda3 > 0 & is.null(coeff_group_numbers[[j]])){
          # Penalized kmeans
          gamma_int_p_new <- DPmeans(nu_int_p_new, lambda3)
        } else {
          if(coeff_group_numbers[[j]] < K){
            # Regular kmeans
            gamma_int_p_new <- kmeans(nu_int_p_new, coeff_group_numbers[[j]])$centers %>%
              as.numeric() %>% sort()
          } else if (coeff_group_numbers[[j]] == K){
            gamma_int_p_new <- sort(nu_int_p_new )
          } else {
            stop("# of groups cannot be larger than # of studies.")
          }
        }

      new_f_val <- pen_obj_by_p(nu_int_p_new, gamma_int_p_new,
                                beta_int_p, Lamb_int_p,
                                kappa, lambda2, lambda3)
      diff <- f_val - new_f_val
      f_val <- new_f_val

      nu_int_p <- nu_int_p_new
      gamma_int_p <- gamma_int_p_new

      # cat("i = ", i, '\n')
      # cat("nu = ", nu_int_p, "\n")
      # cat("gamma = ", gamma_int_p, "\n")
      # cat("f value=", f_val, "\n")
      # cat("diff = ", round(diff, 2), "\n")
      i <- i+1
    }

    nu_est <- cbind(nu_est, nu_int_p)
    gamma_est[[j]] <- gamma_int_p
  }

  list(nu = as.numeric(t(nu_est)),
       gamma = gamma_est)
}

#' Objective function of nu and gamma
#'
pen_obj_by_p <- function(nu_int_p, gamma_int_p,
                         beta_int_p, Lamb_int_p,
                         kappa, lambda2, lambda3){

  min_l2 <- lapply(1:length(gamma_int_p), function(g){
    (nu_int_p - gamma_int_p[g])^2
  }) %>%
    do.call(cbind, .) %>%
    apply(1, min)

  l1_term <- lambda2 * sum(min_l2)
  l2_term <- sum((beta_int_p - nu_int_p + Lamb_int_p/kappa)^2) * kappa/2
  l3_term <- lambda2*lambda3*sum(length(unique(gamma_int_p)))
  l1_term + l2_term + l3_term
}

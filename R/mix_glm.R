#' Fit generalized linear models for integrated data
#'
#' @param y_var response variable
#' @param x_vars predictors with heterogeneous covariate effects
#' @param w_vars predictors with homogeneous covariate effects
#' @param dat_list a list of data sets by studies
#' @param family a description of the error distribution and link function to be used in the model
#' @param kappa parameter for argumented Lagrangian function
#' @param lambda1 penalty parameter for L1 norm
#' @param lambda2 penalty parameter for multi-directional L2 norm
#' @param lambda3 parameter for the penalized K-means algorithm, leading to lambda3/lambda2 penalty paramter on the subgroup sizes
#' @param coeff_group_numbers a vector specifying numbers of subgroups within each predictor. If NULL, penalized K-means is applied.
#' @param beta_int initial values of heterogeneous effects
#' @param alpha_int initial values of homogeneous effects
#' @param sigma_int initial values of scale parameters
#' @param sigma_fixed logic
#' @param penalty_method penalty method
#' @export
mix_glm <- function(y_var, x_vars, w_vars = NULL, dat_list, family = c("gaussian", "binomial"),
                    kappa = 1, lambda1 = 0, lambda2 = 0, lambda3 = 0, weighted = 0,
                    penalty.factor = 1,
                    coeff_group_numbers = NULL,
                    beta_int = NULL, alpha_int = NULL,
                    tol1 = 1e-5, tol2 = 1e-5, iter_max = 500, check = 0){

  family <- match.arg(family)
  if(family == "binomial"){ stop("To be implemented.")}

  p_x <- length(x_vars)
  p_w <- length(w_vars)
  K <- length(dat_list)
  sample_size <- sapply(dat_list, nrow)
  y_list <- lapply(dat_list, function(dat){ select(dat, one_of(y_var)) %>% as.matrix() })
  X_list <- lapply(dat_list, function(dat){ select(dat, one_of(x_vars)) %>% as.matrix() })
  W_list <- lapply(dat_list, function(dat){ select(dat, one_of(w_vars)) %>% as.matrix() })
  y_all <- do.call(rbind, y_list) %>% as.matrix()
  X_all <- do.call(bdiag, X_list) %>% as.matrix()
  W_all <- do.call(rbind, W_list) %>% as.matrix()

  if(is.null(beta_int)){ beta_int <- matrix(0, nrow = K, ncol = p_x) } # K*p_x
  if(is.null(alpha_int)){ alpha_int <- rep(0, p_w) }
  if(is.null(w_vars)){ alpha_int <- numeric() }
  if(lambda3 == 0){ # no penalty on # of groups
    if(is.null(coeff_group_numbers)){
      stop("The coeff_group_numbers is missing without default value.")
    }
  } #else {
    #if(is.null(coeff_group_numbers)){
    #  coeff_group_numbers <- rep(K, p_x)
    #}
  #}

  gamma_int <- lapply(1:p_x, function(i){
    if(!is.null(coeff_group_numbers[[i]])){
      if(coeff_group_numbers[[i]] < K){
        kmeans(beta_int[, i], centers = coeff_group_numbers[[i]])$centers %>%
          as.numeric() %>% sort()
      } else if(coeff_group_numbers[[i]] == K){
        as.numeric(beta_int[, i]) %>% sort()
      } else {
        stop("Number of desired groups cannot be larger than number of studies.")
      }
    } else { # if # of group is not specified
      as.numeric(beta_int[, i]) %>% sort()
    }
  })

  beta_int <- as.vector(t(beta_int)) # beta_int concatenated by rows
  Lamb_int <- rep(0, K*p_x)
  nu_int <- beta_int
  eta_int <- rep(0, K*p_x)

  if(weighted > 0){
    l1_weights <- (abs(c(beta_int, alpha_int))^{-weighted})
    l1_weights[l1_weights > 1e2] <- 1e2
  } else {
    l1_weights <- rep(1, length(c(beta_int, alpha_int)))
  }
  l1_weights <- l1_weights * penalty.factor
  # cat("l1_weights = ", round(l1_weights,2), "\n")

  f_val <- get_obj_val(beta_int, alpha_int, gamma_int,
                       y_all, X_all, W_all, K, p_x,
                       lambda1, lambda2, lambda3, l1_weights)
  # cat("Initial objective function value:", f_val, "\n")

  i <- 0
  fn_vals <- c()
  criter1 <- criter2 <- 1
  while((criter1 > tol1 | criter2 > tol2) & i < iter_max){

      theta_new <- try(update_theta(beta_int, alpha_int, y_all, X_all, W_all,
                                    Lamb_int, nu_int, kappa, lambda1, l1_weights))

    if(class(theta_new) == "try-error"){
      print(theta_new)
      stop("Updating theta failed...")
    }
    beta_new <- theta_new$beta
    alpha_new <- theta_new$alpha

    # MDSP method: multi-directional separation penalty
    res <- try(update_nu_gamma(nu_int, gamma_int, beta_new, Lamb_int,
                               lambda2, lambda3, kappa, p_x, K,
                               coeff_group_numbers), silent = T)
    if(class(res) == "try-error"){
      print(res)
      stop("Updating gamma and nu failed ...")
    }
    nu_new <- res$nu
    gamma_new <- res$gamma
    eta_new <- beta_new - nu_new
    Lamb_new <- Lamb_int + kappa * eta_new

    criter1 <- sqrt(sum((beta_new - beta_int)^2))/(K * p_x) +
      sqrt(sum((alpha_new - alpha_int)^2))/(p_w + as.numeric(p_w == 0))
    criter2 <- sqrt(sum((eta_new - eta_int)^2))

    # update parameters
    beta_int <- beta_new
    alpha_int <- alpha_new
    eta_int <- beta_int - nu_int
    Lamb_int <- Lamb_new
    nu_int <- nu_new
    gamma_int <- gamma_new

    # calculate objective function
    new_f_val <- get_obj_val(beta_int, alpha_int, gamma_int,
                             y_all, X_all, W_all, K, p_x,
                             lambda1, lambda2, lambda3, l1_weights)
    fn_vals <- c(fn_vals, new_f_val)
    f_change <- f_val - new_f_val
    f_val <- new_f_val

    i <- i + 1

    if(check == 1){
      cat("##### Iteration ", i, "\n")
      cat("alpha = ", round(alpha_int, 2), "\n")
      cat("beta = ", "\n")
      print(matrix(round(beta_int, 2), ncol = p_x, nrow = K, byrow = T))
      cat("gamma = ", round(unlist(gamma_int), 2), "\n")
      cat("criter1 = ", criter1, "\n")
      cat("criter2 = ", criter2, "\n")
      cat("fn change = ", f_change, "\n")
      cat("objective fn value = ", f_val, "\n")
    }
  }

  if(criter1 > tol1 | criter2 > tol2){
    return(NULL)
  } else {
    cat("Final objective fn value = ", f_val, "\n")
    list(beta = matrix(beta_int, ncol = p_x, nrow = K, byrow = TRUE),
         alpha = alpha_int,
         nu = matrix(nu_int, ncol = p_x, nrow = K, byrow = TRUE),
         gamma = lapply(gamma_int, function(x){unique(x)}),
         penalty_par = c(lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3),
         fn_vals = fn_vals)
  }
}


group_by_gamma <- function(md_res, x_names){
  md_est <- as.data.frame(md_res$beta)
  names(md_est) <- x_names
  selected <- which(colSums(md_est != 0) > 0)

  res <- data.frame(matrix(ncol=length(selected),
    nrow = nrow(md_est),
    dimnames=list(NULL, names(selected))))
  for(i in 1:length(selected)){
    res[,i] <- sapply(md_est[, selected[i]], function(x){
      which.min(abs(x-
        c(0, setdiff(md_res$gamma[[selected[i]]], 0))
      )) %>% as.factor()
    })
  }
  res
}

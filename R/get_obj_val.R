#' Calculate objective function (gaussian only)
#'
get_obj_val <- function(beta_int, alpha_int, gamma_int,
                        y_all, X_all, W_all, K, p_x,
                        lambda1, lambda2, lambda3, l1_weights){

  err_all <- y_all - X_all %*% beta_int - W_all %*% alpha_int

  beta_mat <- matrix(beta_int, nrow = K, ncol = p_x, byrow = T)
  md_min_l2_val <- lapply(1:p_x, function(t){
    lapply(gamma_int[[t]], function(tt){
      (beta_mat[, t] - tt)^2
    }) %>% do.call(cbind, .) %>% apply(1, min)
  }) %>% do.call(cbind, .)
  l1_penalty <- sum(c(abs(beta_int), abs(alpha_int)) * l1_weights) * lambda1
  l2_penalty <- lambda2 * sum(md_min_l2_val)
  l3_penalty <- lambda3 * sum(sapply(gamma_int, function(x){length(unique(x))}))

  mean(err_all^2)/2 + l1_penalty + l2_penalty + l3_penalty
}

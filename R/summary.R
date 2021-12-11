#' Predict response
#'
#' @export
predict_mixglm <- function(res, dat_list, y_var, x_vars, w_vars, family){
  K <- length(dat_list)
  all_pred <- c()

  for(i in 1:K){
    dat <- dat_list[[i]]
    y <- as.matrix(dat[y_var])
    X <- as.matrix(dat[x_vars])
    if(is.null(w_vars)){
      mu <- X %*% as.matrix(res$beta[i, ])
    } else {
      W <- as.matrix(dat[w_vars])
      mu <- X %*% as.matrix(res$beta[i, ]) + W %*% as.matrix(res$alpha)
    }
    pred_i <- data.frame(y = dat[y_var],
                         mu = mu,
                         study = ifelse(is.null(names(dat_list)[i]), i, names(dat_list)[i]))
    all_pred <- rbind(all_pred, pred_i)
  }

  all_pred
}


#' Summarize selection results
#'
#' @export
selection_result <- function(..., methods, p_x, p_w){
  all_res <- list(...)

  res <- list()
  p <- p_x + p_w
  for(j in 1:p){
    df <- data.frame(prob=numeric(), methods = character(), idx = integer())
    for(i in 1:length(all_res)){
      if(methods[i] != "MDSP"){
        prob <- lapply(all_res[[i]], function(x){abs(x[, j]) > 1e-8}) %>% do.call(rbind, .) %>%
          apply(2, mean)
      } else {
        if(p <= p_w){
          prob <- NA
        } else {
          prob <- lapply(all_res[[i]][[2]], function(x){
            x[j-p_w, ]
          }) %>% do.call(rbind, .) %>% apply(2, mean)
        }
      }
      df_i <- data.frame(prob = prob, methods = methods[i], idx = 1:length(prob))
      df <- rbind(df, df_i)
    }
    res[[j]] <- df
  }
  res
}


get_sqrt_mses <- function(..., methods, p, thetas){
  all_res <- list(...)
  stopifnot(length(all_res) == length(methods))

  res <- list()
  for(j in 1:p){
    df <- data.frame(value = numeric(), methods = character(), idx = integer())
    for(i in 1:length(all_res)){
      values_i <- lapply(all_res[[i]], function(x){
        (x[, j] - thetas[, j])^2 %>% as.vector()
      }) %>% do.call(rbind, .) %>% apply(2, mean) %>% sqrt()

      df_i <- data.frame(value = values_i, methods = methods[i],
                         idx = 1:length(values_i))
      df <- rbind(df, df_i)
    }
    res[[j]] <- df
  }
  res
}


get_abs_biase <- function(..., methods, p, thetas){
  all_res <- list(...)
  stopifnot(length(all_res) == length(methods))

  res <- list()
  for(j in 1:p){
    df <- data.frame(value = numeric(), methods = character(), idx = integer())
    for(i in 1:length(all_res)){
      values_i <- lapply(all_res[[i]], function(x){
        (x[, j] - thetas[, j]) %>% as.vector()
      }) %>% do.call(rbind, .) %>% apply(2, mean) %>% abs()

      df_i <- data.frame(value = values_i, methods = methods[i],
                         idx = 1:length(values_i))
      df <- rbind(df, df_i)
    }
    res[[j]] <- df
  }
  res
}

get_mses <- function(..., methods, p, thetas){
  all_res <- list(...)
  stopifnot(length(all_res) == length(methods))

  res <- list()
  for(j in 1:p){
    df <- data.frame(value = numeric(), methods = character(), idx = integer())
    for(i in 1:length(all_res)){
      values_i <- lapply(all_res[[i]], function(x){
        (x[, j] - thetas[, j])^2 %>% as.vector()
      }) %>% do.call(rbind, .) %>% apply(2, mean)

      df_i <- data.frame(value = values_i, methods = methods[i],
                         idx = 1:length(values_i))
      df <- rbind(df, df_i)
    }
    res[[j]] <- df
  }
  res
}

get_bias <- function(..., methods, p, thetas){
  all_res <- list(...)
  stopifnot(length(all_res) == length(methods))

  res <- list()
  for(j in 1:p){
    df <- data.frame(value = numeric(), methods = character(), idx = integer())
    for(i in 1:length(all_res)){
      values_i <- lapply(all_res[[i]], function(x){
        (x[, j] - thetas[, j]) %>% as.vector()
      }) %>% do.call(rbind, .) %>% apply(2, mean)

      df_i <- data.frame(value = values_i, methods = methods[i],
                         idx = 1:length(values_i))
      df <- rbind(df, df_i)
    }
    res[[j]] <- df
  }
  res
}

get_interval <- function(..., methods, p, thetas){
  all_res <- list(...)
  stopifnot(length(all_res) == length(methods))

  res <- list()
  for(j in 1:p){
    df <- data.frame(value = numeric(), methods = character(), idx = integer())
    for(i in 1:length(all_res)){
      values_i <- lapply(all_res[[i]], function(x){
         as.vector(x[, j])
      }) %>% do.call(rbind, .)

      df_i <- data.frame(mu = apply(values_i, 2, mean),
                         low_b =  apply(values_i, 2, function(x){quantile(x, 0.025)}),
                         upp_b =  apply(values_i, 2, function(x){quantile(x, 0.975)}),
                         methods = methods[i],
                         idx = 1:ncol(values_i))
      df <- rbind(df, df_i)
    }
    res[[j]] <- df
  }
  res
}


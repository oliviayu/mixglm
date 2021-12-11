#'
#' @importFrom penalized cvl
cv_fused_lasso <- function(y, X, lambda1s, lambda2s, fold = 5, penalty.factor = 1){
  ss <- expand.grid(lambda1s, lambda2s)
  cvls <- c()
  for(k in 1:nrow(ss)){
    cv_res <- penalized::cvl(y, X, lambda1 = ss[k, 1]*penalty.factor, lambda2 = ss[k, 2],
                             fusedl = TRUE, model = "linear", fold = fold,
                             epsilon = 1e-5, maxiter = 200)
    cvls <- c(cvls, cv_res$cvl)
    cat("case ", k, " is done. \n")
  }
  res <- data.frame(cbind(ss, cvls))
  names(res) <- c("lambda1", "lambda2", "cvl")
  res
}

#'
cv_mix_glm <- function(fold = 10, lambda1_range = 1e-3, lambda2_range = 3, lambda3_range = .1,
                       y_var, x_vars, w_vars, dat_list, family = family, kappa = 1,
                       weighted = 1, coeff_group_numbers = rep(K, length(x_vars)),
                       beta_int =  as.matrix(nv_coef[, 2:(p_x+1)]),
                       alpha_int = mean(nv_coef[, 1]),
                       tol1 = 1e-3, tol2 = 1e-3, iter_max = 500, check = 1, seed = NULL){

  lambdas <- expand.grid(lambda1_range, lambda2_range, lambda3_range)

  if(!is.null(seed)){ set.seed(seed) }

  subg_dat_list <- lapply(dat_list, function(dat){
    split_data_list(dat, fold = fold, by_id = FALSE)
  })

  res <- data.frame(lambdas, cv_mse = NA)
  names(res)[1:3] <- paste0("lambda", 1:3)

  all_res <- c()

  for(i in 1:nrow(lambdas)){

    mses <- rep(NA, fold)
    for(j in 1:fold){
      training_dat_list <- lapply(subg_dat_list, function(sub_list){
        sub_list[[j]] <- c()
        dd <- do.call(rbind, sub_list)
        rownames(dd) <- c()
        dd
      })
      testing_dat_list <- lapply(subg_dat_list, function(sub_list){
        sub_list[[j]]
      })
      md_i <- try(mix_glm(y_var, x_vars, w_vars, training_dat_list, family,
                          kappa, lambdas[i, 1], lambdas[i, 2], lambdas[i, 3], weighted,
                          coeff_group_numbers=coeff_group_numbers,
                          beta_int=beta_int, alpha_int=alpha_int,
                          tol1=tol1, tol2=tol2, iter_max=iter_max, check=check), silent = TRUE)
      if(is.null(md_i)|class(md_i) == "try-error"){
        break
      } else {
        pred_res <- predict_mixglm(md_i, testing_dat_list, y_var, x_vars, w_vars, family)
        mses[j] <- mean((pred_res[, y_var] - pred_res$mu)^2)

        df_appx <- sum(abs(md_i$alpha) > 0) +
          sum(apply(md_i$beta, 2, function(x){
            unique(round(x[abs(x) > 0], 3)) %>% length()
          }))

        all_res <- rbind(all_res, as.numeric(c(lambdas[i, ], j, mses[j], df_appx)))
        cat("MSE of case=", i, ", fold=", j, ":", mses[j], ",", df_appx, "\n")
      }
    }
    res$cv_mse[i] <- mean(mses)
  }

  all_res <- as.data.frame(all_res)
  names(all_res) <- c(paste0("lambda", 1:3), "fold", "mse", "p")

  list(res, all_res)
}

#'
gcv_mix_glm <- function(lambda1_range = 1e-3, lambda2_range = 3, lambda3_range = .1,
                        y_var, x_vars, w_vars, dat_list, family = family, kappa = 1,
                        weighted = 1,
                        penalty.factor = 1,
                        coeff_group_numbers = rep(K, length(x_vars)),
                        beta_int =  as.matrix(nv_coef[, 2:(p_x+1)]),
                        alpha_int = mean(nv_coef[, 1]),
                        tol1 = 1e-3, tol2 = 1e-3, iter_max = 500, check = 1, seed = NULL){

  lambdas <- expand.grid(lambda1_range, lambda2_range, lambda3_range)
  N <- sapply(dat_list, nrow) %>% sum()

  if(!is.null(seed)){ set.seed(seed) }

  res <- data.frame(lambdas, cv_mse = NA, bic = NA, aic=NA, aicc=NA)
  names(res)[1:3] <- paste0("lambda", 1:3)
  b_N <- 1

  for(i in 1:nrow(lambdas)){

    md_i <- try(mix_glm(y_var, x_vars, w_vars, dat_list, family,
                        kappa, lambdas[i, 1], lambdas[i, 2], lambdas[i, 3],
                        weighted, penalty.factor,
                        coeff_group_numbers, beta_int, alpha_int,
                        tol1, tol2, iter_max, check), silent = TRUE)

    if(is.null(md_i)|class(md_i) == "try-error"){
      next
    } else {
      # df2 <- sum(abs(md_i$alpha) > 0) + sum(sapply(md_i$gamma, function(x){sum(unique(x) > 0)}))
      # num_gamma <- sum(sapply(md_i$gamma, length))

      pred_res <- predict_mixglm(md_i, dat_list, y_var, x_vars, w_vars, family)
      df_appx <- sum(abs(md_i$alpha) > 0) +
        sum(apply(md_i$beta, 2, function(x){
          unique(round(x[abs(x) > 0], 3)) %>% length()
        }))
      cat("appx df = ", df_appx, "\n")

      res$cv_mse[i] <- mean((pred_res[, y_var] - pred_res$mu)^2)/(1 - df_appx/N)^2
      res$aic[i] <- mean((pred_res[, y_var] - pred_res$mu)^2) + 2/N * df_appx
      res$aicc[i] <- mean((pred_res[, y_var] - pred_res$mu)^2) + 2/N * df_appx +
        2*df_appx*(df_appx+1)/(N - df_appx - 1)/N
      res$bic[i] <- mean((pred_res[, y_var] - pred_res$mu)^2) + log(N)/N * df_appx
      # res$mbic[i] <- log(mean((pred_res[, y_var] - pred_res$mu)^2)) + log(N)/N * df_appx
      cat("Case=", i, ": cvMSE =", round(res$cv_mse[i], 2),
          ", BIC =", round(res$bic[i],2),
          ", AIC =", round(res$aic[i],2),
          ", AICC =", round(res$aicc[i],2), "\n")
    }
  }

  res
}

split_data_list <- function(dat, id_var = NULL, fold = 10, by_id = TRUE){
  if(by_id == TRUE){
    unique_id <- unique(dat[, id_var])
    grouped_id <- split(sample(unique_id, length(unique_id)), 1:fold)
    res <- lapply(grouped_id, function(ids){
      dd <- subset(dat, dat[, id_var] %in% ids) %>% arrange_(id_var)
      rownames(dd) <- c()
      dd
    })
    names(res) <- grouped_id
  } else {
    idx <- 1:nrow(dat)
    grouped_idx <- split(sample(idx, length(idx)), 1:fold)
    res <- lapply(grouped_idx, function(i){
      dd <- dat[i, ]
      rownames(dd) <- c()
      dd
    })
  }
  res
}

cv_gflasso <- function(y_all, X_all, fold, edges, s1s, s2s){
  res <- c()
  n_f <- floor(length(y_all)/fold)
  row_split <- split(sample(1:length(y_all), length(y_all)), rep(1:fold, each = n_f))

  for(s1 in s1s){
    for(s2 in s2s){
      mses <- rep(0, fold)

      # for(k in 1:fold){
      #   test_rows <- row_split[[k]]
      #   gmd_k <- try(gflasso(X_all[-test_rows, ],
      #                        y_all[-test_rows],
      #                        edges, s1=s1, s2=s2), silent = TRUE)
      #   mses[k] <- mean((y_all[test_rows] - X_all[test_rows, ]%*%gmd_k$weight)^2)
      # }

      # BIC
      gmd <- try(gflasso(X_all, y_all, edges, s1=s1, s2=s2), silent = TRUE)
      df_appx <- unique(round(gmd$weight[abs(gmd$weight) > 1e-8], 3)) %>% length()
      cat(c(s1, s2,
            mean((y_all - X_all%*%gmd$weight)^2) + log(length(y_all))/length(y_all) * df_appx,
            log(mean((y_all - X_all%*%gmd$weight)^2)) + log(length(y_all))/length(y_all) * df_appx,
            mean(mses)), "\n")
      res <- rbind(res,
                   c(s1, s2,
                     mean((y_all - X_all%*%gmd$weight)^2) + log(length(y_all))/length(y_all) * df_appx,
                     log(mean((y_all - X_all%*%gmd$weight)^2)) + log(length(y_all))/length(y_all) * df_appx,
                     mean(mses)))
      cat("s1=",s1, "s2=", s2, "is done. \n")
      # cat("res=", res[nrow(res),], "\n")
    }
  }
  res <- as.data.frame(res)
  names(res) <- c("s1", "s2", "bic", "mbic", "mse")
  res
}

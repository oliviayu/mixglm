#'
#' @importFrom glmnet glmnet
update_theta <- function(beta_int, alpha_int, y_all, X_all, W_all,
                              Lamb_int, nu_int, kappa, lambda1, l1_weights){

  # Use L1 norm for variable selection
  n <- length(y_all)
  y_str <- rbind(y_all, matrix(sqrt(kappa*n)*(nu_int - Lamb_int/kappa), ncol = 1))
  X_str <- rbind(cbind(X_all, W_all),
                 cbind(sqrt(kappa*n)*diag(1, length(beta_int), length(beta_int)),
                       matrix(0, nrow = length(beta_int), ncol = length(alpha_int))))/sqrt(1 + n*kappa)

  pmd <- glmnet(X_str, y_str, family = "gaussian", alpha = 1, lambda = lambda1/sqrt(1 + n*kappa),
                intercept = FALSE, standardize = FALSE, penalty.factor = l1_weights)
  hat_beta <- pmd$beta / sqrt(1 + n*kappa)

  list(beta = hat_beta[1:length(beta_int)],
       alpha = hat_beta[-c(1:length(beta_int))])

  ## Use SCAD for variable selection
  ## Warnings: SCAD is not a convex function which may not be applicable to ADMM algorithm!
  ## Note: ncpen seems not perform well even for penalty = "lasso".

  # n <- length(y_all)
  # y_str <- rbind(y_all, matrix(sqrt(kappa*n)*(nu_int - Lamb_int/kappa), ncol = 1))
  # X_str <- rbind(cbind(X_all, W_all),
  #                cbind(sqrt(kappa*n)*diag(1, length(beta_int), length(beta_int)),
  #                      matrix(0, nrow = length(beta_int), ncol = length(alpha_int))))/sqrt(1 + n*kappa)
  #
  # pmd <- ncpen::ncpen(y_str, X_str, family = "gaussian",
  #                     penalty = "lasso", x.standardize = FALSE, intercept = FALSE, tau = 1,
  #                     alpha = 1, lambda = lambda1/sqrt(1 + n*kappa))
  # betas <- pmd$beta[, ncol(pmd$beta)]/sqrt(1 + n*kappa)
  # list(beta = as.numeric(betas[1:length(beta_int)]),
  #      alpha = as.numeric(betas[-c(1:length(beta_int))]))
}

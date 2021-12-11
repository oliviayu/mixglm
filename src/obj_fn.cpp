#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <string>

using namespace Rcpp;
using namespace std;
using namespace Eigen;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
double logdet(Eigen::MatrixXd M) {
  LLT<MatrixXd> chol(M);
  double log_det = 2 * chol.matrixL().toDenseMatrix().diagonal().array().log().sum();

  return log_det;
}

// [[Rcpp::export]]
double log_like_glm_cpp(Eigen::VectorXd y_k, Eigen::MatrixXd X_k, Eigen::VectorXd beta_k,
                        Eigen::MatrixXd W_k, Eigen::VectorXd alpha,
                        Eigen::VectorXd sigma2, std::string family) {

  Eigen::VectorXd mu_k = X_k * beta_k;
  if(alpha.size() > 0){
    Eigen::VectorXd mu2_k = W_k * alpha;
    mu_k += mu2_k;
  }

  double log_like_k = 0;
  if(family.compare("gaussian") == 0){
    Eigen::MatrixXd R_k = sigma2.asDiagonal();
    Eigen::MatrixXd invSigma_k = R_k.inverse();
    Eigen::VectorXd err_k = y_k - mu_k;
    log_like_k = - 0.5 * err_k.transpose() * invSigma_k * err_k - 0.5 * logdet(R_k);
  }
  if(family.compare("binomial") == 0){
    Eigen::VectorXd prob_k = mu_k.array().exp() + Eigen::VectorXd::Ones(mu_k.size()).array();
    //std::cout << mu_k << " ";

    double log_like_k_1 = - prob_k.array().log().sum();
    double log_like_k_2 = y_k.transpose() * mu_k;
    log_like_k = log_like_k_1 + log_like_k_2;
  }

  return log_like_k;
}

// [[Rcpp::export]]
double complete_log_like_glm_cpp(List y, List X, List W, List beta, Eigen::VectorXd alpha,
                                 std::string family, Eigen::VectorXd sigma2, bool sigma_fixed,
                                 Eigen::VectorXi sample_size, int K){

  double res = 0;

  for(int k = 0; k < K; ++k){

    int n_k = sample_size(k);
    Eigen::VectorXd beta_k = beta[k];

    Eigen::VectorXd y_k = y[k];
    Eigen::MatrixXd X_k = X[k];
    Eigen::MatrixXd W_k = W[k];

    Eigen::VectorXd sigma2_vec;
    if(sigma_fixed){
      sigma2_vec = Eigen::VectorXd::Ones(n_k) * sigma2;
    } else {
      sigma2_vec = Eigen::VectorXd::Ones(n_k) * sigma2(k);
    }
    double res_k = log_like_glm_cpp(y_k, X_k, beta_k, W_k, alpha, sigma2_vec, family);
    //std::cout << res_k << " ";

    res += res_k;
  }

  return res;
}

// [[Rcpp::export]]
List gr_theta_k_binomial_cpp(Eigen::VectorXd y_k, Eigen::MatrixXd X_k, Eigen::MatrixXd W_k,
                             Eigen::VectorXd beta_k, Eigen::VectorXd alpha) {

  Eigen::VectorXd mu_k = X_k * beta_k;
  if(alpha.size() > 0){
    Eigen::VectorXd mu2_k = W_k * alpha;
    mu_k += mu2_k;
  }

  Eigen::VectorXd exp_mu_k = mu_k.array().exp();
  Eigen::VectorXd temp = (exp_mu_k + Eigen::VectorXd::Ones(mu_k.size())).array().inverse();
  Eigen::VectorXd prob_k = exp_mu_k.array() * temp.array();
  Eigen::VectorXd gr_beta_k =  X_k.transpose()*(prob_k - y_k);

  Eigen::VectorXd gr_alpha_k = Eigen::VectorXd::Zero(alpha.size());
  if(alpha.size() > 0){
    Eigen::VectorXd gr_alpha_k_2 = W_k.transpose()*(prob_k - y_k);
    gr_alpha_k += gr_alpha_k_2;
  }

  List res = List::create(Named("gr_beta") = gr_beta_k,
                          Named("gr_alpha") = gr_alpha_k,
                          Named("gr_sigma") = 0);

  return res;
}

// [[Rcpp::export]]
List gr_theta_k_gaussian_cpp(Eigen::VectorXd y_k, Eigen::MatrixXd X_k, Eigen::MatrixXd W_k,
                             Eigen::VectorXd beta_k, Eigen::VectorXd alpha,
                             Eigen::MatrixXd R_k, double sigma_k) {

  Eigen::VectorXd mu_k = X_k * beta_k;
  if(alpha.size() > 0){
    Eigen::VectorXd mu2_k = W_k * alpha;
    mu_k += mu2_k;
  }

  Eigen::MatrixXd invSigma_k = R_k.inverse();
  Eigen::VectorXd std_err = invSigma_k * (y_k - mu_k);

  Eigen::VectorXd gr_alpha_k = Eigen::VectorXd::Zero(alpha.size());
  if(alpha.size() > 0){
    Eigen::VectorXd gr_alpha_k_2 = - W_k.transpose() * std_err;
    gr_alpha_k += gr_alpha_k_2;
  }
  Eigen::VectorXd gr_beta_k = - X_k.transpose() * std_err;

  double gr_sigma = (invSigma_k.trace() - std_err.transpose() * std_err) * sigma_k;

  List res = List::create(Named("gr_beta") = gr_beta_k,
                          Named("gr_alpha") = gr_alpha_k,
                          Named("gr_sigma") = gr_sigma);

  return res;
}

// [[Rcpp::export]]
List gr_theta_cpp(List y, List X, List W, std::string family,
                  List beta, Eigen::VectorXd alpha, Eigen::VectorXd sigma, bool sigma_fixed,
                  Eigen::VectorXi sample_size, int K, int p_x, int p_w,
                  List nu, List Lamb, double kappa){

  List gr_beta = List::create(); // initiate list
  Eigen::VectorXd gr_alpha = Eigen::VectorXd::Zero(p_w);
  Eigen::VectorXd gr_sigma = Eigen::VectorXd::Zero(K);

  for(int k = 0; k < K; ++k){

    int n_k = sample_size(k);
    Eigen::VectorXd beta_k = beta[k];
    Eigen::VectorXd nu_k = nu[k];
    Eigen::VectorXd Lamb_k = Lamb[k];

    Eigen::VectorXd y_k = y[k];
    Eigen::MatrixXd X_k = X[k];
    Eigen::MatrixXd W_k = W[k];

    List res_k = List::create();

    if(family.compare("gaussian") == 0){
      // double sigma_k = sigma_fixed ? sigma(0) : sigma(k);
      double sigma_k = 0;
      if(sigma_fixed){
        sigma_k = sigma(0);
      } else {
        sigma_k = sigma(k);
      }
      Eigen::VectorXd sigma2_vec = Eigen::VectorXd::Ones(n_k) * pow(sigma_k, 2);
      Eigen::MatrixXd R_k = sigma2_vec.asDiagonal();
      res_k = gr_theta_k_gaussian_cpp(y_k, X_k, W_k, beta_k, alpha, R_k, sigma_k);
    }
    if(family.compare("binomial") == 0){
      res_k = gr_theta_k_binomial_cpp(y_k, X_k, W_k, beta_k, alpha);
    }

    if(p_w > 0){
      Eigen::VectorXd gr_alpha_k = res_k["gr_alpha"];
      gr_alpha += gr_alpha_k / sample_size.sum();
    }
    Eigen::VectorXd gr_beta_k = res_k["gr_beta"];
    Eigen::VectorXd gr_beta_k_final = gr_beta_k / sample_size.sum() + kappa * (beta_k - nu_k) + Lamb_k;
    gr_beta.push_back(gr_beta_k_final);

    double gr_sigma_k = res_k["gr_sigma"];
    gr_sigma(k) = gr_sigma_k / sample_size.sum();
  }

  List res = List::create(Named("gr_beta") = gr_beta,
                          Named("gr_alpha") = gr_alpha,
                          Named("gr_sigma") = gr_sigma);

  return res;
}

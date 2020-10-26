//i = 1:m indexes the 43 PUMAs each of size N_i
//j = 1:10,000 indexes response y_ij where unit j is in PUMA_i
data {
  int<lower = 1> n_obs;   //this is N. (We'll sample 10,000 items)
  int<lower = 1> n_cov;   //this is p  (ncol(Xmat))
  int<lower = 1> n_cty;   //this is m, the number of PUMAs/counties
  int<lower=0, upper=1> y[n_obs];        //bernoulli responses
  int<lower=0, upper=n_cty> cty[n_obs];  //PUMA number for each observation
  matrix[n_obs, n_cov] x; //covariates: X_i is the n_i by p matrix with rows x'_ij
  real kappa;      //fixed = 5
  real sigma_beta; //fixed = sqrt(10)
}
parameters {
  vector[n_cov] beta;       //p-dimensional vector of regression coeff
  real<lower = 0> sigma_u;  //truncate to half-cauchy
  vector[n_cty] u;          //logit intercept prior
}

model {
  y ~ bernoulli_logit(x*beta+u[cty]); 
  u ~ normal(0, sigma_u);
  sigma_u ~ cauchy(0, kappa);
  beta ~ multi_normal(rep_vector(0,n_cov),  
                      diag_matrix(rep_vector(sigma_beta,n_cov)));  
                      //diag_matrix() is discouraged in documentation
}

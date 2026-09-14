data {
  int N; // the number of observations
  int K; // the number of columns in the model matrix
  int P; // the number of columns in the distance matrix (P should match N for multivariate model)
  matrix[N, K] X; // the model matrix
  matrix[P, P] distMat; // the distance matrix for covariance (ensure dimensions match P)
  vector[N] y; // response (using continuous values)
}

parameters {
  vector[K] betas; // the regression parameters
  real<lower=0> ssq; // the variance parameter for sigma^2
  real<lower=0> phi; // the range parameter for the distance matrix for Matern covariance function
}

transformed parameters {
  vector[N] linpred;
  vector[N] mu; // the expected values
  
  // Linear predictor calculation
  linpred = X * betas;
  // No link function applied, so mu = linpred
  mu = linpred; 
}

model {
  matrix[P, P] Sigma = ssq * exp(-distMat / phi); // Covariance matrix
  
  // Priors
  betas ~ normal(0, 1); // prior for betas
  ssq ~ inv_gamma(2, 2); // prior for variance parameter
  phi ~ uniform(0, sqrt(2)); // prior for range parameter
  
  // Likelihood: Multivariate Normal distribution for the response y
  y ~ multi_normal(mu, Sigma); 
}

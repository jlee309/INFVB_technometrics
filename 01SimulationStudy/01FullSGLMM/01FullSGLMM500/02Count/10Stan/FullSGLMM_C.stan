data {
  int N; // the number of observations
  int K; // the number of columns in the model matrix
  int P; // the number of columns in the Identity matrix
  matrix[N, K] X; // the model matrix
  matrix[N, P] M; // the Identity matrix
  matrix[P, P] distMat; // the Distance Matrix for covariance (ensure dimensions match P)
  int<lower=0> y[N]; // response (integer number of successes)
}

parameters {
  vector[K] betas; // the regression parameters
  vector[P] omegas; // the reparameterized random effects
  real<lower=0> ssq; // the variance parameter for sigma^2 
  real<lower=0> phi; // the range parameter for the distance matrix for Matern Covariance function
}

transformed parameters {
  vector[N] linpred;
  vector[N] mu; // the expected values (probabilities of success)
  
  // Linear predictor calculation
  linpred = X * betas + M * omegas;
  
  // Applying the log link function to get expected counts
  mu = exp(linpred); 
}

model {
  matrix[P, P] Sigma = ssq * exp(-distMat / phi); // Covariance matrix
  
  // Priors
  betas ~ normal(0, 10); // prior for betas
  omegas ~ multi_normal(rep_vector(0, P), Sigma); // prior for the random effects
  ssq ~ inv_gamma(2, 2); // prior for variance parameter
  phi ~ uniform(0, sqrt(2)); // prior for range parameter
  
  // Likelihood: Poisson distribution for the response y
  y ~ poisson(mu);
}

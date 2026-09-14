data {
  int N; // the number of observations
  int K; // the number of columns in the model matrix
  int P; // the number of columns in the basis matrix
  matrix[N, K] X; // the model matrix
  matrix[N, P] M; // the basis matrix
  int<lower=0> y[N]; // response (integer number of successes)
}

parameters {
  vector[K] betas; // the regression parameters
  vector[P] deltas; // the reparameterized random effects
  real<lower=0> ssq; // the variance parameter for sigma^2 
}

transformed parameters {
  vector[N] linpred;
  vector[N] mu; // the expected values (probabilities of success)
  
  // Linear predictor calculation
  linpred = X * betas + M * deltas;
  
  // Applying the logit link function to get probabilities
  mu = exp(linpred); 
}

model {  
  betas ~ normal(0, 10); // prior for betas
  ssq ~ inv_gamma(0.1, 0.1); // prior for variance parameter
  deltas ~ normal(0, sqrt(ssq)); // prior for the random effects
  
  // Likelihood: Bernoulli distribution for the response y
  y ~ poisson (mu); 
}

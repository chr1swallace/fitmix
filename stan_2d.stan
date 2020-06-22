data {
 int<lower = 0> N;
  int<lower=0> K;
  matrix[N,2] obsData;
  matrix[N,2] obsVars;
}

parameters {
  real<lower=0> sigma[2];
  real<lower=-1,upper=1> rho[K];
  simplex[K] theta;
}

transformed parameters {
  cov_matrix[2] S[K];
// convert sigmas to usable Variance-Covariance matrix
  S[1] = [[ sigma[1], sigma[1]*rho[1] ],
          [ sigma[1]*rho[1], sigma[1] ]];
  S[2] = [[ sigma[2], sqrt(sigma[1]*sigma[2])*rho[2] ],
          [ sqrt(sigma[1]*sigma[2])*rho[2], sigma[1] ]];
  S[3] = [[ sigma[1], sqrt(sigma[1]*sigma[2])*rho[3] ],
          [ sqrt(sigma[1]*sigma[2])*rho[3], sigma[2] ]];
  S[4] = [[ sigma[2], sigma[2]*rho[4] ],
          [ sigma[2]*rho[4], sigma[2] ]];
}

model {
  sigma[1] ~ lognormal(-7, 0.1);
  sigma[2] ~ lognormal(-3, 0.1);
  for(k in 1:K) {
    rho[k] ~ uniform(-1,1);
  }
  theta ~ dirichlet(rep_vector(1.0,K)); 
  /* print(add_diag(S[1], obsVars[1,])) */
  /* print(add_diag(S[2], obsVars[1,])) */
  /* print(add_diag(S[3], obsVars[1,])) */
  /* print(add_diag(S[4], obsVars[1,])) */
  print("log density before =", target());
  for (n in 1:N) {
    vector[K] lps; 
    for(k in 1:K) {
      lps[k] = log(theta[k]) +
        multi_normal_lpdf(obsData[n,] | [0,0]', add_diag(S[k], obsVars[n,]));
      /* print("log density after n=",n," k=",k," is ",target()); */
}
    target += log_sum_exp(lps);
  }
}

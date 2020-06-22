data {
 int<lower = 0> N;
 vector[N] obsData;
 vector[N] obsVars;
}

parameters {
  real<lower=0> sigma[2];
  real<lower=0, upper=1> theta;
}

model {
  sigma ~ lognormal(log(0.0225), 1);
  theta ~ beta(5, 1); // prob of mixture 1 - deliberately skewed to right
 for (n in 1:N)
   target += log_mix(theta,
                     normal_lpdf(obsData[n] | 0, sigma[1] + obsVars[n]),
                     normal_lpdf(obsData[n] | 0, sigma[2] + obsVars[n]));
}

functions {
  real log_pochhammer(real a, int m) {
    return lgamma(a + m) - lgamma(a);
  }

  real py_partition_log(int[] n_j, int K, real alpha, real sigma) {
    int N = 0;
    real out = 0.0;

    for (k in 1:K)
      N += n_j[k];

    if (alpha <= 0 || sigma <= 0 || sigma >= 1)
      return negative_infinity();

    for (i in 1:(K - 1))
      out += log(alpha + i * sigma);

    out -= log_pochhammer(alpha + 1, N - 1);
    out += K * log(sigma);

    for (j in 1:K)
      out += log_pochhammer(1 - sigma, n_j[j] - 1);

    return out;
  }
}

data {
  int<lower=1> K;
  int<lower=1> n_j[K];
}

parameters {
  real alpha_raw;
  real sigma_raw;
}

transformed parameters {
  real<lower=0> alpha = exp(alpha_raw);
  real<lower=0, upper=1> sigma = inv_logit(sigma_raw);
}

model {
  alpha_raw ~ normal(log(5), 0.5);
  sigma_raw ~ normal(logit(0.6), 1);

  target += py_partition_log(n_j, K, alpha, sigma);
}

generated quantities {
  real alpha_out = alpha;
  real sigma_out = sigma;
}

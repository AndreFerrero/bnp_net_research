// file: py_sampler_sigma.stan

functions {
  // rising Pochhammer: (a)_m = Gamma(a+m)/Gamma(a)
  real log_pochhammer(real a, int m) {
    return lgamma(a + m) - lgamma(a);
  }

  // Pitman–Yor partition *log-likelihood* (not an _lpdf function)
  real eppf_lp(int[] n_j, int K, real alpha, real sigma) {
    int N = 0;
    real out = 0.0;

    // Total sample size
    for (k in 1:K)
      N += n_j[k];

    // Enforce domain
    if (alpha <= 0 || sigma <= 0 || sigma >= 1)
      return negative_infinity();

    // 1) sum_{i=1}^{K-1} log(alpha + i*sigma)
    for (i in 1:(K-1))
      out += log(alpha + i * sigma);

    // 2) - log_pochhammer(alpha+1, N-1)
    out -= log_pochhammer(alpha + 1, N - 1);

    // 3) sum_j log_pochhammer(1-sigma, n_j[j]-1)
    for (j in 1:K)
      out += log_pochhammer(1 - sigma, n_j[j] - 1);

    return out;
  }
}

data {
  int<lower=1> K;          // number of clusters
  int<lower=1> n_j[K];     // cluster-size counts
  real<lower=0> a_sigma;   // Prior for sigma
  real<lower=0> b_sigma;   // Prior for sigma
  real<lower=0> alpha_fixed;  // Fixed value for alpha
}

parameters {
  real<lower=0, upper=1> sigma;  // discount ∈ (0,1)
}

model {

  // Prior for sigma
  sigma ~ beta(a_sigma, b_sigma);

  // Add the Pitman–Yor partition log-likelihood with fixed alpha
  target += eppf_lp(n_j, K, alpha_fixed, sigma);
}

// file: py_sampler_sigma_transf.stan

functions {
  // rising Pochhammer: (a)_m = Γ(a+m) / Γ(a)
  real log_pochhammer(real a, int m) {
    return lgamma(a + m) - lgamma(a);
  }

  // Pitman–Yor partition *log-likelihood* (EPPF)
  real eppf_lp(int[] n_j, int K, real alpha, real sigma) {
    int N = 0;
    real out = 0.0;
    // total N
    for (k in 1:K) N += n_j[k];

    // support
    if (alpha <= 0 || sigma <= 0 || sigma >= 1)
      return negative_infinity();

    // 1) ∑_{i=1}^{K-1} log(α + i σ)
    for (i in 1:(K-1))
      out += log(alpha + i * sigma);

    // 2) − log_pochhammer(α+1, N−1)
    out -= log_pochhammer(alpha + 1, N - 1);

    // 3) ∑_j log_pochhammer(1−σ, n_j[j]−1)
    for (j in 1:K)
      out += log_pochhammer(1 - sigma, n_j[j] - 1);

    return out;
  }
}

data {
  int<lower=1>   K;             // # clusters
  int<lower=1>   n_j[K];        // cluster counts
  real<lower=0>  mu;
  real<lower=0>  tau;
  real<lower=0>  alpha_fixed;   // fixed Pitman–Yor α
}

parameters {
  real eta;                     // unconstrained parameter
}

transformed parameters {
  real<lower=0,upper=1> sigma = inv_logit(eta);
}

model {
  // prior on transformed parameter
  eta ~ logistic(mu, tau);

  // --- Pitman–Yor EPPF log-likelihood ---
  target += eppf_lp(n_j, K, alpha_fixed, sigma);
}

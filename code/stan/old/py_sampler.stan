// file: py_sampler_sigma_transf.stan

functions {
  real eppf_lp(int[] n_j, int K, real alpha, real sigma) {
    int N = sum(n_j);

    if (alpha <= 0 || sigma <= 0 || sigma >= 1)
      return negative_infinity();

    // Term 1: sum(log(alpha + i * sigma)) for i = 1 to K - 1
    vector[K - 1] i_vec = linspaced_vector(K - 1, 1, K - 1);
    real term1 = sum(log(alpha + sigma * i_vec));

    // Term 2: lgamma(alpha + N) - lgamma(alpha + 1)
    real term2 = lgamma(alpha + N) - lgamma(alpha + 1);

    // Term 3: sum of lgamma(n_j - sigma) - K * lgamma(1 - sigma)
    // Use to_vector to convert int[] to vector, then use element-wise ops
    vector[K] n_j_vec = to_vector(n_j);
    vector[K] gamma_terms = lgamma(n_j_vec - sigma);  // vectorized lgamma
    real term3 = sum(gamma_terms) - K * lgamma(1 - sigma);

    return term1 - term2 + term3;
  }
}


data {
  int<lower=1>   K;             // # clusters
  int<lower=1>   n_j[K];        // cluster counts
  vector<lower=0>[2]  prior_alpha;
  vector<lower=0>[2]  prior_sigma;
}

parameters {
  real <lower=0,upper=1> sigma ;
  real<lower=0> alpha;
}

model {
  // prior on alpha
  alpha ~ gamma(prior_alpha[1],prior_alpha[2]);
  // prior on sigma
  sigma ~ beta(prior_sigma[1], prior_sigma[2]);
  

  // --- Pitman–Yor EPPF log-likelihood ---
  target += eppf_lp(n_j, K, alpha, sigma);
}

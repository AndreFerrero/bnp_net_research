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
  
  void py_sample_rng(int N,
                     real alpha,
                     real sigma,
                     int[] x,      // <- unsized array arg
                     int[] counts, // <- unsized array arg
                     int K_out) {
    int active = 1;
    // initialize counts & first assignment
    for (i in 1:N) counts[i] = 0;
    x[1]       = 1;
    counts[1]  = 1;

    // sequential seating
    for (n in 2:N) {
      vector[active] cnt_vec = to_vector(counts[1:active]);
      vector[active + 1] prob;
      prob[1:active]   = cnt_vec - sigma;
      prob[active + 1] = alpha + sigma * active;

      int draw = categorical_rng(prob);
      x[n] = draw;
      if (draw == active + 1)
        active += 1;
      counts[draw] += 1;
    }

    K_out = active;
  }
}

data {
  int<lower=1>   K_A;             // nodes
  int<lower=1>   K_B;
  int<lower=1>   n_A[K_A];        // node counts
  int<lower=1>   n_B[K_B];
  vector<lower=0>[2]  prior_alpha_A;
  vector<lower=0>[2]  prior_alpha_B;
  vector<lower=0>[2]  prior_sigma_A;
  vector<lower=0>[2]  prior_sigma_B;
}

parameters {
  real<lower = 0> alpha_A;
  real<lower = 0> alpha_B;
  real<lower = 0, upper = 1> sigma_A;
  real<lower = 0, upper = 1> sigma_B;
  int<lower=1> N_edges;
}

model {
  alpha_A ~ gamma(prior_alpha_A[1],prior_alpha_A[2]);
  alpha_B ~ gamma(prior_alpha_B[1],prior_alpha_B[2]);
  
  sigma_A ~ beta(prior_sigma_A[1], prior_sigma_A[2]);
  sigma_B ~ beta(prior_sigma_B[1], prior_sigma_B[2]);
  
  target += eppf_lp(n_A, K_A, alpha_A, sigma_A) +
  eppf_lp(n_B, K_B, alpha_B, sigma_B);
}

generated quantities {
  //––– Simulated PY draws –––––––––––––––––––––––––––––––
  int xA_sim[N_edges];
  int countsA[N_edges];
  int K_A_sim;

  int xB_sim[N_edges];
  int countsB[N_edges];
  int K_B_sim;

  py_sample_rng(N_edges, alpha_A, sigma_A,
                xA_sim, countsA, K_A_sim);

  py_sample_rng(N_edges, alpha_B, sigma_B,
                xB_sim, countsB, K_B_sim);

  //––– Count unique (A,B) pairs –––––––––––––––––––––––––
  int pair_counts[N_edges * N_edges];
  for (i in 1:(N_edges * N_edges))
    pair_counts[i] = 0;

  for (n in 1:N_edges) {
    int idx = (xA_sim[n] - 1) * K_B_sim + xB_sim[n];
    pair_counts[idx] += 1;
  }

  int unique_pairs = 0;
  for (i in 1:(K_A_sim * K_B_sim)) {
    if (pair_counts[i] > 0)
      unique_pairs += 1;
  }

  //––– Posterior predictive density –––––––––––––––––––––
  real density_sim = unique_pairs / ((real) K_A_sim * K_B_sim);
}

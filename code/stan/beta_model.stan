functions {
  // Log EPPF for Pitman–Yor process
  real eppf_lp(array[] int n_j, int K, real alpha, real sigma) {
    int N = sum(n_j);
    if (alpha <= 0 || sigma < 0 || sigma >= 1) 
      return negative_infinity();
    
    // Term 1: \sum_{i=1}^{K-1} log(alpha + i * sigma)
    vector[K - 1] i_vec = linspaced_vector(K - 1, 1, K - 1);
    real term1 = sum(log(alpha + sigma * i_vec));
    
    // Term 2: lgamma(alpha + N) - lgamma(alpha + 1)
    real term2 = lgamma(alpha + N) - lgamma(alpha + 1);
    
    // Term 3: \sum_j lgamma(n_j - sigma) - K * lgamma(1 - sigma)
    vector[K] n_j_vec = to_vector(n_j);
    vector[K] gamma_terms = lgamma(n_j_vec - sigma);
    real term3 = sum(gamma_terms) - K * lgamma(1 - sigma);
    
    return term1 - term2 + term3;
  }
}
data {
  int<lower=1> K_A; // number of A-side nodes
  int<lower=1> K_B; // number of B-side nodes
  array[K_A] int<lower=1> n_A; // counts for A
  array[K_B] int<lower=1> n_B; // counts for B
  
  vector<lower=0>[2] prior_alpha_A; // Gamma(shape, rate)
  vector<lower=0>[2] prior_alpha_B;
  
  real<lower=0> eps; // small epsilon for spike~Beta(eps,1)
}
parameters {
  real<lower=0> alpha_A;
  real<lower=0> alpha_B;
  real<lower=0, upper=1> sigma_A;
  real<lower=0, upper=1> sigma_B;
}
model {
  // Priors on alpha
  alpha_A ~ gamma(prior_alpha_A[1], prior_alpha_A[2]);
  alpha_B ~ gamma(prior_alpha_B[1], prior_alpha_B[2]);
  
  // Approximate spike-and-slab prior for sigma
  sigma_A ~ beta(eps, 1);
  sigma_B ~ beta(eps, 1);
  
  // Likelihood contributions via EPPF
  target += eppf_lp(n_A, K_A, alpha_A, sigma_A)
            + eppf_lp(n_B, K_B, alpha_B, sigma_B);
}

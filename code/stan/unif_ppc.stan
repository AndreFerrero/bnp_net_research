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
  vector<lower=0>[2] prior_alpha_A;
  vector<lower=0>[2] prior_alpha_B;
  vector<lower=0>[2] prior_sigma_A;
  vector<lower=0>[2] prior_sigma_B;
  
  int<lower=1> e_obs; // observed edges
}
parameters {
  real<lower=0> alpha_A;
  real<lower=0> alpha_B;
  real<lower=0, upper=1> sigma_A;
  real<lower=0, upper=1> sigma_B;
}
model {
  alpha_A ~ gamma(prior_alpha_A[1], prior_alpha_A[2]);
  alpha_B ~ gamma(prior_alpha_B[1], prior_alpha_B[2]);
  
  sigma_A ~ beta(prior_sigma_A[1], prior_sigma_A[2]);
  sigma_B ~ beta(prior_sigma_B[1], prior_sigma_B[2]);
  
  target += eppf_lp(n_A, K_A, alpha_A, sigma_A)
            + eppf_lp(n_B, K_B, alpha_B, sigma_B);
}
generated quantities {
  real density_ppc; // save only density
  
  {
    // local scope to discard the sampled network 
    int maxN = e_obs; // upper limit on number of draws
    array[maxN] int countsA;
    array[maxN] int countsB;
    array[maxN, maxN] int seen;
    int unique_edges;
    int KactA;
    int KactB;
    int currA;
    int currB;
    
    // Create zero valued objects
    for (i in 1 : maxN) {
      countsA[i] = 0;
      countsB[i] = 0;
      for (j in 1 : maxN) 
        seen[i, j] = 0;
    }
    
    // Initialize first edge = (1,1)
    countsA[1] = 1;
    countsB[1] = 1;
    seen[1, 1] = 1;
    unique_edges = 1;
    KactA = 1;
    KactB = 1;
    
    for (n in 2 : e_obs) {
      // --- A‐side ---
      {
        real denomA = alpha_A + n - 1;
        vector[KactA] cA = to_vector(countsA)[1 : KactA];
        vector[KactA + 1] pA;
        pA[1 : KactA] = (cA - sigma_A) / denomA;
        pA[KactA + 1] = (alpha_A + sigma_A * KactA) / denomA;
        currA = categorical_rng(pA);
        countsA[currA] += 1;
        if (currA == KactA + 1) 
          KactA += 1;
      }
      
      // --- B‐side ---
      {
        real denomB = alpha_B + n - 1;
        vector[KactB] cB = to_vector(countsB)[1 : KactB];
        vector[KactB + 1] pB;
        pB[1 : KactB] = (cB - sigma_B) / denomB;
        pB[KactB + 1] = (alpha_B + sigma_B * KactB) / denomB;
        currB = categorical_rng(pB);
        countsB[currB] += 1;
        if (currB == KactB + 1) 
          KactB += 1;
      }
      
      // Track if new edge is unique
      if (seen[currA, currB] == 0) {
        seen[currA, currB] = 1;
        unique_edges += 1;
      }
    }
    
    // Compute density
    density_ppc = unique_edges / (1.0 * KactA * KactB);
  }
}

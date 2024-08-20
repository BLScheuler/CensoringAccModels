functions { 
  //real lba_cumhaz(real rt, real t0, real bMinusA, real A, real v) { 
  //  real muT;
  //  real b;
  //  array[4] real terms;
  //  b = bMinusA + A;

  //  if(!(rt-t0 >=0))
  //    return(0);
  //  if(!(b >=0))
  //    reject("lba_cumhaz(): b must be positive; found b=", b);
  //  if(!(A >=0))
  //    reject("lba_cumhaz(): A must be positive; found b=", A);
  //  if(!(b-A >=0))
  //    reject("lba_cumhaz(): b-A must be positive; found b-A=", b-A);

  //  muT = b - (rt-t0)*v;
  //    terms[1] = log(A-muT) - log(A) + std_normal_lcdf((muT-A)/(rt-t0));
  //    terms[2] = log(muT) - log(A) + std_normal_lcdf(muT/(rt-t0));
  //    terms[3] = log(t0-rt) - log(A) + std_normal_lpdf((muT-A)/(rt-t0));
  //    terms[4] = log(rt-t0) - log(A) + std_normal_lpdf(muT/(rt-t0));
  //  return(-log_sum_exp(terms));
  //}
  real lba_cumhaz(real rt, real t0, real bMinusA, real A, real v) {
    real muT;
    real b;
    array[4] real terms;
    b = bMinusA + A;

    if(!(rt-t0 >=0))
      return(0);
    if(!(b >=0))
      reject("lba_cumhaz(): b must be positive; found b=", b);
    if(!(A >=0))
      reject("lba_cumhaz(): A must be positive; found b=", A);
    if(!(b-A >=0))
      reject("lba_cumhaz(): b-A must be positive; found b-A=", b-A);

    muT = b - (rt-t0)*v;
    terms[1] = (A-muT)/A * std_normal_cdf((muT-A)/(rt-t0));
    terms[2] = muT/A * std_normal_cdf(muT/(rt-t0));
    terms[3] = -(rt-t0)/A * exp(std_normal_lpdf((muT-A)/(rt-t0)));
    terms[4] = (rt-t0)/A * exp(std_normal_lpdf(muT/(rt-t0)));
    return(-log(sum(terms)));
  }
}
data {
  int<lower=1> N;
  array[N] real x;
  vector[N] y;
}
transformed data {
  matrix[N,N] Kbm;
  for (n in 1:N) {
    for (m in 1:N) {
      Kbm[n, m] = min([x[n], x[m]]);
    }
  }
}
parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
  real<lower=0> sigma_diff;
  real<lower=0> bMinusA;
  real<lower=0> A;
  real<lower=0> t0;
  real<lower=0> mean_v;
}
transformed parameters {
  real<lower=0> b;
  vector[N] mu;
  for (s in 1:N) {
    mu[s] = lba_cumhaz(x[s], t0, bMinusA, A, mean_v);
  }
  b = bMinusA + A;
}
model {
  matrix[N, N] L_K;
  // matrix[N, N] K = gp_exp_quad_cov(x, alpha, rho);
  matrix[N, N] K = Kbm;
  real sq_sigma = square(sigma);
  real sq_sigma_diff = square(sigma_diff);

  // diagonal elements 
  K = K * sq_sigma_diff;
  for (n in 1:N)
    K[n, n] = K[n, n] + sq_sigma;

  L_K = cholesky_decompose(K);

  //rho ~ inv_gamma(5, 5);
  //alpha ~ std_normal();
  sigma ~ std_normal();

  bMinusA ~ normal(.4, .4);
  A ~ normal(1, .5);
  t0 ~ normal(.3, .3);
  mean_v ~ normal(1, .1);

  y ~ multi_normal_cholesky(mu, L_K);
}

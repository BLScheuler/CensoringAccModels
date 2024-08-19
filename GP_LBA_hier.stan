functions { 
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
  int<lower=1> Nt;
  int<lower=1> Nsubj;
  array[Nt] real x;
  matrix[Nt,Nsubj] y;
}
parameters {
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
  real<lower=0> group_bMinusA;
  real<lower=0> group_A;
  real<lower=0> group_t0;
  real<lower=0> group_mean_v;
  vector<lower=0>[Nsubj] bMinusA;
  vector<lower=0>[Nsubj] A;
  vector<lower=0>[Nsubj] t0;
  vector<lower=0>[Nsubj] mean_v;
}
transformed parameters {
  vector<lower=0>[Nsubj] b;
  matrix[Nt,Nsubj] mu;
  for (sx in 1:Nsubj) {
    for (tx in 1:Nt) {
      mu[tx,sx] = lba_cumhaz(x[tx], t0[sx], bMinusA[sx], A[sx], mean_v[sx]);
    }
    b[sx] = bMinusA[sx] + A[sx];
  }
}
model {
  matrix[Nt, Nt] K = gp_exp_quad_cov(x, alpha, rho);
  matrix[Nt, Nt] L_K;
  real sq_sigma = square(sigma);

  // diagonal elements 
  for (n in 1:Nt)
    K[n, n] = K[n, n] + sq_sigma * x[n];
  
  L_K = cholesky_decompose(K);

  rho ~ inv_gamma(5, 5);
  alpha ~ normal(0, .1);
  sigma ~ normal(0, .1);

  group_bMinusA ~ normal(10, 3);
  group_A ~ normal(1, .5);
  group_t0 ~ normal(250, 10);
  group_mean_v ~ normal(.03, .01);

  bMinusA ~ normal(group_bMinusA, 3);
  A ~ normal(group_A, .5);
  t0 ~ normal(group_t0, 50);
  mean_v ~ normal(group_mean_v, .01);

  for (sj in 1:Nsubj) {
    y[,sj] ~ multi_normal_cholesky(mu[,sj], L_K);
  }
}

functions{
  
  real lba_pdf(real t, real b, real A, real v, real s){
    //PDF of the LBA model
    
    real b_A_tv_ts;
    real b_tv_ts;
    real term_1;
    real term_2;
    real term_3;
    real term_4;
    real pdf;
    
    b_A_tv_ts = (b - A - t*v)/(t*s);
    b_tv_ts = (b - t*v)/(t*s);
    term_1 = v*Phi(b_A_tv_ts);
    term_2 = s*exp(normal_lpdf(b_A_tv_ts|0,1)); 
    term_3 = v*Phi(b_tv_ts);
    term_4 = s*exp(normal_lpdf(b_tv_ts|0,1)); 
    pdf = (1/A)*(-term_1 + term_2 + term_3 - term_4);
    
    return pdf;
  }
  
  real lba_cdfx(real t, real b, real A, real v, real s){
    //CDF of the LBA model
    
    real b_A_tv;
    real b_tv;
    real ts;
    real term_1;
    real term_2;
    real term_3;
    real term_4;
    real cdf;	
    
    b_A_tv = b - A - t*v;
    b_tv = b - t*v;
    ts = t*s;
    term_1 = b_A_tv/A * Phi(b_A_tv/ts);	
    term_2 = b_tv/A   * Phi(b_tv/ts);
    term_3 = ts/A     * exp(normal_lpdf(b_A_tv/ts|0,1)); 
    term_4 = ts/A     * exp(normal_lpdf(b_tv/ts|0,1)); 
    cdf = 1 + term_1 - term_2 + term_3 - term_4;
    return cdf;
  }
  
  real lba_lpdf(matrix RT, real bMinusA, real A, vector v, real s, real tau){
  
    real t;
    real b;
    real cdf;
    real pdf;		
    vector[rows(RT)] prob;
    real out;
    real prob_neg;

    b = A + bMinusA;
    for (i in 1:rows(RT)){	
      t = RT[i,1] - tau;
      if(t > 0){			
        cdf = 1;
        
        for(j in 1:num_elements(v)){
          if(RT[i,2] == j){
            pdf = lba_pdf(t, b, A, v[j], s);
          }else{	
            cdf = (1-lba_cdfx(t, b, A, v[j], s)) * cdf;
          }
        }
        prob_neg = 1;
        for(j in 1:num_elements(v)){
          prob_neg = Phi(-v[j]/s) * prob_neg;    
        }
        prob[i] = pdf*cdf;		
        prob[i] = prob[i]/(1-prob_neg);	
        if(prob[i] < 1e-10){
          prob[i] = 1e-10;				
        }
        
      }else{
        prob[i] = 1e-10;			
      }		
    }
    out = sum(log(prob));
    return out;		
  }
  
  vector lba_rng(real bMinusA, real A, vector v, real s, real tau){
  
    int get_pos_drift;	
    int no_pos_drift;
    int get_first_pos;
    vector[num_elements(v)] drift;
    int max_iter;
    int iter;
    array[num_elements(v)] real start;
    array[num_elements(v)] real ttf;
    array[num_elements(v)] int resp;
    real rt;
    vector[2] pred;
    real b;
    
    //try to get a positive drift rate
    get_pos_drift = 1;
    no_pos_drift = 0;
    max_iter = 1000;
    iter = 0;
    while(get_pos_drift){
      for(j in 1:num_elements(v)){
        drift[j] = normal_rng(v[j],s);
        if(drift[j] > 0){
          get_pos_drift = 0;
        }
      }
      iter = iter + 1;
      if(iter > max_iter){
        get_pos_drift = 0;
        no_pos_drift = 1;
      }	
    }
    //if both drift rates are <= 0
    //return an infinite response time
    if(no_pos_drift){
      pred[1] = -1;
      pred[2] = -1;
    }else{
      b = A + bMinusA;
      for(i in 1:num_elements(v)){
        //start time of each accumulator	
        start[i] = uniform_rng(0,A);
        //finish times
        ttf[i] = (b-start[i])/drift[i];
      }
      //rt is the fastest accumulator finish time	
      //if one is negative get the positive drift
      resp = sort_indices_asc(ttf);
      ttf = sort_asc(ttf);
      get_first_pos = 1;
      iter = 1;
      while(get_first_pos){
        if(ttf[iter] > 0){
          pred[1] = ttf[iter] + tau;
          pred[2] = resp[iter]; 
          get_first_pos = 0;
        }
        iter = iter + 1;
      }
    }
    return pred;	
  }
}

data{
  int N;
  matrix[N,2] RT;
  int NUM_CHOICES;
}

parameters {
  real<lower=0> bMinusA;
  real<lower=0> A;
  real<lower=0> tau;
  vector<lower=0>[NUM_CHOICES] v;
}

transformed parameters {
  real s;
  s = 1;
}

model {
  bMinusA ~ normal(.4,.4)T[0,];
  A ~ normal(1,.5)T[0,];
  tau ~ normal(.3,.3)T[0,];
  for(n in 1:NUM_CHOICES) { 
    v[n] ~ normal(1,.1)T[0,];
  }
  RT ~ lba(bMinusA,A,v,s,tau);
}

generated quantities {
  vector[2] pred;
  pred = lba_rng(bMinusA,A,v,s,tau);
}

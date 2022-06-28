/**
 * Bradley-Terry Model
 *
 * full model
 *
 * hyp controls the hyper_priors
 * scale is the scale parameter of the hyper_priors
 * use_davidson implements the Davidson extension of the BT model
 * use_log_lik to compute the log likelyhood
 * ties   the number of ties for each comparison
 */
data {
  int<lower = 0, upper = 3> hyp;
  real<lower = 0.0>       scale;
  int<lower = 0, upper = 1> use_davidson; // use Davidson model?
  int<lower = 0, upper = 1> use_log_lik;  // compute log likehood?
  int<lower = 1> K;                     // players
  int<lower = 1> N;                     // pairs
  array[N] int<lower=1, upper = K> player1;   // player 1 for pairs n
  array[N] int<lower=1, upper = K> player2;   // player 2 for pairs n
  array[N] int<lower = 0> win1;     // number of wins for player 1
  array[N] int<lower = 0> win2;     // number of wins for player 2
  array[N] int<lower = 0> ties;
}

transformed data{
  array[N] int<lower = 0> n;
  array[use_davidson ? N : 0] int<lower = 0> nn;

   for (i in 1:N) {
      n[i] = win1[i] + win2[i];
      if (use_davidson) nn[i] = n[i]+ties[i];
   }
}

parameters {
  real<lower = 0.001> sigma;                // scale of ability variation
  real<lower = 0.001> sigmanu;
  real<lower = 0> nu;                  // Davidson
  vector[K] beta;                      // ability for player n
}

transformed parameters{
  vector[K] w;
  vector[N] pwin;
  vector[use_davidson ? N : 0] ptie;

  w = exp(beta);
  for (i in 1:N){
    real aux;
    if (use_davidson) {
      aux = exp(nu+(beta[player1[i]]+beta[player2[i]])/2.0);
      ptie[i] =  aux/(w[player1[i]]+w[player2[i]]+aux);
      if (is_nan(ptie[i])) ptie[i] = 0.0 ; // get some nan errors sometimes
    } else {
      aux = 0.0;
    }
    pwin[i] =  w[player1[i]]/(w[player1[i]]+w[player2[i]]+aux);
    if (is_nan(pwin[i])) pwin[i] = 0.0 ; // get some nan errors sometimes

  }
}

model {

  if (hyp==0)      sigma ~ lognormal(0, 0.5) ;
  else if (hyp==1) sigma ~ lognormal(0, scale);
  else if (hyp==2) sigma ~ cauchy(0,scale);
  else if (hyp==3) sigma ~ normal(0,scale);

  sigmanu ~ lognormal(0, 0.5);

  nu ~ normal(0, sigmanu);
  beta ~ normal(0, sigma);

  win1 ~ binomial(n, pwin);
  if (use_davidson) ties ~ binomial(nn, ptie);
}

generated quantities {
  array[N] int win1_rep;
  array[use_davidson ? N : 0] int tie_rep;
  array[use_log_lik ? N : 0] real log_lik;

  win1_rep = binomial_rng(n, pwin) ;

  if (use_davidson) tie_rep = binomial_rng(nn, ptie) ;

  if (use_log_lik){
    for (i in 1:N) {
      if (use_davidson) {
        log_lik[i] = binomial_lpmf(ties[i] | nn[i], ptie[i]) + binomial_lpmf( win1[i] | n[i] , pwin[i]);
      } else {
        log_lik[i] = binomial_lpmf( win1[i] | n[i] , pwin[i]);
    }}}

}



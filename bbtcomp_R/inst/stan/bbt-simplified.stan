/**
 * Bradley-Terry Model
 *
 * simplified version
 * performs the basic BT model (not davidson)
 * hyper_prior fixed
 * generate the win1_rep for ppc
 * uses binomial_logit
 *
 *
 */
 
data {
  int<lower=1> K;                           // players
  int<lower=1> N;                           // pairs          
  array[N] int<lower=1, upper=K> player1;   // player 1 for pairs n
  array[N] int<lower=1, upper=K> player2;   // player 2 for pairs n
  array[N] int<lower=0> win1;               // number of wins for player 1
  array[N] int<lower=0> win2;               // number of wins for player 2
}

transformed data{
  
  array[N] int<lower=0> n;
  for (i in 1:N) n[i] = win1[i] + win2[i];
}

parameters {
  
  real<lower = 0.001> sigma;                // scale of ability variation
  vector[K] beta;                           // ability for player n
}

model {
  
  sigma ~ lognormal(0, 0.5) ;
  beta ~ normal(0, sigma);                  // hierarchical
  win1 ~ binomial_logit(n, beta[player1] - beta[player2]);
}

generated quantities {
  array[N] int win1_rep;
  
  win1_rep = binomial_rng(n,inv_logit(beta[player1]-beta[player2])) ;
}


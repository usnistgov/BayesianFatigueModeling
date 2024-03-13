functions {
  real logistic(real b0, real b1, real xval) {
    return 1/(1+exp(-(b0+b1*(xval))));
  }
}


data {
  int<lower=0> N; // total data size
  int<lower=0> N_cens; // number of runs > 1e7
  int<lower=0> N_obs; // number of runouts
  int<lower=0> N_pred; // number of stres locations for predictions
  real<lower=0> x_obs[N_obs]; // runout stress values
  real<lower=0> x_cens[N_cens]; // non-runout stress values
  real x_pred[N_pred]; // stress locations for predictions
  real all_x[N]; // all stress values
  real all_y[N]; // all lifetime values
  vector[N_obs] y_obs; // all runout lifetimes values
  
  int is_7_obs[N_obs];
  int is_7_cens[N_cens];
  
  // prior parameters
  real beta0_prior_loc; 
  real beta1_prior_loc;
  real<lower=0> sigma_prior_scale;
  real logistic_b0_loc;
  real logistic_b0_scale;
  real logistic_b1_loc;
  real logistic_b1_scale;
}


parameters {
  real<lower=0> beta0;
  real beta1;
  real<lower=0> sigma;
  real logistic_b0;
  real<upper=0> logistic_b1;
}


model {
  
  // priors for linear model 
  beta0 ~ normal(beta0_prior_loc,abs(beta0_prior_loc));
  beta1 ~ normal(beta1_prior_loc,abs(beta1_prior_loc));
  sigma ~ normal(0,sigma_prior_scale);
  
  // priors for logistic model
  logistic_b0 ~ normal(logistic_b0_loc,logistic_b0_scale);
  logistic_b1 ~ normal(logistic_b1_loc,logistic_b1_scale);
  
  // likelihood for observed data
  for(ii in 1:N_obs) {
    is_7_obs[ii] ~ bernoulli(logistic(logistic_b0,logistic_b1,x_obs[ii]));
    y_obs[ii] ~ normal(beta0 + beta1*x_obs[ii],sigma);
  }
  
  // likelihood for censored data
  for(ii in 1:N_cens) {
    is_7_cens[ii] ~ bernoulli(logistic(logistic_b0,logistic_b1,x_cens[ii]));
  }
  
}

generated quantities {
  // generated predictions to compute prediction intervals
  
  real y_rep[N_pred];
  real logistic_probs[N_pred];
  int runout_rep[N_pred];
  real mean_func[N_pred];
  
  for(ii in 1:N_pred) {
    logistic_probs[ii] = logistic(logistic_b0,logistic_b1,x_pred[ii]);
    runout_rep[ii] = bernoulli_rng(logistic_probs[ii]);
    y_rep[ii] = (runout_rep[ii])*(-.01) + (1-runout_rep[ii])*normal_rng(beta0 + beta1*x_pred[ii],sigma);
    mean_func[ii] = logistic_probs[ii]*(-.01) + (1-logistic_probs[ii])*(beta0+beta1*x_pred[ii]);
  }
  
}
#-------------------------------------------------------------------------------
#                           AUGMENTED GAUSSIAN MODEL
#-------------------------------------------------------------------------------


Run_Aug_Gaussian_Mod <- function(dataList, modelName) {
  
  write("
    data {
      int<lower=1> nSubj;
      int<lower=1> nStim;
      real xs[nStim];
      real<lower=0,upper=100> responses[nSubj,nStim];
    }
    
    parameters {
      real<lower=-.5,upper=.5> M[nSubj];
      real<lower=0,upper=1> SDPlus[nSubj];
      real<lower=0,upper=1> SDMinus[nSubj];
      real<lower=0,upper=100> height[nSubj];
      real<lower=0> noise;
      real predR[nSubj,nStim];
    }
    
    transformed parameters {
      matrix[nSubj,nStim] simFuncMinus;
      matrix[nSubj,nStim] simFuncPlus;
      matrix[nSubj,nStim] simFunc;
    
      for (subj in 1:nSubj) {
        for (stim in 1:nStim) {
          if (xs[stim] > M[subj])
             simFunc[subj,stim] = height[subj] * exp(1)^-(square(xs[stim]-M[subj])/(2*square(SDPlus[subj])));
          else
            simFunc[subj,stim] = height[subj] * exp(1)^-(square(xs[stim]-M[subj])/(2*square(SDMinus[subj])));
        }
      }
    }
    
    model {
    
      // likelihood
      for (subj in 1:nSubj) {
        for (stim in 1:nStim) {
        responses[subj,stim] ~ normal(simFunc[subj,stim], noise);
        predR[subj,stim] ~ normal(simFunc[subj,stim], noise);
      }
      
      // priors
      M[subj] ~ normal(0,.1)T[-.5,.5];
      SDPlus[subj] ~ gamma(1.5, 1.5)T[0,1]; // normal(0,.1)T[0,1];
      SDMinus[subj] ~ gamma(1.5, 1.5)T[0,1]; // normal(0,.1)T[0,1];
      height[subj] ~ normal(75, 10)T[1,100];
      noise ~ normal(0,5)T[0,]; // gamma(1,1);
      }
    }
    
    generated quantities {
      real log_lik[nSubj,nStim];
      for (subj in 1:nSubj) {
        for (stim in 1:nStim) {
          log_lik[subj,stim] = normal_lpdf(responses[subj,stim] | simFunc[subj,stim], noise);
        }
      }
    }
    ", file = "models/gausA.stan")
  
  stanfit <- stan(file = "models/gausA.stan",
                  data = dataList, 
                  pars = c("M", "SDPlus", "SDMinus", "height", "noise", "predR", "log_lik"),
                  iter = n_iter, 
                  warmup = n_burnin, 
                  thin = n_thin, 
                  chains = n_chains, 
                  init = "random",
                  algorithm = "NUTS",
                  cores = 1)
  
  diag <- rstan::get_sampler_params(stanfit, inc_warmup = FALSE)
  samples <- rstan::extract(stanfit)
  
  summary <- rstan::summary(stanfit, probs = c(0.025, 0.50, 0.975))$summary
  write.csv(summary, file = paste0(file_name_root, modelName, "-summary.csv"), row.names = TRUE)
  
  # calculate waic
  loglik <- loo::extract_log_lik(stanfit)
  waic <- loo::waic(loglik)
  
  # output
  out <- list(stanfit, diag, samples, summary, waic)
  names(out) <- c("stanfit", "diag", "samples", "summary", "waic")
  return(out)
}

#-------------------------------------------------------------------------------
#                     LOGISTIC AUGMENTED GAUSSIAN MODEL
#-------------------------------------------------------------------------------


Run_Aug_Gaussian_Logis_Mod <- function(dataList, modelName) {
  
  write("
    data {
      int<lower=1> nSubj;
      int<lower=1> nStim;
      int<lower=1> nTotal;
      real xs[nStim];
      int<lower=0,upper=100> subj[nTotal];
      int<lower=0,upper=100> responses[nTotal];
      int<lower=0,upper=100> stim[nTotal];
      int<lower=0,upper=100> trial[nTotal];
    }
    
    parameters {
      real<lower=-.5,upper=.5> M[nSubj];
      real<lower=0,upper=1> SDPlus[nSubj];
      real<lower=0,upper=1> SDMinus[nSubj];
      real<lower=0,upper=100> height[nSubj];
      real<lower=0> noise;
      real predR[nTotal];
      real mu[nTotal];
    }
    
    transformed parameters {
      matrix[nSubj,nStim] simFuncMinus;
      matrix[nSubj,nStim] simFuncPlus;
      matrix[nSubj,nStim] simFunc;
    
      for (s in 1:nSubj) {
        for (dim in 1:nStim) {
          if (xs[dim] > M[s])
             simFunc[s,dim] = height[s] * exp(1)^-(square(xs[dim]-M[s])/(2*square(SDPlus[s])));
          else
            simFunc[s,dim] = height[s] * exp(1)^-(square(xs[dim]-M[s])/(2*square(SDMinus[s])));
        }
      }
    }
    
    model {
    
      // likelihood
      for (i in 1:nTotal)
        mu[i] ~ normal(simFunc[subj[i],stim[i]], noise));
        responses[i] ~ bernoulli_logit(mu[i]);
        predR[i] ~ normal(simFunc[subj[i],stim[i]], noise);
      }
      
      // priors
      for (s in 1:nSubj) {
        M[s] ~ normal(0,.1)T[-.5,.5];
        SDPlus[s] ~ gamma(1.5, 1.5)T[0,1]; // normal(0,.1)T[0,1];
        SDMinus[s] ~ gamma(1.5, 1.5)T[0,1]; // normal(0,.1)T[0,1];
        height[s] ~ normal(75, 10)T[1,100];
        noise ~ normal(0,5)T[0,]; // gamma(1,1);
      }
    }
    
    generated quantities {
      real log_lik[nTotal];
      for (i in 1:nTotal) {
        log_lik[i] = bernouilli_logit_lpdf(normal_lpdf(responses[subj[i],stim[i]] | simFunc[subj[i],stim[i]], noise));
      }
    }
    ", file = "models/gausAlogis.stan")
  
  stanfit <- stan(file = "models/gausAlogis.stan",
                  data = dataList, 
                  pars = c("M", "SDPlus", "SDMinus", "height", "noise", "predR", "log_lik"),
                  iter = n_iter, 
                  warmup = n_burnin, 
                  thin = n_thin, 
                  chains = n_chains, 
                  init = "random",
                  algorithm = "NUTS",
                  cores = 1)
  
  diag <- rstan::get_sampler_params(stanfit, inc_warmup = FALSE)
  samples <- rstan::extract(stanfit)
  
  summary <- rstan::summary(stanfit, probs = c(0.025, 0.50, 0.975))$summary
  write.csv(summary, file = paste0(file_name_root, modelName, "-summary.csv"), row.names = TRUE)
  
  # calculate waic
  loglik <- loo::extract_log_lik(stanfit)
  waic <- loo::waic(loglik)
  
  # output
  out <- list(stanfit, diag, samples, summary, waic)
  names(out) <- c("stanfit", "diag", "samples", "summary", "waic")
  return(out)
}
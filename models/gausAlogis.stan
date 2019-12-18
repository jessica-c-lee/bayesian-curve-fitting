
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
      for (i in 1:nTotal) {
        responses[i] ~ bernoulli_logit(mu[i]);
        mu[i] ~ normal(simFunc[subj[i],stim[i]], noise));
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
    

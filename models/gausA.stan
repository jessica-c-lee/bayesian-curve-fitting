
    data {
      int<lower=1> nSubj;
      int<lower=1> nStim;
      int<lower=1> nTotal;
      real xs[nStim];
      real<lower=0,upper=100> subj[nTotal];
      real<lower=0,upper=100> responses[nTotal];
      real<lower=0,upper=100> stim[nTotal];
      real<lower=0,upper=100> trial[nTotal];
    }
    
    parameters {
      real<lower=-.5,upper=.5> M[nSubj];
      real<lower=0,upper=1> SDPlus[nSubj];
      real<lower=0,upper=1> SDMinus[nSubj];
      real<lower=0,upper=100> height[nSubj];
      real<lower=0> noise;
      real predR[nTotal];
      real theta[nTotal];
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
    
    model {
    
      // likelihood
      for (i in 1:nTotal) {
        responses[i] ~ bernouilli_logit(theta[i]);
        theta[i] ~ normal(simFunc[subj[i],stim[i]], noise);
        predR[i] ~ normal(simFunc[subj[i],stim[i]], noise);
      }
      
      // priors
      for (subj in 1:nSubj) {
        M[subj] ~ normal(0,.1)T[-.5,.5];
        SDPlus[subj] ~ gamma(1.5, 1.5)T[0,1]; // normal(0,.1)T[0,1];
        SDMinus[subj] ~ gamma(1.5, 1.5)T[0,1]; // normal(0,.1)T[0,1];
        height[subj] ~ normal(75, 10)T[1,100];
        noise ~ normal(0,5)T[0,]; // gamma(1,1);
      }
    }
    
    generated quantities {
      real log_lik[nTotal];
      for (i in 1:nTotal) {
        log_lik[i] = normal_lpdf(responses[subj[i],stim[i]] | simFunc[subj[i],stim[i]], noise);
      }
    }
    

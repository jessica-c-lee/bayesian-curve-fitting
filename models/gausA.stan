
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
        // simFunc[subj,:] = append_col(simFuncMinus[subj,1:6], simFuncPlus[subj,7:11]);
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
    

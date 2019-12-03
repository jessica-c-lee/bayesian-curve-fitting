
    data {
      int<lower=1> nSubj;
      int<lower=1> nStim;
      real xs[nStim];
      real<lower=0,upper=100> responses[nSubj,nStim];
    }

    parameters {
      real<lower=-.5,upper=.5> M[nSubj];
      real<lower=0,upper=1> SD[nSubj];
      real<lower=0,upper=100> height[nSubj];
      real<lower=0> noise;
      real predR[nSubj,nStim];
    }

    transformed parameters {
      real simFunc[nSubj,nStim];

      for (subj in 1:nSubj) {
        for (stim in 1:nStim) {
          simFunc[subj,stim] = height[subj] * exp(1)^-(square(xs[stim]-M[subj])/(2*square(SD[subj])));
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
        SD[subj] ~ normal(0,.1)T[0,1];
        height[subj] ~ normal(75, 10)T[1,100];
        noise ~ gamma(1,1);
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
  

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
      real<lower=1,upper=100> height[nSubj];
      real<lower=0> noise;
      real predR[nSubj,nStim];
      real<lower=-.5,upper=.5> M_group;
      real<lower=0,upper=1> SDPlus_group;
      real<lower=0,upper=1> SDMinus_group;
      real<lower=1,upper=100> height_group;
}

transformed parameters {
   matrix[nSubj,nStim] gausFunc;
    
    for (subj in 1:nSubj) {
      for (stim in 1:nStim) {
        if (xs[stim] > M[subj])
          gausFunc[subj,stim] = height[subj] * exp(1)^-(square(xs[stim]-M[subj])/(2*square(SDPlus[subj])));
        else
          gausFunc[subj,stim] = height[subj] * exp(1)^-(square(xs[stim]-M[subj])/(2*square(SDMinus[subj])));
        }
      }
}

model {
   // likelihood
      for (subj in 1:nSubj) {
        for (stim in 1:nStim) {
        responses[subj,stim] ~ normal(gausFunc[subj,stim], noise);
        predR[subj,stim] ~ normal(gausFunc[subj,stim], noise);
        M[subj] ~ normal(M_group, .5);
        SDPlus[subj] ~ normal(SDPlus_group, .5);
        SDMinus[subj] ~ normal(SDMinus_group, .5);
        height[subj] ~ normal(height_group, 10);
      }
   // group priors
      M_group ~ normal(0,.25)T[-.5,.5];
      SDPlus_group ~ gamma(1.5, 1.5)T[0,1]; 
      SDMinus_group ~ gamma(1.5, 1.5)T[0,1];
      height_group ~ normal(75, 10)T[1,100];
      noise ~ normal(0,5)T[0,];
      }
}

generated quantities {
   real log_lik[nSubj,nStim];
    for (subj in 1:nSubj) {
      for (stim in 1:nStim) {
        log_lik[subj,stim] = normal_lpdf(responses[subj,stim] | gausFunc[subj,stim], noise);
      }
    }
}

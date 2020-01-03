#-------------------------------------------------------------------------------
#                                BUILD MODEL STRING
#-------------------------------------------------------------------------------
mod_data <- "int<lower=1> nSubj;
      int<lower=1> nStim;
      real xs[nStim];
      real<lower=0,upper=100> responses[nSubj,nStim];"

mod_params <- "real<lower=-.5,upper=.5> M[nSubj];
      real<lower=0,upper=1> SDPlus[nSubj];
      real<lower=0,upper=1> SDMinus[nSubj];
      real<lower=1,upper=100> height[nSubj];
      real<lower=0> noise;
      real predR[nSubj,nStim];"

mod_group_params <- "real<lower=-.5,upper=.5> M_group;
      real<lower=0,upper=1> SDPlus_group;
      real<lower=0,upper=1> SDMinus_group;
      real<lower=1,upper=100> height_group;"

mod_tParams <- "matrix[nSubj,nStim] gausFuncMinus;
      matrix[nSubj,nStim] gausFuncPlus;
      matrix[nSubj,nStim] gausFunc;
    
      for (subj in 1:nSubj) {
        for (stim in 1:nStim) {
          if (xs[stim] > M[subj])
             gausFunc[subj,stim] = height[subj] * exp(1)^-(square(xs[stim]-M[subj])/(2*square(SDPlus[subj])));
          else
            gausFunc[subj,stim] = height[subj] * exp(1)^-(square(xs[stim]-M[subj])/(2*square(SDMinus[subj])));
        }
      }"

mod_likelihood <- "// likelihood
      for (subj in 1:nSubj) {
        for (stim in 1:nStim) {
        responses[subj,stim] ~ normal(gausFunc[subj,stim], noise);
        predR[subj,stim] ~ normal(gausFunc[subj,stim], noise);
      }"

mod_likelihood_hier <- "// likelihood
      for (subj in 1:nSubj) {
        for (stim in 1:nStim) {
        responses[subj,stim] ~ normal(gausFunc[subj,stim], noise);
        predR[subj,stim] ~ normal(gausFunc[subj,stim], noise);
        M[subj] ~ normal(M_group, .5);
        SDPlus[subj] ~ normal(SDPlus_group, .5);
        SDMinus[subj] ~ normal(SDMinus_group, .5);
        height[subj] ~ normal(height_group, 10);
      }"

mod_prior <- "// priors
      M[subj] ~ normal(0,.25)T[-.5,.5];
      SDPlus[subj] ~ gamma(1.5, 1.5)T[0,1]; 
      SDMinus[subj] ~ gamma(1.5, 1.5)T[0,1];
      height[subj] ~ normal(75, 10)T[1,100];
      noise ~ normal(0,5)T[0,];
      }"

mod_flat_prior <- "// priors 
      M[subj] ~ uniform(-.5,.5);
      SDPlus[subj] ~ uniform(0,1); 
      SDMinus[subj] ~ uniform(0,1); 
      height[subj] ~ uniform(1,100);
      noise ~ normal(0,5)T[0,];
      }"

mod_prior_hier <- "// group priors
      M_group ~ normal(0,.25)T[-.5,.5];
      SDPlus_group ~ gamma(1.5, 1.5)T[0,1]; 
      SDMinus_group ~ gamma(1.5, 1.5)T[0,1];
      height_group ~ normal(75, 10)T[1,100];
      noise ~ normal(0,5)T[0,];
      }"

mod_flat_prior_hier <- "// group priors
      M_group ~ uniform(-.5,.5);
      SDPlus_group ~ uniform(0,1); 
      SDMinus_group ~ uniform(0,1); 
      height_group ~ uniform(1,100);
      noise ~ normal(0,5)T[0,];
      }"

mod_quants <- "real log_lik[nSubj,nStim];
      for (subj in 1:nSubj) {
        for (stim in 1:nStim) {
          log_lik[subj,stim] = normal_lpdf(responses[subj,stim] | gausFunc[subj,stim], noise);
        }
      }"
#_______________________________________________________________________________

Build_Model <- function(funcs, data, params, tParams, likelihood, prior, quants) {
  # This function inserts strings in each block in the stan skeleton model
  
  modString <- paste0(read_lines("./models/skeleton.stan"), collapse = "\n")
  modString <- str_replace_all(modString, "INS_DATA", data)
  modString <- str_replace_all(modString, "INS_PARAMS", params)
  modString <- str_replace_all(modString, "INS_T_PARAMS", tParams)
  modString <- str_replace_all(modString, "INS_LIKELIHOOD", likelihood)
  modString <- str_replace_all(modString, "INS_PRIOR", prior)
  modString <- str_replace_all(modString, "INS_QUANTS", quants)
}
#_______________________________________________________________________________
# write augmented gaussian model (non-hierarchical, estimating subject params only)
# write(Build_Model(data = mod_data, 
#                   params = mod_params, 
#                   tParams = mod_tParams, 
#                   likelihood = mod_likelihood,
#                   prior = mod_prior, 
#                   quants = mod_quants), file = paste0("models/", "gausA", ".stan"))
#_______________________________________________________________________________
# write augmented gaussian model (hierarchical, estimating subject and group params)
write(Build_Model(data = mod_data, 
                  params = paste0(mod_params, "\n", mod_group_params), 
                  tParams = mod_tParams, 
                  likelihood = mod_likelihood_hier,
                  prior = mod_prior_hier, 
                  quants = mod_quants), file = paste0("models/", "gausAhier", ".stan"))

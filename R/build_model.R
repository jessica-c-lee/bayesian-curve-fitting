#-------------------------------------------------------------------------------
#                               BUILD MODEL STRING
#-------------------------------------------------------------------------------

# function that builds stan model string

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
# model strings

# data
mod_data <- "int<lower=1> nSubj;
      int<lower=1> nStim;
      real xs[nStim];
      real<lower=0,upper=100> responses[nSubj,nStim];"

# parameters for normal/augmented Gaussians
mod_params_norm <- "real<lower=-.5,upper=.5> M[nSubj];
      real<lower=0,upper=1> SD[nSubj];
      real<lower=1,upper=100> height[nSubj];
      real<lower=0> noise;
      real predR[nSubj,nStim];
      real<lower=-.5,upper=.5> M_group;
      real<lower=0,upper=1> SD_group;
      real<lower=1,upper=100> height_group;"

mod_params_aug <- "real<lower=-.5,upper=.5> M[nSubj];
      real<lower=0,upper=1> SDPlus[nSubj];
      real<lower=0,upper=1> SDMinus[nSubj];
      real<lower=1,upper=100> height[nSubj];
      real<lower=0> noise;
      real predR[nSubj,nStim];
      real<lower=-.5,upper=.5> M_group;
      real<lower=0,upper=1> SDPlus_group;
      real<lower=0,upper=1> SDMinus_group;
      real<lower=1,upper=100> height_group;"

# transformed parameters for normal/augmented Gaussians
mod_tParams_aug <- "matrix[nSubj,nStim] gausFunc;
    
    for (subj in 1:nSubj) {
      for (stim in 1:nStim) {
        if (xs[stim] > M[subj])
          gausFunc[subj,stim] = height[subj] * exp(1)^-(square(xs[stim]-M[subj])/(2*square(SDPlus[subj])));
        else
          gausFunc[subj,stim] = height[subj] * exp(1)^-(square(xs[stim]-M[subj])/(2*square(SDMinus[subj])));
        }
      }"

mod_tParams_norm <- "matrix[nSubj,nStim] gausFunc;
  
    for (subj in 1:nSubj) {
      for (stim in 1:nStim) {
        gausFunc[subj,stim] = height[subj] * exp(1)^-(square(xs[stim]-M[subj])/(2*square(SD[subj])));
      }
    }"

# likelihood function
mod_likelihood_aug <- "// likelihood
      for (subj in 1:nSubj) {
        for (stim in 1:nStim) {
        responses[subj,stim] ~ normal(gausFunc[subj,stim], noise);
        predR[subj,stim] ~ normal(gausFunc[subj,stim], noise);
        M[subj] ~ normal(M_group, .5);
        SDPlus[subj] ~ normal(SDPlus_group, .5);
        SDMinus[subj] ~ normal(SDMinus_group, .5);
        height[subj] ~ normal(height_group, 10);
      }"

mod_likelihood_norm <- "// likelihood
    for (subj in 1:nSubj) {
      for (stim in 1:nStim) {
        responses[subj,stim] ~ normal(gausFunc[subj,stim], noise);
        predR[subj,stim] ~ normal(gausFunc[subj,stim], noise);
        M[subj] ~ normal(M_group, .5);
        SD[subj] ~ normal(SD_group, .5);
        height[subj] ~ normal(height_group, 10);
    }"

# priors
mod_prior_aug <- "// group priors
      M_group ~ normal(0,.25)T[-.5,.5];
      SDPlus_group ~ gamma(1.5, 1.5)T[0,1]; 
      SDMinus_group ~ gamma(1.5, 1.5)T[0,1];
      height_group ~ normal(75, 10)T[1,100];
      noise ~ normal(0,5)T[0,];
      }"

mod_prior_norm <- "// group priors
      M_group ~ normal(0,.25)T[-.5,.5];
      SD_group ~ gamma(1.5, 1.5)T[0,1]; 
      height_group ~ normal(75, 10)T[1,100];
      noise ~ normal(0,5)T[0,];
      }"

mod_flat_prior_aug <- "// group priors
      M_group ~ uniform(-.5,.5);
      SDPlus_group ~ uniform(0,1); 
      SDMinus_group ~ uniform(0,1); 
      height_group ~ uniform(1,100);
      noise ~ normal(0,5)T[0,];
      }"

mod_flat_prior_norm <- "// group priors
      M_group ~ uniform(-.5,.5);
      SD_group ~ uniform(0,1); 
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
# build augmented gaussian (hierarchical, estimating subject and group params)

write(Build_Model(data = mod_data, 
                  params = mod_params_aug, 
                  tParams = mod_tParams_aug, 
                  likelihood = mod_likelihood_aug,
                  prior = mod_prior_aug, 
                  quants = mod_quants), file = paste0("models/", "AugGaus", ".stan"))

#_______________________________________________________________________________
# build augmented gaussian (hierarchical, estimating subject and group params)

write(Build_Model(data = mod_data, 
                  params = mod_params_norm, 
                  tParams = mod_tParams_norm, 
                  likelihood = mod_likelihood_norm,
                  prior = mod_prior_norm, 
                  quants = mod_quants), file = paste0("models/", "Gaus", ".stan"))
#_______________________________________________________________________________
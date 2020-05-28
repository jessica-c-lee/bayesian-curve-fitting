# DEMO 1
nSubj_toy <- 25

# generate gradients for two groups of subjects differing in mean, width+ and height
toy_params_1 <- list(
  list(M = rnorm(nSubj_toy, 0, .1), H = rnorm(nSubj_toy, 70, 1), 
       WM = rnorm(nSubj_toy, .2, .1) + .05, WP = rnorm(nSubj_toy, .2, .1) + .05, 
       noise = rep(2, nSubj_toy)
  ),
  list(M = rnorm(nSubj_toy, .1, .1), H = rnorm(nSubj_toy, 90, 1), 
       WM = rnorm(nSubj_toy, .2, .1) + .05, WP = rnorm(nSubj_toy, .4, .1) + .05, 
       noise = rep(2, nSubj_toy)
  )
)

toy_data_1 <- Simulate_Data(toy_params_1, nSubj = nSubj_toy, fName = "demo1")

# run analysis to recover parameters
file_name_root <- paste0("output/", "demo1", "-")
toyout1 <- Fit_Aug_Gaussian(fileName = "data/demo1.csv", 
                            modelFile = "models/AugGaus.stan", 
                            graphName = "", groupName1 = "group1", groupName2 = "group2",
                            dimVals = dim_vals, params = augG_params, groupParams = augG_group_params, 
                            ropeLow = rope_low, ropeHigh = rope_high, 
                            ropeLowDiffs = rope_low_diffs, ropeHighDiffs = rope_high_diffs,
                            hdiLim = hdi_limit, nRow = c(5,5), figMult = 2.5)

# plot simulated vs. recovered values
Plot_Param_Recovery(params = toy_params_1, samples1 = toyout1$mcmc_out_1$samples, 
                    samples2 = toyout1$mcmc_out_2$samples, nSubj = nSubj_toy, fName = "demo1")

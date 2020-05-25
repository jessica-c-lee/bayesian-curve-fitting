Fit_Norm_Gaussian <- function(fileName, modelFile,  dimVals, groupName1, groupName2, 
                              graphName, params, nRow, figMult, labels = FALSE) {
  
  # Function that reads data, runs the analysis, and saves the output
  # Note that groupName1 and groupName2 must match those in the data files
  
  # 1. read data
  out <- Read_Gen_Data(fileName, dimVals, groupName1, groupName2)
  data_list_1 <- out[[1]][[1]]
  data_list_2 <- out[[1]][[2]]
  
  # 2. fit models for each group
  mcmc_out_1 <- Run_Model(data_list_1, modelName = groupName1, modelFile, params = params)
  samples_1 <- mcmc_out_1[["samples"]]
  mcmc_out_2 <- Run_Model(data_list_2, modelName = groupName2, modelFile, params = params)
  samples_2 <- mcmc_out_2[["samples"]]
  
  # 3. plot posterior predictives
  Posterior_Preds(samples = samples_1, responses = as.vector(t(data_list_1$responses)), 
                  modelName = groupName1, nSubj = data_list_1$nSubj, 
                  nStim = data_list_1$nStim, summary = as.data.frame(mcmc_out_1$summary), 
                  nRow = nRow[1], figMult = figMult, labels = labels)
  Posterior_Preds(samples = samples_2, responses = as.vector(t(data_list_2$responses)), 
                  modelName = groupName2, nSubj = data_list_2$nSubj, 
                  nStim = data_list_2$nStim, summary = as.data.frame(mcmc_out_2$summary), 
                  nRow = nRow[2], figMult = figMult, labels = labels)
  
  # 4. write WAIC
  loo_1 <- mcmc_out_1[["loo"]]
  loo_2 <- mcmc_out_2[["loo"]]
  waic_1 <- mcmc_out_1[["waic"]]
  waic_2 <- mcmc_out_2[["waic"]]
  temp <- data.frame(group = c(groupName1, groupName2),
                     elpd_waic = c(waic_1$estimates[1,1], waic_2$estimates[1,1]),
                     p_waic = c(waic_1$estimates[2,1], waic_2$estimates[2,1]),
                     waic = c(waic_1$estimates[3,1], waic_2$estimates[3,1]))
  write_csv(temp, paste0(file_name_root, "waics.csv"))
  
  # return output list
  out <- list(data_list_1, data_list_2, mcmc_out_1, mcmc_out_2, waic_1, waic_2,
              loo_1, loo_2)
  names(out) <- c("data_list_1", "data_list_2", "mcmc_out_1", "mcmc_out_2", 
                  "waic_1", "waic_2", "loo_1", "loo_2")
  return(out)
}

# Fit normal Gaussian to NSW02 data
file_name_root <- paste0("output/", "NSW02-norm", "-")
nsw02out_normG <- Fit_Norm_Gaussian(fileName = "data/NSW02-Data.csv", modelFile = "models/Gaus.stan", 
                                    graphName = "simlin-", groupName1 = "similarity", groupName2 = "linear",
                                    params = c("M", "SD", "height", "noise", "predR", "log_lik",
                                               "M_group", "SD_group", "height_group"),
                                    dimVals = dim_vals, nRow = c(5,4), figMult = 2)

# Compare standard Gaussian to augmented Gaussian
# Group 1: similiarity subgroup
group1_comp <- Compare_Gaussians(mod1_loo = nsw02out_normG[["loo_1"]], mod2_loo = nsw02out[["loo_1"]], 
                                 mod1_waic = nsw02out_normG[["waic_1"]], mod2_waic = nsw02out[["waic_1"]])
# Group 2: linear subgroup
group2_comp <- Compare_Gaussians(mod1_loo = nsw02out_normG[["loo_2"]], mod2_loo = nsw02out[["loo_2"]], 
                                 mod1_waic = nsw02out_normG[["waic_2"]], mod2_waic = nsw02out[["waic_2"]])

print(group1_comp)
print(group2_comp)

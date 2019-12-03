#-------------------------------------------------------------------------------
#                     FITTING GENERALIZATION GRADIENTS
#-------------------------------------------------------------------------------

seed_num <- 1000
set.seed(seed_num)
# mc.cores = parallel::detectCores()
file_name_root <- "output/"

# load packages 
library(tidyverse)
library(gridExtra)
library(rstan)

# mcmc parameters
n_chains <- 4
n_iter <- 10000
n_burnin <- 1000
n_thin <- 1
n_samp <- 50
param_names <- c("Mean", "WidthMinus", "WidthPlus", "Height")

# figures
graph_file_type = ".jpeg"
fig_cols <- c("black", "red")
scat_shape <- 16
scat_size <- 1
scat_col <- alpha(fig_cols[2], .1)
gg_height <- 10
gg_width <- 10
dpi = 600
n_row <- 5
density_cols <- c("darkblue", "orange") 
dim_vals <- seq(-.5, +.5, .1)

# source
source("R/functions.R")
source("models/models.R")

# ------------------------------ GROUP ANALYSIS --------------------------------

# single cue vs. differential (Lovibond, Lee, & Hayes, 2019)
Run_Analysis(fileName = "data/NSW19-Data.csv", dimVals = dim_vals, nRow = c(10,10),
             figMult = 4, graphName = "SingDiff-", paramNames = param_names,
             groupName1 = "single", groupName2 = "differential")

# # single cue vs. differential (Lee, Hayes, & Lovibond, 2018)   
# Run_Analysis(fileName = "data/NSW02-Data.csv", dimVals = dim_vals, nRow = c(6,6),
#              figMult = 2.5, graphName = "SingDiff-", paramNames = param_names,
#              groupName1 = "single", groupName2 = "differential")

# similarity vs. linear (Lee, Hayes, & Lovibond, 2018)
Run_Analysis(fileName = "data/NSW02-Data-Diff.csv", dimVals = dim_vals, nRow = c(5,5),
             figMult = 2, graphName = "SimLin-", paramNames = param_names,
             groupName1 = "similarity", groupName2 = "linear")

# # similarity: single vs. diff (Lee, Hayes, & Lovibond, 2018)
# Run_Analysis(fileName = "data/NSW02-Data-Sim.csv", dimVals = dim_vals, nRow = c(4,4),
#              figMult = 2.5, graphName = "SingDiff-", paramNames = param_names,
#              groupName1 = "single cue", groupName2 = "differential")

# single vs. distant neg (Lee, Lovibond, Hayes, & Navarro, 2019)
Run_Analysis(fileName = "data/NSW09-Data.csv", dimVals = dim_vals, nRow = c(7,7),
             figMult = 2.5, graphName = "SingDist-", paramNames = param_names,
             groupName1 = "single pos", groupName2 = "distant neg")


# ------------------------------ SIMULATED DATA --------------------------------

# simulate <- Simulate_Data(M = rep(c(-.2, -.1, 0, .1, .2), times = 5), 
#                           SD = rep(c(.05, .1, .15, .25, .5), each = 5), 
#                           simHeight = rep(75, 25), dimVals = dim_vals,
#                           noise = 5, nSubj = 25)
# simulate$fig
# demo_data <- simulate$data
# write.csv(demo_data, file = paste0("data/","demo_data.csv"), row.names = FALSE)
# 
# # 1. read data
# demo <- Read_Demo(fileName = paste0("data/","demo_data.csv"), dimVals = dim_vals)
# data_list <- demo[[1]]
# 
# # 2. fit models 
# demo_output <- Run_GaussianA_Mod(data_list, modelName = "demo")
# demo_samples <- demo_output[["samples"]]
# 
# # 3. posterior predictives
# Gaus_Post_Preds <- Posterior_Preds(samples = demo_samples, responses = demo[[4]],
#                                    modelName = "demo", nSubj = demo[[2]], 
#                                    nStim = demo[[3]], nRow = n_row,
#                                    summary = as.data.frame(demo_output$summary))

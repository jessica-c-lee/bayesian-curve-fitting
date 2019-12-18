#-------------------------------------------------------------------------------
#                     FITTING GENERALIZATION GRADIENTS
#-------------------------------------------------------------------------------

seed_num <- 1000
set.seed(seed_num)
file_name_root <- "output/"

# load packages 
library(tidyverse)
library(gridExtra)
library(rstan)
library(sjstats)

# mcmc parameters
n_chains <- 2
n_iter <- 10000
n_burnin <- 1000
n_thin <- 1
n_samp <- 50
param_names <- c("M", "SDMinus", "SDPlus", "height")
hdi_limit <- c(.9, .95)

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
out <- Read_Trial_Gen_Data("data/LiveseyMcLaren2019.csv", dim_vals, "variable", "fixed")
data_list_1 <- out[[1]][[1]]
data_list_2 <- out[[1]][[2]]

# 2. fit models for each group
mcmc_out_1 <- Run_Aug_Gaussian_Logis_Mod(data_list_1, modelName = "variable")
samples_1 <- mcmc_out_1[["samples"]]
mcmc_out_2 <- Run_Aug_Gaussian_Logis_Mod(data_list_2, modelName = "fixed")
samples_2 <- mcmc_out_2[["samples"]]



# single cue vs. differential (Lovibond, Lee, & Hayes, 2019)
nsw19out <- Run_Analysis(fileName = "data/NSW19-Data.csv", dimVals = dim_vals, nRow = c(10,10),
                         figMult = 4, graphName = "SingDiff-", paramNames = param_names,
                         groupName1 = "single", groupName2 = "differential")

# similarity vs. linear (Lee, Hayes, & Lovibond, 2018)
nsw18out <- Run_Analysis(fileName = "data/NSW02-Data-Diff.csv", dimVals = dim_vals, nRow = c(5,5),
                         figMult = 2, graphName = "SimLin-", paramNames = param_names,
                         groupName1 = "similarity", groupName2 = "linear")

# single vs. distant neg (Lee, Lovibond, Hayes, & Navarro, 2019)
nsw09out <- Run_Analysis(fileName = "data/NSW09-Data.csv", dimVals = dim_vals, nRow = c(7,7),
                         figMult = 2.5, graphName = "SingDist-", paramNames = param_names,
                         groupName1 = "single pos", groupName2 = "distant neg")
#_______________________________________________________________________________
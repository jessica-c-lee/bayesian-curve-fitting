#-------------------------------------------------------------------------------
#                     FITTING GENERALIZATION GRADIENTS
#-------------------------------------------------------------------------------

seed_num <- 1000
set.seed(seed_num)
file_name_root <- "output/"
dir.create(file_name_root)

# load packages 
library(tidyverse)
library(gridExtra)
library(rstan)
# library(sjstats)

# mcmc parameters
n_chains <- 4
n_iter <- 10000
n_burnin <- 1000
n_thin <- 1
n_samp <- 50
param_names <- c("M", "SDMinus", "SDPlus", "height")
hdi_limit <- c(.95) 
params <- c("M", "SDPlus", "SDMinus", "height", "noise", "predR", "log_lik",
            "M_group", "SDPlus_group", "SDMinus_group", "height_group")

# figures
graph_file_type <- ".jpeg"
fig_cols <- c("black", "red")
scat_shape <- 16
scat_size <- 1
scat_col <- alpha(fig_cols[2], .1)
gg_height <- 10
gg_width <- 10
dpi <- 600
n_row <- 5
density_cols <- c("darkblue", "orange") 
dim_vals <- seq(-.5, +.5, .1)

# source
source("R/build_model.R")
source("R/functions.R")

# --------------------------------- ANALYSES -----------------------------------

# similarity vs. linear (Lee, Hayes, & Lovibond, 2018)
file_name_root <- paste0("output/", "NSW02", "-")
nsw02out <- Run_Analysis(fileName = "data/NSW02-Data-Diff.csv", 
                         modelFile = "models/gausAhier.stan", params = params,
                         dimVals = dim_vals, nRow = c(5,5), figMult = 2, 
                         graphName = "SimLin-", paramNames = param_names,
                         groupName1 = "similarity", groupName2 = "linear")

# single vs. distant neg (Lee, Lovibond, Hayes, & Navarro, 2019)
file_name_root <- paste0("output/", "NSW09", "-")
nsw09out <- Run_Analysis(fileName = "data/NSW09-Data.csv", 
                         modelFile = "models/gausAhier.stan", params = params,
                         dimVals = dim_vals, nRow = c(7,7), figMult = 2.5, 
                         graphName = "SingDist-", paramNames = param_names,
                         groupName1 = "single pos", groupName2 = "distant neg")

# single cue vs. differential (Lovibond, Lee, & Hayes, 2019)
file_name_root <- paste0("output/", "NSW19", "-")
nsw19out <- Run_Analysis(fileName = "data/NSW19-Data.csv", 
                         modelFile = "models/gausAhier.stan", params = params,
                         dimVals = dim_vals, nRow = c(10,10), figMult = 4, 
                         graphName = "SingDiff-", paramNames = param_names, 
                         groupName1 = "single", groupName2 = "differential")
#_______________________________________________________________________________
#-------------------------------------------------------------------------------
#          MODELLING GENERALIZATION GRADIENTS AS AUGMENTED GAUSSIANS
#-------------------------------------------------------------------------------

seed_num <- 1000
set.seed(seed_num)
file_name_root <- "output/"
dir.create(file_name_root)

# load packages 
library(tidyverse)
library(patchwork)
library(rstan)

# experiment-specific parameters
dim_vals <- seq(-.5, +.5, .1) # dimension values, must range between -.5 and +.5, CS+ should be at 0

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

# figure parameters
graph_file_type <- ".jpeg"
fig_cols <- c("black", "red")
scat_shape <- 16
scat_size <- 1
scat_col <- alpha(fig_cols[2], .1)
gg_height <- 10
gg_width <- 10
dpi <- 600
n_row <- 5
density_cols <- c("darkblue", "orange") # group1, group2

# build augmented Gaussian model string
source("R/build_model.R")

# source functions
source("R/functions.R")

# --------------------------- SIMULATION & RECOVERY ----------------------------
# DEMO 1
nSubj_toy <- 25
# generate gradients for two groups of subjects differing in mean, width+ and height
toy_params_1 <- list(
  list(M = rnorm(nSubj_toy, 0, .1), H = rnorm(nSubj_toy, 70, 1), 
       WM = rnorm(nSubj_toy, .2, .1), WP = rnorm(nSubj_toy, .2, .1), 
       noise = rep(2, nSubj_toy)
  ),
  list(M = rnorm(nSubj_toy, .1, .1), H = rnorm(nSubj_toy, 90, 1), 
       WM = rnorm(nSubj_toy, .2, .1), WP = rnorm(nSubj_toy, .4, .1), 
       noise = rep(2, nSubj_toy)
  )
)

toy_data_1 <- Simulate_Data(toy_params_1, nSubj = nSubj_toy, fName = "demo1")

# run analysis to recover parameters
file_name_root <- paste0("output/", "demo1", "-")
toyout1 <- Run_Analysis(fileName = "data/demo1.csv", 
                       modelFile = "models/gausAhier.stan", params = params,
                       dimVals = dim_vals, nRow = c(4,4), figMult = 2, 
                       graphName = "", paramNames = param_names,
                       groupName1 = "group1", groupName2 = "group2")

# plot simulated vs. recovered values
Plot_Param_Recovery(params = toy_params_1, samples1 = toyout1$samples_1, 
                    samples2 = toyout1$samples_2, nSubj = nSubj_toy, fName = "demo1")

# --------------------------------- ANALYSES -----------------------------------
# Fit augmented Gaussian functions to individual gradients and estimate mean,
# width-, width+, and height parameters for each subject and each group


# Re-analysis of Lee, Hayes, & Lovibond (2018, Exp 2) --------------------------
# This analysis compares gradients between a subgroup of participants who 
# reported generalising on the basis of similarity to the CS+ (similarity) 
# and a subgroup who reported generalising on he basis of the linear 
# relationship between the CS+ and CS- (linear) following intradimensional 
# training
file_name_root <- paste0("output/", "NSW02", "-")
nsw02out <- Run_Analysis(fileName = "data/NSW02-Data.csv", 
                         modelFile = "models/gausAhier.stan", params = params,
                         dimVals = dim_vals, nRow = c(5,5), figMult = 2, 
                         graphName = "SimLin-", paramNames = param_names,
                         groupName1 = "similarity", groupName2 = "linear")

# Re-analysis of Lee, Lovibond, Hayes, & Navarro (2019, Exp 1) -----------------
# This analysis compares gradients between a group given single cue training 
# (single pos group) and a group given interdimensional discrimination 
# training (distant neg group)
file_name_root <- paste0("output/", "NSW09", "-")
nsw09out <- Run_Analysis(fileName = "data/NSW09-Data.csv", 
                         modelFile = "models/gausAhier.stan", params = params,
                         dimVals = dim_vals, nRow = c(7,7), figMult = 2.5, 
                         graphName = "SingDist-", paramNames = param_names,
                         groupName1 = "single pos", groupName2 = "distant neg")

# Re-analysis of Lovibond, Lee, & Hayes (2020, Exp 1) --------------------------
# This analysis compares gradients between a group given single cue training
# (single) and a group given intradimensional discrimination training 
# (differential)
file_name_root <- paste0("output/", "NSW19", "-")
nsw19out <- Run_Analysis(fileName = "data/NSW19-Data.csv", 
                         modelFile = "models/gausAhier.stan", params = params,
                         dimVals = dim_vals, nRow = c(10,10), figMult = 4, 
                         graphName = "SingDiff-", paramNames = param_names, 
                         groupName1 = "single", groupName2 = "differential")
#_______________________________________________________________________________
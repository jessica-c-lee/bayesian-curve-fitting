#-------------------------------------------------------------------------------
#          MODELLING GENERALISATION GRADIENTS AS AUGMENTED GAUSSIANS
#-------------------------------------------------------------------------------

seed_num <- 1000
set.seed(seed_num)
file_name_root <- "output/"
dir.create(file_name_root)

# load packages 
library(tidyverse)
library(patchwork)
library(rstan)
library(bayestestR)
library(loo)

# user-defined parameters ------------------------------------------------------

# experiment-specific parameters
dim_vals <- seq(-.5, +.5, .1) # dimension values, CS+ should be at 0

# mcmc parameters
n_chains <- 4
n_iter <- 20000
n_burnin <- 1000
n_thin <- 1
n_samp <- 50
hdi_limit <- .95
rope_low <- c(-0.05, .1, .1, 70) # ROPE limits for raw parameters (M, W-, W+, H)
rope_high <- c(+0.05, .2, .2, 80)
rope_low_diffs <- c(-0.05, -.05, -.05, -5) # ROPE limits for group diffs (M, W-, W+, H)
rope_high_diffs <- abs(rope_low_diffs)
augG_params <- c("M", "SDPlus", "SDMinus", "height", "noise", "predR", "log_lik",
                 "M_group", "SDPlus_group", "SDMinus_group", "height_group")
augG_group_params <- c("M_group", "SDPlus_group", "SDMinus_group", "height_group")

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

# ------------------------------- RE-ANALYSES ----------------------------------
# Fit augmented Gaussian functions to individual gradients and estimate mean,
# width-, width+, and height parameters for each subject and group


# Re-analysis of Lee, Hayes, & Lovibond (2018, Exp 2) --------------------------
# This analysis compares gradients between a subgroup of participants who 
# reported generalising on the basis of similarity to the CS+ (similarity) 
# and a subgroup who reported generalising on he basis of the linear 
# relationship between the CS+ and CS- (linear) following intradimensional 
# training

file_name_root <- paste0("output/", "NSW02", "-")
nsw02out <- Fit_Aug_Gaussian(fileName = "data/NSW02-Data.csv", 
                             modelFile = "models/AugGaus.stan", 
                             groupName1 = "similarity", groupName2 = "linear", 
                             graphName = "simlin-", dimVals = dim_vals, 
                             params = augG_params, groupParams = augG_group_params,
                             ropeLow = rope_low, ropeHigh = rope_high, 
                             ropeLowDiffs = rope_low_diffs, ropeHighDiffs = rope_high_diffs,
                             hdiLim = hdi_limit, nRow = c(5,4), figMult = 2)

# Re-analysis of Lee, Lovibond, Hayes, & Navarro (2019, Exp 1) -----------------
# This analysis compares gradients between a group given single cue training 
# (single pos group) and a group given interdimensional discrimination 
# training (distant neg group)

file_name_root <- paste0("output/", "NSW09", "-")
nsw09out <- Fit_Aug_Gaussian(fileName = "data/NSW09-Data.csv",
                             modelFile = "models/AugGaus.stan", 
                             groupName1 = "single pos", groupName2 = "distant neg", 
                             graphName = "singdist-", dimVals = dim_vals, 
                             params = augG_params, groupParams = augG_group_params,
                             ropeLow = rope_low, ropeHigh = rope_high, 
                             ropeLowDiffs = rope_low_diffs, ropeHighDiffs = rope_high_diffs,
                             hdiLim = hdi_limit, nRow = c(7,7), figMult = 4)

# Re-analysis of Lovibond, Lee, & Hayes (2020, Exp 1) --------------------------
# This analysis compares gradients between a group given single cue training
# (single) and a group given intradimensional discrimination training 
# (differential)

file_name_root <- paste0("output/", "NSW19", "-")
nsw19out <- Fit_Aug_Gaussian(fileName = "data/NSW19-Data.csv",
                             modelFile = "models/AugGaus.stan", 
                             groupName1 = "single", groupName2 = "differential", 
                             graphName = "singdiff-", dimVals = dim_vals, 
                             params = augG_params, groupParams = augG_group_params,
                             ropeLow = rope_low, ropeHigh = rope_high, 
                             ropeLowDiffs = rope_low_diffs, ropeHighDiffs = rope_high_diffs,
                             hdiLim = hdi_limit, nRow = c(10,10), figMult = 5)

# --------------------------- SUPPLEMENTAL ANALYSES ----------------------------
# simulate gradients for parameter recovery exercise (see Supplemental Materials)
source("R/param_recovery.R")

# fit normal Gaussians to Lee et al. (2018)
source("R/norm_gaussian.R")

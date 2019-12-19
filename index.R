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
cat_dim_vals <- seq(0, .6, .1)

# source
source("models/models.R")
source("R/functions.R")

# # set this to analysis name
# analysis_name <- "Livesey2019Spat"
# file_name_root <- paste0("output/", analysis_name, "-")

# --------------------------------- ANALYSES -----------------------------------

# single cue vs. differential (Lovibond, Lee, & Hayes, 2019)
nsw19out <- Run_Analysis(fileName = "data/NSW19-Data.csv", analysisName = "NSW19",
                         dimVals = dim_vals, nRow = c(10,10), figMult = 4, 
                         graphName = "SingDiff-", paramNames = param_names, 
                         groupName1 = "single", groupName2 = "differential")

# similarity vs. linear (Lee, Hayes, & Lovibond, 2018)
nsw02out <- Run_Analysis(fileName = "data/NSW02-Data-Diff.csv", analysisName = "NSW02", 
                         dimVals = dim_vals, nRow = c(5,5), figMult = 2, 
                         graphName = "SimLin-", paramNames = param_names,
                         groupName1 = "similarity", groupName2 = "linear")

# single vs. distant neg (Lee, Lovibond, Hayes, & Navarro, 2019)
nsw09out <- Run_Analysis(fileName = "data/NSW09-Data.csv", analysisName = "NSW09", 
                         dimVals = dim_vals, nRow = c(7,7), figMult = 2.5, 
                         graphName = "SingDist-", paramNames = param_names,
                         groupName1 = "single pos", groupName2 = "distant neg")

# spatially variable vs. fixed stimuli (Livesey & McLaren 2019)
liv2019spat <- Run_Analysis(fileName = "data/Livesey2019Spat.csv", analysisName = "Liv2019Spat", 
                            dimVals = cat_dim_vals, nRow = c(7,7), figMult = 2.5, 
                            graphName = "SpatVarFixed-", paramNames = param_names,
                            groupName1 = "variable", groupName2 = "fixed")

# frequency variable vs. fixed stimuli (Livesey & McLaren 2019)
liv2019freq <- Run_Analysis(fileName = "data/Livesey2019Freq.csv", analysisName = "Liv2019Freq",
                            dimVals = cat_dim_vals, nRow = c(7,7), figMult = 2.5, 
                            graphName = "FreqVarFixed-", paramNames = param_names,
                            groupName1 = "variable", groupName2 = "fixed")
#_______________________________________________________________________________
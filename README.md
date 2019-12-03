# bayesian-curve-fitting
 R code to fit augmented Gaussians to individual gradients in a hierarchical Bayesian model
 
## index.R
Code that fits Gaussian and logistic functions to generalisation gradients in separate hierarchical Bayesian models.

##  R/functions.R
See R script for details on each function.

## models/models.R
Code for the Gaussian model and the Logistic model, implemented in stan (rstan R package).

## data/demo.csv
Demo data from an experiment with partial reinforcement (75%) differential training, and generalization testing along the blue-green dimension (see Lee, Hayes, & Lovibond, 2018).

## output
* Data.jpeg - mean generalisation gradient.
* WAICs.csv - Widely Applicable Information Criterions for the Gaussian and Logistic models.
* BFs.csv - Bayes Factors comparing Gaussian and Logistic models calculated using the "bridgesampling" R package.
* Gaus_Summary.csv - descriptive statistics for the posteriors from the Gaussian model.
* Logis_Summary.csv - descriptive statistics for the posteriors from the Logistic model.
* Guassian-PostPreds.jpeg - posterior predictives (100 samples from the posterior) for each participant for the Gaussian model.
* Logistic-PostPreds.jpeg - posterior predictives (100 samples from the posterior) for each participant for the Logistic model.

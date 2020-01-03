# bayesian-curve-fitting
R code to fit augmented Gaussians to individual generalization gradients in a hierarchical Bayesian model. The model estimates the posteriors for 4 parameters of the augmented Gaussian:
* mean: the location of the gradient peak
* width-: the width (SD) of the left side of the gradient
* width+: the width (SD) of the right side of the gradient
* height: the height (value at the peak) of the gradient

Accompanies the paper Lee, Mills, & Livesey (under review). Provides a re-analysis of results from 3 studies:
* Experiment 1 of Lovibond, Lee, & Hayes (2019): single cue vs. differential training groups
* Experiment 2 of Lee, Hayes, & Lovibond (2018): similarity vs. linear rule subgroups
* Experiment 2 of Lee, Lovibond, Hayes, & Navarro (2019): single pos vs. distant neg training groups.
 
## How to use
Run the index.R file. Note that the code is set up to compare gradients between two groups.

## Output generated
* summary.csv: 
* gradients.jpeg: plots the empirical gradients facetted by subject
* postpreds.jpeg: plots the empirical gradients facetted by subject with posterior predictives overlayed
* HDIs.csv: 95% Highest Density Intervals (HDIs) for each parameter for each group
* group_diff_HDIs.csv: 95% HDIs for the group difference for each parameter

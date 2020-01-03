# bayesian-curve-fitting

R code to fit augmented Gaussians to individual generalization gradients in a hierarchical Bayesian model. The augmented Gaussian has 4 parameters that allow it to fit asymmetrical gradients: 
* mean: the location of the gradient peak
* width-: the width (SD) of the left side of the gradient
* width+: the width (SD) of the right side of the gradient
* height: the height (peak) of the gradient

Accompanies the paper Lee, Mills, & Livesey (under review). Provides a re-analysis of results from 3 studies:
* Experiment 1 of Lovibond, Lee, & Hayes (2019): single cue vs. differential training groups
* Experiment 2 of Lee, Hayes, & Lovibond (2018): similarity vs. linear rule subgroups
* Experiment 2 of Lee, Lovibond, Hayes, & Navarro (2019): single pos vs. distant neg training groups.

References:
Lee, J. C., Hayes, B. K., & Lovibond, P. F. (2018). Peak shift and rules in human generalization. *Journal of Experimental Psychology: Learning, Memory, and Cognition, 44*, 1955-1970.
Lee, J. C., Lovibond, P. F., Hayes, B. K., & Navarro, D. J. (2019). Negative evidence and inductive reasoning in generalization of associative learning. *Journal of Experimental Psychology: General, 148*, 289-303.
Lovibond, P. F., Lee, J. C., & Hayes, B. K. (2019). Stimulus discriminability and induction as independent components of generalization. *Journal of Experimental Psychology: Learning, Memory, and Cognition.*

## Author note
Still in development! For any issues please email jessica.lee@unsw.edu.au.

## How to use
Run the index.R file. Note that the code is set up to compare gradients between two groups.

## Output generated
* summary.csv: summary statistics of the posterior samples
* gradients.jpeg: plots the empirical gradients facetted by subject
* postpreds.jpeg: plots the empirical gradients facetted by subject with posterior predictives overlayed
* HDIs.csv: 95% Highest Density Intervals (HDIs) for each parameter for each group
* group_diff_HDIs.csv: 95% HDIs for the group difference for each parameter

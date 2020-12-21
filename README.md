# bayesian-curve-fitting

Accompanies the paper:
Lee, J. C., Mills, L., Hayes, B. K., & Livesey, E. J. (2021). Modelling generalisation gradients as augmented Gaussian functions. Quarterly Journal of Experimental Psychology, 74(1), 106-121. https://doi.org/10.1177/1747021820949470 

R code to fit augmented Gaussians to individual generalization gradients in a hierarchical Bayesian model. Uses the [rstan](https://mc-stan.org/users/interfaces/rstan) and [bayestestR](https://github.com/easystats/bayestestR) packages.

The augmented Gaussian has 4 parameters that allow it to fit asymmetrical gradients:
* mean: the location of the gradient peak
* width-: the width (SD) of the left side of the gradient
* width+: the width (SD) of the right side of the gradient
* height: the height (peak) of the gradient

Provides a re-analysis of results from 3 studies:
* Experiment 1 of Lovibond, Lee, & Hayes (2019): single cue vs. differential training groups
* Experiment 2 of Lee, Hayes, & Lovibond (2018): similarity vs. linear rule subgroups
* Experiment 2 of Lee, Lovibond, Hayes, & Navarro (2019): single pos vs. distant neg training groups.

## References
Lee, J. C., Hayes, B. K., & Lovibond, P. F. (2018). Peak shift and rules in human generalization. *Journal of Experimental Psychology: Learning, Memory, and Cognition, 44(12)*, 1955-1970.

Lee, J. C., Lovibond, P. F., Hayes, B. K., & Navarro, D. J. (2019). Negative evidence and inductive reasoning in generalization of associative learning. *Journal of Experimental Psychology: General, 148(2)*, 289-303.

Lovibond, P. F., Lee, J. C., & Hayes, B. K. (2020). Stimulus discriminability and induction as independent components of generalization. *Journal of Experimental Psychology: Learning, Memory, and Cognition, 46(6)*, 1106-1120.

## How to use
* Open the .Rproj file and run the index.R file. Note that the code is set up to compare gradients between two groups 
* The data file must be in long format, with the grouping variable labelled as "group", the stimulus dimension variable labelled as "x" and the dependent variable labelled as "y" (see demo1.csv)

## Note
* The code is set up to fit gradients across 11 test stimuli, with the CS+ at the midpoint of the dimension (coded as 0). All functions use the fixed stimulus values -0.5:0.1:+0.5. These values are arbitrary, but the model specification is dependent on these values.
* The code is set up to model responses ranging from 0-100. The model must be re-specified to work with a different range.

## Output generated
* summary.csv: summary statistics for the posterior samples for each parameter
* gradients.jpeg: plot of the empirical gradients facetted by subject
* postpreds.jpeg: plot of the empirical gradients facetted by subject with posterior predictives overlayed
* density.jpeg: 5 panelled figure with the mean generalisation gradients (a) and posterior density plots of the 4 group-level augmented Gaussian parameters (b-e)
* groupdiff-density.jpeg: posterior density plots for the group difference for each of the 4 augmented Gaussian parameters
* waics.csv: Widely Applicable Information Criterions for both groups, computed with the [loo package](https://cran.r-project.org/web/packages/loo/index.html)
* HDIs.csv
  - HDI lim: Highest Density Interval calculated (%)
  - HDI low: Highest Density Interval lower limit
  - HDI high: Highest Density Interval upper limit
  - p(direction): proportion of the posterior that is positive or negative (whatever is most probable). Note that this is meaningless for the width and height parameters since they can only be positive
* group_diff_HDIs.csv: 
  - HDI lim: Highest Density Interval calculated (%)
  - HDI low: Highest Density Interval lower limit
  - HDI high: Highest Density Interval upper limit
  - p(direction): proportion of the posterior that is positive or negative (whatever is most probable)
  - ROPE low: Region of Practical Equivalence lower limit 
  - ROPE high: Region of Practical Equivalence upper limit 
  - p(ROPE): proportion of the posterior that lies within the ROPE

## Contact
Contact jessica.lee@unsw.edu.au. 

# gsusie-package

The `gsusie` package implements the Generalized Sum of Single-Effects (GSuSiE) model. 
Developed from the Sum of Single-Effects (SuSiE) model [Wang et al. (2020)](https://doi.org/10.1111/rssb.12388) 
that performs variable selection in linear regression models, 
the GSuSiE model expands the usage by allowing to select variables 
in generalized linear models (GLMs) such as logistic and Poisson regressions. 
In other words, the `gsusie` function can model binary and count type responses. 
Likewise, this method is particularly compatible for high-dimensional 
settings where a high correlation exists among explanatory variables $\mathbf{X}$, 
and the number of true effect variables are much smaller than the total number 
of explanatory variables (i.e., sparse variable selection). 
Yet, it can also be effective in more scenarios. 

This method provides summaries such as posterior inclusion probabilities (PIPs) 
and credible sets to make inferences on which variables should be selected. 
Additionally, robust estimations are provided to avoid the potential effect of 
outliers on model fitting. The adapted iterative Bayesian Stepwise Selection 
(IBSS) algorithm is combined with the iterative reweighted least squared (IRLS) 
approach to address GLMs and robust estimation. 



## Quick Start

To install this R package directly from Github: 
```
# install.packages("devtools")
devtools::install_github("ymyuan98/gsusie")
```


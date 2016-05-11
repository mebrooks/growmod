# growmod
###An R package for estimating growth using mark-recapture data of multiple individuals
This package fits growth models to mark-recapture data of individual sizes. Growth from one time to the next is modeled as a Markov process. Covariates such as environmental variables or individual variables can affect both the intercept and autoregressive parts of the model. Models can include random effects of individual, time, and birth cohort. Time must be in discrete steps. 
Observations of individual sizes can be missing, but covariates cannot be missing. Individual could represent an individual organism, population, or other unit of replication.

Estimation is done in a state-space framework using maximum marginal likelihood estimation via [Template Model Builder](https://github.com/kaskr/adcomp).

This package is part of a manuscript in preparation for publication.

## Installation
To install `growmod`, first install `TMB` from CRAN, then install `growmod` via
```
devtools::install_github("mebrooks/growmod")
```
(you may need to install the `devtools` package, compiler tools, and other dependencies first).
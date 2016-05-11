# growmod
##An R package for estimating growth using mark-recapture data of multiple individuals
This package is useful for fitting growth models to mark-recapture data of individual sizes. Growth is modeled as a Markov process. Covariates such as environmental variables or individual variables can affect both the intercept and autoregressive parts of the model. Models can include random effects of individual, time, and birth cohort. Time must be in discrete steps. 
Observations of individual sizes can be missing, but covariates cannot be missing. Individual could represent an individual organism, population, or other unit of replication.

Models are fit using maximum marginal likelihood estimation via the Template Model Builder (TMB) package https://github.com/kaskr/adcomp

This package is part of a manuscript in preparation for publication.

To install `growmod`, first install `TMB` from CRAN, then install `growmod` via
```
devtools::install_github("mebrooks/growmod")
```
(you may need to install the devtools package, compiler tools, and other dependencies first).
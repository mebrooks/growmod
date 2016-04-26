# growmod
##An R package for fitting state-space models to repeated measures of multiple individuals
Covariates such as environmental variables or individaul variables can affect both the intercept and autoregressive parts of the model. Models can include random effects of individual, time, and birth cohort. Time must be in discrete steps. 
Observations of indivdiaul sizes can be missing, but covariates cannot be missing. Individual could represent an individaul organism or population.

Models are fit using maximum marginal likelihood estimation via the Template Model Builder (TMB) package https://github.com/kaskr/adcomp


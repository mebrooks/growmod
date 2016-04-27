library(growmod)

true=c("(Intercept)"=1, age=0,  "I(age^2)"=0, size=.6, "size:age"=0, sigma_proc=.01, time_growth_sd=.1, indiv_growth_sd=.1, indiv_age_growth_sd=0, indiv_cor=0, size0_mu=2.5, size0_sd=1)

sg=simobs(true, nind=5000, ntime=40)

m1=growmod(~1, ~1, data=sg, REcohort=FALSE)
summary(m1)
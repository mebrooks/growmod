library(growmod)
library(plyr)
library(lme4)
true=c("(Intercept)"=1, age=.05,  "I(age^2)"=0, prevsize=.6, "prevsize:age"=-.05, sigma_proc=.01, time_growth_sd=.1, indiv_growth_sd=.1, indiv_age_growth_sd=0, indiv_cor=0, size0_mu=2.5, size0_sd=1)
set.seed(111)
sim=simobs(true, nind=1000, ntime=20)
m1=growmod(~age, ~age, data=sim)

sim2=ddply(organize_data(sim, sigma_obs = .01)$Ldat$Predictors, ~ID, mutate,
	prevsize= size[match(age, age+1)]
)

m2=lmer(size~age+prevsize+prevsize:age  +(1|ID)+(1|t), sim2)
summary(m1)
summary(m2)

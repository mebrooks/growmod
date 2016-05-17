##' Organize data so that it can be read into the TMB model
##'
##' @param data a data frame of observed sizes and predictors in long format. See details.
##' @param sigma_obs the standard deviation of the observation error with the same units as size (e.g. log kg). To be used in a \code{growmod} model when \code{estobserr=FALSE}.
##' @return a list of objects in the correct format to be read in by the TMB model defined in \code{growmod.cpp} with the exception of beta, eta, M, and X which depend on the model formulas and are defined in the growmod function.
##' @details data must contain at least the following columns: \code{size, t, t0, and ID}. \code{size} is the size of an individual on the log scale. \code{t} is the point in time when an observation was made. \code{t} must be in discrete units (i.e. integers) with increments of 1. \code{t0} is the time when an individual was born. \code{ID} is a unique identifier for each individaul. It must contain at least one row per combination of individual and time. Multiple rows may be needed if an individual was measured more than once at the same timepoint. These multiple measurements are very informative for the estimation of observation error. If multiple rows for the same combination of individual and time are given, then only the covariates from the first instance is used. If a size observation was missed for one combination of individual and time, then that row should have an \code{NA} in the size column, but must contian any covariates. There must be at least one non-missing size for each individual.
##' @export
##' @importFrom reshape cast
##' @importFrom plyr ddply
organize_data=function(data, sigma_obs)
{
if(!all(c("size", "t", "t0", "ID") %in% names(data))){stop("data must contain size, t, t0, and ID.")}

#organize data
data$fID=factor(data$ID)
data=transform(data, it=t-min(t),
					iID=as.numeric(fID)-min(as.numeric(fID))
			)
obs=na.omit(data[,c("size", "iID", "it")])
nind=length(unique(obs$iID))
ntimes=length(unique(data$it))
mintime=min(data$t)
if(any(is.na(data$t0))){stop("Each individual must have a birth time (t0).")}

indivdat=ddply(data, ~iID, summarize, t0=t0[1], T=max(t), ID=ID[1],
					bad=any(!(t0[1]:max(t) %in% t)) | !any(!is.na(size))
					)
if(any(indivdat$bad)){stop("Each individual must have a row of data for every timepoint from t0 (for that individual) to the last timepoint for that individual, even when size observations are missing. Those rows can have NA in the size column, but must contain any predictors needed for formulaX and formulaM. Also, each individual must have at least one non-missing size. \n Check the following individuals: ", subset(indivdat, bad==TRUE)$ID)}

#take care of repeated measures at the same time point
noreps=ddply(data, ~iID+it, summarize, size=mean(size))
preddat=join(noreps, data[,-grep("size", names(data))], match="first", by=c('iID', 'it')) 
noreps$size[is.na(noreps$size)]=mean(noreps$size, na.rm=TRUE)
tmp= preddat[, c('it', 'iID')]#no rep sizes
tmp$value=0:(nrow(tmp)-1)
tmp2=cast(tmp, iID~it, fill= NAint)[,-1]
lookup =as.matrix(tmp2)

Ldat=list(obs=obs$size , #logsize observed
		indiv_obs=obs$iID, #indiv corresponding to observed size
		time_obs=obs$it, #time corresponding to observed size
		t0= indivdat$t0-mintime,
		T= indivdat$T-mintime,
		indiv_size= preddat $ iID,
		time_size= preddat $it,
		age_size=preddat$t-preddat$t0,
		lookup= lookup,
		sigma_obs=sigma_obs,
		Predictors= preddat,
		NAint=NAint
	)
Lpin=list(size= noreps$size, #rep(0, length(Minit)), 
		log_time_growth_sd=0,
		log_indiv_growth_sd=0,
		log_cohort_growth_sd=0,
		log_resid_growth_sd=0,
		log_size0_sd=0,
		size0_mu=mean(subset(data, t==t0)$size, na.rm=TRUE),
		indiv_growth_dev=rep(0, nind),
		time_growth_dev=rep(0, ntimes),
		log_sigma_proc=0,
		cohort_growth_dev=rep(0, ntimes),
		scale_indiv_cor=0,
		log_indiv_age_growth_sd=0,
		indiv_age_growth_dev=rep(0, nind)
		)
	return(list(Ldat= Ldat, Lpin= Lpin))
}
###########################################
##' A function for fitting an autoregressive model including fixed effect covariates and random effects. 
##'
##' @param formulaX a formula for the covariates that affect the intercept of the AR(1) model
##' @param formulaM a formula for the covariates that affect the autoregressive coefficient of the AR(1) model
##' @param data a data frame of observed sizes and covariates in long format. See details.
##' @param estobserr logical - should observation error 
##' @param sigma_obs the standard deviation of the observation error with the same units as size (e.g. log kg). Ignored if \code{estobserr=TRUE}.
##' @param predfirstsize NULL if there are no predictors on first size, or a design matrix of predictors with one row per individual in the order as the individuals first apear in data
##' @param DLL the name of the compiled TMB/c++ file
##' @param silent logical - Disable all tracing information in maximum likelihood estimation?
##' @param selecting logical - when true, save time by not estimating confidence intervals. Only return AIC and convergence.
##' @param REtime logical - estimate a random intercept for each time point
##' @param REID logical - estimate a random intercept for each individual
##' @param REcohort logical - estimate a random intercept for each cohort with the same \code{t0} (i.e. birth cohort)
##' @param REageID logical - estimate a random slope on age for each individual (only use in models that include age in formulaX)
##'
##' @return a list of objects in the correct format to be read in by the TMB model defined in \code{growmod.cpp} with the exception of \code{beta}, \code{eta}, \code{M}, and \code{X} which depend on the model formulas and are defined in the growmod function.
##' @details \code{data} must at least the following columns: \code{size, t, t0}, and \code{ID}. 
##' \code{size} is the size of an individual on the log scale. 
##' \code{t} is the point in time when an observation was made. 
##' \code{t} must be in discrete units (i.e. integers) with increments of 1. 
##' \code{t0} is the time when an individual was born. 
##' \code{ID} is a unique identifier for each individaul. 
##' \code{data} must contain at least one row per combination of individual and time. Multiple rows may be needed if an individual was measured more than once at the same timepoint. These multiple measurements inform the estimation of observation error. If multiple rows for the same combination of individual and time are given, then only the predictors from the first instance is used. If a size observation was missed for one combination of individual and time, then that row should have an \code{NA} in the size column, but must contian any covariates. There must be at least one non-missing size for each individual.
##' @export
##' @import TMB
##' @useDynLib growmod
##' @examples
##' true=c("(Intercept)"=1, age=0,  "I(age^2)"=0, size=.6, "size:age"=0, sigma_proc=.01, time_growth_sd=.1, indiv_growth_sd=.1, indiv_age_growth_sd=0, indiv_cor=0, size0_mu=2.5, size0_sd=1)
##' sg=simobs(true, nind=1000, ntime=20)
##' m1=growmod(~1, ~1, data=sg, REcohort=FALSE)
##' summary(m1)
##' 

growmod=function(formulaX=~1, formulaM=~1, data, estobserr=TRUE, sigma_obs = NULL, predfirstsize =NULL, DLL="growmod", silent=TRUE, selecting=FALSE, REtime=TRUE, REID=TRUE, REcohort=FALSE, REageID=FALSE,...)
{
	ogd=organize_data(data, sigma_obs)
	Ldat= ogd$Ldat
	Lpin= ogd$Lpin
	
	#Design matricies 
	Ldat$X=model.matrix(formulaX, data=Ldat$Predictors)
	Lpin$beta=rep(0, ncol(Ldat$X))
	
	Ldat$M=model.matrix(formulaM, data=Ldat$Predictors)
	Lpin$eta=c(.1, rep(0, ncol(Ldat$M)-1))
	##################################################
	#Which random effects are not being used
	map=list()
	if(!REcohort)
	{
		map$log_cohort_growth_sd=as.factor(NA) 
		map$cohort_growth_dev=as.factor(rep(NA, length(Lpin$cohort_growth_dev)))
	}
	if(!REtime)
	{
		map$log_time_growth_sd=as.factor(NA) 
		map$time_growth_dev=as.factor(rep(NA, length(Lpin$time_growth_dev)))
	}
	if(!REID)
	{
		if(REageID){stop("Models with a random slope of ID, but no intercept is not implemented.")}
		map$log_indiv_growth_sd=as.factor(NA)
		map$indiv_growth_dev=as.factor(rep(NA, length(Lpin$indiv_growth_dev)))
	}
	if(!REageID)
	{
		map$scale_indiv_cor=as.factor(NA)
		map$log_indiv_age_growth_sd=as.factor(NA) 
		map$indiv_age_growth_dev=as.factor(rep(NA, length(Lpin$indiv_age_growth_dev)))
	}
	random=c("indiv_growth_dev", "indiv_age_growth_dev", "time_growth_dev","cohort_growth_dev", "size")[c(REID, REageID, REtime, REcohort, TRUE)]
	##################################################
	#What to do about observation error
	if(!estobserr)
	{
		map$log_sigma_obs=as.factor(NA)
		if(is.null(Ldat$sigma_obs)){stop("If you do not want growmod to estimate observation error, then you must give a value for sigma_obs.")}
		if(is.null(Lpin$log_sigma_obs)){Lpin$log_sigma_obs=log(mean(Ldat$obs)/20)}#start at 5% error
	}
	if(estobserr)
	{
		Lpin$log_sigma_obs=0		
		Ldat$sigma_obs=0 #just in case it is NULL	
	}
	##################################################
	#What to do about first size
	if(is.null(predfirstsize))
	{
		Ldat$B=matrix(1, nrow=length(Ldat$t0), ncol=1)#just using t0 to get length right
		colnames(Ldat$B)="size0_mu"
		Lpin$theta=1
	}	
	if(!is.null(predfirstsize))
	{
		Ldat$B=predfirstsize
		Lpin$theta=rep(0, ncol(predfirstsize))
	}
	##################################################
	#TMB fit
	obj=MakeADFun(data=Ldat, parameters=Lpin, random=random, DLL=DLL, silent=silent, map=map,...)
	fit = nlminb(obj $par, obj $fn, obj $gr, control=list(eval.max=10000, iter.max=10000))
	if(selecting){
		mod=c(AIC=TMBAIC(fit, n=length(Ldat$obs)), convergence=fit$convergence)
		class(mod)="selecting_growmod"
		return(mod)
	}
	#else
	sdr=sdreport(obj)
	if(!any(is.na(summary(sdr, "fixed")))){
		mod=list(parList=obj$env$parList(), 
				fit=fit, 
				sdr=sdr, 
				Xnames=colnames(Ldat$X), Mnames=colnames(Ldat$M), Bnames=colnames(Ldat$B),
				formulaX= formulaX, formulaM= formulaM, 
				AIC=TMBAIC(fit, n=length(Ldat$obs)),
				nobs=length(Ldat$obs),
				nmis=length(Lpin$size)-length(Ldat$obs),
				recap=length(Ldat$obs)/length(Lpin$size)
				)
		class(mod)="growmod"
		return(mod)
	}
	#else
	warning("some values of sdreport were NA")
	mod=list(parList=obj$env$parList(), 
				fit=fit, 
				sdr=sdr, 
				Xnames=colnames(Ldat$X), Mnames=colnames(Ldat$M), Bnames=colnames(Ldat$B),
				formulaX= formulaX, formulaM= formulaM, 
				AIC=NA,
				nobs=length(Ldat$obs),
				nmis=length(Lpin$size)-length(Ldat$obs),
				recap=length(Ldat$obs)/length(Lpin$size) 
				)
	class(mod)="growmod"
	return(mod)
}
###########################################
##' Calculate the Akaike Information Criteria corrected for small sample size
##' 
##' @param opt an optimized TMB model
##' @param n number of observations
##' @param correction logical - Correct for small sample size? (recommended as it converges to other formula for large samples)
##' @export
TMBAIC=function(opt, n=NULL, correction=TRUE)
{
	if(correction & is.null(n)){stop("Sample size n needed to calculate AICc.")}
	l=opt[["objective"]]
	k=length(opt[["par"]])
	if(correction){return(2*k - 2*l + 2*k*(k+1)/(n-k-1))}
	return(2*k - 2*l)
}
###########################################
##' Extract coefficients from a fitted model and give the coefficients more informative names based on the specified formulas.
##'
##' @param mod an object of type \code{growmod} that was fit using the \code{growmod()} function. See \code{?growmod} for details.
##' @param CV logical - Should estimates of the coefficients of variance (variance/intercept) be reported?
##' @param size0 logical - Should estimates of the distribution of initial sizes be reported? 
##' @export
extract_coefs=function(mod, CV=FALSE, size0=FALSE)
{
#takes  results of growmod() and returns a data frame of est, sd, z.value, and names of coefficients
	cf=summary(mod$sdr, "fixed")
	nf=rownames(cf)
	take0=which(nf=="beta")
	take0=c(take0, which(nf=="eta"))
	if(size0) {take0=c(take0, which(nf=="theta"))}
	nf[nf=="beta"]=mod$Xnames
	nf[nf=="eta"]=paste0("size:",mod$Mnames)
	nf[nf=="theta"]=paste0("size0 model ",mod$Bnames)

	nf[nf=="size:(Intercept)"]="size"
	rownames(cf)=nf
	cf2=cf[take0, ]

	cr=summary(mod$sdr, "report")
	nr=rownames(cr)
	take=grep("growth_sd",nr)
	take=c(take, grep("sigma_obs",nr))
	take=c(take, grep("sigma_proc",nr))
	take=c(take, grep("indiv_cor",nr))
	if(CV) { take = c(take, grep("CV",nr))}
	if(size0) { take = c(take, grep("size0_sd",nr))}
	cr2=cr[take, ]
	c2=as.data.frame(rbind(cf2,cr2))
  	c2$var=rownames(c2)
  	colnames(c2)[1]="est"
	colnames(c2)[2]="se"
	c3=transform(c2, "t value"=est/se)
	return(c3)
}	
###########################################
##' Output estimates and summary statistics from a sucessfully fitted model
##'
##' return summary statistics for a model fitted by the \code{growmod()} function.
##' @param mod an object of type \code{growmod} that was fit using the \code{growmod()} function. See \code{?growmod} for details.
##' @export
summary.growmod=function(mod)
{
	cat("State-space growth model fit by maximum marginal likelihood estimation \n")
	cat("Formula X: ", as.character(mod$formulaX), "\n")
	cat("Formula M: ", as.character(mod$formulaM), "\n\n")
	cat("Convergence: ", as.character(mod$fit$message), "\n")
	cat("AIC: ", as.character(mod$AIC), "\n\n")
	cat("number of observations: ", as.character(mod$nobs), "\n")
	cat("number of missing observations: ", as.character(mod$nmis), "\n")
	cat("recapture rate: ", as.character(round(mod$recap, 3)), "\n")
	cat("number of individuals: ", length(mod$parList$indiv_growth_dev), "\n")
	cat("number of time points: ", length(mod$parList$time_growth_dev), "\n")
	cat("Fixed and Random effects: \n")
	print(extract_coefs(mod)[,-3])
}	
###########################################
##' Extract estimates of the latent variable
##'
##' @param mod an object of type \code{growmod} that was fit to \code{data} using the \code{growmod()} function. See \code{?growmod} for details.
##' @param data the same data set that was used for fitting \code{mod}
##' @return predictions from a growmod model on the same scale as the sizes provided (e.g. log kg)
##' @export
extract_pred=function(mod, data){
	od=organize_data(data, 0)
	rr=summary(mod$sdr, "random")
	nr=rownames(rr)
	size=rr[nr=="size",]
	pred=data.frame(est=size[,1], se=size[,2], age=od$Ldat$age, ID=d$Ldat$Predictors$ID, lamb=od$Ldat$Predictors$lamb, wnao=od$Ldat$Predictors$wnao, pop=od$Ldat$Predictors$prevpop)
	return(pred)
}	
###########################################
##' Extract estimates of individual deviates from the average pattern.
##'
##' @param mod an object of type \code{growmod} that was fit to \code{data} using the \code{growmod()} function. See \code{?growmod} for details.
##' @param data the same data set that was used for fitting \code{mod}
##' @return a data frame with the estimated individual deviations from the average intercept added to the original observations and covariates. The data fram has repeated observations of the same individual at the same point averaged out.
##' @export
extract_ID_dev=function(mod, data){
	od=organize_data(data, 0)
	Ldat=od$Ldat
	i=mod$parList$indiv_growth_dev
	dat=Ldat$Predictors
	dat$dev=i[Ldat$Predictors$iID+1]
	return(dat)
}	
###########################################
##' Make predictions of size at the next point in time.
##'
##' @param mod an object of type \code{growmod} that was fit to \code{data} using the \code{growmod()} function. See \code{?growmod} for details.
##' @param newdata a data frame containing a column for each of the covariates in the formulas usind for fitting \code{mod} and a column \code{size} which is the size in the previous time point.
##' @param exp logical - Should the predictions be exponentiated to put them back on a natural scale? See details.
##' @return predicted sizes
##' @details It is assumed that size is measured on the log scale (e.g. log kg). Therefore, when \code{exp=FALSE}, predictions are made on the log scale. When \code{exp=TRUE}, predictions are made using the formula \code{exp(m + 0.5*v)} where \code{m} is the prediction on the log scale, and \code{v} is the sum of estimated variances including random effects and process error.
##' @export
predict.growmod=function(mod, newdata, exp=FALSE){

	X=model.matrix(mod$formulaX, newdata)
	M=model.matrix(mod$formulaM, newdata)
	if(!exp){return(X%*% mod$parList$beta + M%*% mod$parList$eta*newdata$size)}
	
	r=summary(mod$sdr, "report")
	if("indiv_age_growth_sd" %in%rownames(r))
		stop("pred_expsize is only implemented for models without random slopes")
	sumsigmasq=sum(r[grep("growth_sd",rownames(r)) , "Estimate"]^2)+r["sigma_proc","Estimate"]^2
	return(exp(X%*% mod$parList$beta + M%*% mod$parList$eta*newdata$size + 0.5*sumsigmasq))

}	
###########################################
##' Simulate growth
##'
##' @details Covariates other than age are not implemented. Birth times are uniformly distributed across the time series. Observations of an individual stop either when it reaches the end of its lifespan, or when time reaches \code{ntime}. 
##' @param pars a numeric vector contianing the following named terms representing regression coefficients: \code{'(Intercept)', 'age',  'I(age^2)', 'size', 'size:age', 'sigma_proc', 'time_growth_sd', 'indiv_growth_sd', 'indiv_age_growth_sd', 'indiv_cor', 'size0_mu', 'size0_sd'} 
##' @param nind the number of individuals to simulate
##' @param ntime the number of times spanned bt the dataset
##' @param maxage the maximum age that all individuals reach (not used if lifespans is specified)
##' @param lifespans a set of intergers to be resampled from. It should reflect the frequency of each lifespan.
##' @export
##' @examples
##' true=c("(Intercept)"=1, age=0,  "I(age^2)"=0, size=.6, "size:age"=0, sigma_proc=.01, time_growth_sd=.1, indiv_growth_sd=.1, indiv_age_growth_sd=0, indiv_cor=0, size0_mu=2.5, size0_sd=1)
##' sg=simgrow(true, nind=10, ntime=20)
simgrow=function(pars, nind=150, ntime=29, maxage=12, lifespans=NULL){

	time_growth_dev=rnorm(ntime, sd= pars['time_growth_sd'])#time random effects

	t0=sample(x=ntime, size=nind, replace=TRUE)#birth time for each indiv
	if(is.null(lifespans)){	T=pmin(t0+ maxage, ntime)
		}else{
		T=pmin(t0+ sample(x=lifespans, size=nind, replace=TRUE), ntime)
	}
	if(is.na(pars['indiv_age_growth_sd'])){
		pars['indiv_age_growth_sd']=0
		pars['indiv_cor']=0
	}
	Sigma=matrix(data=c(pars['indiv_growth_sd']^2, 
				pars['indiv_growth_sd']*pars['indiv_age_growth_sd']*pars['indiv_cor'],
				pars['indiv_growth_sd']*pars['indiv_age_growth_sd']*pars['indiv_cor'],
				pars['indiv_age_growth_sd']^2),
			nrow=2, ncol=2)
	idev=MASS::mvrnorm(n=nind, mu=c(0,0), Sigma= Sigma)

	dat=data.frame(NULL)
	count=1
	for(i in 1:nind)
	{
		b=t0[i]
		d=T[i]
		temp=data.frame(size=rep(NA, d-b+1), t=b:d, ID=i, age=0:(d-b))
		temp$size[1]=rnorm(1, pars['size0_mu'], pars['size0_sd'])#first size
		temp$t0=b
		temp$T=d
		if(d>b)
		{
			for(y in b:(d-1))#process model
			{
				temp$size[y-b+2] = rnorm(1, pars['(Intercept)']+idev[i,1] +
									(pars['age']+idev[i,2])*temp$age[y-b+1] +
									pars['I(age^2)']*temp$age[y-b+1]^2 + 
									(pars['size'] + pars['size:age']*temp$age[y-b+1])*temp$size[y-b+1] +
									time_growth_dev[y], 
									pars['sigma_proc'])		
				count=count+1
			}
		}
		dat=rbind(dat, temp)
	}	
	
	return(dat)
}	
###########################################
##' Simulate observed sizes.
##'
##' @param pars a numeric vector contianing the following named terms representing regression coefficients: \code{'(Intercept)', 'age',  'I(age^2)', 'size', 'size:age', 'sigma_proc', 'time_growth_sd', 'indiv_growth_sd', 'indiv_age_growth_sd', 'indiv_cor', 'size0_mu', 'size0_sd'} 
##' @param nind the number of individuals to simulate
##' @param ntime the number of times spanned by the dataset
##' @param maxage the maximum age that all individuals reach (not used if lifespans is specified)
##' @param lifespans a set of intergers to be resampled from. It should reflect the frequency of each lifespan. The units are the same as for times.
##' @param recap the probability of captuing an individual in any time when it was alive.
##' @param sigma_obs the standard deviation of the observation error with the same units as size (e.g. log kg).
##'
##' @details Covariates other than age are not implemented yet. This function calls \code{simgrow()}.
##'
##' @return A data frame of the format needed for input to \code{growmod}. Size observations that are missing due to imperfect recapture are entered as NA because rows of predictors are still needed for those combinations of individual and time. 
##' @export
##' @examples
##' true=c("(Intercept)"=1, age=0,  "I(age^2)"=0, size=.6, "size:age"=0, sigma_proc=.01, time_growth_sd=.1, indiv_growth_sd=.1, indiv_age_growth_sd=0, indiv_cor=0, size0_mu=2.5, size0_sd=1)
##' sg=simobs(true, nind=10, ntime=20)
simobs=function(pars, nind=150, ntime=29, maxage=12, lifespans=NULL, recap=.5, sigma_obs=0.03554119)
{
	dat= simgrow(pars=pars, nind=nind, ntime=ntime, maxage= maxage, lifespans= lifespans)

	#make sure there's at least one obser per indiv.
	sdat=split(dat, dat$ID)
	sd2=lapply(sdat, function(x){
		n=nrow(x)
		x[sample(1:n, max(round(n* recap),1)),]
	})
	dat2=do.call(rbind, sd2)	
	dat2 $size=rnorm(n=nrow(dat2), mean=dat2$size, sd=sigma_obs) #observed sizes
	
	obs=join(dat2, dat[,-1], type="full")			
	return(obs)
}
##########################################
##' Likelihood Ratio Test
##' @param full an object of type \code{growmod} that was fit using the \code{growmod()} function.
##' @param restricted an object of type \code{growmod} that was fit using the \code{growmod()} function with a subset of the covariates used for \code{full}
##' @export
LRtest=function(full,restricted){
    statistic <- 2*(restricted$fit$objective - full$fit$objective)
    df <- length(full$fit$par) - length(restricted$fit$par)
    p.value <- 1-pchisq(statistic,df=df)
    data.frame(statistic,df,p.value)
}
###########################################
NAint=-999999

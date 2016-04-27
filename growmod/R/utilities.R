NAint=-999999
###########################################
organize_data=function(data, sigma_obs){

if(!all(c("size", "t", "t0", "ID") %in% names(data))){stop("data must contain size, t, t0, and ID.")}
if(!min(data$t0)>=min(data$t)){stop("The minimum t0 must be >= to minimum t so that predictors are available from the beginning.")}

#organize data
data$fID=factor(data$ID)
data=transform(data, it=t-min(t),
					iID=as.numeric(fID)-min(as.numeric(fID))
			)
obs=na.omit(data[,c("size", "iID", "it")])
nind=length(unique(obs$iID))
ntimes=length(unique(data$it))
mintime=min(data$t)

indivdat=ddply(data, ~iID, summarize, t0=t0[1], T=max(t), 
					bad=any(!(t0[1]:max(t) %in% t)))
if(any(indivdat$bad)){stop("Each individual must have a row of data for every timepoint from t0 to max(t), even when size observations are missing. Those rows can have NA in the size column, but must contain any predictors needed for formulaX and formulaM.")}
if(any(is.na(indivdat$t0))){stop("Each individual must have a birth time (t0).")}

#take care of repeated measures at the same time point
noreps=ddply(data, ~iID+it, summarize, msize=mean(size))
preddat=join(noreps, data, match="first") 
preddat$size[is.na(preddat$size)]=mean(preddat$size, na.rm=TRUE)
tmp= preddat[, c('it', 'iID')]#no rep sizes
tmp$value=0:(nrow(tmp)-1)
tmp2=reshape::cast(tmp, iID~it, fill= NAint)[,-1]
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
Lpin=list(size=preddat$size, #rep(0, length(Minit)), 
		log_time_growth_sd=0,
		log_indiv_growth_sd=0,
		log_cohort_growth_sd=0,
		log_resid_growth_sd=0,
		log_size0_sd=0,
		size0_mu=mean(subset(data, t==t0)$size),
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
growmod=function(formulaX, formulaM, data, estobserr=FALSE, sigma_obs = 0.03554119, predfirstsize =NULL, DLL="growmod", silent=TRUE, selecting=FALSE, REtime=TRUE, REID=TRUE, REcohort=TRUE, REageID=FALSE,...)
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
		if(is.null(Lpin$log_sigma_obs)){Lpin$log_sigma_obs=log(Ldat$sigma_obs)}
	}
	if(estobserr)
	{
		Lpin$log_sigma_obs=0			
	}
	##################################################
	#What to do about first size
	if(is.null(predfirstsize))
	{
		Ldat$B=matrix(1, nrow=length(Ldat$t0), ncol=1)#just using t0 to get length right
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
		mod=c(AIC=TMBAIC(fit), convergence=fit$convergence)
		class(mod)="selecting_growmod"
		return(mod)
	}
	#else
	sdr=sdreport(obj)
	if(!any(is.na(summary(sdr, "fixed")))){
		mod=list(parList=obj$env$parList(), 
				fit=fit, 
				sdr=sdr, 
				Xnames=colnames(Ldat$X), Mnames=colnames(Ldat$M), 
				formulaX= formulaX, formulaM= formulaM, 
				AIC=TMBAIC(fit),
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
				Xnames=colnames(Ldat$X), Mnames=colnames(Ldat$M), 
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
TMBAIC=function(opt){2*length(opt[["par"]])+2*opt[["objective"]]}
###########################################
extract_coefs=function(x, CV=FALSE, size0=FALSE){
#takes  results of growmod() and returns a data frame of est, sd, z.value, and names of coefficients
	n0=1+length(x$Xnames)+length(x$Mnames)
	cf=summary(x$sdr, "fixed")[(1:n0),]
	nf=rownames(cf)
	nf[nf=="beta"]=x$Xnames
	nf[nf=="eta"]=paste0("size:",x$Mnames)
	nf[nf=="size:(Intercept)"]="size"
	rownames(cf)=nf
	cr=summary(x$sdr, "report")
	take=grep("growth_sd",rownames(cr))
	take=c(take, grep("sigma_obs",rownames(cr)))
	take=c(take, grep("sigma_proc",rownames(cr)))
	take=c(take, grep("indiv_cor",rownames(cr)))
	if(CV) { take = c(take, grep("CV",rownames(cr)))}
	if(size0) { take = c(take, grep("size0",rownames(cr)))}
	cr2=cr[take, ]
	c2=as.data.frame(rbind(cf,cr2))
  	c2$var=rownames(c2)
  	colnames(c2)[1]="est"
	colnames(c2)[2]="se"
	c3=transform(c2, "t value"=est/se)
	return(c3)
}	
###########################################
summary.growmod=function(mod){
	cat("State-space growth model fit by maximum marginal likelihood estimation \n")
	cat("Formula X: ", as.character(mod$formulaX), "\n")
	cat("Formula M: ", as.character(mod$formulaM), "\n\n")
	cat("Convergence: ", as.character(mod$fit$message), "\n")
	cat("AIC: ", as.character(mod$AIC), "\n\n")
	cat("number of observations: ", as.character(mod$nobs), "\n")
	cat("number of missing observations: ", as.character(mod$nmis), "\n")
	cat("recapture rate: ", as.character(round(mod$recap, 3)), "\n")

	cat("Fixed and Random effects: \n")
	print(extract_coefs(mod)[,-3])
}	
###########################################
extract_pred=function(mod, data){
	od=organize_data(data, 0)
	rr=summary(mod$sdr, "random")
	nr=rownames(rr)
	size=rr[nr=="size",]
	pred=data.frame(est=size[,1], se=size[,2], age=od$Ldat$age, ID=d$Ldat$Predictors$ID, lamb=od$Ldat$Predictors$lamb, wnao=od$Ldat$Predictors$wnao, pop=od$Ldat$Predictors$prevpop)
	return(pred)
}	
###########################################
extract_indiv_dev=function(mod, data){
	od=organize_data(data, 0)
	Ldat=od$Ldat
	i=mod$parList$indiv_growth_dev
	dat=Ldat$Predictors
	dat$dev=i[Ldat$Predictors$iID+1]
	return(dat)
}	
###########################################
pred_age=function(m1, newdata, i=FALSE){
#there can only be one indiv ontogeny in newdata	
	X=model.matrix(m1$formulaX, newdata)
	M=model.matrix(m1$formulaM, newdata)
	newdata$size[newdata$age==0]=m1$parList$size0_mu
	if(i) 
	{
		newdata$sizelo[newdata$age==0]=m1$parList$size0_mu-2*exp(m1$parList$log_size0_sd)
		newdata$sizehi[newdata$age==0]=m1$parList$size0_mu+2*exp(m1$parList$log_size0_sd)
	}	
	beta= m1$parList$beta
	eta=m1$parList$eta
	for(a in 2:length(newdata$age))
	{
		newdata$size[a]=X[a,]%*%beta + M[a,]%*%eta*newdata$size[a-1]
		if(i)
		{
			if(is.null(m1$parList$log_indiv_age_growth_sd))
			{
				newdata$sizelo[a]=X[a,]%*%beta-2*exp(m1$parList$log_indiv_growth_sd)+ M[a,]%*%eta*newdata$sizelo[a-1]
				newdata$sizehi[a]=X[a,]%*%beta+2*exp(m1$parList$log_indiv_growth_sd)+ M[a,]%*%eta*newdata$sizehi[a-1]
			}else{
				
			}	
			
		}	
	}
	newdata			
}
###########################################
pred_size=function(m1, newdata){

	X=model.matrix(m1$formulaX, newdata)
	M=model.matrix(m1$formulaM, newdata)
	return(X%*%m1$parList$beta + M%*%m1$parList$eta*newdata$size)
}	
###########################################
pred_expsize=function(m1, newdata){

	X=model.matrix(m1$formulaX, newdata)
	M=model.matrix(m1$formulaM, newdata)
	r=summary(m1$sdr, "report")
	if("indiv_age_growth_sd" %in%rownames(r))
		stop("pred_expsize is only implemented for models without random slopes")
	sumsigmasq=sum(r[grep("growth_sd",rownames(r)) , "Estimate"]^2)+r["sigma_proc","Estimate"]^2
	return(exp(X%*%m1$parList$beta + M%*%m1$parList$eta*newdata$size + 0.5*sumsigmasq))
}	

###########################################
##' Simulate the growth of individuals without any environmental predictors
##' param pars is a vector contianing the following named terms: (Intercept), age,  I(age^2), size, size:age, sigma_proc, time_growth_sd, indiv_growth_sd, indiv_age_growth_sd, indiv_cor, size0_mu, size0_sd, 
##' param nind is the number of individuals to simulate
##' param ntime is the number of times spanned bt the dataset
##' param maxage is the maximum age that all individuals reach (not used if lifespans is specified)
##' param lifespans is a set of intergers to be resampled from. It should reflect the frequency of each lifespan.
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
##' param pars is a vector contianing the following named terms: (Intercept), age,  I(age^2), size, size:age, sigma_proc, time_growth_sd, indiv_growth_sd, indiv_age_growth_sd, indiv_cor, size0_mu, size0_sd, 
##' param nind is the number of individuals to simulate
##' param ntime is the number of times spanned bt the dataset
##' param maxage is the maximum age that all individuals reach (not used if lifespans is specified)
##' param lifespans is a set of intergers to be resampled from. It should reflect the frequency of each lifespan.
##' param recap is the probability of captuing an individual in any time when it was alive.
##' param sigma_obs is the standard deviation of the observation error with the same units as size (e.g. log kg).
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
LRtest=function(full,restricted){
    statistic <- 2*(restricted$fit$objective - full$fit$objective)
    df <- length(full$fit$par) - length(restricted$fit$par)
    p.value <- 1-pchisq(statistic,df=df)
    data.frame(statistic,df,p.value)
}
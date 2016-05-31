#include <TMB.hpp>                                
template<class Type>
Type objective_function<Type>::operator() ()
{
	/***************INITIALIZATION SECTION*******************************/
	//Read in data
	DATA_VECTOR(obs); //observed sizes [nobs]
	DATA_FACTOR(indiv_obs);//indiv associated with each elem of obs [nobs]
	DATA_FACTOR(time_obs); //time associated with each elem of obs [nobs]
	DATA_MATRIX(X); //covariates [nlatent by npred]
	DATA_MATRIX(M); //covariates & intercept to interact with size [nlatent by npred]
	DATA_MATRIX(B); //covariates for first size [nind by npred]
	DATA_FACTOR(t0); //birth time of all indiv [nind]
	DATA_FACTOR(T); //last time of all indiv (death time or 2014 if still alive) [nind]
	DATA_FACTOR(indiv_size);//indiv associated with each elem of size & row of X [nlatent]
	DATA_FACTOR(time_size); //time associated with each elem of size & row of X [nlatent]
	DATA_VECTOR(age_size); //age associated with each elem of size & row of X [nlatent]
	DATA_IMATRIX(lookup);//elem of size or row of X for indiv i (row) in time y (col) [nind by ntimes]
	DATA_INTEGER(NAint);
	
	//Define parameters to be estimated
	PARAMETER_VECTOR(size); //latent size
//	PARAMETER(size0_mu);
	PARAMETER_VECTOR(beta); //coeffs for growth model except previous size
	PARAMETER_VECTOR(eta); //coeffs for terms interacting with previous size
	PARAMETER_VECTOR(theta); //coeffs for effects on first size
	PARAMETER_VECTOR(indiv_growth_dev); //individual growth deviates from intercept
	PARAMETER_VECTOR(indiv_age_growth_dev); //individual growth deviates from age slope
	PARAMETER_VECTOR(time_growth_dev); //time growth deviates from intercept
	PARAMETER_VECTOR(cohort_growth_dev); //cohort growth deviates
	
	PARAMETER(log_time_growth_sd); //SD of time RE intercept
	Type time_growth_sd=exp(log_time_growth_sd);
	PARAMETER(log_indiv_growth_sd); //SD of individual RE intercept
	Type indiv_growth_sd=exp(log_indiv_growth_sd);
	PARAMETER(log_indiv_age_growth_sd); //SD of individual RE age-slope
	Type indiv_age_growth_sd=exp(log_indiv_age_growth_sd);
	PARAMETER(log_cohort_growth_sd);
	Type cohort_growth_sd=exp(log_cohort_growth_sd);
	PARAMETER(log_size0_sd);	//SD of sizes in first timepoint of life
	Type size0_sd=exp(log_size0_sd);
	PARAMETER(log_sigma_proc); //SD of process error
	Type sigma_proc=exp(log_sigma_proc);
	PARAMETER(log_sigma_obs);
	Type sigma_obs=exp(log_sigma_obs);
	PARAMETER(scale_indiv_cor);
	//correlation is on a scale that bounds -1 < indiv_cor < 1
	Type indiv_cor=Type(-1.0)+Type(2)/(Type(1.0)+exp(-scale_indiv_cor));
	
	Type pred; //to store predicted size
	Type jnll=0; //the joint negative log-likelihood is the objective function
	
	/***************CALCULATION SECTION*******************************/
	//Do matrix multiplication on the linear predictors to use in growth process model
	vector<Type> Xbeta= X*beta;
	vector<Type> Meta= M*eta;
	vector<Type> Btheta= B*theta;
	
	//PROCESS MODEL
	for(int i=0; i<indiv_growth_dev.size(); i++)//for each individual
	{
		jnll-= dnorm(size(lookup(i,t0(i))), Btheta(i), size0_sd, true); //lamb size
		for(int t=t0(i); t<T(i); t++)//for each time of individual i's life
		{
			if(lookup(i,t+1)==NAint)
				std::cerr << "\n Accessing an invalid (individual, time) combination at ("<<i<<", "<<t+1<<"). There was a problem with data organization.";
			//predict the next size based on the linear mixed model
			pred=Meta(lookup(i,t+1))*size(lookup(i,t)) +
				Xbeta(lookup(i,t+1)) + 
				indiv_growth_dev(i) + 
				indiv_age_growth_dev(i)*age_size(lookup(i,t+1)) + 
				time_growth_dev(t+1) +
				cohort_growth_dev(t0(i));
				
			//add the process to the joint negative log-likelihood
			jnll -= dnorm(size(lookup(i,t+1)), pred, sigma_proc, true);	
		}
	}
	//OBSERVATIONS
	for(int o=0; o<obs.size(); o++)//for each observation
		jnll-= dnorm(obs(o), size(lookup(indiv_obs(o), time_obs(o))), sigma_obs, true); 
	//RANDOM EFFECTS
	if(CppAD::Variable(log_indiv_growth_sd))//indiv random intercept
	{
		if(CppAD::Variable(log_indiv_age_growth_sd))//indiv random slope
		{
			if(CppAD::Variable(scale_indiv_cor))//indiv random intercept and slope are correlated
			{
				using namespace density; //a special namespace for multivariate normal densities
				matrix<Type> Sigma(2,2); //to store the covariance matrix
				//fill the covariace matrix
				Sigma(0,0)=indiv_growth_sd*indiv_growth_sd;
				Sigma(1,1)=indiv_age_growth_sd*indiv_age_growth_sd;
				Sigma(0,1)=indiv_growth_sd*indiv_age_growth_sd*indiv_cor;
				Sigma(1,0)=indiv_growth_sd*indiv_age_growth_sd*indiv_cor;
				MVNORM_t<Type> N_0_Sigma(Sigma);//an MVRNORM object returns the negative log-likelihood
				vector<Type> x(2); //to store the individual intercept and slope deviates
				for(int i=0; i<indiv_growth_dev.size(); i++)
				{
					x(0)=indiv_growth_dev(i);//intercept
					x(1)=indiv_age_growth_dev(i);//slope
					jnll+= N_0_Sigma(x);//add the nll of this individual's deviates to the joint-nll
				}
			}
			else
			{//with indiv, with indept age slope, no cor
				for(int i=0; i<indiv_growth_dev.size(); i++)
				{
					jnll-= dnorm(indiv_growth_dev(i), Type(0), indiv_growth_sd, true);
					jnll-= dnorm(indiv_age_growth_dev(i), Type(0), indiv_age_growth_sd, true);
				}
			}
		}
		else
		{//with indiv, no age slope, no cor
			for(int i=0; i<indiv_growth_dev.size(); i++)
			{
				jnll-= dnorm(indiv_growth_dev(i), Type(0), indiv_growth_sd, true);
			}
		}
	}
	//add random  effects of time into the joint negative log-likelihood
	for(int t=0; t<time_growth_dev.size(); t++)
	{
		if(CppAD::Variable(time_growth_dev(t)))
			jnll-= dnorm(time_growth_dev(t), Type(0), time_growth_sd, true);
		
		if(CppAD::Variable(cohort_growth_dev(t)))
			jnll-= dnorm(cohort_growth_dev(t), Type(0), cohort_growth_sd, true);
	}	
	/***************POSTHOC OUTPUT SECTION*******************************/
	//Get CI of derived parameters
	//COEFFICIENTS OF VARIATION
	if (CppAD::Variable(log_indiv_growth_sd))
	{
		Type indivCV=indiv_growth_sd/beta[0];
		ADREPORT(indivCV);
		ADREPORT(indiv_growth_sd);
		if(CppAD::Variable(log_indiv_age_growth_sd))
		{
			ADREPORT(indiv_age_growth_sd);
			if (CppAD::Variable(scale_indiv_cor))
				ADREPORT(indiv_cor);
		}
	}
	if(CppAD::Variable(log_time_growth_sd))
	{
		Type timeCV=time_growth_sd/beta[0];
		ADREPORT(timeCV);
		ADREPORT(time_growth_sd);
	}
	if(CppAD::Variable(log_cohort_growth_sd))
	{
		Type cohortCV=cohort_growth_sd/beta[0];
		ADREPORT(cohortCV);
		ADREPORT(cohort_growth_sd);
	}
	Type procCV=sigma_proc/beta[0];
	ADREPORT(procCV);
	ADREPORT(sigma_proc);
	ADREPORT(size0_sd);
	ADREPORT(sigma_obs);

	return jnll; //return the objective function
}


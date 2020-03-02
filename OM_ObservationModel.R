ObsModel<-function(L, survey, L.age, survey.age, cv.L, cv.survey, n.L, n.survey){

  ## Purpose: Add lognormal observation error to landings and survey, and sampling error to age comps
  ##
  ## INPUT DATA:
  ##   L, survey = error free time series of landings and survey
  ##   L.age, survey.age = error free composition data from which to draw samples
  ## INPUT PARAMETERS:
  ##   cv.L, cv.survey = CVs
  ##   n.L, n.survey = annual sample sizes of age comps 
  ##
  ## OUTPUT:
  ##   time series with lognormal error
  ##   comp data with sampling error

#####################################################################################
#Lognormal error for time series    
  nobs.L=length(L)
  nobs.survey=length(survey)
  
  #SD in log space, given CV in arithmetic space
  sd.L=sqrt(log(1+cv.L^2))
  sd.survey=sqrt(log(1+cv.survey^2))
  
  #log error
  ln.L=rnorm(nobs.L, mean=0, sd=sd.L)
  ln.survey=rnorm(nobs.survey, mean=0, sd=sd.survey)
  
  #multiplicative lognormal error
  L.obs=L*exp(ln.L)
  survey.obs=survey*exp(ln.survey)

#####################################################################################  
#Sampling of age comps (rows=years, columns=ages)
  nages=ncol(L.age)
  ages=1:nages

  nyr.L.age<-nrow(L.age)
  L.age.obs<-matrix(0,nrow=nyr.L.age, ncol=nages)
  for (i in 1:nyr.L.age){
    probs=L.age[i,]/sum(L.age[i,])
    L.age.obs[i,]=rmultinom(n=1, size=n.L, prob=probs)/n.L
  }
  
  nyr.survey.age=nrow(survey.age)
  survey.age.obs=matrix(0,nrow=nyr.survey.age, ncol=nages)
  for (i in 1:nyr.survey.age){
    probs=survey.age[i,]/sum(survey.age[i,])
    survey.age.obs[i,]=rmultinom(n=1, size=n.survey, prob=probs)/n.survey
  }  

  #####################################################################################    
  return(list(L.obs=L.obs, survey.obs=survey.obs, L.age.obs=L.age.obs, survey.age.obs=survey.age.obs, q=sim1.q))
  
} # end ObsModel  
#####################################################################################
#Run the simulation code, write data input files 
#Written: 11 Oct 2019 by KWS
#Last updated: 12 Dec 2019 by BL

#Simulation details: single fleet, single survey, logistic selectivity for both,
#von Bert growth, Bev-Hol recruitment
#####################################################################################
library(readxl)

source("OM_PopSim.R")  #Population simulator
source("OM_ObservationModel.R") #Add observation error to the simulation output (operating model)
source("OM_msy.calcs.R") #Compute MSY-related benchmarks

#set the RN seed for repeatability if generating stochastic results
set.seed(9924)
#####################################################################################
#Define number of years and time series of annual F values
nyr=30
#F= seq(0.01,0.4, length=nyr)*exp(rnorm(nyr, 0, 0.2)) #sequence of annual F values 

cv.L=cv.L      #CV of landings
cv.survey=cv.survey  #CV of survey
n.L<-n.L #annual sample size (nfish) of age comp samples
n.survey<-n.survey #annual sample size (nfish) of age comp samples
logR.sd=logR.sd  #Standard deviation of log recruitment

#####################################################################################
#Most life-history params from Siegfried et al. (data prioritization)
ages=1:12           #Age structure of the popn
nages=length(ages)  #Number ages modeled

R0=1000000  #Average annual unfished recruitment (scales the popn) 
h=0.75	    #Steepness of the Beverton-Holt spawner-recruit relationship. 
M=0.2       #Age-invariant natural mortality

Linf=800	  #Asymptotic average length
K=0.18     	#Growth coefficient 
a0=-1.36    #Theoretical age at size 0 
a.lw=0.000000025  #Length-weight coefficient
b.lw=3.0     	    #Length-weight exponent 
A50.mat=2.25      #Age at 50% maturity
slope.mat=3       #Slope of maturity ogive  

#Define selectivity
A50.sel=2 #Age at 50% selection for landings
slope.sel=1 #Slope of selectivity for landings
A50.sel.survey=1.5 #Age at 50% selection for survey
slope.sel.survey=2 #Slope of selectivity for survey
#####################################################################################
#Create and fill vectors to be used in the population model 
len=Linf*(1-exp(-K*(ages-a0))) #von Bertalanffy growth 
W.kg=a.lw*len^b.lw  #Weight-length relationship, assumed to be in kg 
W.mt=W.kg/1000 #Weight in mt
M.age=rep(M,nages) #natural mortality at age
mat.age=logistic(ages, slope.mat, A50.mat) #maturity at age
proportion.female=rep(0.5, nages) #proportion female at age
selex=logistic(ages,slope.sel,A50.sel) #selectivity
selex.survey=logistic(ages,slope.sel.survey,A50.sel.survey)
#####################################################################################
#Compute the number of spawners per recruit of an unfished population (Phi.0)
N.pr0=rep(1,nages) #Number of spawners per recruit at age
for (a in 1:(nages-1)) 
  {N.pr0[a+1]=N.pr0[a]*exp(-M.age[a])}		
  N.pr0[nages]=N.pr0[nages]/(1-exp(-M.age[nages]))  #Plus group
Phi.0=sum(N.pr0*proportion.female*mat.age*W.mt)     #Spawners per recruit based on mature female biomass 
#logR.resid=rec_dev_matrix[om_sim,]
#logR.resid=rnorm(nyr, mean=0, sd=logR.sd) #generate recruitment residuals to be passed to sim module

par.sim1 <- list(nyr=nyr, ages=ages, cv.L=cv.L, cv.survey=cv.survey, n.L=n.L, n.survey=n.survey, logR.sd=logR.sd, R0=R0, h=h, M=M, Linf=Linf, K=K, a0=a0, a.lw=a.lw, b.lw=b.lw, A50.mat=A50.mat, slope.mat=slope.mat, A50.sel=A50.sel, slope.sel.survey=slope.sel.survey, len=len, W.kg=W.kg, W.mt=W.mt, M.age=M.age, mat.age=mat.age, proportion.female=proportion.female, selex=selex, selex.survey=selex.survey, N.pr0=N.pr0, Phi.0=Phi.0, logR.resid=logR.resid)
#####################################################################################
#Simulate the stock dynamics
input1<-list(nyr=nyr, F=F, ages=ages, nages=nages, R0=R0, h=h, Phi.0=Phi.0, M.age=M.age, W.mt=W.mt, mat.age=mat.age, prop.f=proportion.female, selex=selex, logR.resid=logR.resid)
sim1=popsim(x=input1)

#####################################################################################
#Observation model; add noise to "true" data
survey.sim1.age=sim1$N.age%*%diag(selex.survey)
survey.sim1.raw=rowSums(survey.sim1.age)
survey.sim1=survey.sim1.raw/mean(survey.sim1.raw)
sim1.q=1/mean(survey.sim1.raw)
dat.sim1=ObsModel(L=sim1$L.mt, survey=survey.sim1,
                  L.age=sim1$L.age, survey.age=survey.sim1.age,
                  cv.L=cv.L, cv.survey=cv.survey, n.L=n.L, n.survey=n.survey)

#save.image(file=file.path(mainwd, "op_data.RData"))




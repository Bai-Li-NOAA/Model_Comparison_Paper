#####################################################################################
#Run the simulation code, write data input files 
#Written: 11 Oct 2019 by KWS
#Last updated: 12 Dec 2019 by BL

#Simulation details: single fleet, single survey, logistic selectivity for both,
#von Bert growth, Bev-Hol recruitment
#####################################################################################

#### Set up OM input values ####
yr <- 1:30
nyr <- length(yr)
fleet_num <- 1
survey_num <- 1

# #load recruitment deviations and f deviations
# logR.resid <- logR.resid
# logf.resid <- logf.resid
# f <- f

#CV of landings 
cv.L <- list()
cv.L$fleet1 <- 0.005

input.cv.L <- list()
input.cv.L$fleet1 <- 0.01

#CV of surveys
cv.survey <- list()
cv.survey$survey1 <- 0.1 
#cv.survey$survey2 <- 0.02

input.cv.survey <- list()
input.cv.survey$survey1 <- 0.2
#input.cv.survey$survey2 <- 0.2

#annual sample size (nfish) of age comp samples
n.L <- list()
n.L$fleet1 <- 200

#annual sample size (nfish) of age comp samples
n.survey <- list()
n.survey$survey1 <- 200
#n.survey$survey2 <- 150

#define fleet selectivity
sel_fleet <- list()

sel_fleet$fleet1$A50.sel <- 2
sel_fleet$fleet1$slope.sel <- 1

#define survey selectivity
sel_survey <- list()
sel_survey$survey1$A50.sel <- 1.5
sel_survey$survey1$slope.sel <- 2

#sel_survey$survey2$A50.sel <- 1.5
#sel_survey$survey2$slope.sel <- 1.5
  
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

#Create and fill vectors to be used in the population model 
len=Linf*(1-exp(-K*(ages-a0))) #von Bertalanffy growth 
W.kg=a.lw*len^b.lw  #Weight-length relationship, assumed to be in kg 
W.mt=W.kg/1000 #Weight in mt
M.age=rep(M,nages) #natural mortality at age
mat.age=logistic(ages, slope.mat, A50.mat) #maturity at age
proportion.female=rep(0.5, nages) #proportion female at age

selex_fleet <- vector(mode="list", length=fleet_num)
names(selex_fleet) <- paste("fleet", 1:fleet_num, sep="")
invisible(sapply(1:length(selex_fleet), function(x) selex_fleet[[x]] <<- logistic(ages, sel_fleet[[x]]$slope.sel, sel_fleet[[x]]$A50.sel)))

selex_survey <- vector(mode="list", length=survey_num)
names(selex_survey) <- paste("survey", 1:survey_num, sep="")
invisible(sapply(1:length(selex_survey), function(x) selex_survey[[x]] <<- logistic(ages, sel_survey[[x]]$slope.sel, sel_survey[[x]]$A50.sel)))

#Compute the number of spawners per recruit of an unfished population (Phi.0)
N.pr0=rep(1,nages) #Number of spawners per recruit at age
for (a in 1:(nages-1)) 
  {N.pr0[a+1]=N.pr0[a]*exp(-M.age[a])}		
  N.pr0[nages]=N.pr0[nages]/(1-exp(-M.age[nages]))  #Plus group
Phi.0=sum(N.pr0*proportion.female*mat.age*W.mt)     #Spawners per recruit based on mature female biomass 

om_input <- list(fleet_num=fleet_num, survey_num=survey_num, nyr=nyr, yr=yr, ages=ages, nages=nages, cv.L=cv.L, cv.survey=cv.survey, n.L=n.L, n.survey=n.survey, logR_sd=logR_sd, logf_sd=logf_sd, R0=R0, h=h, M=M, Linf=Linf, K=K, a0=a0, a.lw=a.lw, b.lw=b.lw, A50.mat=A50.mat, slope.mat=slope.mat, sel_fleet=sel_fleet, sel_survey=sel_survey, len=len, W.kg=W.kg, W.mt=W.mt, M.age=M.age, mat.age=mat.age, proportion.female=proportion.female, selex_fleet=selex_fleet, selex_survey=selex_survey, N.pr0=N.pr0, Phi.0=Phi.0, logR.resid=logR.resid, logf.resid=logf.resid, f=f)

#### Simulate the stock dynamics ####
input1<-list(nyr=nyr, f=f, ages=ages, nages=nages, R0=R0, h=h, Phi.0=Phi.0, M.age=M.age, W.mt=W.mt, mat.age=mat.age, prop.f=proportion.female, selex_fleet=selex_fleet, logR.resid=logR.resid)
om_output <- popsim(x=input1)

### Simulate the survey index ####
survey_age_comp <- vector(mode="list", length=survey_num)
names(survey_age_comp) <- paste("survey", 1:survey_num, sep="")
survey_index <- vector(mode="list", length=survey_num)
names(survey_index) <- paste("survey", 1:survey_num, sep="")
survey_q <- vector(mode="list", length=survey_num)
names(survey_q) <- paste("survey", 1:survey_num, sep="")

invisible(sapply(1:length(survey_age_comp), function(x) {
  survey_age_comp[[x]] <<- om_output$N.age%*%diag(selex_survey[[x]])
  survey_annual_sum <- rowSums(survey_age_comp[[x]])
  survey_index[[x]] <<- survey_annual_sum/mean(survey_annual_sum)
  survey_q[[x]] <<- 1/mean(survey_annual_sum)
}))

om_output$survey_age_comp <- survey_age_comp
om_output$survey_index <- survey_index
om_output$survey_q <- survey_q






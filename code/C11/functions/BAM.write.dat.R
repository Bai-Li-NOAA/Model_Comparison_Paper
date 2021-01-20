#Function to write dat file for BAM

BAM.write.dat<-function(fname, nyr, nages, dat.survey, dat.L, parms,
                     a.lw, b.lw, prop.f, mat.age, M.age, BC){

  sink(fname)

  ##########GENERAL HEADER INPUT############################################################
  cat(as.character(paste("#######GENERAL HEADER INPUT#######",collapse="\t")), sep="\n")
  #Starting and ending year of the model (year data starts)
  cat(as.character(paste(1,collapse="\t")), sep="\n")
  cat(as.character(paste(nyr,collapse="\t")), sep="\n", append=TRUE)

  #Starting and ending year to estimate recruitment deviation from S-R curve
  cat(as.character(paste(1,collapse="\t")), sep="\n", append=TRUE)
  cat(as.character(paste(nyr,collapse="\t")), sep="\n", append=TRUE)

  #Possible 3 phases of constraints on recruitment deviations
  cat(as.character(paste(1,collapse="\t")), sep="\n", append=TRUE)
  cat(as.character(paste((nyr),collapse="\t")), sep="\n", append=TRUE)

  #ending year for selectivity blocks (not used in this sim, but architecture left in place)
  cat(as.character(paste((nyr),collapse="\t")), sep="\n", append=TRUE)

  #nages and vector of ages in the population
  cat(as.character(paste(nages,collapse="\t")), sep="\n", append=TRUE)
  cat(as.character(paste(1:nages,collapse="\t")), sep="\n", append=TRUE)

  #nages and vector of ages in the age comps
  cat(as.character(paste(nages,collapse="\t")), sep="\n", append=TRUE)
  cat(as.character(paste(1:nages,collapse="\t")), sep="\n", append=TRUE)

  #Max F used in spr and msy calcs
  cat(as.character(paste(1,collapse="\t")), sep="\n", append=TRUE)
  #Total number of iterations for spr calcs
  cat(as.character(paste(1001,collapse="\t")), sep="\n", append=TRUE)
  #Number years at end of time series over which to average sector F's, for weighted selectivities
  cat(as.character(paste(1,collapse="\t")), sep="\n", append=TRUE)
  #Bias correction (set to 1.0 for no bias correction or a negative value to compute from rec variance)
  if (em_bias_cor==FALSE){
    cat(as.character(paste(1,collapse="\t")), sep="\n", append=TRUE)
    }else{
    cat(as.character(paste(-1,collapse="\t")), sep="\n", append=TRUE)
  }


  ##########SURVEY DATA############################################################
  cat(as.character(paste("#######SURVEY DATA#######",collapse="\t")), sep="\n")
  #number of years
  cat(as.character(paste(dat.survey$nyr,collapse="\t")), sep="\n", append=TRUE)
  #survey years
  cat(as.character(paste(dat.survey$yrs,collapse="\t")), sep="\n", append=TRUE)
  #observed index
  cat(as.character(paste(dat.survey$vals.obs,collapse="\t")), sep="\n", append=TRUE)
  #cv's of observed index
  cat(as.character(paste(dat.survey$cv,collapse="\t")), sep="\n", append=TRUE)

  #number of years of age comps
  cat(as.character(paste(dat.survey$nyr.ages,collapse="\t")), sep="\n", append=TRUE)
  #years of age comps
  cat(as.character(paste(dat.survey$yrs.age,collapse="\t")), sep="\n", append=TRUE)
  #sample size of age comps
  cat(as.character(paste(dat.survey$nsamp,collapse="\t")), sep="\n", append=TRUE)
  #nfish of age comps
  cat(as.character(paste(dat.survey$nfish,collapse="\t")), sep="\n", append=TRUE)
  #age comps
  for (i in 1:dat.survey$nyr.ages){
    cat(as.character(paste(dat.survey$acomp[i,],collapse="\t")), sep="\n", append=TRUE)
  }

  ##########LANDINGS DATA############################################################
  cat(as.character(paste("#######LANDINGS DATA#######",collapse="\t")), sep="\n")
  #start year of landings
  cat(as.character(paste(dat.L$styr,collapse="\t")), sep="\n", append=TRUE)
  #end year of landings
  cat(as.character(paste(dat.L$endyr,collapse="\t")), sep="\n", append=TRUE)
  #observed L
  cat(as.character(paste(dat.L$vals.obs,collapse="\t")), sep="\n", append=TRUE)
  #cv's of observed L
  cat(as.character(paste(dat.L$cv,collapse="\t")), sep="\n", append=TRUE)

  #number of years of age comps
  cat(as.character(paste(dat.L$nyr.ages,collapse="\t")), sep="\n", append=TRUE)
  #years of age comps
  cat(as.character(paste(dat.L$yrs.age,collapse="\t")), sep="\n", append=TRUE)
  #sample size of age comps
  cat(as.character(paste(dat.L$nsamp,collapse="\t")), sep="\n", append=TRUE)
  #nfish of age comps
  cat(as.character(paste(dat.L$nfish,collapse="\t")), sep="\n", append=TRUE)
  #age comps
  for (i in 1:dat.L$nyr.ages){
    cat(as.character(paste(dat.L$acomp[i,],collapse="\t")), sep="\n", append=TRUE)
  }

  ##########PARAMETER SECTION############################################################
  cat(as.character(paste("#######PARAMETER SECTION#######",collapse="\t")), sep="\n")
  nparm=nrow(parms)
  for (i in 1:nparm){
    cat(as.character(paste(parms[i,1:7],collapse="\t")), sep="\n", append=TRUE)
  }

  ##########DEV CONSTRAINTS and DEVS############################################################
  cat(as.character(paste("#######DEV CONSTRAINTS and DEVS#######",collapse="\t")), sep="\n")
  #Fleet1 F devs (lower bound, upper bound, phase)
  cat(as.character(paste(c(-15,15,2),collapse="\t")), sep="\n", append=TRUE)
  #Recruitment devs (lower bound, upper bound, phase)
  cat(as.character(paste(c(-15,15,2),collapse="\t")), sep="\n", append=TRUE)
  #Nage devs (lower bound, upper bound, phase)
  cat(as.character(paste(c(-15,15,2),collapse="\t")), sep="\n", append=TRUE)

  #Fleet1 F devs (initial guesses)
  cat(as.character(paste(rep(0, nyr),collapse="\t")), sep="\n", append=TRUE)
  #Recruitment devs (initial guesses)
  cat(as.character(paste(rep(0, nyr),collapse="\t")), sep="\n", append=TRUE)
  #Nage devs (initial guesses)
  cat(as.character(paste(rep(0, (nages-1)),collapse="\t")), sep="\n", append=TRUE)

  ##########LIKELIHOOD WEIGHTS############################################################
  cat(as.character(paste("#######LIKELIHOOD WEIGHTS#######",collapse="\t")), sep="\n")
  #Landings, survey
  cat(as.character(paste(1.0,collapse="\t")), sep="\n", append=TRUE)
  cat(as.character(paste(1.0,collapse="\t")), sep="\n", append=TRUE)
  #Age comps for landings, survey
  cat(as.character(paste(1.0,collapse="\t")), sep="\n", append=TRUE)
  cat(as.character(paste(1.0,collapse="\t")), sep="\n", append=TRUE)
  #Additional weights: 1.Nage_init, 2.rec devs, 3.rec devs early, 4.rec devs late, 5.Ftune
  cat(as.character(paste(0.0,collapse="\t")), sep="\n", append=TRUE)
  cat(as.character(paste(1.0,collapse="\t")), sep="\n", append=TRUE)
  cat(as.character(paste(0.0,collapse="\t")), sep="\n", append=TRUE)
  cat(as.character(paste(0.0,collapse="\t")), sep="\n", append=TRUE)
  cat(as.character(paste(0.0,collapse="\t")), sep="\n", append=TRUE)

  ##########MISCELLANEOUS STUFF SECTION############################################################
  cat(as.character(paste("#######MISCELLANEOUS STUFF SECTION#######",collapse="\t")), sep="\n")
  #length-weight (TL-whole wgt) coefficients a and b, W=aL^b (W assumed in kg)
  cat(as.character(paste(a.lw,collapse="\t")), sep="\n", append=TRUE)
  cat(as.character(paste(b.lw,collapse="\t")), sep="\n", append=TRUE)
  #prop female, maturity, and M at age
  cat(as.character(paste(prop.f,collapse="\t")), sep="\n", append=TRUE)
  cat(as.character(paste(mat.age,collapse="\t")), sep="\n", append=TRUE)
  cat(as.character(paste(M.age,collapse="\t")), sep="\n", append=TRUE)
  #spawn time fraction
  cat(as.character(paste(0.0,collapse="\t")), sep="\n", append=TRUE)
  #SR switch, tuning year, tuning F,min SS for age comps
  cat(as.character(paste(1,collapse="\t")), sep="\n", append=TRUE)
  cat(as.character(paste(10,collapse="\t")), sep="\n", append=TRUE)
  cat(as.character(paste(0.2,collapse="\t")), sep="\n", append=TRUE)
  cat(as.character(paste(1,collapse="\t")), sep="\n", append=TRUE)
  #ageing error (identity matrix for none)
  I=diag(nages)
  for (i in 1:nages){
    cat(as.character(paste(I[i,],collapse="\t")), sep="\n", append=TRUE)
  }
  #projections: endyr, styr regs, Fproj type, Fproj mult
  cat(as.character(paste((nyr+10),collapse="\t")), sep="\n", append=TRUE)
  cat(as.character(paste((nyr+1),collapse="\t")), sep="\n", append=TRUE)
  cat(as.character(paste(2,collapse="\t")), sep="\n", append=TRUE)
  cat(as.character(paste(1,collapse="\t")), sep="\n", append=TRUE)

  ##########END OF DATA FILE############################################################
  cat(as.character(paste("#######END OF DATA FILE#######",collapse="\t")), sep="\n")
  cat(as.character(paste(999,collapse="\t")), sep="\n", append=TRUE, fill=TRUE)
  sink()
  closeAllConnections()
}


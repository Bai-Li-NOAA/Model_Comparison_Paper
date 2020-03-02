OM_Run <- function(maindir=maindir, subdir="OM", input.cv.L=0.0001, input.cv.survey=0.0001, cv.L=0, cv.survey=0, n.L=200, n.survey=200, logR.sd=0.1, nyr=30, om_sim_num=1, seed_num=seed_num, rec_dev_matrix=rec_dev_matrix, F=F){
  setwd(paste(maindir))
  
  for(om_sim in 1:om_sim_num){
    logR.resid <<- rec_dev_matrix[om_sim,]
    
    source(file= file.path(maindir, "OM_100RunSim.R")) 
    
    set.seed(seed_num[om_sim])
    
    dat.input=ObsModel(L=sim1$L.mt, survey=survey.sim1,
                       L.age=sim1$L.age, survey.age=survey.sim1.age,
                       cv.L=cv.L, cv.survey=cv.survey, n.L=n.L, n.survey=n.survey)
    dat.input$cv.L <- input.cv.L
    dat.input$cv.survey <- input.cv.survey
    
    
    save(sim1, dat.input, dat.sim1, par.sim1, survey.sim1, survey.sim1.age, file=file.path(maindir, subdir, paste("OM", om_sim, ".RData", sep="")))
  }
  graphics.off()
  rm(list=setdiff(ls(), c("maindir", "om_sim_num", lsf.str())))
}
AMAK_Run <- function(maindir=maindir, subdir="AMAK", om_sim_num=1){
  setwd(file.path(maindir, subdir))
  sapply(1:om_sim_num, function(x) dir.create(file.path(maindir, subdir, paste("s", x, sep=""))))
  
  file.copy(file.path(maindir, "amak.dat"), file.path(maindir, subdir, "amak.dat"), overwrite = T)
  file.copy(file.path(maindir, "amak_data.dat"), file.path(maindir, subdir, "amak_data.dat"), overwrite = T)
  file.copy(file.path(maindir, "amak.exe"), file.path(maindir, subdir, "amak.exe"), overwrite = T)
  
  modify_input = "partial"
  for (om_sim in 1:om_sim_num){
    load(file=file.path(maindir, "OM", paste("OM", om_sim, ".RData", sep="")))
    
    if(modify_input == "all") {
      
    } 
    if(modify_input == "partial") {
      ctlf <- file.path(maindir, subdir, "amak.dat")
      
      char.lines <- readLines(ctlf)
      char.lines[grep("#SigmaR", char.lines)+1] <- gsub(gsub(" .*$", "", char.lines[grep("#SigmaR", char.lines)+1]), par.sim1$logR.sd, char.lines[grep("#SigmaR", char.lines)+1]) #Replace the value before 1st space
      char.lines[grep("#catchability", char.lines)+1] <- gsub(gsub(" .*$", "", char.lines[grep("#catchability", char.lines)+1]), dat.input$q, char.lines[grep("#catchability", char.lines)+1])
      writeLines(char.lines, con=file.path(maindir, subdir, paste("s", om_sim, sep=""), "amak.dat"))
      
      datf <- file.path(maindir, subdir, "amak_data.dat")
      char.lines <- readLines(datf)
      char.lines[grep("#catch\t", char.lines)+1] <- as.character(paste(dat.input$L.obs, collapse="\t"))
      char.lines[grep("#catch_cv", char.lines)+1] <- as.character(paste(rep(dat.input$cv.L, par.sim1$nyr), collapse="\t"))
      for (i in 1:par.sim1$nyr){
        char.lines[grep("#page_fsh", char.lines)+i]<-as.character(paste(dat.input$L.age.obs[i,],collapse="\t"))
      }
      char.lines[grep("#sample_ages_fsh", char.lines)+1] <- as.character(paste(rep(par.sim1$n.L, par.sim1$nyr), collapse="\t"))
      
      char.lines[grep("#biom_ind\t", char.lines)+1] <- as.character(paste(dat.input$survey.obs, collapse="\t"))
      char.lines[grep("#biom_cv", char.lines)+1] <- as.character(paste(dat.input$cv.survey*dat.input$survey.obs, collapse="\t"))
      for (i in 1:par.sim1$nyr){
        char.lines[grep("#page_ind", char.lines)+i]<-as.character(paste(dat.input$survey.age.obs[i,],collapse="\t"))
      }
      char.lines[grep("#sample_ages_ind", char.lines)+1] <- as.character(paste(rep(par.sim1$n.survey, par.sim1$nyr), collapse="\t"))
      
      writeLines(char.lines, con=file.path(maindir, subdir, paste("s", om_sim, sep=""),  "amak_data.dat"))
    }
    
    file.copy(file.path(maindir, subdir, "amak.exe"), file.path(maindir, subdir, paste("s", om_sim, sep=""), "amak.exe"), overwrite = T)
  }
  
  for (om_sim in 1:om_sim_num){
    setwd(file.path(maindir, subdir, paste("s", om_sim, sep="")))
    system(paste("amak.exe amak.dat", sep = ""), invisible = T)
  }
}
ASAP_Run <- function(maindir=maindir, subdir="ASAP", om_sim_num=1){
  setwd(file.path(maindir, subdir))
  sapply(1:om_sim_num, function(x) dir.create(file.path(maindir, subdir, paste("s", x, sep=""))))
  
  file.copy(file.path(maindir, "asap3.DAT"), file.path(maindir, subdir, "asap3.DAT"), overwrite = T)
  file.copy(file.path(maindir, "ASAP3.exe"), file.path(maindir, subdir, "ASAP3.exe"), overwrite = T)
  
  modify_input = "partial"
  for (om_sim in 1:om_sim_num){
    asap_input <- ReadASAP3DatFile(file.path(maindir, subdir, "asap3.DAT"))
    load(file=file.path(maindir, "OM", paste("OM", om_sim, ".RData", sep="")))
    if(modify_input == "all") {
      asap_input$dat$n_years <- par.sim1$nyr
      asap_input$dat$year1 <- sim1$yr[1]
      asap_input$dat$n_ages <- length(par.sim1$ages)
      asap_input$dat$n_fleets <- 1 #needs modification in om
      asap_input$dat$n_fleet_sel_blocks <- 1 #needs modification in om
      asap_input$dat$n_indices <- 1 #needs modification in om
      asap_input$dat$M <- matrix(rep(par.sim1$M.age, times=par.sim1$nyr), nrow=par.sim1$nyr, byrow = T)
      asap_input$dat$fec_opt <- 0 #needs modifcation in asap file
      asap_input$dat$fracyr_spawn <- 0 #needs modifcation in om/asap file (month-1)/12
      asap_input$dat$maturity <- matrix(rep(par.sim1$mat.age*par.sim1$proportion.female, times=par.sim1$nyr), nrow=par.sim1$nyr, byrow = T)
      asap_input$dat$n_WAA_mats <- 1 #needs modification in asap file
      asap_input$dat$WAA_mats[[1]] <- matrix(rep(par.sim1$W.kg, times=par.sim1$nyr), nrow=par.sim1$nyr, byrow = T) #needs modification in asap file (list index)
      asap_input$dat$WAA_pointers <- matrix(c(1,1,1,1,1,1), nrow=6, byrow = T) #needs modification in asap file 
    } #incomplement
    if(modify_input == "partial") {
      asap_input$dat$CAA_mats[[1]] <- cbind(dat.input$L.age.obs, dat.input$L.obs)
      asap_input$dat$catch_cv <- cbind(rep(dat.input$cv.L, times=par.sim1$nyr))
      asap_input$dat$IAA_mats[[1]] <- cbind(sim1$yr, dat.input$survey.obs, rep(dat.input$cv.survey, times=par.sim1$nyr), dat.input$survey.age.obs, rep(par.sim1$n.survey, times=par.sim1$nyr))
      asap_input$dat$catch_Neff <- cbind(rep(par.sim1$n.L, par.sim1$nyr))
      asap_input$dat$F1_ini <- sim1$F[1]
      asap_input$dat$q_ini <- dat.input$q*1000
      asap_input$dat$N1_ini <- par.sim1$N.pr0*1000
      asap_input$dat$recruit_cv <- cbind(rep(sqrt(exp(par.sim1$logR.sd^2)-1), par.sim1$nyr))
    }
    WriteASAP3DatFile(fname = file.path(maindir, subdir, paste("s", om_sim, sep=""), "asap3.DAT"), dat.object=asap_input, header.text = "")
    file.copy(file.path(maindir, subdir, "ASAP3.exe"), file.path(maindir, subdir, paste("s", om_sim, sep=""),"ASAP3.exe"), overwrite = T)
  }
  
  for (om_sim in 1:om_sim_num){
    setwd(file.path(maindir, subdir, paste("s", om_sim, sep="")))
    system(paste("ASAP3.exe asap3.DAT", sep = ""), invisible = T)
  }
}
BAM_Run <- function(maindir=maindir, subdir="BAM", om_sim_num=1){
  setwd(file.path(maindir, subdir))
  sapply(1:om_sim_num, function(x) dir.create(file.path(maindir, subdir, paste("s", x, sep=""))))
  
  file.copy(file.path(maindir, "BAM.DatInput.Parms.xlsx"), file.path(maindir, subdir, "BAM.DatInput.Parms.xlsx"), overwrite = T)
  file.copy(file.path(maindir, "bam-sim.exe"), file.path(maindir, subdir, "bam-sim.exe"), overwrite = T)
  
  modify_input = "partial"
  for (om_sim in 1:om_sim_num){
    load(file=file.path(maindir, "OM", paste("OM", om_sim, ".RData", sep="")))
    
    if(modify_input == "all") {
      
    } #incomplete
    if (modify_input == "partial"){
      bam_parms=read_xlsx(path=file.path(maindir, subdir, "BAM.DatInput.Parms.xlsx"))
      bam_parms[bam_parms$notes=="sd of recruitment in log space",c(1,5)]=par.sim1$logR.sd
      bam_parms[bam_parms$notes=="catchability of survey",c(1,5)]=log(dat.sim1$q)
      bam_parms[bam_parms$notes=="log mean F",c(1,5)]=log(mean(sim1$F))
      
      dat1.survey=list(nyr=par.sim1$nyr, yrs=sim1$yr, vals.obs=dat.input$survey.obs, cv=rep(dat.input$cv.survey, par.sim1$nyr), nyr.ages=par.sim1$nyr, yrs.age=sim1$yr, nsamp=rep(par.sim1$n.survey, par.sim1$nyr), nfish=rep(par.sim1$n.survey, par.sim1$nyr), acomp=dat.input$survey.age.obs)
      
      dat1.L=list(styr=sim1$yr[1],endyr=sim1$yr[length(sim1$yr)], vals.obs=dat.input$L.obs, cv=rep(dat.input$cv.L, par.sim1$nyr), nyr.ages=par.sim1$nyr, yrs.age=sim1$yr, nsamp=rep(par.sim1$n.L, par.sim1$nyr), nfish=rep(par.sim1$n.L, par.sim1$nyr), acomp=dat.input$L.age.obs)
      
      BAM.write.dat(fname=file.path(maindir, subdir, paste("s", om_sim, sep=""), "BAM-Sim.dat"), nyr=par.sim1$nyr, nages=length(par.sim1$ages), dat.survey=dat1.survey, dat.L=dat1.L, parms=bam_parms, a.lw=par.sim1$a.lw, b.lw=par.sim1$b.lw, prop.f=par.sim1$proportion.female, par.sim1$mat.age, par.sim1$M.age)
    }
    file.copy(file.path(maindir, subdir, "bam-sim.exe"), file.path(maindir, subdir, paste("s", om_sim, sep=""), "bam-sim.exe"), overwrite = T)
  }
  
  for (om_sim in 1:om_sim_num){
    setwd(file.path(maindir, subdir, paste("s", om_sim, sep="")))
    system(paste("bam-sim.exe BAM-Sim.DAT", sep = ""), invisible = T)
  }
}
SS_Run <- function(maindir=maindir, subdir="SS", om_sim_num=1){
  setwd(file.path(maindir, subdir))
  sapply(1:om_sim_num, function(x) dir.create(file.path(maindir, subdir, paste("s", x, sep=""))))
  
  file.copy(file.path(maindir, "starter.ss"), file.path(maindir, subdir, "starter.ss"), overwrite = T)
  file.copy(file.path(maindir, "forecast.ss"), file.path(maindir, subdir, "forecast.ss"), overwrite = T)
  file.copy(file.path(maindir, "data.ss"), file.path(maindir, subdir, "data.ss"), overwrite = T)
  file.copy(file.path(maindir, "control.ss"), file.path(maindir, subdir, "control.ss"), overwrite = T)
  file.copy(file.path(maindir, "wtatage.ss"), file.path(maindir, subdir, "wtatage.ss"), overwrite = T)
  file.copy(file.path(maindir, "ss.exe"), file.path(maindir, subdir, "ss.exe"), overwrite = T)
  
  modify_input = "partial"
  for (om_sim in 1:om_sim_num){
    load(file=file.path(maindir, "OM", paste("OM", om_sim, ".RData", sep="")))
    
    ss_data <- SS_readdat_3.30(file.path(maindir, subdir, "data.ss"), verbose=FALSE)
    
    ss_ctl <- SS_readctl_3.30(file.path(maindir, subdir, "control.ss"), verbose=FALSE, nseas=1, N_areas=1, Nages=length(par.sim1$ages), Ngenders=1, Npopbins=NA, Nfleet=1, Nsurveys=1, Do_AgeKey=FALSE, N_tag_groups=NA, N_CPUE_obs=c(0,0,9,12))
    
    if(modify_input == "all") {
      
    } #incomplement
    if(modify_input == "partial") {
      ss_data$catch <- as.data.frame(cbind(sim1$yr, rep(1, times=par.sim1$nyr), rep(1, times=par.sim1$nyr), dat.input$L.obs, rep(sqrt(log(1+dat.input$cv.L^2)), times=par.sim1$nyr)))
      names(ss_data$catch) <- c("year", "seas", "fleet", "catch", "catch_se")
      
      ss_data$CPUE <- as.data.frame(cbind(sim1$yr, rep(1, times=par.sim1$nyr), rep(2, times=par.sim1$nyr), dat.input$survey.obs, rep(sqrt(log(1+dat.input$cv.survey^2)), times=par.sim1$nyr)))
      names(ss_data$CPUE) <- c("year", "seas", "index", "obs", "se_log")
      ss_data$agecomp[, 9:ncol(ss_data$agecomp)] <- rbind(cbind(rep(par.sim1$n.L, par.sim1$nyr), dat.input$L.age.obs), cbind(rep(par.sim1$n.survey, par.sim1$nyr), dat.input$survey.age.obs))
      
      ss_ctl$Q_parms$INIT <- log(dat.input$q*1000)
      ss_ctl$Q_parms$PRIOR <- log(dat.input$q*1000)
      #ctl$SRparm["SR_sigmaR", "INIT"]
      ss_ctl$SRparm$INIT[which(rownames(ss_ctl$SRparm)=="SR_sigmaR")] <- par.sim1$logR.sd
      ss_ctl$SRparm$PRIOR[which(rownames(ss_ctl$SRparm)=="SR_sigmaR")] <- par.sim1$logR.sd
      ss_ctl$F_setup[1] <- sim1$F[1]
    }
    SS_writedat_3.30(datlist = ss_data, outfile = file.path(maindir, subdir, paste("s", om_sim, sep=""), "data.ss"), overwrite = T, verbose = FALSE)
    SS_writectl_3.30(ctllist = ss_ctl, outfile = file.path(maindir, subdir, paste("s", om_sim, sep=""), "control.ss"), overwrite = T, verbose = F)
    
    file.copy(file.path(maindir, subdir, "starter.ss"), file.path(maindir, subdir, paste("s", om_sim, sep=""), "starter.ss"), overwrite = T)
    file.copy(file.path(maindir, subdir, "forecast.ss"), file.path(maindir, subdir, paste("s", om_sim, sep=""), "forecast.ss"), overwrite = T)
    file.copy(file.path(maindir, subdir, "wtatage.ss"), file.path(maindir, subdir, paste("s", om_sim, sep=""), "wtatage.ss"), overwrite = T)
    file.copy(file.path(maindir, subdir, "ss.exe"), file.path(maindir, subdir, paste("s", om_sim, sep=""), "ss.exe"), overwrite = T)
  }
  
  for (om_sim in 1:om_sim_num){
    setwd(file.path(maindir, subdir, paste("s", om_sim, sep="")))
    system(paste("ss.exe data.ss", sep = ""), invisible = T)
  }
}

recruit_dev_case <- function(rec_dev_change=T, nyr=30, logR.sd=0.1, om_sim_num=om_sim_num){
  if (rec_dev_change==F){
    rec_dev_matrix <- matrix(NA, nrow=om_sim_num, ncol=nyr)
    for(om_sim in 1:om_sim_num){
      set.seed(9924)
      rec_dev_matrix[om_sim,]<-rnorm(nyr, mean=0, sd=logR.sd) #generate recruitment residuals to be passed to sim module
    }
  }
  if (rec_dev_change==T){
    rec_dev_matrix <- matrix(NA, nrow=om_sim_num, ncol=nyr)
    for(om_sim in 1:om_sim_num){
      set.seed(seed_num[om_sim])
      rec_dev_matrix[om_sim,]<-rnorm(nyr, mean=0, sd=logR.sd) #generate recruitment residuals to be passed to sim module
    }
  }
  return(rec_dev_matrix)
}
f_case <- function(f_pattern=1, f_min=0.01, f_max=0.4, f_end=NULL, f_mean=NULL, f_sd=0.2, nyr=30){
  ## f_pattern: 1. increase
  ##            2. decrease
  ##            3. increase first, then decrease
  ##            4. decrease first, then increase
  ##            5. fluctuate around a value 
  if(f_pattern == 1) F=seq(f_min, f_max, length=nyr)*exp(rnorm(nyr, 0, f_sd))
  if(f_pattern == 2) F=seq(f_max, f_min, length=nyr)*exp(rnorm(nyr, 0, f_sd))
  if(f_pattern == 3) F=c(seq(f_min, f_max, length=round(nyr/5*4)), seq(f_max, f_end, length=nyr-round(nyr/5*4)))*exp(rnorm(nyr, 0, f_sd))
  if(f_pattern == 4) F=c(seq(f_max, f_min, length=round(nyr/5*4)), seq(f_min, f_end, length=nyr-round(nyr/5*4)))*exp(rnorm(nyr, 0, f_sd))
  if(f_pattern == 5) F=seq(f_mean, f_mean, length=nyr)*exp(rnorm(nyr, 0, f_sd))
  return(F)
}

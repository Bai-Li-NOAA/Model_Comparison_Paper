BAM_Run <- function(maindir=maindir, subdir="BAM", om_sim_num=1){
  setwd(file.path(maindir, "output", subdir))
  sapply(1:om_sim_num, function(x) dir.create(file.path(maindir, "output", subdir, paste("s", x, sep=""))))
  
  bam_parms=read_xlsx(path=file.path(maindir, "em_input", "BAM.DatInput.Parms.xlsx"))
  
  modify_input = "partial"
  for (om_sim in 1:om_sim_num){
    load(file=file.path(maindir, "output", "OM", paste("OM", om_sim, ".RData", sep="")))
    
    if(modify_input == "all") {
      
    } #incomplete
    if (modify_input == "partial"){
      bam_parms[bam_parms$notes=="sd of recruitment in log space",c(1,5)]=om_input$logR_sd
      bam_parms[bam_parms$notes=="catchability of survey",c(1,5)]=log(em_input$survey_q$survey1)
      bam_parms[bam_parms$notes=="log mean F",c(1,5)]=log(mean(om_output$f))
      
      dat1.survey=list(nyr=length(em_input$survey.obs$survey1), yrs=om_input$yr, vals.obs=em_input$survey.obs$survey1, cv=rep(em_input$cv.survey$survey1, length(em_input$survey.obs$survey1)), nyr.ages=length(em_input$survey.obs$survey1), yrs.age=om_input$yr, nsamp=rep(em_input$n.survey$survey1, length(em_input$survey.obs$survey1)), nfish=rep(em_input$n.survey$survey1, length(em_input$survey.obs$survey1)), acomp=em_input$survey.age.obs$survey1)
      
      dat1.L=list(styr=om_input$yr[1],endyr=om_input$yr[length(om_input$yr)], vals.obs=em_input$L.obs$fleet1, cv=rep(em_input$cv.L$fleet1, length(em_input$L.obs$fleet1)), nyr.ages=length(em_input$L.obs$fleet1), yrs.age=om_input$yr, nsamp=rep(em_input$n.L$fleet1, length(em_input$L.obs$fleet1)), nfish=rep(em_input$n.L$fleet1, length(em_input$L.obs$fleet1)), acomp=em_input$L.age.obs$fleet1)
      
      BAM.write.dat(fname=file.path(maindir, "output", subdir, paste("s", om_sim, sep=""), "BAM-Sim.dat"), nyr=om_input$nyr, nages=om_input$nages, dat.survey=dat1.survey, dat.L=dat1.L, parms=bam_parms, a.lw=om_input$a.lw, b.lw=om_input$b.lw, prop.f=om_input$proportion.female, om_input$mat.age,om_input$M.age)
    }
    file.copy(file.path(maindir, "em_input", "bam-sim.exe"), file.path(maindir,"output", subdir, paste("s", om_sim, sep=""), "bam-sim.exe"), overwrite = T)
  }
  
  for (om_sim in 1:om_sim_num){
    setwd(file.path(maindir, "output", subdir, paste("s", om_sim, sep="")))
    system(paste("bam-sim.exe BAM-Sim.DAT", sep = ""), show.output.on.console = F)
  }
}
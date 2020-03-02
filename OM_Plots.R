################################################################
## Plot outputs from 100 runs of the operating model
## Written: 30 Dec 2019 by BL
## Last updated: 
################################################################
library(PBSadmb)

maindir = "C:/Users/bai.li/Desktop/sim100/om1_1"
model_names = c("OM", "AMAK", "ASAP", "BAM", "SS")
om_sim_num <- 1

subdir = "r4ss_R"
setwd(file.path(maindir, subdir))
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)

timeseries_line_plot <- function(time, ydata, xlab, ylab){
  x_data <- time
  xlim <- range(time)
  y_data <- ydata
  ylim <- range(y_data)
  median_y <- apply(y_data, 1, median)
  
  plot(1, type="n", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, panel.first=grid(lty=1))
  sapply(1:om_sim_num, function(x) lines(x_data, y_data[, x], type="l", lty=1, col="gray70"))
  lines(x_data, median_y, type="o", pch=19, lty=1, col="deepskyblue3", cex=0.5)
}
timeseries_residuals_plot <- function(time, om_data, est_data, xlab, ylab, legend, type){
  ## time = years
  ## om_data = output data from operating model
  ## est_data = output data from assessment model
  ## xlab = x labels
  ## ylab = y labels
  ## legend = title of each subfigure
  ## type = "line" or "bar" or "boxplot"
  x_data <- time
  xlim <- range(time)
  y0_data <- om_data
  y1_data <- est_data
  #residual_data <- (y1_data - y0_data)
  residual_data <- (y1_data - y0_data)/y0_data
  #ylim <- range (residual_data)
  ylim <- c(-0.2, 0.4)
  median_y <- apply(residual_data, 1, median)
  
  if(type=="line") {
    temp<-plot(1, type="n", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim)
    sapply(1:om_sim_num, function(x) lines(x_data, residual_data[, x], type="l", lty=1, col="gray70"))
    abline(h=0, lty=2)
    lines(x_data, median_y, type="o", pch=19, lty=1, col="deepskyblue3", cex=0.5)
    
    legend("topright", legend, bty="n")
  }
  
  if(type=="bar") {
    tmp <- barplot(median_y, xlab=xlab, ylab=ylab, col="gray70")
    box()
    axis(1, at=tmp, labels=x_data)
    legend("topleft", legend, bty="n")
  }
  
  if(type=="boxplot") {
    boxplot(t(residual_data), xlab=xlab, ylab=ylab, col="gray90", notch=F, staplelwd = 1, pch=19, cex=0.5)
    abline(h=0, col="coral3")
    legend("topleft", legend, bty="n")
  }
  
}

#### Aggregate data of biomass, abundance, SSB, recruitment, F (apical F*selectivity), F multiplier, landings, and survey from models to matrix ####

## OM
subdir = "OM"
load(file.path(maindir, subdir, paste("OM", 1, ".RData", sep="")))
om_biomass <- om_abundance <- om_ssb <- om_recruit <- om_Ftot <- om_Fmul <- om_landing <- om_survey <- matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)
om_msy <- om_fmsy <- om_ssbmsy <- c()
om_fratio <- om_ssbratio <- matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)

for(om_sim in 1:om_sim_num){
  load(file.path(maindir, subdir, paste("OM", om_sim, ".RData", sep="")))
  om_biomass[,om_sim] <- sim1$biomass.mt
  om_abundance[,om_sim] <- sim1$abundance/1000
  om_ssb[,om_sim] <- sim1$SSB
  om_recruit[,om_sim] <- sim1$N.age[,1]/1000
  om_Ftot[,om_sim] <- apply(sim1$FAA*par.sim1$selex, 1, max)
  om_Fmul[,om_sim] <- sim1$F
  #om_landing[,om_sim] <- sim1$L.mt
  #om_survey[,om_sim] <- survey.sim1
  om_landing[,om_sim] <- dat.sim1$L.obs
  om_survey[,om_sim] <- dat.sim1$survey.obs
  om_msy[om_sim] <- sim1$msy$msy
  om_fmsy[om_sim] <- sim1$msy$Fmsy
  om_ssbmsy[om_sim] <- sim1$msy$SSBmsy
  om_fratio[, om_sim] <- om_Ftot[, om_sim]/om_fmsy[om_sim]
  om_ssbratio[, om_sim] <- om_ssb[, om_sim]/om_ssbmsy[om_sim]
}
om_list <- list(om_biomass, om_abundance, om_ssb, om_recruit, om_Ftot, om_landing, om_survey, om_msy, om_fmsy, om_ssbmsy, om_fratio, om_ssbratio)
names(om_list) <- c("biomass", "abundance", "ssb", "recruit", "Ftot", "landing", "survey", "msy", "fmsy", "ssbmsy", "fratio", "ssbratio")

## AMAK
amak_biomass <- amak_abundance <- amak_ssb <- amak_recruit <- amak_Ftot <- amak_Fmul <- amak_landing <- amak_survey <- matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)
amak_msy <- amak_fmsy <- amak_ssbmsy <- c()
amak_fratio <- amak_ssbratio <- matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)

subdir = "AMAK"
for(om_sim in 1:om_sim_num){
  setwd(file.path(maindir, subdir, paste("s", om_sim, sep="")))
  amak_output <- readRep("For_R", suffix = ".rep")
  amak_std <- readRep("amak", suffix = ".std")
  amak_biomass[,om_sim] <- amak_output$TotBiom[which(amak_output$TotBiom[,1]>=sim1$yr[1] & amak_output$TotBiom[,1]<=sim1$yr[length(sim1$yr)]),2]
  amak_abundance[,om_sim] <- apply(amak_output$N[,2:ncol(amak_output$N)], 1, sum)/1000
  amak_ssb[,om_sim] <- amak_output$SSB[which(amak_output$SSB[,1]>=sim1$yr[1] &  amak_output$SSB[,1]<= sim1$yr[length(sim1$yr)]),2]
  amak_recruit[,om_sim] <- amak_output$R[,2]/1000
  amak_Ftot[,om_sim] <- apply(amak_output$TotF, 1, max)
  amak_Fmul[,om_sim] <- NA
  amak_landing[,om_sim] <- amak_output$Pred_catch_1
  amak_survey[,om_sim] <- amak_output$Obs_Survey_1[,3]
  amak_msy[om_sim] <- amak_std$value[which(amak_std$name=="MSY")]
  amak_fmsy[om_sim] <- amak_std$value[which(amak_std$name=="Fmsy")]
  amak_ssbmsy[om_sim] <- amak_std$value[which(amak_std$name=="Bmsy")]
  amak_fratio[, om_sim] <- amak_Ftot[, om_sim]/amak_fmsy[om_sim]
  amak_ssbratio[, om_sim] <- amak_ssb[,om_sim]/amak_ssbmsy[om_sim]
}
amak_list <- list(amak_biomass, amak_abundance, amak_ssb, amak_recruit, amak_Ftot, amak_landing, amak_survey, amak_msy, amak_fmsy, amak_ssbmsy, amak_fratio, amak_ssbratio)
names(amak_list) <- c("biomass", "abundance", "ssb", "recruit", "Ftot", "landing", "survey", "msy", "fmsy", "ssbmsy", "fratio", "ssbratio")

## ASAP
asap_biomass <- asap_abundance <- asap_ssb <- asap_recruit <- asap_Ftot <- asap_Fmul <- asap_landing <- asap_survey <- matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)
asap_msy <- asap_fmsy <- asap_ssbmsy <- c()
asap_fratio <- asap_ssbratio <- matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)

subdir = "ASAP"
for(om_sim in 1:om_sim_num){
  asap_output <- dget(file.path(maindir, subdir, paste("s", om_sim, sep=""), "asap3.rdat"))
  setwd(file.path(maindir, subdir, paste("s", om_sim, sep="")))
  asap_std <- readRep("asap3", suffix = ".std")
  asap_biomass[,om_sim] <- asap_output$tot.jan1.B
  asap_abundance[,om_sim] <- apply(asap_output$N.age, 1, sum)
  asap_ssb[,om_sim] <- asap_output$SSB
  asap_recruit[,om_sim] <- asap_output$N.age[,1]
  asap_Ftot[,om_sim] <- apply(asap_output$fleet.FAA$FAA.directed.fleet1, 1, max)
  asap_Fmul[,om_sim] <- asap_output$fleet.Fmult
  asap_landing[,om_sim] <- asap_output$catch.pred
  asap_survey[,om_sim] <- asap_output$index.pred$ind01
  asap_msy[om_sim] <- asap_std$value[which(asap_std$name=="MSY")]
  asap_fmsy[om_sim] <- asap_std$value[which(asap_std$name=="Fmsy_report")]
  asap_ssbmsy[om_sim] <- asap_std$value[which(asap_std$name=="SSBmsy_report")]
  asap_fratio[, om_sim] <- asap_Ftot[, om_sim]/asap_fmsy[om_sim]
  asap_ssbratio[, om_sim] <- asap_ssb[,om_sim]/asap_ssbmsy[om_sim]
}
asap_list <- list(asap_biomass, asap_abundance, asap_ssb, asap_recruit, asap_Ftot, asap_landing, asap_survey, asap_msy, asap_fmsy, asap_ssbmsy, asap_fratio, asap_ssbratio)
names(asap_list) <- c("biomass", "abundance", "ssb", "recruit", "Ftot", "landing", "survey", "msy", "fmsy", "ssbmsy", "fratio", "ssbratio")

## BAM
bam_biomass <- bam_abundance <- bam_ssb <- bam_recruit <- bam_Ftot <- bam_Fmul <- bam_landing <- bam_survey <- matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)
bam_msy <- bam_fmsy <- bam_ssbmsy <- c()
bam_fratio <- bam_ssbratio <- matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)

subdir = "BAM"
for(om_sim in 1:om_sim_num){
  bam_output <- dget(file.path(maindir, subdir, paste("s", om_sim, sep=""), "bam-sim.rdat"))
  
  bam_biomass[,om_sim] <- bam_output$t.series$B[1:sim1$yr[length(sim1$yr)]]
  bam_abundance[,om_sim] <- bam_output$t.series$N[1:sim1$yr[length(sim1$yr)]]/1000
  bam_ssb[,om_sim] <- bam_output$t.series$SSB[1:sim1$yr[length(sim1$yr)]]
  bam_recruit[,om_sim] <- bam_output$t.series$recruits[1:sim1$yr[length(sim1$yr)]]/1000
  bam_Ftot[,om_sim] <- bam_output$t.series$F.full[1:sim1$yr[length(sim1$yr)]]
  bam_Fmul[,om_sim] <- NA
  bam_landing[,om_sim] <- bam_output$t.series$L.fleet1.pr[1:sim1$yr[length(sim1$yr)]]
  bam_survey[,om_sim] <- bam_output$t.series$U.survey1.pr[1:sim1$yr[length(sim1$yr)]]
  bam_msy[om_sim] <- bam_output$parms$msy.mt
  bam_fmsy[om_sim] <- bam_output$parms$Fmsy
  bam_ssbmsy[om_sim] <- bam_output$parms$SSBmsy
  bam_fratio[, om_sim] <- bam_output$t.series$F.Fmsy[1:sim1$yr[length(sim1$yr)]]
  bam_ssbratio[, om_sim] <- bam_output$t.series$SSB.SSBmsy[1:sim1$yr[length(sim1$yr)]]
}
bam_list <- list(bam_biomass, bam_abundance, bam_ssb, bam_recruit, bam_Ftot, bam_landing, bam_survey, bam_msy, bam_fmsy, bam_ssbmsy, bam_fratio, bam_ssbratio)
names(bam_list) <- c("biomass", "abundance", "ssb", "recruit", "Ftot", "landing", "survey", "msy", "fmsy", "ssbmsy", "fratio", "ssbratio")

## SS
ss_biomass <- ss_abundance <- ss_ssb <- ss_recruit <- ss_Ftot <- ss_Fmul <- ss_landing <- ss_survey <- matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)
ss_msy <- ss_fmsy <- ss_ssbmsy <- c()
ss_fratio <- ss_ssbratio <- matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)

subdir = "SS"
for(om_sim in 1:om_sim_num){
  ss_output <- SS_output(dir=file.path(maindir, subdir, paste("s", om_sim, sep="")), ncols = 300)
  
  setwd(file.path(maindir, subdir, paste("s", om_sim, sep="")))
  ss_std <- readRep("ss", suffix = ".std")
  msy_std <- ss_std[ss_std$name=="Mgmt_quant",]
  
  ss_biomass[,om_sim] <- ss_output$timeseries$`SmryBio_SX:1_GP:1`[which(ss_output$timeseries$Yr>=sim1$yr[1] & ss_output$timeseries$Yr<=sim1$yr[length(sim1$yr)])]
  ss_abundance[,om_sim] <- ss_output$timeseries$`SmryNum_SX:1_GP:1`[which(ss_output$timeseries$Yr>=sim1$yr[1] & ss_output$timeseries$Yr<=sim1$yr[length(sim1$yr)])]
  ss_ssb[,om_sim] <- ss_output$timeseries$SpawnBio[which(ss_output$timeseries$Yr>=sim1$yr[1] & ss_output$timeseries$Yr<=sim1$yr[length(sim1$yr)])]
  ss_recruit[,om_sim] <- ss_output$natage_annual_2_with_fishery[which(ss_output$natage_annual_2_with_fishery$Yr>=sim1$yr[1] & ss_output$natage_annual_2_with_fishery$Yr<=sim1$yr[length(sim1$yr)]),5]
  ss_Ftot[,om_sim] <- ss_output$timeseries$`F:_1`[which(ss_output$timeseries$Yr>=sim1$yr[1] & ss_output$timeseries$Yr<=sim1$yr[length(sim1$yr)])]
  ss_Fmul[,om_sim] <- NA
  ss_landing[,om_sim] <- ss_output$timeseries$`sel(B):_1`[which(ss_output$timeseries$Yr>=sim1$yr[1] & ss_output$timeseries$Yr<=sim1$yr[length(sim1$yr)])]
  ss_survey[,om_sim] <- ss_output$cpue$Exp
  ss_msy[om_sim] <- msy_std[15, "value"]
  ss_fmsy[om_sim] <- msy_std[14, "value"]
  ss_ssbmsy[om_sim] <- msy_std[12, "value"]
  
  ss_fratio[, om_sim] <- ss_output$Kobe$F.Fmsy[which(ss_output$Kobe$Yr>=sim1$yr[1] & ss_output$Kobe$Yr<=sim1$yr[length(sim1$yr)])]
  ss_ssbratio[, om_sim] <- ss_output$Kobe$B.Bmsy[which(ss_output$Kobe$Yr>=sim1$yr[1] & ss_output$Kobe$Yr<=sim1$yr[length(sim1$yr)])]
}
ss_list <- list(ss_biomass, ss_abundance, ss_ssb, ss_recruit, ss_Ftot, ss_landing, ss_survey, ss_msy, ss_fmsy, ss_ssbmsy, ss_fratio, ss_ssbratio)
names(ss_list) <- c("biomass", "abundance", "ssb", "recruit", "Ftot", "landing", "survey", "msy", "fmsy", "ssbmsy", "fratio", "ssbratio")

save(om_list, amak_list, asap_list, bam_list, ss_list, file=file.path(maindir, "output_details.RData"))
#### Plot biomass, abundance, SSB, recruitment, F (apical F*selectivity), landings, and survey over time ####
subdir = "Figures"

jpeg(file=file.path(maindir, subdir, "om_timeseries.jpg"), width=170, height=180, units="mm", res=300)
mat=matrix(1:8, ncol=2, nrow=4, byrow = T)
layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
par(las=1, mar=c(4.5,4.5,1,0.5))

ylab = c("Biomass (mt)", "Abundane (1000 fish)", "SSB (mt)", "Recruitment (1000 fish)", "F", "Landings (mt)", "Survey Index (scaled)")
xlab = rep("Year", times=length(ylab))

for(i in 1:length(ylab)){
  timeseries_line_plot(time=sim1$yr, ydata=om_list[[i]], ylab=ylab[i], xlab=xlab[i])
}
dev.off()

jpeg(file=file.path(maindir, subdir, "om_fbratio_timeseries.jpg"), width=170, height=180, units="mm", res=300)
mat=matrix(1:2, ncol=1, nrow=2, byrow = T)
layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
par(las=1, mar=c(4.5,4.5,1,0.5))

ylab = c("F/FMSY", "SSB/SSBMSY")
xlab = rep("Year", times=length(ylab))

for(i in 1:length(ylab)){
  timeseries_line_plot(time=sim1$yr, ydata=om_list[[10+i]], ylab=ylab[i], xlab=xlab[i])
  abline(h=1, col="coral3")
}
dev.off()
#### Plot residuals over time####
subdir = "Figures"

## AMAK
jpeg(file=file.path(maindir, subdir, "amak_residuals_timeseries_lines.jpg"), width=170, height=160, units="mm", res=300)
mat=matrix(1:9, ncol=3, nrow=3, byrow = T)
layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
par(las=1, mar=c(4.5,4.5,1,0.5))

legend = c("Biomass (mt)", "Abundane (1000 fish)", "SSB (mt)", "Recruitment (1000 fish)", "F", "Landings (mt)", "Survey Index (scaled)", "F/FMSY", "SSB/SSBMSY")
xlab = rep("Year", times=length(legend))
ylab = rep("Residuals", times=length(legend))

plot_data_id <- c(1:7, 11:12)
for(i in 1:length(plot_data_id)){
  timeseries_residuals_plot(time=sim1$yr, om_data=om_list[[plot_data_id[i]]], est_data=amak_list[[plot_data_id[i]]], ylab="Residuals", xlab=xlab[i], legend=legend[i], type="line")
}
dev.off()

jpeg(file=file.path(maindir, subdir, "amak_residuals_timeseries_bars.jpg"), width=170, height=160, units="mm", res=300)
mat=matrix(1:9, ncol=3, nrow=3, byrow = T)
layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
par(las=1, mar=c(4.5,4.5,1,0.5))

legend = c("Biomass (mt)", "Abundane (1000 fish)", "SSB (mt)", "Recruitment (1000 fish)", "F", "Landings (mt)", "Survey Index (scaled)", "F/FMSY", "SSB/SSBMSY")
xlab = rep("Year", times=length(legend))
ylab = rep("Residuals", times=length(legend))

plot_data_id <- c(1:7, 11:12)
for(i in 1:length(plot_data_id)){
  timeseries_residuals_plot(time=sim1$yr, om_data=om_list[[plot_data_id[i]]], est_data=amak_list[[plot_data_id[i]]], ylab="Residuals", xlab=xlab[i], legend=legend[i], type="bar")
}
dev.off()

jpeg(file=file.path(maindir, subdir, "amak_residuals_timeseries_boxplot.jpg"), width=170, height=160, units="mm", res=300)
mat=matrix(1:9, ncol=3, nrow=3, byrow = T)
layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
par(las=1, mar=c(4.5,4.5,1,0.5))

legend = c("Biomass (mt)", "Abundane (1000 fish)", "SSB (mt)", "Recruitment (1000 fish)", "F", "Landings (mt)", "Survey Index (scaled)", "F/FMSY", "SSB/SSBMSY")
xlab = rep("Year", times=length(legend))
ylab = rep("Residuals", times=length(legend))

plot_data_id <- c(1:7, 11:12)
for(i in 1:length(plot_data_id)){
  timeseries_residuals_plot(time=sim1$yr, om_data=om_list[[plot_data_id[i]]], est_data=amak_list[[plot_data_id[i]]], ylab="Residuals", xlab=xlab[i], legend=legend[i], type="boxplot")
}
dev.off()

## ASAP
jpeg(file=file.path(maindir, subdir, "asap_residuals_timeseries_lines.jpg"), width=170, height=160, units="mm", res=300)
mat=matrix(1:9, ncol=3, nrow=3, byrow = T)
layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
par(las=1, mar=c(4.5,4.5,1,0.5))

legend = c("Biomass (mt)", "Abundane (1000 fish)", "SSB (mt)", "Recruitment (1000 fish)", "F", "Landings (mt)", "Survey Index (scaled)", "F/FMSY", "SSB/SSBMSY")
xlab = rep("Year", times=length(legend))
ylab = rep("Residuals", times=length(legend))

plot_data_id <- c(1:7, 11:12)
for(i in 1:length(plot_data_id)){
  timeseries_residuals_plot(time=sim1$yr, om_data=om_list[[plot_data_id[i]]], est_data=asap_list[[plot_data_id[i]]], ylab="Residuals", xlab=xlab[i], legend=legend[i], type="line")
}
dev.off()

jpeg(file=file.path(maindir, subdir, "asap_residuals_timeseries_bars.jpg"), width=170, height=160, units="mm", res=300)
mat=matrix(1:9, ncol=3, nrow=3, byrow = T)
layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
par(las=1, mar=c(4.5,4.5,1,0.5))

legend = c("Biomass (mt)", "Abundane (1000 fish)", "SSB (mt)", "Recruitment (1000 fish)", "F", "Landings (mt)", "Survey Index (scaled)", "F/FMSY", "SSB/SSBMSY")
xlab = rep("Year", times=length(legend))
ylab = rep("Residuals", times=length(legend))

plot_data_id <- c(1:7, 11:12)
for(i in 1:length(plot_data_id)){
  timeseries_residuals_plot(time=sim1$yr, om_data=om_list[[plot_data_id[i]]], est_data=asap_list[[plot_data_id[i]]], ylab="Residuals", xlab=xlab[i], legend=legend[i], type="bar")
}
dev.off()

jpeg(file=file.path(maindir, subdir, "asap_residuals_timeseries_boxplot.jpg"), width=170, height=160, units="mm", res=300)
mat=matrix(1:9, ncol=3, nrow=3, byrow = T)
layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
par(las=1, mar=c(4.5,4.5,1,0.5))

legend = c("Biomass (mt)", "Abundane (1000 fish)", "SSB (mt)", "Recruitment (1000 fish)", "F", "Landings (mt)", "Survey Index (scaled)", "F/FMSY", "SSB/SSBMSY")
xlab = rep("Year", times=length(legend))
ylab = rep("Residuals", times=length(legend))

plot_data_id <- c(1:7, 11:12)
for(i in 1:length(plot_data_id)){
  timeseries_residuals_plot(time=sim1$yr, om_data=om_list[[plot_data_id[i]]], est_data=asap_list[[plot_data_id[i]]], ylab="Residuals", xlab=xlab[i], legend=legend[i], type="boxplot")
}
dev.off()

## BAM
jpeg(file=file.path(maindir, subdir, "bam_residuals_timeseries_lines.jpg"), width=170, height=160, units="mm", res=300)
mat=matrix(1:9, ncol=3, nrow=3, byrow = T)
layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
par(las=1, mar=c(4.5,4.5,1,0.5))

legend = c("Biomass (mt)", "Abundane (1000 fish)", "SSB (mt)", "Recruitment (1000 fish)", "F", "Landings (mt)", "Survey Index (scaled)", "F/FMSY", "SSB/SSBMSY")
xlab = rep("Year", times=length(legend))
ylab = rep("Residuals", times=length(legend))

plot_data_id <- c(1:7, 11:12)
for(i in 1:length(plot_data_id)){
  timeseries_residuals_plot(time=sim1$yr, om_data=om_list[[plot_data_id[i]]], est_data=bam_list[[plot_data_id[i]]], ylab="Residuals", xlab=xlab[i], legend=legend[i], type="line")
}
dev.off()

jpeg(file=file.path(maindir, subdir, "bam_residuals_timeseries_bars.jpg"), width=170, height=160, units="mm", res=300)
mat=matrix(1:9, ncol=3, nrow=3, byrow = T)
layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
par(las=1, mar=c(4.5,4.5,1,0.5))

legend = c("Biomass (mt)", "Abundane (1000 fish)", "SSB (mt)", "Recruitment (1000 fish)", "F", "Landings (mt)", "Survey Index (scaled)", "F/FMSY", "SSB/SSBMSY")
xlab = rep("Year", times=length(legend))
ylab = rep("Residuals", times=length(legend))

plot_data_id <- c(1:7, 11:12)
for(i in 1:length(plot_data_id)){
  timeseries_residuals_plot(time=sim1$yr, om_data=om_list[[plot_data_id[i]]], est_data=bam_list[[plot_data_id[i]]], ylab="Residuals", xlab=xlab[i], legend=legend[i], type="bar")
}
dev.off()

jpeg(file=file.path(maindir, subdir, "bam_residuals_timeseries_boxplot.jpg"), width=170, height=160, units="mm", res=300)
mat=matrix(1:9, ncol=3, nrow=3, byrow = T)
layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
par(las=1, mar=c(4.5,4.5,1,0.5))

legend = c("Biomass (mt)", "Abundane (1000 fish)", "SSB (mt)", "Recruitment (1000 fish)", "F", "Landings (mt)", "Survey Index (scaled)", "F/FMSY", "SSB/SSBMSY")
xlab = rep("Year", times=length(legend))
ylab = rep("Residuals", times=length(legend))

plot_data_id <- c(1:7, 11:12)
for(i in 1:length(plot_data_id)){
  timeseries_residuals_plot(time=sim1$yr, om_data=om_list[[plot_data_id[i]]], est_data=bam_list[[plot_data_id[i]]], ylab="Residuals", xlab=xlab[i], legend=legend[i], type="boxplot")
}
dev.off()

## SS
jpeg(file=file.path(maindir, subdir, "ss_residuals_timeseries_lines.jpg"), width=170, height=160, units="mm", res=300)
mat=matrix(1:9, ncol=3, nrow=3, byrow = T)
layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
par(las=1, mar=c(4.5,4.5,1,0.5))

legend = c("Biomass (mt)", "Abundane (1000 fish)", "SSB (mt)", "Recruitment (1000 fish)", "F", "Landings (mt)", "Survey Index (scaled)", "F/FMSY", "SSB/SSBMSY")
xlab = rep("Year", times=length(legend))
ylab = rep("Residuals", times=length(legend))

plot_data_id <- c(1:7, 11:12)
for(i in 1:length(plot_data_id)){
  timeseries_residuals_plot(time=sim1$yr, om_data=om_list[[plot_data_id[i]]], est_data=ss_list[[plot_data_id[i]]], ylab="Residuals", xlab=xlab[i], legend=legend[i], type="line")
}
dev.off()

jpeg(file=file.path(maindir, subdir, "ss_residuals_timeseries_bars.jpg"), width=170, height=160, units="mm", res=300)
mat=matrix(1:9, ncol=3, nrow=3, byrow = T)
layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
par(las=1, mar=c(4.5,4.5,1,0.5))

legend = c("Biomass (mt)", "Abundane (1000 fish)", "SSB (mt)", "Recruitment (1000 fish)", "F", "Landings (mt)", "Survey Index (scaled)", "F/FMSY", "SSB/SSBMSY")
xlab = rep("Year", times=length(legend))
ylab = rep("Residuals", times=length(legend))

plot_data_id <- c(1:7, 11:12)
for(i in 1:length(plot_data_id)){
  timeseries_residuals_plot(time=sim1$yr, om_data=om_list[[plot_data_id[i]]], est_data=ss_list[[plot_data_id[i]]], ylab="Residuals", xlab=xlab[i], legend=legend[i], type="bar")
}
dev.off()

jpeg(file=file.path(maindir, subdir, "ss_residuals_timeseries_boxplot.jpg"), width=170, height=160, units="mm", res=300)
mat=matrix(1:9, ncol=3, nrow=3, byrow = T)
layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
par(las=1, mar=c(4.5,4.5,1,0.5))

legend = c("Biomass (mt)", "Abundane (1000 fish)", "SSB (mt)", "Recruitment (1000 fish)", "F", "Landings (mt)", "Survey Index (scaled)", "F/FMSY", "SSB/SSBMSY")
xlab = rep("Year", times=length(legend))
ylab = rep("Residuals", times=length(legend))

plot_data_id <- c(1:7, 11:12)
for(i in 1:length(plot_data_id)){
  timeseries_residuals_plot(time=sim1$yr, om_data=om_list[[plot_data_id[i]]], est_data=ss_list[[plot_data_id[i]]], ylab="Residuals", xlab=xlab[i], legend=legend[i], type="boxplot")
}
dev.off()
#### Plot relative error for MSY, FMSY, SSBMSY (boxplot and histgram?) ####
## AMAK
amak_msy_re <- (amak_list$msy - om_list$msy)/om_list$msy
amak_fmsy_re <- (round(amak_list$fmsy, digits = 2) - om_list$fmsy)/om_list$fmsy
amak_ssbmsy_re <- (amak_list$ssbmsy - om_list$ssbmsy)/om_list$ssbmsy

##ASAP
asap_msy_re <- (asap_list$msy - om_list$msy)/om_list$msy
asap_fmsy_re <- (round(asap_list$fmsy, digits = 2) - om_list$fmsy)/om_list$fmsy
asap_ssbmsy_re <- (asap_list$ssbmsy - om_list$ssbmsy)/om_list$ssbmsy

## BAM
bam_msy_re <- (bam_list$msy - om_list$msy)/om_list$msy
bam_fmsy_re <- (round(bam_list$fmsy, digits = 2) - om_list$fmsy)/om_list$fmsy
bam_ssbmsy_re <- (bam_list$ssbmsy - om_list$ssbmsy)/om_list$ssbmsy

## SS
ss_msy_re <- (ss_list$msy - om_list$msy)/om_list$msy
ss_fmsy_re <- (round(ss_list$fmsy, digits = 2) - om_list$fmsy)/om_list$fmsy
ss_ssbmsy_re <- (ss_list$ssbmsy - om_list$ssbmsy)/om_list$ssbmsy

jpeg(file=file.path(maindir, subdir, "msy_fmsy_ssbmsy_boxplot.jpg"), width=170, height=160, units="mm", res=300)
ylim1=c(-0.15, 0.15)
ylim2=c(-1, 1)
ylim3=c(-0.1, 0.1)
mat=matrix(1:12, ncol=3, nrow=4, byrow = T)
layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
par(las=1, mar=c(4.5,4.5,1,0.5))

boxplot(amak_msy_re, ylim=ylim1)
abline(h=0, col="coral3", lty=2)
boxplot(amak_fmsy_re, ylim=ylim2)
abline(h=0, col="coral3", lty=2)
boxplot(amak_ssbmsy_re, ylim=ylim3)
abline(h=0, col="coral3", lty=2)

boxplot(asap_msy_re, ylim=ylim1)
abline(h=0, col="coral3", lty=2)
boxplot(asap_fmsy_re, ylim=ylim2)
abline(h=0, col="coral3", lty=2)
boxplot(asap_ssbmsy_re, ylim=ylim3)
abline(h=0, col="coral3", lty=2)

boxplot(bam_msy_re, ylim=ylim1)
abline(h=0, col="coral3", lty=2)
boxplot(bam_fmsy_re, ylim=ylim2)
abline(h=0, col="coral3", lty=2)
boxplot(bam_ssbmsy_re, ylim=ylim3)
abline(h=0, col="coral3", lty=2)


boxplot(ss_msy_re, ylim=ylim1)
abline(h=0, col="coral3", lty=2)
boxplot(ss_fmsy_re, ylim=ylim2)
abline(h=0, col="coral3", lty=2)
boxplot(ss_ssbmsy_re, ylim=ylim3)
abline(h=0, col="coral3", lty=2)
dev.off()

for (j in seq(1, om_sim_num, by=10)){
#for (j in 1){
  jpeg(file=file.path(maindir, subdir, paste("Fig", j, "_true_estimate_timeseries.jpg", sep="")), width=170, height=160, units="mm", res=300)
  mat=matrix(1:9, ncol=3, nrow=3, byrow = T)
  layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
  par(las=1, mar=c(4.5,4.5,1,0.5))
  
  #jpeg(file=file.path(maindir, subdir, paste("Fig", j, "_true_estimate_timeseries.jpg", sep="")), width=170, height=180, units="mm", res=300)
  #col=c("gray30", "deepskyblue3")
  #mat=matrix(1:8, ncol=2, nrow=4, byrow = T)
  #layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
  #par(las=1, mar=c(4.5,4.5,1,0.5))
  
  for(i in c(1:7, 11:12)){
    ylim=range(om_list[[i]][,j], amak_list[[i]][,j], asap_list[[i]][,j], bam_list[[i]][,j], ss_list[[i]][,j])
    #ylim=range(om_list[[i]][,j], amak_list[[i]][,j])
    plot(sim1$yr, om_list[[i]][,j], pch=19, col=1, cex=0.5, ylim=ylim)
    lines(sim1$yr, amak_list[[i]][,j], typ="l", col=2)
    lines(sim1$yr, asap_list[[i]][,j], typ="l", col=3)
    lines(sim1$yr, bam_list[[i]][,j], typ="l", col=4)
    lines(sim1$yr, ss_list[[i]][,j], typ="l", col=5)
    legend("topright", c("OP", "AMAK", "ASAP", "BAM", "SS"), pch=c(19, NA, NA, NA, NA), lty=c(NA,1,1,1,1), col=1:5, bty="n", cex=0.7)
  }
  dev.off()
}

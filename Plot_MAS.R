#### Check and install missing packages and library all packages ####
list_of_packages <- c("rstudioapi", "gdata", "PBSadmb", "PBSmodelling", "stringr", "matrixcalc", "r4ss", "ASAPplots", "future", "readxl", "scales", "corrplot", "future", "glue", "jsonlite")
missing_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages)) install.packages(missing_packages)
if("ASAPplots" %in% missing_packages) devtools::install_github("cmlegault/ASAPplots", build_vignettes = TRUE)

if(!("PBSmodelling" %in% installed.packages()[,"Package"])) {
  packageurl = "https://cran.r-project.org/src/contrib/PBSmodelling_2.68.8.tar.gz"
  install.packages(packageurl, repos=NULL, type="source")
}

if(!("PBSadmb" %in% installed.packages()[,"Package"])) {
  packageurl = "https://cran.r-project.org/src/contrib/PBSadmb_1.1.4.tar.gz"
  install.packages(packageurl, repos=NULL, type="source")
}

invisible(lapply(list_of_packages, library, character.only = TRUE))

current_path <- getActiveDocumentContext()$path
maindir <- dirname(current_path)

setwd(maindir)
folder_names <- c("OM", "AMAK", "ASAP", "BAM", "MAS", "Figures", "Outputs")

om_sim_num <- 160
keep_sim_num <- 100

load(file.path(maindir, "seed_num.RData"))
seed_num <- seed_num[1:om_sim_num]

input.cv.L=0.01
input.cv.survey=0.2
cv.L=0.005
cv.survey=0.1
n.L=200
n.survey=200
logR.sd=0.2
nyr=30
#### Create main folders ####
sapply(folder_names, function(x) {
  if (!file.exists(file.path(maindir, x))) dir.create(file.path(maindir, x))
})

#### Run Read_Write_Inputs_Fun ####
subdir = "Read_Write_Inputs_Fun"
setwd(file.path(maindir, subdir))
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)

setwd(file.path(maindir, subdir, "r4ss"))
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)
#### Check convergence of iterations ####
if (T){
  positive_hessian = matrix(NA, ncol=4, nrow=om_sim_num)
  gradient_0.0001 = matrix(NA, ncol=4, nrow=om_sim_num)
  gradient_0.001 = matrix(NA, ncol=4, nrow=om_sim_num)
  gradient_0.01 = matrix(NA, ncol=4, nrow=om_sim_num)
  gradient = matrix(NA, ncol=4, nrow=om_sim_num)
  convergence_measures <- list(positive_hessian=positive_hessian, gradient_0.0001=gradient_0.0001, gradient_0.001=gradient_0.0001, gradient_0.01=gradient_0.0001, gradient=gradient)
  
  for (om_sim in 1:om_sim_num){
    subdir="AMAK"
    setwd(file.path(maindir, subdir, paste("s", om_sim, sep="")))
    temp <- getADMBHessian(File=file.path(maindir, subdir, paste("s", om_sim, sep="")), FileName="admodel.hes")$hes
    convergence_measures$positive_hessian[om_sim, 1] <- ifelse(is.positive.definite((temp+t(temp))/2), 1, 0)
    temp <- as.numeric(scan("amak.par", what='', n=16, quiet=TRUE)[c(6,11,16)])[3]
    convergence_measures$gradient_0.0001[om_sim, 1]<- ifelse(temp>0.0001, 0, 1)
    convergence_measures$gradient_0.001[om_sim, 1]<- ifelse(temp>0.001, 0, 1)
    convergence_measures$gradient_0.01[om_sim, 1]<- ifelse(temp>0.01, 0, 1)
    convergence_measures$gradient[om_sim, 1]<- temp
    
    subdir="ASAP"
    setwd(file.path(maindir, subdir, paste("s", om_sim, sep="")))
    temp <- getADMBHessian(File=file.path(maindir, subdir, paste("s", om_sim, sep="")), FileName="admodel.hes")$hes
    convergence_measures$positive_hessian[om_sim, 2] <- ifelse(is.positive.definite((temp+t(temp))/2), 1, 0)
    temp <- as.numeric(scan("asap3.par", what='', n=16, quiet=TRUE)[c(6,11,16)])[3]
    convergence_measures$gradient_0.0001[om_sim, 2]<- ifelse(temp>0.0001, 0, 1)
    convergence_measures$gradient_0.001[om_sim, 2]<- ifelse(temp>0.001, 0, 1)
    convergence_measures$gradient_0.01[om_sim, 2]<- ifelse(temp>0.01, 0, 1)
    convergence_measures$gradient[om_sim, 2]<- temp
    
    subdir="BAM"
    setwd(file.path(maindir, subdir, paste("s", om_sim, sep="")))
    temp <- getADMBHessian(File=file.path(maindir, subdir, paste("s", om_sim, sep="")), FileName="admodel.hes")$hes
    convergence_measures$positive_hessian[om_sim, 3] <- ifelse(is.positive.definite((temp+t(temp))/2), 1, 0)
    temp <- as.numeric(scan("bam-sim.par", what='', n=16, quiet=TRUE)[c(6,11,16)])[3]
    convergence_measures$gradient_0.0001[om_sim, 3]<- ifelse(temp>0.0001, 0, 1)
    convergence_measures$gradient_0.001[om_sim, 3]<- ifelse(temp>0.001, 0, 1)
    convergence_measures$gradient_0.01[om_sim, 3]<- ifelse(temp>0.01, 0, 1)
    convergence_measures$gradient[om_sim, 3]<- temp
    
    subdir="BAM"
    setwd(file.path(maindir, subdir, paste("s", om_sim, sep="")))
    temp <- getADMBHessian(File=file.path(maindir, subdir, paste("s", om_sim, sep="")), FileName="admodel.hes")$hes #Hessian matrix
    convergence_measures$positive_hessian[om_sim, 4] <- ifelse(is.positive.definite((temp+t(temp))/2), 1, 0)
    temp <- as.numeric(scan("bam-sim.par", what='', n=16, quiet=TRUE)[c(6,11,16)])[3] #gradient value
    convergence_measures$gradient_0.0001[om_sim, 4]<- ifelse(temp>0.0001, 0, 1)
    convergence_measures$gradient_0.001[om_sim, 4]<- ifelse(temp>0.001, 0, 1)
    convergence_measures$gradient_0.01[om_sim, 4]<- ifelse(temp>0.01, 0, 1)
    convergence_measures$gradient[om_sim, 4]<- temp
  }
  save(convergence_measures, file=file.path(maindir, "Outputs", "convergence_measures.RData"))
  load(file.path(maindir, "Outputs", "convergence_measures.RData"))
  
  not_positive_hessian <- unique(unlist(sapply(1:ncol(convergence_measures$positive_hessian), function(x) which(convergence_measures$positive_hessian[,x]==0))))
  write.csv(not_positive_hessian, file=file.path(maindir, "Outputs", "Not_positive_hessian.csv"))
  
  col= c("black", "orange", "green", "red", "deepskyblue3")
  if(max(convergence_measures$gradient)<0.1){
    jpeg(file=file.path(maindir, "Figures", "Gradient.jpg"), width=170, height=150, units="mm", res=300)
    par(mfrow=c(2,2), mar=c(4,4,1,1))
    xlim = c(0, max(convergence_measures$gradient))
    bins <- seq(0, max(convergence_measures$gradient)*1.05, by=0.0005)
    hist(convergence_measures$gradient[,1], xlim=xlim, xlab = "Gradient", main="", col=col[2], breaks = bins)
    legend("topright", "AMAK", bty="n")
    box()
    
    hist(convergence_measures$gradient[,2], xlim=xlim, xlab = "Gradient", main="", breaks = bins, col=col[3])
    legend("topright", "ASAP", bty="n")
    box()
    
    hist(convergence_measures$gradient[,3], xlim=xlim, xlab = "Gradient", main="", breaks = bins, col=col[4])
    legend("topright", "BAM", bty="n")
    box()
    
    hist(convergence_measures$gradient[,4], xlim=xlim, xlab = "Gradient", main="", breaks = bins, col=col[5])
    legend("topright", "MAS", bty="n")
    box()
    dev.off()
  }
  
  
  good_gradient_percentage <- matrix(NA, nrow=1, ncol=4)
  good_gradient_percentage[,1] <- length(which(convergence_measures$gradient[,1]<=0.001))/om_sim_num*100
  good_gradient_percentage[,2] <- length(which(convergence_measures$gradient[,2]<=0.001))/om_sim_num*100
  good_gradient_percentage[,3] <- length(which(convergence_measures$gradient[,3]<=0.001))/om_sim_num*100
  good_gradient_percentage[,4] <- length(which(convergence_measures$gradient[,4]<=0.001))/om_sim_num*100
  colnames(good_gradient_percentage) <- c("AMAK", "ASAP", "BAM", "MAS")
  rownames(good_gradient_percentage) <- "Percentage"
  write.csv(good_gradient_percentage, file=file.path(maindir, "Outputs", "good_gradient_ratio.csv"))
  
  keep_sim_id <- c(1:om_sim_num)[-unique(c(which(convergence_measures$gradient[,1] %in% boxplot.stats(convergence_measures$gradient[,1])$out), which(convergence_measures$gradient[,2] %in% boxplot.stats(convergence_measures$gradient[,2])$out), which(convergence_measures$gradient[,3] %in% boxplot.stats(convergence_measures$gradient[,3])$out), which(convergence_measures$gradient[,4] %in% boxplot.stats(convergence_measures$gradient[,4])$out), unique(unlist(sapply(1:ncol(convergence_measures$positive_hessian), function(x) which(convergence_measures$positive_hessian[,x]==0))))))][1:keep_sim_num]
  
  om_sim_num <- length(keep_sim_id)
}

#### Plot functions ####
model_names = c("OM", "AMAK", "ASAP", "BAM", "MAS")
figure_number <- 10
if(om_sim_num <= figure_number) {
  figure_id <- 1:om_sim_num
} else{
  figure_id <- seq(1, om_sim_num, by=round((om_sim_num-1)/(figure_number-1)))
}

real_figure_id <- keep_sim_id[figure_id]
write.csv(real_figure_id, file=file.path(maindir, "Outputs", "Real_Figure_ID.csv"))

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
om_msy <- om_fmsy <- om_ssbmsy <- matrix(NA, nrow=1, ncol=om_sim_num)
om_fratio <- om_ssbratio <- matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)
om_agecomp <- list()
om_landing_err <- om_survey_err <- matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)

for(om_sim in 1:om_sim_num){
  load(file.path(maindir, subdir, paste("OM", keep_sim_id[om_sim], ".RData", sep="")))
  om_biomass[,om_sim] <- sim1$biomass.mt
  om_abundance[,om_sim] <- sim1$abundance/1000
  om_ssb[,om_sim] <- sim1$SSB
  om_recruit[,om_sim] <- sim1$N.age[,1]/1000
  om_Ftot[,om_sim] <- apply(sim1$FAA, 1, max)
  om_Fmul[,om_sim] <- sim1$F
  om_landing[,om_sim] <- sim1$L.mt
  om_survey[,om_sim] <- survey.sim1
  om_msy[, om_sim] <- sim1$msy$msy
  om_fmsy[, om_sim] <- sim1$msy$Fmsy
  om_ssbmsy[, om_sim] <- sim1$msy$SSBmsy
  om_fratio[, om_sim] <- om_Ftot[, om_sim]/om_fmsy[om_sim]
  om_ssbratio[, om_sim] <- om_ssb[, om_sim]/om_ssbmsy[om_sim]
  om_agecomp[[om_sim]] <- apply(sim1$N.age/1000, 1, function(x) x/sum(x))
  om_landing_err[,om_sim] <- dat.input$L.obs
  om_survey_err[,om_sim] <- dat.input$survey.obs
}
om_list <- list(om_biomass, om_abundance, om_ssb, om_recruit, om_Ftot, om_landing, om_survey, om_msy, om_fmsy, om_ssbmsy, om_fratio, om_ssbratio, om_agecomp, om_landing_err, om_survey_err)
names(om_list) <- c("biomass", "abundance", "ssb", "recruit", "Ftot", "landing", "survey", "msy", "fmsy", "ssbmsy", "fratio", "ssbratio", "agecomp", "om_landing_err", "om_survey_err")

## AMAK
amak_biomass <- amak_abundance <- amak_ssb <- amak_recruit <- amak_Ftot <- amak_Fmul <- amak_landing <- amak_survey <- matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)
amak_msy <- amak_fmsy <- amak_ssbmsy <- matrix(NA, nrow=1, ncol=om_sim_num)
amak_fratio <- amak_ssbratio <- matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)
amak_agecomp <- list()

subdir = "AMAK"
for(om_sim in 1:om_sim_num){
  setwd(file.path(maindir, subdir, paste("s", keep_sim_id[om_sim], sep="")))
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
  amak_msy[, om_sim] <- amak_std$value[which(amak_std$name=="MSY")]
  amak_fmsy[, om_sim] <- round(amak_std$value[which(amak_std$name=="Fmsy")], digits = 2)
  amak_ssbmsy[, om_sim] <- amak_std$value[which(amak_std$name=="Bmsy")]
  amak_fratio[, om_sim] <- amak_Ftot[, om_sim]/amak_fmsy[om_sim]
  amak_ssbratio[, om_sim] <- amak_ssb[,om_sim]/amak_ssbmsy[om_sim]
  amak_agecomp[[om_sim]] <- apply(amak_output$N[,2:ncol(amak_output$N)]/1000, 1, function(x) x/sum(x))
  
}
amak_list <- list(amak_biomass, amak_abundance, amak_ssb, amak_recruit, amak_Ftot, amak_landing, amak_survey, amak_msy, amak_fmsy, amak_ssbmsy, amak_fratio, amak_ssbratio, amak_agecomp)
names(amak_list) <- c("biomass", "abundance", "ssb", "recruit", "Ftot", "landing", "survey", "msy", "fmsy", "ssbmsy", "fratio", "ssbratio", "agecomp")

## ASAP
asap_biomass <- asap_abundance <- asap_ssb <- asap_recruit <- asap_Ftot <- asap_Fmul <- asap_landing <- asap_survey <- matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)
asap_msy <- asap_fmsy <- asap_ssbmsy <- matrix(NA, nrow=1, ncol=om_sim_num)
asap_fratio <- asap_ssbratio <- matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)
asap_agecomp <- list()

subdir = "ASAP"
for(om_sim in 1:om_sim_num){
  asap_output <- dget(file.path(maindir, subdir, paste("s", keep_sim_id[om_sim], sep=""), "asap3.rdat"))
  setwd(file.path(maindir, subdir, paste("s", keep_sim_id[om_sim], sep="")))
  asap_std <- readRep("asap3", suffix = ".std")
  asap_biomass[,om_sim] <- asap_output$tot.jan1.B
  asap_abundance[,om_sim] <- apply(asap_output$N.age, 1, sum)
  asap_ssb[,om_sim] <- asap_output$SSB
  asap_recruit[,om_sim] <- asap_output$N.age[,1]
  asap_Ftot[,om_sim] <- apply(asap_output$fleet.FAA$FAA.directed.fleet1, 1, max)
  asap_Fmul[,om_sim] <- asap_output$fleet.Fmult
  asap_landing[,om_sim] <- asap_output$catch.pred
  asap_survey[,om_sim] <- asap_output$index.pred$ind01
  asap_msy[, om_sim] <- asap_std$value[which(asap_std$name=="MSY")]
  asap_fmsy[, om_sim] <- round(asap_std$value[which(asap_std$name=="Fmsy_report")], digits = 2)
  asap_ssbmsy[, om_sim] <- asap_std$value[which(asap_std$name=="SSBmsy_report")]
  asap_fratio[, om_sim] <- asap_Ftot[, om_sim]/asap_fmsy[om_sim]
  asap_ssbratio[, om_sim] <- asap_ssb[,om_sim]/asap_ssbmsy[om_sim]
  asap_agecomp[[om_sim]] <- apply(asap_output$N.age, 1, function(x) x/sum(x))
}
asap_list <- list(asap_biomass, asap_abundance, asap_ssb, asap_recruit, asap_Ftot, asap_landing, asap_survey, asap_msy, asap_fmsy, asap_ssbmsy, asap_fratio, asap_ssbratio, asap_agecomp)
names(asap_list) <- c("biomass", "abundance", "ssb", "recruit", "Ftot", "landing", "survey", "msy", "fmsy", "ssbmsy", "fratio", "ssbratio", "agecomp")

## BAM
bam_biomass <- bam_abundance <- bam_ssb <- bam_recruit <- bam_Ftot <- bam_Fmul <- bam_landing <- bam_survey <- matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)
bam_msy <- bam_fmsy <- bam_ssbmsy <- matrix(NA, nrow=1, ncol=om_sim_num)
bam_fratio <- bam_ssbratio <- matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)
bam_agecomp <- list()

subdir = "BAM"
for(om_sim in 1:om_sim_num){
  bam_output <- dget(file.path(maindir, subdir, paste("s", keep_sim_id[om_sim], sep=""), "bam-sim.rdat"))
  
  bam_biomass[,om_sim] <- bam_output$t.series$B[1:sim1$yr[length(sim1$yr)]]
  bam_abundance[,om_sim] <- bam_output$t.series$N[1:sim1$yr[length(sim1$yr)]]/1000
  bam_ssb[,om_sim] <- bam_output$t.series$SSB[1:sim1$yr[length(sim1$yr)]]
  bam_recruit[,om_sim] <- bam_output$t.series$recruits[1:sim1$yr[length(sim1$yr)]]/1000
  bam_Ftot[,om_sim] <- bam_output$t.series$F.full[1:sim1$yr[length(sim1$yr)]]
  bam_Fmul[,om_sim] <- NA
  bam_landing[,om_sim] <- bam_output$t.series$L.fleet1.pr[1:sim1$yr[length(sim1$yr)]]
  bam_survey[,om_sim] <- bam_output$t.series$U.survey1.pr[1:sim1$yr[length(sim1$yr)]]
  bam_msy[, om_sim] <- bam_output$parms$msy.mt
  bam_fmsy[, om_sim] <- round(bam_output$parms$Fmsy, digits = 2)
  bam_ssbmsy[, om_sim] <- bam_output$parms$SSBmsy
  bam_fratio[, om_sim] <- bam_output$t.series$F.Fmsy[1:sim1$yr[length(sim1$yr)]]
  bam_ssbratio[, om_sim] <- bam_output$t.series$SSB.SSBmsy[1:sim1$yr[length(sim1$yr)]]
  bam_agecomp[[om_sim]] <- apply(bam_output$N.age[1:par.sim1$nyr,]/1000, 1, function(x) x/sum(x))
}
bam_list <- list(bam_biomass, bam_abundance, bam_ssb, bam_recruit, bam_Ftot, bam_landing, bam_survey, bam_msy, bam_fmsy, bam_ssbmsy, bam_fratio, bam_ssbratio, bam_agecomp)
names(bam_list) <- c("biomass", "abundance", "ssb", "recruit", "Ftot", "landing", "survey", "msy", "fmsy", "ssbmsy", "fratio", "ssbratio", "agecomp")

## MAS
mas_biomass <- mas_abundance <- mas_ssb <- mas_recruit <- mas_Ftot <- mas_Fmul <- mas_landing <- mas_survey <- matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)
mas_msy <- mas_fmsy <- mas_ssbmsy <- matrix(NA, nrow=1, ncol=om_sim_num)
mas_fratio <- mas_ssbratio <- matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)
mas_agecomp <- list()

subdir = "MAS"
for(om_sim in 1:om_sim_num){
  output_file <- file.path(maindir, subdir, paste("mas_s1_em", keep_sim_id[om_sim], ".json", sep=""))
  args <- toString(output_file)
  json_output <- read_json(args[1])
  popdy<-json_output$population_dynamics
  years <- popdy$nyears
  seasons<-popdy$nseasons
  areas<-popdy$nareas
  ages_ <- popdy$nages
  
  #population
  pop<-popdy$populations[[1]]
  flt<-popdy$fleets[[1]]
  srvy<-popdy$surveys[[1]]
  area<-popdy$areas[[1]]
  pop_naa<-pop$undifferentiated$numbers_at_age
  pop_recruits<-pop$undifferentiated$recruits
  pop_abundance<-pop$undifferentiated$abundance
  pop_abundance_f<-pop$females$abundance
  pop_abundance_m<-pop$males$abundance
  pop_ssb<-pop$undifferentiated$spawning_stock_biomass
  pop_biomass<-pop$undifferentiated$biomass
  pop_F<-pop$undifferentiated$fishing_mortality
  flt_landings<-flt$undifferentiated$catch_biomass
  srvy_landings<-srvy$undifferentiated$survey_biomass
  srvy_props<-srvy$undifferentiated$survey_proportion_at_age
  
  
  mas_biomass[,om_sim] <- unlist(pop_biomass$values)
  mas_abundance[,om_sim] <- unlist(pop_abundance$values)
  mas_ssb[,om_sim] <- unlist(pop_ssb$values)
  mas_recruit[,om_sim] <- unlist(pop_recruits$values)
  mas_Ftot[,om_sim] <- unlist(pop_F$values)
  mas_Fmul[,om_sim] <- NA
  mas_landing[,om_sim] <- unlist(flt_landings$values)
  mas_survey[,om_sim] <- unlist(srvy_landings$values)
  mas_msy[, om_sim] <- sim1$msy$msy
  mas_fmsy[, om_sim] <- sim1$msy$Fmsy
  mas_ssbmsy[, om_sim] <- sim1$msy$SSBmsy
  mas_fratio[, om_sim] <- om_Ftot[, om_sim]/om_fmsy[om_sim]
  mas_ssbratio[, om_sim] <- om_ssb[, om_sim]/om_ssbmsy[om_sim]
  mas_agecomp[[om_sim]] <- apply(sim1$N.age/1000, 1, function(x) x/sum(x))

}
mas_list <- list(mas_biomass, mas_abundance, mas_ssb, mas_recruit, mas_Ftot, mas_landing, mas_survey, mas_msy, mas_fmsy, mas_ssbmsy, mas_fratio, mas_ssbratio, mas_agecomp)
names(mas_list) <- c("biomass", "abundance", "ssb", "recruit", "Ftot", "landing", "survey", "msy", "fmsy", "ssbmsy", "fratio", "ssbratio", "agecomp")

save(om_list, amak_list, asap_list, bam_list, mas_list, file=file.path(maindir, "Outputs", "om_em_output.RData"))
load(file.path(maindir, "Outputs", "om_em_output.RData"))
#### RE, ARE, RMSE, CC Table ####
performance_measures <- function(TrueVal, EstVal1, EstVal2, EstVal3, EstVal4){
  model_names <- c("AMAK", "ASAP", "BAM", "MAS")
  temp <- list(EstVal1, EstVal2, EstVal3, EstVal4)
  
  re_matrix <- matrix(NA, ncol=length(model_names), nrow = length(TrueVal))
  if (length(TrueVal)==1) {
    re_matrix <- sapply(1:length(model_names), function(x) re_matrix[,x] <- (temp[[x]]-TrueVal)/TrueVal)
    re_matrix <- as.data.frame(t(re_matrix))
    colnames(re_matrix) <- model_names
  }
  else{
    re_matrix <- as.data.frame(sapply(1:length(model_names), function(x) re_matrix[,x] <- (temp[[x]]-TrueVal)/TrueVal))
    colnames(re_matrix) <- model_names
  }
  
  are_matrix <- abs(re_matrix)
  
  cc_matrix <- matrix(NA, ncol=length(model_names), nrow=1)
  cc_matrix <- sapply(1:length(model_names), function(x) cc_matrix[,x] <- sum((TrueVal - mean(TrueVal))*(temp[[x]]-mean(temp[[x]])))/sqrt(sum((TrueVal - mean(TrueVal))^2)*sum((temp[[x]]-mean(temp[[x]]))^2)))
  cc_matrix <- as.data.frame(t(cc_matrix))
  colnames(cc_matrix) <- model_names
  
  rmse_matrix <- matrix(NA, ncol=length(model_names), nrow=1)
  rmse_matrix <- sapply(1:length(model_names), function(x) rmse_matrix[,x] <- sqrt(sum((temp[[x]]-TrueVal)^2)/length(TrueVal)))
  rmse_matrix <- as.data.frame(t(rmse_matrix))
  colnames(rmse_matrix) <- model_names
  
  return(list(re=re_matrix, are=are_matrix, cc=cc_matrix, rmse=rmse_matrix))
}

data_id <- which(names(om_list) %in% c("biomass", "abundance", "ssb", "recruit", "Ftot", "landing", "survey", "msy", "fmsy", "ssbmsy", "fratio", "ssbratio"))
pm_matrix_list <- list()
pm_list <- list()
for(j in 1:om_sim_num){
  for(i in 1:length(data_id)){
    pm_matrix_list[[i]] <- performance_measures(TrueVal=om_list[[data_id[i]]][,j], EstVal1=amak_list[[data_id[i]]][,j], EstVal2=asap_list[[data_id[i]]][,j], EstVal3=bam_list[[data_id[i]]][,j], EstVal4=mas_list[[data_id[i]]][,j])
  }
  names(pm_matrix_list) <- names(om_list)[data_id]
  pm_list[[j]] <- pm_matrix_list
}

#### AE (age composition) ####
ae <- function(TrueVal, EstVal1, EstVal2, EstVal3, EstVal4){
  model_names <- c("AMAK", "ASAP", "BAM", "MAS")
  temp <- list(EstVal1, EstVal2, EstVal3, EstVal4)
  
  ae_matrix <- matrix(NA, ncol=length(model_names), nrow = ncol(TrueVal))
  ae_matrix <- sapply(1:length(model_names), function(x) ae_matrix[,x] <- apply(abs(temp[[x]]-TrueVal), 2, sum)/length(par.sim1$ages))
  colnames(ae_matrix) <- model_names
  
  return(list(ae=ae_matrix))
}

ae_list <- list()
for(j in 1:om_sim_num){
  ae_list[[j]] <- ae(TrueVal=om_list$agecomp[[j]], amak_list$agecomp[[j]], asap_list$agecomp[[j]], bam_list$agecomp[[j]], mas_list$agecomp[[j]])
}

save(pm_list, ae_list, file=file.path(maindir, "Outputs", "performance_measures_output.RData"))
load(file.path(maindir, "Outputs", "performance_measures_output.RData"))
#### Plot biomass, abundance, SSB, recruitment, F (apical F*selectivity), landings, and survey over time ####
subdir = "Figures"

ylab = c("Biomass (mt)", "Abundance (1000 fish)", "SSB (mt)", "Recruitment (1000 fish)", "F", "Landings (mt)", "Survey Index (scaled)")
xlab = rep("Year", times=length(ylab))
plot_data_id <- 1:length(ylab)
col= c("black", "orange", "green", "red", "deepskyblue3")
nrow=2
ncol=4

for (j in figure_id){
  jpeg(file=file.path(maindir, subdir, paste("Fig", j, "_true_estimate_timeseries.jpg", sep="")), width=205, height=75, units="mm", res=300)
  
  par(mar=c(0.7, 4, 0.2, 0.2), mfrow=c(nrow,ncol), oma = c(4, 4, 0.2, 0.2))
  
  for(i in 1:length(plot_data_id)){
    ylim=range(om_list[[plot_data_id[i]]][,j], amak_list[[plot_data_id[i]]][,j], asap_list[[plot_data_id[i]]][,j], bam_list[[plot_data_id[i]]][,j], mas_list[[plot_data_id[i]]][,j])
    plot(sim1$yr, om_list[[plot_data_id[i]]][,j], pch=19, col=col[1], cex=0.7, ylim=ylim, axes=F, xlab="", ylab="")
    lines(sim1$yr, amak_list[[plot_data_id[i]]][,j], typ="l", col=col[2], lwd=1.5)
    lines(sim1$yr, asap_list[[plot_data_id[i]]][,j], typ="l", col=col[3], lwd=1.5)
    lines(sim1$yr, bam_list[[plot_data_id[i]]][,j], typ="l", col=col[4], lwd=1.5)
    lines(sim1$yr, mas_list[[plot_data_id[i]]][,j], typ="l", col=col[5], lwd=1.5)
    if(i < ncol) axis(1, labels = F)
    else {
      axis(1)
      mtext(side=1, text=xlab, line=2, cex=0.6, font=2, col="blue")
    }
    mtext(side=2, text=ylab[i], line=3.1, cex=0.6, font=2, col="blue")
    axis(2, las=2)
    box()
  }
  plot.new()
  legend("bottom", c("OM", "AMAK", "ASAP", "BAM", "MAS"), pch=c(19, NA, NA, NA, NA), lty=c(NA,1,1,1,1), lwd=c(NA, 1.5, 1.5, 1.5, 1,5), col=col, bty="n", cex=0.8)
  dev.off()
}

#### Plot F/FMSY and SSB/SSBMSY over time ####
ylab = c("F/FMSY", "SSB/SSBMSY")
xlab = rep("Year", times=length(ylab))

plot_data_id <- which(names(om_list) %in% c("fratio", "ssbratio"))
col= c("black", "orange", "green", "red", "deepskyblue3")
legend = c("AMAK", "ASAP", "BAM", "MAS")

for (j in figure_id){
  jpeg(file=file.path(maindir, subdir, paste("Fig", j, "_true_estimate_f_ssbratio_overtime.jpg", sep="")), width=160, height=40, units="mm", res=300)
  
  par(mar=c(0.7, 4, 0.2, 0.2), mfrow=c(1,3), oma = c(4, 4, 0.2, 0.2))
  
  for(i in 1:length(plot_data_id)){
    ylim=range(om_list[[plot_data_id[i]]][,j], amak_list[[plot_data_id[i]]][,j], asap_list[[plot_data_id[i]]][,j], bam_list[[plot_data_id[i]]][,j], mas_list[[plot_data_id[i]]][,j])
    plot(sim1$yr, om_list[[plot_data_id[i]]][,j], pch=19, col=col[1], cex=0.5, ylim=ylim, axes=F, xlab="", ylab="")
    lines(sim1$yr, amak_list[[plot_data_id[i]]][,j], typ="l", col=col[2], lwd=1)
    lines(sim1$yr, asap_list[[plot_data_id[i]]][,j], typ="l", col=col[3], lwd=1)
    lines(sim1$yr, bam_list[[plot_data_id[i]]][,j], typ="l", col=col[4], lwd=1)
    lines(sim1$yr, mas_list[[plot_data_id[i]]][,j], typ="l", col=col[5], lwd=1)
    
    axis(1)
    mtext(side=1, text=xlab, line=2, cex=0.6, font=2, col="blue")
    mtext(side=2, text=ylab[i], line=2.5, cex=0.6, font=2, col="blue")
    axis(2, las=2)
    box()
  }
  plot.new()
  legend("topleft", c("OM", "AMAK", "ASAP", "BAM", "MAS"), pch=c(19, NA, NA, NA, NA), lty=c(NA,1,1,1,1), lwd=c(NA, 1.5, 1.5, 1.5, 1,5), col=col, bty="n", cex=0.8)
  dev.off()
}

#### Correlation matrix ####

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

title = paste("OM_", c("Biomass", "Abundance", "SSB", "Recruitment", "F", "Landings", "Survey Index", "Fratio", "SSB_ratio"), sep="")
plot_data_id <- c(1:7, 11:12)
cor_matrix <- matrix (NA, ncol=length(title), nrow=4)
colnames(cor_matrix) <- title
rownames(cor_matrix) <- c("AMAK", "ASAP", "BAM", "MAS")
p_matrix <- matrix (NA, ncol=length(title), nrow=4)
cor_list = list()
p_list = list()
for (j in figure_id){
  for (i in 1:length(plot_data_id)){
    cor_matrix[1,i] <- cor.test(om_list[[plot_data_id[i]]][,j], amak_list[[plot_data_id[i]]][,j])$estimate
    p_matrix[1,i] <- cor.test(om_list[[plot_data_id[i]]][,j], amak_list[[plot_data_id[i]]][,j])$p.value
    
    cor_matrix[2,i] <- cor.test(om_list[[plot_data_id[i]]][,j], asap_list[[plot_data_id[i]]][,j])$estimate
    p_matrix[2,i] <- cor.test(om_list[[plot_data_id[i]]][,j], asap_list[[plot_data_id[i]]][,j])$p.value
    
    cor_matrix[3,i] <- cor.test(om_list[[plot_data_id[i]]][,j], bam_list[[plot_data_id[i]]][,j])$estimate
    p_matrix[3,i] <- cor.test(om_list[[plot_data_id[i]]][,j], bam_list[[plot_data_id[i]]][,j])$p.value
    
    cor_matrix[4,i] <- cor.test(om_list[[plot_data_id[i]]][,j], mas_list[[plot_data_id[i]]][,j])$estimate
    p_matrix[4,i] <- cor.test(om_list[[plot_data_id[i]]][,j], mas_list[[plot_data_id[i]]][,j])$p.value
  }
  cor_list[[j]] <- cor_matrix
  p_list[[j]] <- p_matrix
  jpeg(file=file.path(maindir, subdir, paste("Fig", j, "_cormatrix.jpg", sep="")), width=150, height=100, units="mm", res=300)
  corrplot(cor_list[[j]],
           method="shade", # visualisation method
           shade.col=NA, # colour of shade line
           tl.col="black", # colour of text label
           tl.srt=45, # text label rotation
           #col=col(200), # colour of glyphs
           addCoef.col="white", # colour of coefficients
           order="original", # ordering method
           p.mat = p_list[[j]], sig.level = 0.01, insig = "blank",
           addgrid.col = "darkgray"
  )
  dev.off()
}


#### Relative error of point estimates ####
plot_data_id <- which(names(om_list) %in% c("abundance", "Ftot", "msy", "fmsy", "ssbmsy", "fratio", "ssbratio"))

plot_data <- as.data.frame(matrix(NA, ncol=4, nrow=length(plot_data_id)))
colnames(plot_data) <- c("AMAK", "ASAP", "BAM", "MAS")
rownames(plot_data) <- c("Mean Abundance", "Mean F", "MSY", "FMSY", "SSBMSY", "Fcur/FMSY", "SSBcur/SSBMSY") 
for(j in figure_id){
  for (i in 1:length(plot_data_id)){
    if(i %in% which(rownames(plot_data) %in% c("Fcur/FMSY", "SSBcur/SSBMSY"))){
      plot_data[i, ] <- pm_list[[j]][[plot_data_id[i]]]$re[sim1$yr[length(sim1$yr)], ]
    }
    else{
      if(i %in% which(rownames(plot_data) %in% c("Abundance", "F"))){
        plot_data[i, ] <- colMeans(pm_list[[j]][[plot_data_id[i]]]$re) 
      }
      else plot_data[i, ] <- pm_list[[j]][[plot_data_id[i]]]$re
    }
  }
  jpeg(file=file.path(maindir, subdir, paste("Fig", j, "_re_msy.jpg", sep="")), width=200, height=120, units="mm", res=300)
  par(mar=c(4,4,1,1))
  col= c("orange", "green", "red", "deepskyblue3")
  re_barplot <- barplot(t(plot_data), beside = T, ylab="Relative Error", col=col, width=0.5, cex.names = 0.7)
  box()
  legend("topleft", pch=c(15, 15, 15, 15), legend=colnames(plot_data), col=col, bty="n")
  dev.off()
}

#### Relative error of Fratio and SSB ratio over time ####
plot_data_id <- which(names(om_list) %in% c("fratio", "ssbratio"))
legend = c("AMAK", "ASAP", "BAM", "MAS")
for(j in figure_id){
  jpeg(file=file.path(maindir, subdir, paste("Fig", j, "_re_fratio.jpg", sep="")), width=170, height=110, units="mm", res=300)
  par(mar=c(0.7, 4, 0.2, 0.2), mfrow=c(2,2), oma = c(4, 4, 0.2, 0.2))
  for (i in 1:ncol(pm_list[[j]]$fratio$re)){
    re_bar <- barplot(pm_list[[j]]$fratio$re[,i], col=col[i], ylim=range(pm_list[[j]]$fratio$re))
    axis(side=1, at=re_bar[seq(1, par.sim1$nyr, by=4)], labels=F)
    box()
    if (i %in% c(3,4)) {
      axis(side=1, at=re_bar[seq(1, par.sim1$nyr, by=4)], labels = sim1$yr[seq(1, par.sim1$nyr, by=4)])
      mtext("Year", side=1, line=2, col="blue")
    }
    mtext("Relative Error", side=2, line=2, col="blue")
    legend("topleft", legend=legend[i], bty="n")
  }
  dev.off()
}

for(j in figure_id){
  jpeg(file=file.path(maindir, subdir, paste("Fig", j, "_re_ssbratio.jpg", sep="")), width=170, height=110, units="mm", res=300)
  par(mar=c(0.7, 4, 0.2, 0.2), mfrow=c(2,2), oma = c(4, 4, 0.2, 0.2))
  for (i in 1:ncol(pm_list[[j]]$ssbratio$re)){
    re_bar <- barplot(pm_list[[j]]$ssbratio$re[,i], col=col[i], ylim=range(pm_list[[j]]$ssbratio$re))
    axis(side=1, at=re_bar[seq(1, par.sim1$nyr, by=4)], labels=F)
    box()
    if (i %in% c(3,4)) {
      axis(side=1, at=re_bar[seq(1, par.sim1$nyr, by=4)], labels = sim1$yr[seq(1, par.sim1$nyr, by=4)])
      mtext("Year", side=1, line=2, col="blue")
    }
    mtext("Relative Error", side=2, line=2, col="blue")
    legend("topright", legend[i], bty="n")
  }
  dev.off()
}


#### Age composition ####
plot_data_id <- which(names(om_list) %in% c("agecomp"))
resid_data <- list()
for(j in figure_id){
  resid_data[[1]] <- (amak_list[[plot_data_id]][[j]] - om_list[[plot_data_id]][[j]])
  colnames(resid_data[[1]]) <- sim1$yr
  rownames(resid_data[[1]]) <- par.sim1$ages
  resid_data[[2]] <- (asap_list[[plot_data_id]][[j]] - om_list[[plot_data_id]][[j]])
  colnames(resid_data[[2]]) <- sim1$yr
  rownames(resid_data[[2]]) <- par.sim1$ages
  resid_data[[3]] <- (bam_list[[plot_data_id]][[j]] - om_list[[plot_data_id]][[j]])
  colnames(resid_data[[3]]) <- sim1$yr
  rownames(resid_data[[3]]) <- par.sim1$ages
  resid_data[[4]] <- (mas_list[[plot_data_id]][[j]] - om_list[[plot_data_id]][[j]])
  colnames(resid_data[[4]]) <- sim1$yr
  rownames(resid_data[[4]]) <- par.sim1$ages
  jpeg(file=file.path(maindir, subdir, paste("Fig", j, "_agecomp_resid.jpg", sep="")), width=170, height=110, units="mm", res=300)
  #par(mar=c(0.7, 4, 0.2, 0.2), mfrow=c(2,2), oma = c(0.2, 2, 2, 0.2))
  par(mfrow=c(2,2), mar=c(0, 0, 0, 0))
  for(i in 1:4){
    corrplot(resid_data[[i]], "circle", is.corr = F, tl.cex=0.6, tl.col="black", cl.cex = 0.6, cl.lim = range(resid_data[[1]], resid_data[[2]], resid_data[[3]], resid_data[[4]]))
  }
  dev.off()
}

agecom_last <- matrix(NA, ncol=5, nrow=length(par.sim1$ages))
col= c("black", "orange", "green", "red", "deepskyblue3")

for(j in figure_id){
  agecom_last[,1] <- om_list[[plot_data_id]][[j]][, sim1$yr[length(sim1$yr)]]
  agecom_last[,2] <- amak_list[[plot_data_id]][[j]][, sim1$yr[length(sim1$yr)]]
  agecom_last[,3] <- asap_list[[plot_data_id]][[j]][, sim1$yr[length(sim1$yr)]]
  agecom_last[,4] <- bam_list[[plot_data_id]][[j]][, sim1$yr[length(sim1$yr)]]
  agecom_last[,5] <- as.matrix(mas_list[[plot_data_id]][[j]][, sim1$yr[length(sim1$yr)]])
  ylim=range(agecom_last)
  jpeg(file=file.path(maindir, subdir, paste("Fig", j, "_agecomp_last.jpg", sep="")), width=120, height=100, units="mm", res=300)
  plot(agecom_last[,1], pch=19, col=col[1], ylim=ylim, ylab="Proportion", xlab="Age")
  sapply(2:ncol(agecom_last), function(x) lines(agecom_last[,x], type="l", col=col[x]))
  legend("topright", c("OM", "AMAK", "ASAP", "BAM", "MAS"), pch=c(19, NA, NA, NA, NA), lty=c(NA,1,1,1,1), lwd=c(NA, 1.5, 1.5, 1.5, 1,5), col=col, bty="n", cex=0.8, title=paste("Year", sim1$yr[length(sim1$yr)]))
  dev.off()
}

for(j in figure_id){
  jpeg(file=file.path(maindir, subdir, paste("Fig", j, "_agecomp_ae.jpg", sep="")), width=120, height=100, units="mm", res=300)
  ylim=range(ae_list[[j]]$ae)
  xlim=range(sim1$yr)
  plot(1, type="n", xlab="Year", ylab="Sum of Absolute Error", xlim=xlim, ylim=ylim, panel.first=grid(lty=1))
  sapply(1:ncol(ae_list[[j]]$ae), function(x) lines(sim1$yr, ae_list[[j]]$ae[, x], type="l", lty=1, col=col[x+1]))
  legend("top", c("OM", "AMAK", "ASAP", "BAM", "MAS"), pch=c(19, NA, NA, NA, NA), lty=c(NA,1,1,1,1), lwd=c(NA, 1.5, 1.5, 1.5, 1,5), col=col, bty="n", cex=0.8)
  dev.off()
}

#### AMAK residuals (boxplot) ####
jpeg(file=file.path(maindir, subdir, "amak_residuals_timeseries_boxplot.jpg"), width=170, height=160, units="mm", res=300)
mat=matrix(1:9, ncol=3, nrow=3, byrow = T)
layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
par(las=1, mar=c(4.5,4.5,1,0.5))

legend = c("Biomass (mt)", "Abundance (1000 fish)", "SSB (mt)", "Recruitment (1000 fish)", "F", "Landings (mt)", "Survey Index (scaled)", "F/FMSY", "SSB/SSBMSY")
xlab = rep("Year", times=length(legend))
ylab = rep("Residuals", times=length(legend))

plot_data_id <- c(1:7, 11:12)
for(i in 1:length(plot_data_id)){
  timeseries_residuals_plot(time=sim1$yr, om_data=om_list[[plot_data_id[i]]], est_data=amak_list[[plot_data_id[i]]], ylab="Residuals", xlab=xlab[i], legend=legend[i], type="boxplot")
}
dev.off()

#### ASAP residuals (boxplot) ####
jpeg(file=file.path(maindir, subdir, "asap_residuals_timeseries_boxplot.jpg"), width=170, height=160, units="mm", res=300)
mat=matrix(1:9, ncol=3, nrow=3, byrow = T)
layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
par(las=1, mar=c(4.5,4.5,1,0.5))

legend = c("Biomass (mt)", "Abundance (1000 fish)", "SSB (mt)", "Recruitment (1000 fish)", "F", "Landings (mt)", "Survey Index (scaled)", "F/FMSY", "SSB/SSBMSY")
xlab = rep("Year", times=length(legend))
ylab = rep("Residuals", times=length(legend))

plot_data_id <- c(1:7, 11:12)
for(i in 1:length(plot_data_id)){
  timeseries_residuals_plot(time=sim1$yr, om_data=om_list[[plot_data_id[i]]], est_data=asap_list[[plot_data_id[i]]], ylab="Residuals", xlab=xlab[i], legend=legend[i], type="boxplot")
}
dev.off()

#### BAM residuals (boxlplot) ####
jpeg(file=file.path(maindir, subdir, "bam_residuals_timeseries_boxplot.jpg"), width=170, height=160, units="mm", res=300)
mat=matrix(1:9, ncol=3, nrow=3, byrow = T)
layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
par(las=1, mar=c(4.5,4.5,1,0.5))

legend = c("Biomass (mt)", "Abundance (1000 fish)", "SSB (mt)", "Recruitment (1000 fish)", "F", "Landings (mt)", "Survey Index (scaled)", "F/FMSY", "SSB/SSBMSY")
xlab = rep("Year", times=length(legend))
ylab = rep("Residuals", times=length(legend))

plot_data_id <- c(1:7, 11:12)
for(i in 1:length(plot_data_id)){
  timeseries_residuals_plot(time=sim1$yr, om_data=om_list[[plot_data_id[i]]], est_data=bam_list[[plot_data_id[i]]], ylab="Residuals", xlab=xlab[i], legend=legend[i], type="boxplot")
}
dev.off()
#### MAS residuals (boxplot) ####
jpeg(file=file.path(maindir, subdir, "mas_residuals_timeseries_boxplot.jpg"), width=170, height=160, units="mm", res=300)
mat=matrix(1:9, ncol=3, nrow=3, byrow = T)
layout(mat=mat, widths=rep.int(1, ncol(mat)), heights=rep.int(1, nrow(mat)))
par(las=1, mar=c(4.5,4.5,1,0.5))

legend = c("Biomass (mt)", "Abundance (1000 fish)", "SSB (mt)", "Recruitment (1000 fish)", "F", "Landings (mt)", "Survey Index (scaled)", "F/FMSY", "SSB/SSBMSY")
xlab = rep("Year", times=length(legend))
ylab = rep("Residuals", times=length(legend))

plot_data_id <- c(1:7, 11:12)
for(i in 1:length(plot_data_id)){
  timeseries_residuals_plot(time=sim1$yr, om_data=om_list[[plot_data_id[i]]], est_data=mas_list[[plot_data_id[i]]], ylab="Residuals", xlab=xlab[i], legend=legend[i], type="boxplot")
}
dev.off()
#### stock status determination matrix ####
fratio_list <- list(amak_list$fratio, asap_list$fratio, bam_list$fratio, mas_list$fratio)
ssbratio_list <- list(amak_list$ssbratio, asap_list$ssbratio, bam_list$ssbratio, mas_list$ssbratio)

overfishing_temp <- list()
overfishing_performance<-matrix(NA, nrow=par.sim1$nyr, ncol=4)
for (k in 1:4){
  overfishing_temp[[k]]<-matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)
  for(j in 1:om_sim_num){
    for (i in 1:par.sim1$nyr){
      if((fratio_list[[k]][i,j]<1 & om_list$fratio[i,j]<1) | (fratio_list[[k]][i,j]>1 & om_list$fratio[i,j]>1)){
        overfishing_temp[[k]][i,j] = 1
      } else{
        overfishing_temp[[k]][i,j] = 0
      } 
    }
  }
  overfishing_performance[,k] <- apply(overfishing_temp[[k]], 1, sum)/om_sim_num*100
}
colnames(overfishing_performance) <- model_names[2:length(model_names)]
rownames(overfishing_performance) <- sim1$yr
write.csv(overfishing_performance, file=file.path(maindir, "Outputs", "overfishing_performance.csv"))

ssbratio_list <- list(amak_list$ssbratio*2, asap_list$ssbratio*2, bam_list$ssbratio*2, mas_list$ssbratio*2)
ssb_list <- list(amak_list$ssb, asap_list$ssb, bam_list$ssb, mas_list$ssb)

overfished_temp <- list()
overfished_performance<-matrix(NA, nrow=par.sim1$nyr, ncol=4)
for (k in 1:4){
  overfished_temp[[k]]<-matrix(NA, nrow=par.sim1$nyr, ncol=om_sim_num)
  for(j in 1:om_sim_num){
    for (i in 1:par.sim1$nyr){
      if((ssb_list[[k]][i,j]>ssbratio_list[[k]][i,j] & om_list$ssb[i,j]>(om_list$ssbratio[i,j]*2)) | (ssb_list[[k]][i,j]<ssbratio_list[[k]][i,j] & om_list$ssb[i,j]<(om_list$ssbratio[i,j]*2))){
        overfished_temp[[k]][i,j] = 1
      } else{
        overfished_temp[[k]][i,j] = 0
      } 
    }
  }
  overfished_performance[,k] <- apply(overfished_temp[[k]], 1, sum)/om_sim_num*100
}
colnames(overfished_performance) <- model_names[2:length(model_names)]
rownames(overfished_performance) <- sim1$yr
write.csv(overfished_performance, file=file.path(maindir, "Outputs", "overfished_performance.csv"))


jpeg(file=file.path(maindir, subdir, paste("Overfishing_percentage.jpg", sep="")), width=150, height=120, units="mm", res=300)
overfishing_barplot=barplot(overfishing_performance, beside = T, col=rep(col[2:length(col)], each=par.sim1$nyr), ylab="Percentage", xlab="Year", axes=F)
label_position = c(1, 10, 20, 30)
axis(1, at=as.vector(overfishing_barplot[label_position,]), label=rep(c(1, 10, 20, 30), times=4), las=2, cex.axis=0.5)
axis(2, cex.axis=0.6)
box()
dev.off()

jpeg(file=file.path(maindir, subdir, paste("Overfished_percentage.jpg", sep="")), width=150, height=120, units="mm", res=300)
overfished_barplot=barplot(overfished_performance, beside = T, col=rep(col[2:length(col)], each=par.sim1$nyr), ylab="Percentage", xlab="Year", axes = F)
label_position = c(1, 10, 20, 30)
axis(1, at=as.vector(overfished_barplot[label_position,]), label=rep(c(1, 10, 20, 30), times=4), las=2, cex.axis=0.5)
axis(2, cex.axis=0.6)
box()
dev.off()

#### MSY RE ####
msy_re=matrix(NA, ncol=4*3, nrow=om_sim_num)
for(i in 1:om_sim_num){
  msy_re[i,]<-c(as.numeric(as.character(pm_list[[i]]$msy$re)), as.numeric(as.character(pm_list[[i]]$fmsy$re)), as.numeric(as.character(pm_list[[i]]$ssbmsy$re)))
}

jpeg(file=file.path(maindir, subdir, "msy_re_sim.jpg"), width=150, height=120, units="mm", res=300)
par(mar = c(4, 8, 4, 2) + 0.1)
labels <- c("AMAK_MSY", "ASAP_MSY", "BAM_MSY", "MAS_MSY",
            "AMAK_FMSY", "ASAP_FMSY", "BAM_FMSY", "MAS_FMSY", 
            "AMAK_SSBMSY", "ASAP_SSBMSY", "BAM_SSBMSY", "MAS_SSBMSY")
msy_sim_boxplot <-boxplot(msy_re[,seq(dim(msy_re)[2],1)], horizontal=T, col=col[5:2], pch=16, cex=0.5, axes=F, at=c(1,2,3,4, 6,7,8,9, 11,12,13,14), xlab="Relative Error")

abline(v=0, col="coral3", lty=2)
box()
axis(1)
axis(2, at=c(1,2,3,4, 6,7,8,9, 11,12,13,14), labels = rev(labels), las=2)
dev.off()

#### Residuals over time plot by parameter ####
plot_data_id <- which(names(om_list) %in% c("biomass", "abundance", "ssb", "recruit" , "Ftot","landing","survey", "fratio", "ssbratio"))

title = c("Biomass", "Abundance", "SSB", "Recruitment", "F", "Landings", "Survey Index", "Fratio", "SSBratio")
legend = c("AMAK", "ASAP", "BAM", "MAS")
col= c("orange", "green", "red", "deepskyblue3")
light_col = rgb(t(col2rgb(col)/1.4), maxColorValue = 357)


for (i in 1:length(plot_data_id)){
  jpeg(file=file.path(maindir, subdir, paste(title[i], "_re_sim.jpg", sep="")), width=200, height=120, units="mm", res=300)
  par(mar=c(0.7, 4, 0.2, 0.2), mfrow=c(2,2), oma = c(4, 4, 0.2, 0.2))
  for(k in 1:4){
    temp <- matrix(NA, ncol=om_sim_num, nrow=par.sim1$nyr)
    re_quantile <- matrix(NA, ncol=par.sim1$nyr, nrow=5)
    re_quantile_list = list()
    
    for(j in 1:om_sim_num){
      temp[,j] <- pm_list[[j]][[plot_data_id[i]]]$re[,k]
    }
    re_quantile <- sapply(1:par.sim1$nyr, function(x) quantile(temp[x,], c(0.1, 0.25, 0.5, 0.75, 0.9)))
    plot(re_quantile[3,], type="l", col="white", ylim=c(-max(abs(re_quantile), 0.3), max(abs(re_quantile), 0.3)), axes=F, xlab="", ylab="Relative Error", panel.first=grid(lty=1))
    polygon(c(sim1$yr, rev(sim1$yr)), c(re_quantile[1,], rev(re_quantile[5,])), border=NA, col=col[k])
    polygon(c(sim1$yr, rev(sim1$yr)), c(re_quantile[2,], rev(re_quantile[4,])), border=NA, col=light_col[k])
    lines(re_quantile[3,], type="l", col="white")
    
    box()
    axis(side=1, at=seq(5, par.sim1$nyr, by=5), labels=F)
    axis(2)
    if(k %in% c(3,4)){
      axis(side=1, at=seq(5, par.sim1$nyr, by=5), labels = sim1$yr[seq(5, par.sim1$nyr, by=5)])
      mtext("Year", side=1, line=2)
    }
    legend("topleft", legend=paste(legend[k], ":", title[i], sep=""), bty="n")
    
  }
  dev.off()
}

#### Plot true value and estimations over time ####
plot_data_id <- which(names(om_list) %in% c("biomass", "abundance", "ssb", "recruit" , "Ftot","landing","survey", "fratio", "ssbratio"))

title = c("Biomass", "Abundance", "SSB", "Recruitment", "F", "Landings", "Survey Index", "Fratio", "SSBratio")
legend = c("AMAK", "ASAP", "BAM", "MAS")
for (i in 1:length(plot_data_id)){
  jpeg(file=file.path(maindir, subdir, paste(title[i], "_estimation_over_time_sim.jpg", sep="")), width=200, height=120, units="mm", res=300)
  par(mar=c(0.7, 4, 0.2, 0.2), mfrow=c(2,2), oma = c(4, 4, 0.2, 0.2))
  temp_list <- list(amak_list[[plot_data_id[i]]], asap_list[[plot_data_id[i]]], bam_list[[plot_data_id[i]]], mas_list[[plot_data_id[i]]])
  ylim=range(om_list[[plot_data_id[i]]], unlist(temp_list))
  for(k in 1:length(temp_list)){
    plot(om_list[[plot_data_id[i]]][,1], type="l", col=alpha("deepskyblue", 0.2), axes=F, xlab="", ylab=title[i], panel.first=grid(lty=1), ylim=ylim)
    sapply(1:om_sim_num, function(x) lines(temp_list[[k]][,x], col=alpha("gray10", 0.7)))
    sapply(1:om_sim_num, function(x) lines(om_list[[plot_data_id[i]]][,x], col=alpha("deepskyblue", 0.2), type="l", pch=19, cex=0.5))
    #lines(om_list[[plot_data_id[i]]][,1], type="o", col="blue", pch=19, cex=0.5)
    box()
    axis(side=1, at=seq(5, par.sim1$nyr, by=5), labels=F)
    axis(2)
    if(k %in% c(3,4)){
      axis(side=1, at=seq(5, par.sim1$nyr, by=5), labels = sim1$yr[seq(5, par.sim1$nyr, by=5)])
      mtext("Year", side=1, line=2)
    }
    legend("topleft", legend=legend[k], bty="n")
  }
  dev.off()
}

# for (i in 1:length(plot_data_id)){
#   jpeg(file=file.path(maindir, subdir, paste(title[i], "_estimation_over_time_sim.jpg", sep="")), width=200, height=120, units="mm", res=300)
#   par(mar=c(0.7, 4, 0.2, 0.2), mfrow=c(2,2), oma = c(4, 4, 0.2, 0.2))
#   temp_list <- list(amak_list[[plot_data_id[i]]], asap_list[[plot_data_id[i]]], bam_list[[plot_data_id[i]]], mas_list[[plot_data_id[i]]])
#   ylim=range(om_list[[plot_data_id[i]]], unlist(temp_list))
#   for(k in 1:length(temp_list)){
#     plot(om_list[[plot_data_id[i]]][,1], type="l", col="blue", axes=F, xlab="", ylab=title[i], panel.first=grid(lty=1), ylim=ylim)
#     sapply(1:om_sim_num, function(x) lines(temp_list[[k]][,x], col="gray70"))
#     lines(om_list[[plot_data_id[i]]][,1], type="o", col="blue", pch=19, cex=0.5)
#     box()
#     axis(side=1, at=seq(5, par.sim1$nyr, by=5), labels=F)
#     axis(2)
#     if(k %in% c(3,4)){
#       axis(side=1, at=seq(5, par.sim1$nyr, by=5), labels = sim1$yr[seq(5, par.sim1$nyr, by=5)])
#       mtext("Year", side=1, line=2)
#     }
#     legend("topleft", legend=legend[k], bty="n")
#   }
#   dev.off()
# }

#### Estimated parameters ####
#parameters <- c("logR", "Q", "sel_a50_fishery", "sel_slope_fishery", "sel_a50_survey", "sel_slope_survey")
parameters <- c("Rzero", "Q")
parameter_val_list <- list()
temp <- as.data.frame(matrix(NA, nrow=om_sim_num, ncol=2))
colnames(temp) <- parameters
parameter_re_list <- list()

for(k in 1:5){
  for (om_sim in 1:om_sim_num){
    if(k==1){
      load(file.path(maindir, "OM", paste("OM", keep_sim_id[om_sim], ".RData", sep="")))
      temp$Rzero[om_sim] <- par.sim1$R0/1000
      temp$Q[om_sim] <- dat.sim1$q
    }
    if(k==2){
      setwd(file.path(maindir, "AMAK", paste("s", keep_sim_id[om_sim], sep="")))
      amak_std <- readRep("amak", suffix = ".std")
      temp$Rzero[om_sim] <- exp(amak_std$value[which(amak_std$name=="log_Rzero")])/1000
      temp$Q[om_sim] <- exp(amak_std$value[which(amak_std$name=="log_q_ind[1]")])
    }
    if(k==3){
      asap_output <- dget(file.path(maindir, "ASAP", paste("s", keep_sim_id[om_sim], sep=""), "asap3.rdat"))
      temp$Rzero[om_sim] <- asap_output$SR.parms$SR.R0
      temp$Q[om_sim] <- asap_output$q.indices/1000
    }
    if(k==4){
      bam_output <- dget(file.path(maindir, "BAM", paste("s", keep_sim_id[om_sim], sep=""), "bam-sim.rdat"))
      temp$Rzero[om_sim] <- bam_output$parms$R0/1000
      temp$Q[om_sim] <- bam_output$parms$q.survey1
    }
    if(k==5){
      mas_output <- dget(file.path(maindir, "BAM", paste("s", keep_sim_id[om_sim], sep=""), "bam-sim.rdat"))
      temp$Rzero[om_sim] <- bam_output$parms$R0/1000
      temp$Q[om_sim] <- bam_output$parms$q.survey1
    }
  }
  parameter_val_list[[k]] <- temp
}

for(k in 1:4){
  parameter_re_list[[k]] <- (parameter_val_list[[k+1]]-parameter_val_list[[1]])/parameter_val_list[[1]]
}

jpeg(file=file.path(maindir, subdir, "Rzero_sim.jpg"), width=150, height=120, units="mm", res=300)
temp <- matrix(NA, nrow=om_sim_num, ncol=5)
for(k in 1:5){
  temp[,k] <- parameter_val_list[[k]][,1]
}

temp_boxplot <- boxplot(temp[,2:ncol(temp)], beside=T, axes=F, ylab=parameters[1], col=col, pch=16, ylim=range(temp))
box()
axis(1, at=1:4, labels = c("AMAK", "ASAP", "BAM", "MAS"))
axis(2)
abline(h=parameter_val_list[[1]][,1], col="gray70", lty=2)
dev.off()


jpeg(file=file.path(maindir, subdir, "Q_sim.jpg"), width=150, height=120, units="mm", res=300)
temp <- matrix(NA, nrow=om_sim_num, ncol=5)
for(k in 1:5){
  temp[,k] <- parameter_val_list[[k]][,2]
}

temp_boxplot <- boxplot(temp[,2:ncol(temp)], beside=T, axes=F, ylab=parameters[2], col=col, pch=16)
box()
axis(1, at=1:4, labels = c("AMAK", "ASAP", "BAM", "MAS"))
axis(2)
abline(h=parameter_val_list[[1]][,2], col="gray70", lty=2)
dev.off()

#### Last year condition ####
plot_data_id <- which(names(om_list) %in% c("biomass", "abundance", "ssb", "recruit", "Ftot", "fratio", "ssbratio"))

plot_data <- as.data.frame(matrix(NA, ncol=4, nrow=length(plot_data_id)))
colnames(plot_data) <- c("AMAK", "ASAP", "BAM", "MAS")
rownames(plot_data) <- c("Biomass (last year)", "Abundance (last year)", "SSB (last year)", "Recruit (last year)", "F (last year)", "Fratio (last year)", "SSBratio (last year)") 

for (i in 1:length(plot_data_id)){
  jpeg(file=file.path(maindir, subdir, paste(rownames(plot_data)[i], ".jpg", sep="")), width=150, height=120, units="mm", res=300)
  boxplot(amak_list[[plot_data_id[i]]][par.sim1$nyr,], asap_list[[plot_data_id[i]]][par.sim1$nyr,], bam_list[[plot_data_id[i]]][par.sim1$nyr,], mas_list[[plot_data_id[i]]][par.sim1$nyr,], col=col, axes=F, ylab=rownames(plot_data)[i], pch=16)
  box()
  axis(1, at=1:4, labels = c("AMAK", "ASAP", "BAM", "MAS"))
  axis(2)
  abline(h=om_list[[plot_data_id[i]]][par.sim1$nyr,], col="gray70", lty=2)
  dev.off()
}

#### Plot survey and landings with error as true values ####
plot_data_id <- which(names(om_list) %in% c("landing", "survey"))

title = c("Landings", "Survey Index")
legend = c("AMAK", "ASAP", "BAM", "MAS")

for (i in 1:length(plot_data_id)){
  jpeg(file=file.path(maindir, subdir, paste(title[i], "_with_observation_error.jpg", sep="")), width=200, height=120, units="mm", res=300)
  par(mar=c(0.7, 4, 0.2, 0.2), mfrow=c(2,2), oma = c(4, 4, 0.2, 0.2))
  temp_list <- list(amak_list[[plot_data_id[i]]], asap_list[[plot_data_id[i]]], bam_list[[plot_data_id[i]]], mas_list[[plot_data_id[i]]])
  ylim=range(om_list[[plot_data_id[i]+8]], unlist(temp_list))
  for(k in 1:length(temp_list)){
    plot(om_list[[plot_data_id[i]+8]][,1], type="l", col=alpha("deepskyblue", 0.2), axes=F, xlab="", ylab=title[i], panel.first=grid(lty=1), ylim=ylim)
    sapply(1:om_sim_num, function(x) lines(temp_list[[k]][,x], col=alpha("gray10", 0.7)))
    sapply(1:om_sim_num, function(x) lines(om_list[[plot_data_id[i]+8]][,x], col=alpha("deepskyblue", 0.2), type="l", pch=19, cex=0.5))
    #lines(om_list[[plot_data_id[i]+8]][,1], type="o", col="blue", pch=19, cex=0.5)
    box()
    axis(side=1, at=seq(5, par.sim1$nyr, by=5), labels=F)
    axis(2)
    if(k %in% c(3,4)){
      axis(side=1, at=seq(5, par.sim1$nyr, by=5), labels = sim1$yr[seq(5, par.sim1$nyr, by=5)])
      mtext("Year", side=1, line=2)
    }
    legend("topleft", legend=legend[k], bty="n")
  }
  dev.off()
}


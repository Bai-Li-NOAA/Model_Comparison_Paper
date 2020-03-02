library(rstudioapi) #for obtaining current script location
library(gdata) #for write.fwf function from the SS_writectl_3.30 function
library(PBSadmb)
library(stringr)
library(matrixcalc) # for calculate positive definite Hessian

current_path <- getActiveDocumentContext()$path
maindir <- dirname(current_path)

setwd(maindir)
folder_names <- c("OM", "AMAK", "ASAP", "BAM", "SS", "Figures", "Outputs")
om_sim_num <- 160
keep_sim_num <- 100
# om_sim_num <- 160
# keep_sim_num <- 100
# if (om_sim_num==1){
#   seed_num <- 9924
#   save(seed_num, file=file.path(maindir, "seed_num.RData"))
# } else {
#   seed_num <- sample(100000, om_sim_num, replace=F)
#   save(seed_num, file=file.path(maindir, "seed_num.RData"))
# }
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
#### Run OM and EM ####
rec_dev_matrix <- recruit_dev_case(rec_dev_change=T, nyr=nyr, logR.sd=logR.sd, om_sim_num=om_sim_num)  
F <- f_case(f_pattern=1, f_min=0.01, f_max=0.39, f_sd=0.2, nyr=nyr)

OM_Run(maindir=maindir, subdir="OM", input.cv.L=input.cv.L, input.cv.survey=input.cv.survey, cv.L=cv.L, cv.survey=cv.survey, n.L=n.L, n.survey=n.survey, logR.sd=logR.sd, nyr= nyr, om_sim_num=om_sim_num, seed_num=seed_num, rec_dev_matrix=rec_dev_matrix, F=F)
rm(list=setdiff(ls(), c("maindir", "om_sim_num", "keep_sim_num", lsf.str())))

library(future)
plan(multiprocess)
asapjob %<-% ASAP_Run(maindir=maindir, subdir="ASAP", om_sim_num=om_sim_num)
amakjob %<-% AMAK_Run(maindir=maindir, subdir="AMAK", om_sim_num=om_sim_num)
bamjob %<-% BAM_Run(maindir=maindir, subdir="BAM", om_sim_num=om_sim_num)
ssjob %<-% SS_Run(maindir=maindir, subdir="SS", om_sim_num=om_sim_num)

AMAK_Run(maindir=maindir, subdir="AMAK", om_sim_num=om_sim_num)
rm(list=setdiff(ls(), c("maindir", "om_sim_num", "keep_sim_num", lsf.str())))

ASAP_Run(maindir=maindir, subdir="ASAP", om_sim_num=om_sim_num)
rm(list=setdiff(ls(), c("maindir", "om_sim_num", "keep_sim_num", lsf.str())))

BAM_Run(maindir=maindir, subdir="BAM", om_sim_num=om_sim_num)
rm(list=setdiff(ls(), c("maindir", "om_sim_num", "keep_sim_num", lsf.str())))

SS_Run(maindir=maindir, subdir="SS", om_sim_num=om_sim_num)
rm(list=setdiff(ls(), c("maindir", "om_sim_num", "keep_sim_num", lsf.str())))

sapply(1:om_sim_num, function(x) file.remove(file.path(maindir, "AMAK", paste("s", x, sep=""), "amak.exe")))
sapply(1:om_sim_num, function(x) file.remove(file.path(maindir, "ASAP", paste("s", x, sep=""), "ASAP3.exe")))
sapply(1:om_sim_num, function(x) file.remove(file.path(maindir, "BAM", paste("s", x, sep=""), "BAM-Sim.exe")))
sapply(1:om_sim_num, function(x) file.remove(file.path(maindir, "SS", paste("s", x, sep=""), "ss.exe")))
#### Check convergence of iterations ####
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
  
  subdir="SS"
  setwd(file.path(maindir, subdir, paste("s", om_sim, sep="")))
  temp <- getADMBHessian(File=file.path(maindir, subdir, paste("s", om_sim, sep="")), FileName="admodel.hes")$hes
  convergence_measures$positive_hessian[om_sim, 4] <- ifelse(is.positive.definite((temp+t(temp))/2), 1, 0)
  temp <- as.numeric(scan("ss.par", what='', n=16, quiet=TRUE)[c(6,11,16)])[3]
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
  legend("topright", "SS", bty="n")
  box()
  dev.off()
}


good_gradient_percentage <- matrix(NA, nrow=1, ncol=4)
good_gradient_percentage[,1] <- length(which(convergence_measures$gradient[,1]<=0.001))/om_sim_num*100
good_gradient_percentage[,2] <- length(which(convergence_measures$gradient[,2]<=0.001))/om_sim_num*100
good_gradient_percentage[,3] <- length(which(convergence_measures$gradient[,3]<=0.001))/om_sim_num*100
good_gradient_percentage[,4] <- length(which(convergence_measures$gradient[,4]<=0.001))/om_sim_num*100
colnames(good_gradient_percentage) <- c("AMAK", "ASAP", "BAM", "SS")
rownames(good_gradient_percentage) <- "Percentage"
write.csv(good_gradient_percentage, file=file.path(maindir, "Outputs", "good_gradient_ratio.csv"))

keep_sim_id <- c(1:om_sim_num)[-unique(c(which(convergence_measures$gradient[,1] %in% boxplot.stats(convergence_measures$gradient[,1])$out), which(convergence_measures$gradient[,2] %in% boxplot.stats(convergence_measures$gradient[,2])$out), which(convergence_measures$gradient[,3] %in% boxplot.stats(convergence_measures$gradient[,3])$out), which(convergence_measures$gradient[,4] %in% boxplot.stats(convergence_measures$gradient[,4])$out), unique(unlist(sapply(1:ncol(convergence_measures$positive_hessian), function(x) which(convergence_measures$positive_hessian[,x]==0))))))][1:keep_sim_num]

if(max(convergence_measures$gradient[keep_sim_id,])<0.1){
  jpeg(file=file.path(maindir, "Figures", "Gradient_no_outliers.jpg"), width=170, height=150, units="mm", res=300)
  par(mfrow=c(2,2), mar=c(4,4,1,1))
  xlim = c(0, max(convergence_measures$gradient[keep_sim_id,]))
  bins <- seq(0, max(convergence_measures$gradient[keep_sim_id,])*1.1, by=0.0005)
  hist(convergence_measures$gradient[keep_sim_id,1], xlim=xlim, xlab = "Gradient", main="", col=col[2], breaks = bins)
  legend("topright", "AMAK", bty="n")
  box()
  
  hist(convergence_measures$gradient[keep_sim_id,2], xlim=xlim, xlab = "Gradient", main="", breaks = bins, col=col[3])
  legend("topright", "ASAP", bty="n")
  box()
  
  hist(convergence_measures$gradient[keep_sim_id,3], xlim=xlim, xlab = "Gradient", main="", breaks = bins, col=col[4])
  legend("topright", "BAM", bty="n")
  box()
  
  hist(convergence_measures$gradient[keep_sim_id,4], xlim=xlim, xlab = "Gradient", main="", breaks = bins, col=col[5])
  legend("topright", "SS", bty="n")
  box()
  dev.off()
}

om_sim_num <- length(keep_sim_id)


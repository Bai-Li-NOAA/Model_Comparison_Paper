#### Check and install missing packages and library all packages ####
list_of_packages <- c("rstudioapi", "gdata", "PBSadmb", "stringr", "matrixcalc", "r4ss", "ASAPplots", "future", "readxl", "scales", "corrplot", "future", "glue")
missing_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if(length(missing_packages)) install.packages(missing_packages)
invisible(lapply(list_of_packages, library, character.only = TRUE))

#### Set working directory ####
current_path <- getActiveDocumentContext()$path
maindir <- dirname(current_path)
setwd(maindir)

#### Source all functions ####
subdir = "functions"
setwd(file.path(maindir, subdir))
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)

#### Preliminary settings ####
## Set required working folders (output and figure)
invisible(sapply(c("output", "figure"), function(x) {
  if (!file.exists(file.path(maindir, x))) dir.create(file.path(maindir, x))
}))

em_names <- c("AMAK", "ASAP", "BAM", "SS") # Set estimation models
invisible(sapply(c(em_names, "OM"), function(x) {
  if (!file.exists(file.path(maindir, "output", x))) dir.create(file.path(maindir, "output", x))
}))
em_parfile_names <- c("amak.par", "asap3.par", "bam-sim.par", "ss.par")

## Set number of iterations per case
om_sim_num <- 165 # total iterations per case
keep_sim_num <- 100 # number of kept iterations per case
figure_number <- 10 # number of individual iteration to plot

## Set seeds for reproduciable results
set.seed(9924)

## Set recruitment and fishing mortality
nyr = 30
logR_sd=0.2
logf_sd=0.2

f_dev_matrix <- f_dev_case(f_dev_change=FALSE, nyr=nyr, logf_sd=logf_sd, om_sim_num=om_sim_num)

f_matrix <- f_case(f_pattern=6, f_min=NULL, f_max=NULL, f_end=NULL, f_mean=0.19, start_val=0.01, start_year=6, nyr=nyr, om_sim_num=om_sim_num, f_dev_matrix=f_dev_matrix)

r_dev_matrix <- r_dev_case(r_dev_change=TRUE, nyr=nyr, logR_sd=logR_sd, om_sim_num=om_sim_num)

#### Run OM ####
OM_Run(maindir=maindir, subdir="OM", om_sim_num=om_sim_num, r_dev_matrix=r_dev_matrix, f_dev_matrix=f_dev_matrix, f_matrix=f_matrix)
rm(list=setdiff(ls(), c("maindir", "om_sim_num", "keep_sim_num", "figure_number", "em_names", lsf.str())))

#### Run EM ####
# plan(multiprocess)
# asapjob %<-% ASAP_Run(maindir=maindir, subdir="ASAP", om_sim_num=om_sim_num)
# amakjob %<-% AMAK_Run(maindir=maindir, subdir="AMAK", om_sim_num=om_sim_num)
# bamjob %<-% BAM_Run(maindir=maindir, subdir="BAM", om_sim_num=om_sim_num)
# ssjob %<-% SS_Run(maindir=maindir, subdir="SS", om_sim_num=om_sim_num)

AMAK_Run(maindir=maindir, subdir="AMAK", om_sim_num=om_sim_num)
rm(list=setdiff(ls(), c("maindir", "om_sim_num", "keep_sim_num", "figure_number", "em_names", lsf.str())))

ASAP_Run(maindir=maindir, subdir="ASAP", om_sim_num=om_sim_num)
rm(list=setdiff(ls(), c("maindir", "om_sim_num", "keep_sim_num", "figure_number", "em_names", lsf.str())))

BAM_Run(maindir=maindir, subdir="BAM", om_sim_num=om_sim_num)
rm(list=setdiff(ls(), c("maindir", "om_sim_num", "keep_sim_num", "figure_number", "em_names", lsf.str())))

SS_Run(maindir=maindir, subdir="SS", om_sim_num=om_sim_num)
rm(list=setdiff(ls(), c("maindir", "om_sim_num", "keep_sim_num", "figure_number", "em_names", lsf.str())))

#### Remove .exe to save space ####
sapply(1:om_sim_num, function(x) file.remove(file.path(maindir, "output", "AMAK", paste("s", x, sep=""), "amak.exe")))
sapply(1:om_sim_num, function(x) file.remove(file.path(maindir, "output", "ASAP", paste("s", x, sep=""), "ASAP3.exe")))
sapply(1:om_sim_num, function(x) file.remove(file.path(maindir, "output", "BAM", paste("s", x, sep=""), "BAM-Sim.exe")))
sapply(1:om_sim_num, function(x) file.remove(file.path(maindir, "output", "SS", paste("s", x, sep=""), "ss.exe")))

#### Convergence check ####
positive_hessian = matrix(NA, ncol=length(em_names), nrow=om_sim_num)
gradient = matrix(NA, ncol=length(em_names), nrow=om_sim_num)
convergence_measures <- list(positive_hessian=positive_hessian, gradient=gradient)

for (om_sim in 1:om_sim_num){
  for (em_id in 1:length(em_names)){
    subdir=em_names[em_id]
    setwd(file.path(maindir, "output", subdir, paste("s", om_sim, sep="")))
    if (file.exists(file.path(maindir, "output",  subdir, paste("s", om_sim, sep=""), "admodel.hes"))) {
      temp <- getADMBHessian(File=file.path(maindir, "output",  subdir, paste("s", om_sim, sep="")), FileName="admodel.hes")$hes
      convergence_measures$positive_hessian[om_sim, em_id] <- ifelse(is.positive.definite((temp+t(temp))/2), 1, 0)
      temp <- as.numeric(scan(em_parfile_names[em_id], what='', n=16, quiet=TRUE)[c(6,11,16)])[3]
    convergence_measures$gradient[om_sim, em_id] <- temp
    } else {
      convergence_measures$positive_hessian[om_sim, em_id] <- 0
      convergence_measures$gradient[om_sim, em_id] <- 0
    }
  }
}

save(convergence_measures, file=file.path(maindir, "output", "convergence_measures.RData"))
load(file.path(maindir, "output", "convergence_measures.RData"))

not_positive_hessian <- unique(unlist(sapply(1:ncol(convergence_measures$positive_hessian), function(x) which(convergence_measures$positive_hessian[,x]==0))))
write.csv(not_positive_hessian, file=file.path(maindir, "output", "Not_positive_hessian.csv"))

col= c("black", "orange", "green", "red", "deepskyblue3")
if(max(convergence_measures$gradient)<0.1){
  jpeg(file=file.path(maindir, "figure", "Gradient.jpg"), width=170, height=150, units="mm", res=300)
  par(mfrow=c(2,2), mar=c(4,4,1,1))
  xlim = c(0, max(convergence_measures$gradient))
  bins <- seq(0, max(convergence_measures$gradient)*1.05, by=0.0005)

  for (em_id in 1:length(em_names)){
    hist(convergence_measures$gradient[,em_id], xlim=xlim, xlab = "Gradient", main="", col=col[em_id+1], breaks = bins)
  legend("topright", em_names[em_id], bty="n")
  box()
  }
  dev.off()
}


good_gradient_percentage <- matrix(NA, nrow=1, ncol=length(em_names))
sapply(1:length(em_names), function(x) good_gradient_percentage[,x] <<- length(which(convergence_measures$gradient[,x]<=0.002))/om_sim_num*100)
colnames(good_gradient_percentage) <- em_names
rownames(good_gradient_percentage) <- "Percentage"
write.csv(good_gradient_percentage, file=file.path(maindir, "output", "good_gradient_ratio.csv"))

keep_sim_id <- c(1:om_sim_num)[-unique(c(
  which(convergence_measures$gradient[,1] %in% boxplot.stats(convergence_measures$gradient[,1])$out),
  which(convergence_measures$gradient[,2] %in% boxplot.stats(convergence_measures$gradient[,2])$out),
  which(convergence_measures$gradient[,3] %in% boxplot.stats(convergence_measures$gradient[,3])$out),
  which(convergence_measures$gradient[,4] %in% boxplot.stats(convergence_measures$gradient[,4])$out),
  unique(unlist(sapply(1:ncol(convergence_measures$positive_hessian), function(x) which(convergence_measures$positive_hessian[,x]==0))))))][1:keep_sim_num]

if(max(convergence_measures$gradient[keep_sim_id,])<0.1){
  jpeg(file=file.path(maindir, "figure", "Gradient_no_outliers.jpg"), width=170, height=150, units="mm", res=300)
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
save(keep_sim_id, om_sim_num, file=file.path(maindir, "output", "keep_sim_id.RData"))
#### Plot output ####
source(file.path(maindir, "plot_output.R"))

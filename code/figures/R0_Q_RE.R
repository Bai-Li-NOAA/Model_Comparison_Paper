
maindir_list  <- c("C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case0",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case1",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case2",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case3",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case4",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case5",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case6",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case7",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case8",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case9",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case10",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case11",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case0_F0_0.01_1yr"
                   )

em_num <- 4

#### R0 and Q median relative error table ####
R0_RE <- matrix(NA, nrow=length(maindir_list), ncol=em_num)
Q_RE <- matrix(NA, nrow=length(maindir_list), ncol=em_num)
#BC = exp(0.6^2/2)
median_case_id <- 11
mean_case_id <- 12

for (j in 1:length(maindir_list)){
  maindir <- maindir_list[j]
  setwd(file.path(maindir, "output"))
  load("./keep_sim_id.RData")

  parameters <- c("Rzero", "Q")
  parameter_val_list <- list()
  temp <- as.data.frame(matrix(NA, nrow=om_sim_num, ncol=2))
  colnames(temp) <- parameters
  parameter_re_list <- list()

  for(k in 1:(em_num+1)){
    for (om_sim in 1:om_sim_num){
      if(k==1){
        load(file.path(maindir, "output", "OM", paste("OM", keep_sim_id[om_sim], ".RData", sep="")))
        temp$Rzero[om_sim] <- om_input$R0/1000
        temp$Q[om_sim] <- om_output$survey_q$survey1
      }
      if(k==2){
        setwd(file.path(maindir,"output", "AMAK", paste("s", keep_sim_id[om_sim], sep="")))
        amak_output <- readRep("For_R", suffix = ".rep")
        amak_std <- readRep("amak", suffix = ".std")
        if (j == mean_case_id) {
          temp$Rzero[om_sim] <- convertSRparms(R0=exp(amak_std$value[which(amak_std$name=="log_Rzero")]), h=amak_output$Steep[2], phi=amak_output$phizero, sigmaR=om_input$logR_sd, mean2med=FALSE)$R0BC/1000
        } else {
          temp$Rzero[om_sim] <- exp(amak_std$value[which(amak_std$name=="log_Rzero")])/1000
        }
        temp$Q[om_sim] <- exp(amak_std$value[which(amak_std$name=="log_q_ind[1]")])
      }
      if(k==3){
        asap_output <- dget(file.path(maindir, "output", "ASAP", paste("s", keep_sim_id[om_sim], sep=""), "asap3.rdat"))
        if (j == mean_case_id) {
          temp$Rzero[om_sim] <- convertSRparms(R0=asap_output$SR.parms$SR.R0, h=asap_output$SR.parms$SR.steepness, phi=asap_output$SR.parms$SR.SPR0, sigmaR=om_input$logR_sd, mean2med=FALSE)$R0BC
        } else {
          temp$Rzero[om_sim] <- asap_output$SR.parms$SR.R0
        }
        temp$Q[om_sim] <- asap_output$q.indices/1000
      }
      if(k==4){
        bam_output <- dget(file.path(maindir, "output", "BAM", paste("s", keep_sim_id[om_sim], sep=""), "bam-sim.rdat"))
        if (j == mean_case_id) temp$Rzero[om_sim] <- bam_output$parms$R.virgin.bc/1000 else temp$Rzero[om_sim] <- bam_output$parms$R0/1000
        temp$Q[om_sim] <- bam_output$parms$q.survey1
      }
      if(k==5){
        ss_output <- SS_output(dir=file.path(maindir, "output", "SS", paste("s", keep_sim_id[om_sim], sep="")), ncols = 300, verbose=F, printstats=F)
        if (j == median_case_id) {
          temp$Rzero[om_sim] <- convertSRparms(R0=exp(ss_output$estimated_non_dev_parameters$Value[which(rownames(ss_output$estimated_non_dev_parameters)=="SR_LN(R0)")]),
                                               h=ss_output$parameters$Value[which(rownames(ss_output$parameters)=="SR_BH_steep")],
                                               phi=ss_output$sprseries$SSBzero[1]/exp(ss_output$estimated_non_dev_parameters$Value[which(rownames(ss_output$estimated_non_dev_parameters)=="SR_LN(R0)")]),
                                               sigmaR=om_input$logR_sd,
                                               mean2med=TRUE)$R0BC
        } else {
          temp$Rzero[om_sim] <- exp(ss_output$estimated_non_dev_parameters$Value[which(rownames(ss_output$estimated_non_dev_parameters)=="SR_LN(R0)")])
        }
        temp$Q[om_sim] <- exp(ss_output$estimated_non_dev_parameters$Value[which(rownames(ss_output$estimated_non_dev_parameters)=="LnQ_base_SURVEY1(2)")])/1000
      }
    }
    parameter_val_list[[k]] <- temp
  }

  for(k in 1:em_num){
    parameter_re_list[[k]] <- (parameter_val_list[[k+1]]-parameter_val_list[[1]])/parameter_val_list[[1]]
  }

  for(k in 1:em_num){
    R0_RE[j,k] <- median(parameter_re_list[[k]]$Rzero)
    Q_RE[j,k] <- median(parameter_re_list[[k]]$Q)
    # R0_RE[j,k] <- mean(parameter_re_list[[k]]$Rzero)
    # Q_RE[j,k] <- mean(parameter_re_list[[k]]$Q)
  }
}

write.csv(R0_RE, file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/R0_RE.csv", row.names = F)
write.csv(Q_RE, file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/Q_RE.csv", row.names = F)

#### R0 and Q absolute median relative error table ####
R0_RE <- matrix(NA, nrow=length(maindir_list), ncol=em_num)
Q_RE <- matrix(NA, nrow=length(maindir_list), ncol=em_num)
#BC = exp(0.6^2/2)
median_case_id <- 11
mean_case_id <- 12

for (j in 1:length(maindir_list)){
  maindir <- maindir_list[j]
  setwd(file.path(maindir, "output"))
  load("./keep_sim_id.RData")

  parameters <- c("Rzero", "Q")
  parameter_val_list <- list()
  temp <- as.data.frame(matrix(NA, nrow=om_sim_num, ncol=2))
  colnames(temp) <- parameters
  parameter_re_list <- list()

  for(k in 1:(em_num+1)){
    for (om_sim in 1:om_sim_num){
      if(k==1){
        load(file.path(maindir, "output", "OM", paste("OM", keep_sim_id[om_sim], ".RData", sep="")))
        temp$Rzero[om_sim] <- om_input$R0/1000
        temp$Q[om_sim] <- om_output$survey_q$survey1
      }
      if(k==2){
        setwd(file.path(maindir,"output", "AMAK", paste("s", keep_sim_id[om_sim], sep="")))
        amak_output <- readRep("For_R", suffix = ".rep")
        amak_std <- readRep("amak", suffix = ".std")
        if (j == mean_case_id) {
          temp$Rzero[om_sim] <- convertSRparms(R0=exp(amak_std$value[which(amak_std$name=="log_Rzero")]), h=amak_output$Steep[2], phi=amak_output$phizero, sigmaR=om_input$logR_sd, mean2med=FALSE)$R0BC/1000
        } else {
          temp$Rzero[om_sim] <- exp(amak_std$value[which(amak_std$name=="log_Rzero")])/1000
        }
        temp$Q[om_sim] <- exp(amak_std$value[which(amak_std$name=="log_q_ind[1]")])
      }
      if(k==3){
        asap_output <- dget(file.path(maindir, "output", "ASAP", paste("s", keep_sim_id[om_sim], sep=""), "asap3.rdat"))
        if (j == mean_case_id) {
          temp$Rzero[om_sim] <- convertSRparms(R0=asap_output$SR.parms$SR.R0, h=asap_output$SR.parms$SR.steepness, phi=asap_output$SR.parms$SR.SPR0, sigmaR=om_input$logR_sd, mean2med=FALSE)$R0BC
        } else {
          temp$Rzero[om_sim] <- asap_output$SR.parms$SR.R0
        }
        temp$Q[om_sim] <- asap_output$q.indices/1000
      }
      if(k==4){
        bam_output <- dget(file.path(maindir, "output", "BAM", paste("s", keep_sim_id[om_sim], sep=""), "bam-sim.rdat"))
        if (j == mean_case_id) temp$Rzero[om_sim] <- bam_output$parms$R.virgin.bc/1000 else temp$Rzero[om_sim] <- bam_output$parms$R0/1000
        temp$Q[om_sim] <- bam_output$parms$q.survey1
      }
      if(k==5){
        ss_output <- SS_output(dir=file.path(maindir, "output", "SS", paste("s", keep_sim_id[om_sim], sep="")), ncols = 300, verbose=F, printstats=F)
        if (j == median_case_id) {
          temp$Rzero[om_sim] <- convertSRparms(R0=exp(ss_output$estimated_non_dev_parameters$Value[which(rownames(ss_output$estimated_non_dev_parameters)=="SR_LN(R0)")]),
                                               h=ss_output$parameters$Value[which(rownames(ss_output$parameters)=="SR_BH_steep")],
                                               phi=ss_output$sprseries$SSBzero[1]/exp(ss_output$estimated_non_dev_parameters$Value[which(rownames(ss_output$estimated_non_dev_parameters)=="SR_LN(R0)")]),
                                               sigmaR=om_input$logR_sd,
                                               mean2med=TRUE)$R0BC
        } else {
          temp$Rzero[om_sim] <- exp(ss_output$estimated_non_dev_parameters$Value[which(rownames(ss_output$estimated_non_dev_parameters)=="SR_LN(R0)")])
        }
        temp$Q[om_sim] <- exp(ss_output$estimated_non_dev_parameters$Value[which(rownames(ss_output$estimated_non_dev_parameters)=="LnQ_base_SURVEY1(2)")])/1000
      }
    }
    parameter_val_list[[k]] <- temp
  }

  for(k in 1:em_num){
    parameter_re_list[[k]] <- (parameter_val_list[[k+1]]-parameter_val_list[[1]])/parameter_val_list[[1]]
  }

  for(k in 1:em_num){
    R0_RE[j,k] <- median(abs(parameter_re_list[[k]]$Rzero))
    Q_RE[j,k] <- median(abs(parameter_re_list[[k]]$Q))
    # R0_RE[j,k] <- mean(parameter_re_list[[k]]$Rzero)
    # Q_RE[j,k] <- mean(parameter_re_list[[k]]$Q)
  }
}

write.csv(R0_RE, file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/R0_MARE.csv", row.names = F)
write.csv(Q_RE, file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/Q_MARE.csv", row.names = F)

#### R0 and Q absolute median relative error table ####
R0_RE <- matrix(NA, nrow=length(maindir_list), ncol=em_num)
logQ_RE <- matrix(NA, nrow=length(maindir_list), ncol=em_num)
#BC = exp(0.6^2/2)
median_case_id <- 11
mean_case_id <- 12

for (j in 1:length(maindir_list)){
  maindir <- maindir_list[j]
  setwd(file.path(maindir, "output"))
  load("./keep_sim_id.RData")

  parameters <- c("Rzero", "Q")
  parameter_val_list <- list()
  temp <- as.data.frame(matrix(NA, nrow=om_sim_num, ncol=2))
  colnames(temp) <- parameters
  parameter_re_list <- list()

  for(k in 1:(em_num+1)){
    for (om_sim in 1:om_sim_num){
      if(k==1){
        load(file.path(maindir, "output", "OM", paste("OM", keep_sim_id[om_sim], ".RData", sep="")))
        temp$Rzero[om_sim] <- om_input$R0/1000
        temp$Q[om_sim] <- log(om_output$survey_q$survey1)
      }
      if(k==2){
        setwd(file.path(maindir,"output", "AMAK", paste("s", keep_sim_id[om_sim], sep="")))
        amak_output <- readRep("For_R", suffix = ".rep")
        amak_std <- readRep("amak", suffix = ".std")
        if (j == mean_case_id) {
          temp$Rzero[om_sim] <- convertSRparms(R0=exp(amak_std$value[which(amak_std$name=="log_Rzero")]), h=amak_output$Steep[2], phi=amak_output$phizero, sigmaR=om_input$logR_sd, mean2med=FALSE)$R0BC/1000
        } else {
          temp$Rzero[om_sim] <- exp(amak_std$value[which(amak_std$name=="log_Rzero")])/1000
        }
        temp$Q[om_sim] <- amak_std$value[which(amak_std$name=="log_q_ind[1]")]
      }
      if(k==3){
        asap_output <- dget(file.path(maindir, "output", "ASAP", paste("s", keep_sim_id[om_sim], sep=""), "asap3.rdat"))
        if (j == mean_case_id) {
          temp$Rzero[om_sim] <- convertSRparms(R0=asap_output$SR.parms$SR.R0, h=asap_output$SR.parms$SR.steepness, phi=asap_output$SR.parms$SR.SPR0, sigmaR=om_input$logR_sd, mean2med=FALSE)$R0BC
        } else {
          temp$Rzero[om_sim] <- asap_output$SR.parms$SR.R0
        }
        temp$Q[om_sim] <- log(asap_output$q.indices/1000)
      }
      if(k==4){
        bam_output <- dget(file.path(maindir, "output", "BAM", paste("s", keep_sim_id[om_sim], sep=""), "bam-sim.rdat"))
        if (j == mean_case_id) temp$Rzero[om_sim] <- bam_output$parms$R.virgin.bc/1000 else temp$Rzero[om_sim] <- bam_output$parms$R0/1000
        temp$Q[om_sim] <- log(bam_output$parms$q.survey1)
      }
      if(k==5){
        ss_output <- SS_output(dir=file.path(maindir, "output", "SS", paste("s", keep_sim_id[om_sim], sep="")), ncols = 300, verbose=F, printstats=F)
        if (j == median_case_id) {
          temp$Rzero[om_sim] <- convertSRparms(R0=exp(ss_output$estimated_non_dev_parameters$Value[which(rownames(ss_output$estimated_non_dev_parameters)=="SR_LN(R0)")]),
                                               h=ss_output$parameters$Value[which(rownames(ss_output$parameters)=="SR_BH_steep")],
                                               phi=ss_output$sprseries$SSBzero[1]/exp(ss_output$estimated_non_dev_parameters$Value[which(rownames(ss_output$estimated_non_dev_parameters)=="SR_LN(R0)")]),
                                               sigmaR=om_input$logR_sd,
                                               mean2med=TRUE)$R0BC
        } else {
          temp$Rzero[om_sim] <- exp(ss_output$estimated_non_dev_parameters$Value[which(rownames(ss_output$estimated_non_dev_parameters)=="SR_LN(R0)")])
        }
        temp$Q[om_sim] <- log(exp(ss_output$estimated_non_dev_parameters$Value[which(rownames(ss_output$estimated_non_dev_parameters)=="LnQ_base_SURVEY1(2)")])/1000)
      }
    }
    parameter_val_list[[k]] <- temp
  }

  for(k in 1:em_num){
    parameter_re_list[[k]] <- (parameter_val_list[[k+1]]-parameter_val_list[[1]])/parameter_val_list[[1]]
  }

  for(k in 1:em_num){
    R0_RE[j,k] <- median(abs(parameter_re_list[[k]]$Rzero))
    logQ_RE[j,k] <- median(abs(parameter_re_list[[k]]$Q))
    # R0_RE[j,k] <- mean(parameter_re_list[[k]]$Rzero)
    # Q_RE[j,k] <- mean(parameter_re_list[[k]]$Q)
  }
}

#write.csv(R0_RE, file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/R0_MARE.csv", row.names = F)
write.csv(logQ_RE, file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/logQ_MARE.csv", row.names = F)

#### R0 for AMAK and ASAP original ####
R0_RE <- matrix(NA, nrow=2, ncol=2)
Q_RE <- matrix(NA, nrow=2, ncol=2)

for (j in c(mean_case_id, median_case_id)){
  maindir <- maindir_list[j]
  setwd(file.path(maindir, "output"))
  load("./keep_sim_id.RData")

  parameters <- c("Rzero")
  parameter_val_list <- list()
  temp <- as.data.frame(matrix(NA, nrow=om_sim_num, ncol=1))
  colnames(temp) <- parameters
  parameter_re_list <- list()

  for(k in 1:3){
    for (om_sim in 1:om_sim_num){
      if(k==1){
        load(file.path(maindir, "output", "OM", paste("OM", keep_sim_id[om_sim], ".RData", sep="")))
        if (j == (length(maindir_list)-1)) {
          temp$Rzero[om_sim] <- om_input$R0/1000
        } else {
          temp$Rzero[om_sim] <- om_input$R0/1000
        }
      }
      if(k==2){
        setwd(file.path(maindir,"output", "AMAK", paste("s", keep_sim_id[om_sim], sep="")))
        amak_std <- readRep("amak", suffix = ".std")
        temp$Rzero[om_sim] <- exp(amak_std$value[which(amak_std$name=="log_Rzero")])/1000
      }
      if(k==3){
        asap_output <- dget(file.path(maindir, "output", "ASAP", paste("s", keep_sim_id[om_sim], sep=""), "asap3.rdat"))
        temp$Rzero[om_sim] <- asap_output$SR.parms$SR.R0
      }
    }
    parameter_val_list[[k]] <- temp
  }

  for(k in 1:2){
    parameter_re_list[[k]] <- (parameter_val_list[[k+1]]-parameter_val_list[[1]])/parameter_val_list[[1]]
  }

  for(k in 1:2){
    R0_RE[j-9,k] <- median(abs(parameter_re_list[[k]]$Rzero))
  }
}

write.csv(R0_RE, file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/R0_original_AMAK_ASAP.csv", row.names = F)

#### R0 for AMAK exploration ####
maindir_list  <- c("C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case0",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case0_F0_0.01_1yr")
R0_RE <- matrix(NA, nrow=length(maindir_list), ncol=em_num)
Q_RE <- matrix(NA, nrow=length(maindir_list), ncol=em_num)
estimates_cases <- list()
for (j in 1:length(maindir_list)){
  maindir <- maindir_list[j]
  setwd(file.path(maindir, "output"))
  load("./keep_sim_id.RData")

  parameters <- c("Rzero", "Q")
  parameter_val_list <- list()
  temp <- as.data.frame(matrix(NA, nrow=om_sim_num, ncol=2))
  colnames(temp) <- parameters
  parameter_re_list <- list()

  for(k in 1:(em_num+1)){
    for (om_sim in 1:om_sim_num){
      if(k==1){
        load(file.path(maindir, "output", "OM", paste("OM", keep_sim_id[om_sim], ".RData", sep="")))
        temp$Rzero[om_sim] <- om_input$R0/1000
        temp$Q[om_sim] <- om_output$survey_q$survey1
      }
      if(k==2){
        setwd(file.path(maindir,"output", "AMAK", paste("s", keep_sim_id[om_sim], sep="")))
        amak_std <- readRep("amak", suffix = ".std")
        temp$Rzero[om_sim] <- exp(amak_std$value[which(amak_std$name=="log_Rzero")])/1000
        temp$Q[om_sim] <- exp(amak_std$value[which(amak_std$name=="log_q_ind[1]")])
      }
      if(k==3){
        asap_output <- dget(file.path(maindir, "output", "ASAP", paste("s", keep_sim_id[om_sim], sep=""), "asap3.rdat"))
        temp$Rzero[om_sim] <- asap_output$SR.parms$SR.R0
        temp$Q[om_sim] <- asap_output$q.indices/1000
      }
      if(k==4){
        bam_output <- dget(file.path(maindir, "output", "BAM", paste("s", keep_sim_id[om_sim], sep=""), "bam-sim.rdat"))
        temp$Rzero[om_sim] <- bam_output$parms$R0/1000
        temp$Q[om_sim] <- bam_output$parms$q.survey1
      }
      if(k==5){
        ss_output <- SS_output(dir=file.path(maindir, "output", "SS", paste("s", keep_sim_id[om_sim], sep="")), ncols = 300, verbose=F, printstats=F)
        temp$Rzero[om_sim] <- exp(ss_output$estimated_non_dev_parameters$Value[which(rownames(ss_output$estimated_non_dev_parameters)=="SR_LN(R0)")])
        temp$Q[om_sim] <- exp(ss_output$estimated_non_dev_parameters$Value[which(rownames(ss_output$estimated_non_dev_parameters)=="LnQ_base_SURVEY1(2)")])/1000
      }
    }
    parameter_val_list[[k]] <- temp
  }

  for(k in 1:em_num){
    parameter_re_list[[k]] <- (parameter_val_list[[k+1]]-parameter_val_list[[1]])/parameter_val_list[[1]]
  }

  for(k in 1:em_num){
    R0_RE[j,k] <- median(abs(parameter_re_list[[k]]$Rzero))
    Q_RE[j,k] <- median(abs(parameter_re_list[[k]]$Q))
  }

  estimates_cases[[j]] <- parameter_val_list
}

save(estimates_cases, file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/R0_AMAK_exploration.RData")
load(file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/R0_AMAK_exploration.RData")

boxplot_data <- matrix(NA, nrow=100, ncol=em_num*length(maindir_list))
sapply(1:8, function(x) {
  if (x<5) boxplot_data[,x] <<- estimates_cases[[1]][[x+1]]$Rzero
  else boxplot_data[,x] <<- estimates_cases[[2]][[x-3]]$Rzero
})

jpeg(file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/Rzero_AMAK.jpg", width=100, height=100, units="mm", res=1200)
# boxplot(boxplot_data, beside=T, axes=F, ylab="R0 (thousands of fish)", col=rep(c("orange", "green", "red", "deepskyblue"), times=2), pch=16, ylim=range(boxplot_data), cex=0.6, outline=F)
boxplot(boxplot_data,
        beside=T, axes=F,
        ylab="R0 (thousands of fish)",
        ylim=range(boxplot_data),
        col=rep("white", times=8),
        pch=16, cex=0.6, outline=F)
box()
axis(1, at=1:ncol(boxplot_data), labels = rep(c("AMAK", "ASAP", "BAM", "SS"), times=2), las=2)
axis(2, las=2)
abline(h=1000, col="gray60", lty=2)
abline(v=4.5, col="black", lty=1)
text(0.5, max(boxplot_data), "A)", cex=0.8)
text(4.8, max(boxplot_data), "B)", cex=0.8)
dev.off()



#### q for case ####

## q2
maindir_list  <- c("C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case9")
R0_RE <- matrix(NA, nrow=length(maindir_list), ncol=em_num)
Q_RE <- matrix(NA, nrow=length(maindir_list), ncol=em_num)
estimates_cases <- list()
for (j in 1:length(maindir_list)){
  maindir <- maindir_list[j]
  setwd(file.path(maindir, "output"))
  load("./keep_sim_id.RData")

  parameters <- c("Rzero", "Q")
  parameter_val_list <- list()
  temp <- as.data.frame(matrix(NA, nrow=om_sim_num, ncol=2))
  colnames(temp) <- parameters
  parameter_re_list <- list()

  for(k in 1:(em_num+1)){
    for (om_sim in 1:om_sim_num){
      if(k==1){
        load(file.path(maindir, "output", "OM", paste("OM", keep_sim_id[om_sim], ".RData", sep="")))
        temp$Rzero[om_sim] <- om_input$R0/1000
        temp$Q[om_sim] <- om_output$survey_q$survey2
      }
      if(k==2){
        setwd(file.path(maindir,"output", "AMAK", paste("s", keep_sim_id[om_sim], sep="")))
        amak_std <- readRep("amak", suffix = ".std")
        temp$Rzero[om_sim] <- exp(amak_std$value[which(amak_std$name=="log_Rzero")])/1000
        temp$Q[om_sim] <- exp(amak_std$value[which(amak_std$name=="log_q_ind[2]")])
      }
      if(k==3){
        asap_output <- dget(file.path(maindir, "output", "ASAP", paste("s", keep_sim_id[om_sim], sep=""), "asap3.rdat"))
        temp$Rzero[om_sim] <- asap_output$SR.parms$SR.R0
        temp$Q[om_sim] <- asap_output$q.indices[2]/1000
      }
      if(k==4){
        bam_output <- dget(file.path(maindir, "output", "BAM", paste("s", keep_sim_id[om_sim], sep=""), "bam-sim.rdat"))
        temp$Rzero[om_sim] <- bam_output$parms$R0/1000
        temp$Q[om_sim] <- bam_output$parms$q.survey2
      }
      if(k==5){
        ss_output <- SS_output(dir=file.path(maindir, "output", "SS", paste("s", keep_sim_id[om_sim], sep="")), ncols = 300, verbose=F, printstats=F)
        temp$Rzero[om_sim] <- exp(ss_output$estimated_non_dev_parameters$Value[which(rownames(ss_output$estimated_non_dev_parameters)=="SR_LN(R0)")])
        temp$Q[om_sim] <- exp(ss_output$estimated_non_dev_parameters$Value[which(rownames(ss_output$estimated_non_dev_parameters)=="LnQ_base_SURVEY2(3)")])/1000
      }
    }
    parameter_val_list[[k]] <- temp
  }

  for(k in 1:em_num){
    parameter_re_list[[k]] <- (parameter_val_list[[k+1]]-parameter_val_list[[1]])/parameter_val_list[[1]]
  }

  for(k in 1:em_num){
    R0_RE[j,k] <- median(abs(parameter_re_list[[k]]$Rzero))
    Q_RE[j,k] <- median(abs(parameter_re_list[[k]]$Q))
  }

  estimates_cases[[j]] <- parameter_val_list
}


write.csv(Q_RE, file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/Q2_MARE.csv", row.names = F)

##q1
maindir_list  <- c("C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case9")
R0_RE <- matrix(NA, nrow=length(maindir_list), ncol=em_num)
Q_RE <- matrix(NA, nrow=length(maindir_list), ncol=em_num)
estimates_cases <- list()
for (j in 1:length(maindir_list)){
  maindir <- maindir_list[j]
  setwd(file.path(maindir, "output"))
  load("./keep_sim_id.RData")

  parameters <- c("Rzero", "Q")
  parameter_val_list <- list()
  temp <- as.data.frame(matrix(NA, nrow=om_sim_num, ncol=2))
  colnames(temp) <- parameters
  parameter_re_list <- list()

  for(k in 1:(em_num+1)){
    for (om_sim in 1:om_sim_num){
      if(k==1){
        load(file.path(maindir, "output", "OM", paste("OM", keep_sim_id[om_sim], ".RData", sep="")))
        temp$Rzero[om_sim] <- om_input$R0/1000
        temp$Q[om_sim] <- om_output$survey_q$survey1
      }
      if(k==2){
        setwd(file.path(maindir,"output", "AMAK", paste("s", keep_sim_id[om_sim], sep="")))
        amak_std <- readRep("amak", suffix = ".std")
        temp$Rzero[om_sim] <- exp(amak_std$value[which(amak_std$name=="log_Rzero")])/1000
        temp$Q[om_sim] <- exp(amak_std$value[which(amak_std$name=="log_q_ind[1]")])
      }
      if(k==3){
        asap_output <- dget(file.path(maindir, "output", "ASAP", paste("s", keep_sim_id[om_sim], sep=""), "asap3.rdat"))
        temp$Rzero[om_sim] <- asap_output$SR.parms$SR.R0
        temp$Q[om_sim] <- asap_output$q.indices[1]/1000
      }
      if(k==4){
        bam_output <- dget(file.path(maindir, "output", "BAM", paste("s", keep_sim_id[om_sim], sep=""), "bam-sim.rdat"))
        temp$Rzero[om_sim] <- bam_output$parms$R0/1000
        temp$Q[om_sim] <- bam_output$parms$q.survey1
      }
      if(k==5){
        ss_output <- SS_output(dir=file.path(maindir, "output", "SS", paste("s", keep_sim_id[om_sim], sep="")), ncols = 300, verbose=F, printstats=F)
        temp$Rzero[om_sim] <- exp(ss_output$estimated_non_dev_parameters$Value[which(rownames(ss_output$estimated_non_dev_parameters)=="SR_LN(R0)")])
        temp$Q[om_sim] <- exp(ss_output$estimated_non_dev_parameters$Value[which(rownames(ss_output$estimated_non_dev_parameters)=="LnQ_base_SURVEY1(2)")])/1000
      }
    }
    parameter_val_list[[k]] <- temp
  }

  for(k in 1:em_num){
    parameter_re_list[[k]] <- (parameter_val_list[[k+1]]-parameter_val_list[[1]])/parameter_val_list[[1]]
  }

  for(k in 1:em_num){
    R0_RE[j,k] <- median(abs(parameter_re_list[[k]]$Rzero))
    Q_RE[j,k] <- median(abs(parameter_re_list[[k]]$Q))
  }

  estimates_cases[[j]] <- parameter_val_list
}

write.csv(Q_RE, file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/Q1_MARE.csv", row.names = F)

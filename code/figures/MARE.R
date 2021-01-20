
# Median Absolute Relative Error (MARE) --------------------------------------------------------------
## Create a table of MARE in MSY, FMSY, SSBMSY, SSB, R, F, SSBratio, and Fratio


# Cases ----------------------------------------------------------------------------
dir <- "C:/Users/bai.li/Desktop/mcp_results_r/final_cases"
cases <- paste0("case", c(0:11, "0_F0_0.01_1yr"), sep="")
maindir_list <- file.path(dir, cases)

em_num <- 4

# Survey indices MARE and MRE -----------------------------------------------------------------
survey_MARE <- survey_MRE <- matrix(NA, nrow=length(maindir_list)*30, ncol=em_num)
for(j in 1:length(maindir_list)){
  load(file.path(maindir_list[j], "output", "performance_measures_output.RData"))
  if (j==1){
    survey_are <- survey_re <- list()
  }
  for (i in 1:length(pm_list)){
    survey_are[[i]] <- as.matrix(pm_list[[i]]$survey$are)
    survey_re[[i]] <- as.matrix(pm_list[[i]]$survey$re)
  }

  for(k in 1:em_num){
    survey_MARE[((j-1)*30+1):(j*30),k] <- apply(sapply(1:length(pm_list), function(x) survey_are[[x]][,k]), 1, median)
    survey_MRE[((j-1)*30+1):(j*30),k] <- apply(sapply(1:length(pm_list), function(x) survey_re[[x]][,k]), 1, median)
  }
}

survey_MARE=round(survey_MARE*100, digits = 2)
survey_MARE=cbind(survey_MARE, rep(1:30, times=length(maindir_list)), rep(0:12, each=30))

survey_MRE=round(survey_MRE*100, digits = 2)
survey_MRE=cbind(survey_MRE, rep(1:30, times=length(maindir_list)), rep(0:12, each=30))

write.csv(survey_MARE, file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/survey_mare_percentage.csv", row.names = F)
write.csv(survey_MRE, file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/survey_mre_percentage.csv", row.names = F)


# MSY related MARE ----------------------------------------------------------------------------

msy_mare <- matrix(NA, nrow=length(maindir_list), ncol=em_num*3)

for(j in 1:length(maindir_list)){
  load(file.path(maindir_list[j], "output", "performance_measures_output.RData"))
  if (j==1){
    msy_are <- matrix(NA, nrow=length(pm_list), ncol=em_num*length(maindir_list))
    fmsy_are <- matrix(NA, nrow=length(pm_list), ncol=em_num*length(maindir_list))
    ssbmsy_are <- matrix(NA, nrow=length(pm_list), ncol=em_num*length(maindir_list))
  }
  for (i in 1:length(pm_list)){
    msy_are[i,((j-1)*em_num+1):((j-1)*em_num+em_num)] <- as.matrix(pm_list[[i]]$msy$are)
    fmsy_are[i,((j-1)*em_num+1):((j-1)*em_num+em_num)] <- as.matrix(pm_list[[i]]$fmsy$are)
    ssbmsy_are[i,((j-1)*em_num+1):((j-1)*em_num+em_num)] <- as.matrix(pm_list[[i]]$ssbmsy$are)
  }
}

msy_stat <- matrix(NA, nrow=3, ncol=em_num*length(maindir_list))
fmsy_stat <- matrix(NA, nrow=3, ncol=em_num*length(maindir_list))
ssbmsy_stat <- matrix(NA, nrow=3, ncol=em_num*length(maindir_list))

sapply(1:ncol(msy_stat), function(x) {
  msy_stat[,x] <<- boxplot.stats(msy_are[,x])$`stats`[c(1,3,5)]
  fmsy_stat[,x] <<- boxplot.stats(fmsy_are[,x])$`stats`[c(1,3,5)]
  ssbmsy_stat[,x] <<- boxplot.stats(ssbmsy_are[,x])$`stats`[c(1,3,5)]
})

msy_mare[,1:em_num] <- matrix(msy_stat[2,], ncol=em_num, byrow=T)
msy_mare[,(em_num+1):(em_num*2)] <- matrix(fmsy_stat[2,], ncol=em_num, byrow=T)
msy_mare[,(em_num*2+1):(em_num*3)] <- matrix(ssbmsy_stat[2,], ncol=em_num, byrow=T)

write.csv(round(msy_mare*100, digits = 2), file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/msy_mare_percentage.csv", row.names = paste("C", 0:12, sep=""))

# SSB, R, F -----------------------------------------------------------------------------------

ssb_stat <- r_stat <- f_stat <- matrix(NA, nrow=length(maindir_list), ncol=em_num)
for(j in 1:length(maindir_list)){
  load(file.path(maindir_list[j], "output", "performance_measures_output.RData"))
  if (j==1){
    ssb_are <- r_are <- f_are <- list()
  }
  for (i in 1:length(pm_list)){
    ssb_are[[i]] <- as.matrix(pm_list[[i]]$ssb$are)
    r_are[[i]] <- as.matrix(pm_list[[i]]$recruit$are)
    f_are[[i]] <- as.matrix(pm_list[[i]]$Ftot$are)
  }

  for(k in 1:em_num){
    ssb_stat[j,k] <- median(sapply(1:length(pm_list), function(x) ssb_are[[x]][,k]))
    r_stat[j,k] <- median(sapply(1:length(pm_list), function(x) r_are[[x]][,k]))
    f_stat[j,k] <- median(sapply(1:length(pm_list), function(x) f_are[[x]][,k]))
  }
}

ssbrf_mare <- cbind(ssb_stat, r_stat, f_stat)
write.csv(round(ssbrf_mare*100, digits = 2), file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/ssbrf_mare_percentage.csv", row.names = paste("C", 0:12, sep=""))

# SSBratio and Fratio -------------------------------------------------------------------------

ssbratio_stat <- fratio_stat <- matrix(NA, nrow=length(maindir_list), ncol=em_num)
for(j in 1:length(maindir_list)){
  load(file.path(maindir_list[j], "output", "performance_measures_output.RData"))
  if (j==1){
    ssbratio_are <- fratio_are <- list()
  }
  for (i in 1:length(pm_list)){
    ssbratio_are[[i]] <- as.matrix(pm_list[[i]]$ssbratio$are)
    fratio_are[[i]] <- as.matrix(pm_list[[i]]$fratio$are)
  }

  for(k in 1:em_num){
    ssbratio_stat[j,k] <- median(sapply(1:length(pm_list), function(x) ssbratio_are[[x]][,k]))
    fratio_stat[j,k] <- median(sapply(1:length(pm_list), function(x) fratio_are[[x]][,k]))
  }
}

ssbfratio_mare <- cbind(ssbratio_stat, fratio_stat)
write.csv(round(ssbfratio_mare*100, digits = 2), file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/ssbfratio_mare_percentage.csv", row.names = paste("C", 0:12, sep=""))

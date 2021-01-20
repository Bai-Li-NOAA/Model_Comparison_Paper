#### Boxplot with whiskers (C0 and C12) ####
maindir_list  <- c("C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case0",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case12")

em_num <- 4
legend_names <- c("C0", "C12")

load(file.path(maindir_list[1], "output", "performance_measures_output.RData"))
ssb_median <- ssb_low <- ssb_high <-
  r_median <- r_low <- r_high <-
  f_median <- f_low <- f_high <-
  ssbratio_median <- ssbratio_low <- ssbratio_high <-
  fratio_median <- fratio_low <- fratio_high <-
  matrix(NA, nrow=nrow(pm_list[[1]]$ssb$re), ncol=em_num*length(maindir_list))

for(j in 1:length(maindir_list)){
  load(file.path(maindir_list[j], "output", "performance_measures_output.RData"))
  if (j==1){
    ssb_re <- r_re <- f_re <-
    ssbratio_re <- fratio_re <-
    matrix(NA, nrow=length(pm_list), ncol=nrow(pm_list[[1]]$ssb$re))
  }
  for (k in 1:em_num){
    for (i in 1:length(pm_list)){
      ssb_re[i,] <- as.matrix(pm_list[[i]]$ssb$re[,k])
      r_re[i,] <- as.matrix(pm_list[[i]]$recruit$re[,k])
      f_re[i,] <- as.matrix(pm_list[[i]]$Ftot$re[,k])
      ssbratio_re[i,] <- as.matrix(pm_list[[i]]$ssbratio$re[,k])
      fratio_re[i,] <- as.matrix(pm_list[[i]]$fratio$re[,k])
    }
    ssb_median[,((j-1)*4+k)] <- sapply(1:ncol(ssb_re), function(x) boxplot.stats(ssb_re[,x])$`stats`[3])
    ssb_low[,((j-1)*4+k)] <- sapply(1:ncol(ssb_re), function(x) boxplot.stats(ssb_re[,x])$`stats`[1])
    ssb_high[,((j-1)*4+k)] <- sapply(1:ncol(ssb_re), function(x) boxplot.stats(ssb_re[,x])$`stats`[5])

    r_median[,((j-1)*4+k)] <- sapply(1:ncol(r_re), function(x) boxplot.stats(r_re[,x])$`stats`[3])
    r_low[,((j-1)*4+k)] <- sapply(1:ncol(r_re), function(x) boxplot.stats(r_re[,x])$`stats`[1])
    r_high[,((j-1)*4+k)] <- sapply(1:ncol(r_re), function(x) boxplot.stats(r_re[,x])$`stats`[5])

    f_median[,((j-1)*4+k)] <- sapply(1:ncol(f_re), function(x) boxplot.stats(f_re[,x])$`stats`[3])
    f_low[,((j-1)*4+k)] <- sapply(1:ncol(f_re), function(x) boxplot.stats(f_re[,x])$`stats`[1])
    f_high[,((j-1)*4+k)] <- sapply(1:ncol(f_re), function(x) boxplot.stats(f_re[,x])$`stats`[5])

    ssbratio_median[,((j-1)*4+k)] <- sapply(1:ncol(ssbratio_re), function(x) boxplot.stats(ssbratio_re[,x])$`stats`[3])
    ssbratio_low[,((j-1)*4+k)] <- sapply(1:ncol(ssbratio_re), function(x) boxplot.stats(ssbratio_re[,x])$`stats`[1])
    ssbratio_high[,((j-1)*4+k)] <- sapply(1:ncol(ssbratio_re), function(x) boxplot.stats(ssbratio_re[,x])$`stats`[5])

    fratio_median[,((j-1)*4+k)] <- sapply(1:ncol(fratio_re), function(x) boxplot.stats(fratio_re[,x])$`stats`[3])
    fratio_low[,((j-1)*4+k)] <- sapply(1:ncol(fratio_re), function(x) boxplot.stats(fratio_re[,x])$`stats`[1])
    fratio_high[,((j-1)*4+k)] <- sapply(1:ncol(fratio_re), function(x) boxplot.stats(fratio_re[,x])$`stats`[5])
  }
}

jpeg(file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/SSB_R_F_RE_two_surveys_boxplot_whisker.jpg", width=130, height=100, units="mm", res=600)

op<-par(no.readonly=TRUE)
par(op)
par(oma=c(5,5,0,5),mar=c(0,0,0,0.2),mfrow=c(2,5), pch=16)

x_val <- 1:length(ssb_median[,1])
ylim=c(-0.6, 0.6)

for (j in 1:length(maindir_list)){
  plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  axis(2, at=c(-0.4, 0, 0.4), labels=c(-0.4, 0, 0.4), las=2)
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*4+1], rev(ssb_high[,(j-1)*4+1])), border=NA, col=alpha("orange", 0.2))
  #lines(x_val, ssb_median[,1], type="o", lty=1, col="orange", pch=2)

  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*4+2], rev(ssb_high[,(j-1)*4+2])), border=NA, col=alpha("green", 0.2))
  #lines(x_val, ssb_median[,2], type="o", lty=1, col="green", pch=3)

  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*4+3], rev(ssb_high[,(j-1)*4+3])), border=NA, col=alpha("red", 0.2))
  #lines(x_val, ssb_median[,3], type="o", lty=1, col="red", pch=4)

  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*4+4], rev(ssb_high[,(j-1)*4+4])), border=NA, col=alpha("deepskyblue", 0.2))
  #lines(x_val, ssb_median[,4], type="o", lty=1, col="deepskyblue", pch=5)

  lines(x_val, ssb_median[,(j-1)*4+1], type="l", lty=2, col="orange")
  lines(x_val, ssb_median[,(j-1)*4+2], type="l", lty=3, col="green")
  lines(x_val, ssb_median[,(j-1)*4+3], type="l", lty=4, col="red")
  lines(x_val, ssb_median[,(j-1)*4+4], type="l", lty=5, col="deepskyblue")

  legend("topleft", paste(legend_names[j], "_SSB", sep=""), bty="n", cex=0.8)

  plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*4+1], rev(r_high[,(j-1)*4+1])), border=NA, col=alpha("orange", 0.2))
  #lines(x_val, r_median[,1], type="o", lty=1, col="orange", pch=2)

  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*4+2], rev(r_high[,(j-1)*4+2])), border=NA, col=alpha("green", 0.2))
  #lines(x_val, r_median[,2], type="o", lty=1, col="green", pch=3)

  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*4+3], rev(r_high[,(j-1)*4+3])), border=NA, col=alpha("red", 0.2))
  #lines(x_val, r_median[,3], type="o", lty=1, col="red", pch=4)

  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*4+4], rev(r_high[,(j-1)*4+4])), border=NA, col=alpha("deepskyblue", 0.2))
  #lines(x_val, r_median[,4], type="o", lty=1, col="deepskyblue", pch=5)

  lines(x_val, r_median[,(j-1)*4+1], type="l", lty=2, col="orange")
  lines(x_val, r_median[,(j-1)*4+2], type="l", lty=3, col="green")
  lines(x_val, r_median[,(j-1)*4+3], type="l", lty=4, col="red")
  lines(x_val, r_median[,(j-1)*4+4], type="l", lty=5, col="deepskyblue")

  legend("topleft", paste(legend_names[j], "_R", sep=""), bty="n", cex=0.8)

  plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*4+1], rev(f_high[,(j-1)*4+1])), border=NA, col=alpha("orange", 0.2))
  #lines(x_val, f_median[,1], type="o", lty=1, col="orange", pch=2)

  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*4+2], rev(f_high[,(j-1)*4+2])), border=NA, col=alpha("green", 0.2))
  #lines(x_val, f_median[,2], type="o", lty=1, col="green", pch=3)

  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*4+3], rev(f_high[,(j-1)*4+3])), border=NA, col=alpha("red", 0.2))
  #lines(x_val, f_median[,3], type="o", lty=1, col="red", pch=4)

  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*4+4], rev(f_high[,(j-1)*4+4])), border=NA, col=alpha("deepskyblue", 0.2))
  #lines(x_val, f_median[,4], type="o", lty=1, col="deepskyblue", pch=5)

  lines(x_val, f_median[,(j-1)*4+1], type="l", lty=2, col="orange")
  lines(x_val, f_median[,(j-1)*4+2], type="l", lty=3, col="green")
  lines(x_val, f_median[,(j-1)*4+3], type="l", lty=4, col="red")
  lines(x_val, f_median[,(j-1)*4+4], type="l", lty=5, col="deepskyblue")

  legend("topleft", paste(legend_names[j], "_F", sep=""), bty="n", cex=0.8)

  plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*4+1], rev(ssbratio_high[,(j-1)*4+1])), border=NA, col=alpha("orange", 0.2))
  #lines(x_val, ssbratio_median[,1], type="o", lty=1, col="orange", pch=2)

  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*4+2], rev(ssbratio_high[,(j-1)*4+2])), border=NA, col=alpha("green", 0.2))
  #lines(x_val, ssbratio_median[,2], type="o", lty=1, col="green", pch=3)

  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*4+3], rev(ssbratio_high[,(j-1)*4+3])), border=NA, col=alpha("red", 0.2))
  #lines(x_val, ssbratio_median[,3], type="o", lty=1, col="red", pch=4)

  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*4+4], rev(ssbratio_high[,(j-1)*4+4])), border=NA, col=alpha("deepskyblue", 0.2))
  #lines(x_val, ssbratio_median[,4], type="o", lty=1, col="deepskyblue", pch=5)

  lines(x_val, ssbratio_median[,(j-1)*4+1], type="l", lty=2, col="orange")
  lines(x_val, ssbratio_median[,(j-1)*4+2], type="l", lty=3, col="green")
  lines(x_val, ssbratio_median[,(j-1)*4+3], type="l", lty=4, col="red")
  lines(x_val, ssbratio_median[,(j-1)*4+4], type="l", lty=5, col="deepskyblue")

  legend("topleft", paste(legend_names[j], "_SSBratio", sep=""), bty="n", cex=0.8)

  plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*4+1], rev(fratio_high[,(j-1)*4+1])), border=NA, col=alpha("orange", 0.2))
  #lines(x_val, fratio_median[,1], type="o", lty=1, col="orange", pch=2)

  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*4+2], rev(fratio_high[,(j-1)*4+2])), border=NA, col=alpha("green", 0.2))
  #lines(x_val, fratio_median[,2], type="o", lty=1, col="green", pch=3)

  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*4+3], rev(fratio_high[,(j-1)*4+3])), border=NA, col=alpha("red", 0.2))
  #lines(x_val, fratio_median[,3], type="o", lty=1, col="red", pch=4)

  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*4+4], rev(fratio_high[,(j-1)*4+4])), border=NA, col=alpha("deepskyblue", 0.2))
  #lines(x_val, fratio_median[,4], type="o", lty=1, col="deepskyblue", pch=5)

  lines(x_val, fratio_median[,(j-1)*4+1], type="l", lty=2, col="orange")
  lines(x_val, fratio_median[,(j-1)*4+2], type="l", lty=3, col="green")
  lines(x_val, fratio_median[,(j-1)*4+3], type="l", lty=4, col="red")
  lines(x_val, fratio_median[,(j-1)*4+4], type="l", lty=5, col="deepskyblue")

  legend("topleft", paste(legend_names[j], "_Fratio", sep=""), bty="n", cex=0.8)
}

mtext(text="Year",side=1,line=2.5,outer=TRUE)
mtext(text="RE",side=2,line=2.5,outer=TRUE)

legend(x=30.5,y=1,
       legend=c("AMAK", "ASAP", "BAM", "SS"),
       lty=c(2, 3, 4, 5),
       col=c("orange", "green", "red", "deepskyblue3"),
       cex=1,
       bty="n",
       xpd=NA)
dev.off()

#### Boxplot with confidence interval (C0 and C12) ####
maindir_list  <- c("C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case0",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case12")

em_num <- 4
legend_names <- c("C0", "C12")

load(file.path(maindir_list[1], "output", "performance_measures_output.RData"))
ssb_median <- ssb_low <- ssb_high <-
  r_median <- r_low <- r_high <-
  f_median <- f_low <- f_high <-
  ssbratio_median <- ssbratio_low <- ssbratio_high <-
  fratio_median <- fratio_low <- fratio_high <-
  matrix(NA, nrow=nrow(pm_list[[1]]$ssb$re), ncol=em_num*length(maindir_list))

for(j in 1:length(maindir_list)){
  load(file.path(maindir_list[j], "output", "performance_measures_output.RData"))
  if (j==1){
    ssb_re <- r_re <- f_re <-
      ssbratio_re <- fratio_re <-
      matrix(NA, nrow=length(pm_list), ncol=nrow(pm_list[[1]]$ssb$re))
  }
  for (k in 1:em_num){
    for (i in 1:length(pm_list)){
      ssb_re[i,] <- as.matrix(pm_list[[i]]$ssb$re[,k])
      r_re[i,] <- as.matrix(pm_list[[i]]$recruit$re[,k])
      f_re[i,] <- as.matrix(pm_list[[i]]$Ftot$re[,k])
      ssbratio_re[i,] <- as.matrix(pm_list[[i]]$ssbratio$re[,k])
      fratio_re[i,] <- as.matrix(pm_list[[i]]$fratio$re[,k])
    }
    ssb_median[,((j-1)*4+k)] <- sapply(1:ncol(ssb_re), function(x) boxplot.stats(ssb_re[,x])$`stats`[3])
    ssb_low[,((j-1)*4+k)] <- sapply(1:ncol(ssb_re), function(x) boxplot.stats(ssb_re[,x])$`conf`[1])
    ssb_high[,((j-1)*4+k)] <- sapply(1:ncol(ssb_re), function(x) boxplot.stats(ssb_re[,x])$`conf`[2])

    r_median[,((j-1)*4+k)] <- sapply(1:ncol(r_re), function(x) boxplot.stats(r_re[,x])$`stats`[3])
    r_low[,((j-1)*4+k)] <- sapply(1:ncol(r_re), function(x) boxplot.stats(r_re[,x])$`conf`[1])
    r_high[,((j-1)*4+k)] <- sapply(1:ncol(r_re), function(x) boxplot.stats(r_re[,x])$`conf`[2])

    f_median[,((j-1)*4+k)] <- sapply(1:ncol(f_re), function(x) boxplot.stats(f_re[,x])$`stats`[3])
    f_low[,((j-1)*4+k)] <- sapply(1:ncol(f_re), function(x) boxplot.stats(f_re[,x])$`conf`[1])
    f_high[,((j-1)*4+k)] <- sapply(1:ncol(f_re), function(x) boxplot.stats(f_re[,x])$`conf`[2])

    ssbratio_median[,((j-1)*4+k)] <- sapply(1:ncol(ssbratio_re), function(x) boxplot.stats(ssbratio_re[,x])$`stats`[3])
    ssbratio_low[,((j-1)*4+k)] <- sapply(1:ncol(ssbratio_re), function(x) boxplot.stats(ssbratio_re[,x])$`conf`[1])
    ssbratio_high[,((j-1)*4+k)] <- sapply(1:ncol(ssbratio_re), function(x) boxplot.stats(ssbratio_re[,x])$`conf`[2])

    fratio_median[,((j-1)*4+k)] <- sapply(1:ncol(fratio_re), function(x) boxplot.stats(fratio_re[,x])$`stats`[3])
    fratio_low[,((j-1)*4+k)] <- sapply(1:ncol(fratio_re), function(x) boxplot.stats(fratio_re[,x])$`conf`[1])
    fratio_high[,((j-1)*4+k)] <- sapply(1:ncol(fratio_re), function(x) boxplot.stats(fratio_re[,x])$`conf`[2])
  }
}

jpeg(file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/SSB_R_F_RE_two_surveys_boxplot_ci.jpg", width=130, height=100, units="mm", res=600)

op<-par(no.readonly=TRUE)
par(op)
par(oma=c(5,5,0,5),mar=c(0,0,0,0.2),mfrow=c(2,5), pch=16)

x_val <- 1:length(ssb_median[,1])
ylim=c(-0.2, 0.2)

for (j in 1:length(maindir_list)){
  plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  axis(2, at=c(-0.1, 0, 0.1), labels=c(-0.1, 0, 0.1), las=2)
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*4+1], rev(ssb_high[,(j-1)*4+1])), border=NA, col=alpha("orange", 0.2))
  #lines(x_val, ssb_median[,1], type="o", lty=1, col="orange", pch=2)

  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*4+2], rev(ssb_high[,(j-1)*4+2])), border=NA, col=alpha("green", 0.2))
  #lines(x_val, ssb_median[,2], type="o", lty=1, col="green", pch=3)

  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*4+3], rev(ssb_high[,(j-1)*4+3])), border=NA, col=alpha("red", 0.2))
  #lines(x_val, ssb_median[,3], type="o", lty=1, col="red", pch=4)

  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*4+4], rev(ssb_high[,(j-1)*4+4])), border=NA, col=alpha("deepskyblue", 0.2))
  #lines(x_val, ssb_median[,4], type="o", lty=1, col="deepskyblue", pch=5)

  lines(x_val, ssb_median[,(j-1)*4+1], type="l", lty=2, col="orange")
  lines(x_val, ssb_median[,(j-1)*4+2], type="l", lty=3, col="green")
  lines(x_val, ssb_median[,(j-1)*4+3], type="l", lty=4, col="red")
  lines(x_val, ssb_median[,(j-1)*4+4], type="l", lty=5, col="deepskyblue")

  legend("topleft", paste(legend_names[j], "_SSB", sep=""), bty="n", cex=0.8)

  plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*4+1], rev(r_high[,(j-1)*4+1])), border=NA, col=alpha("orange", 0.2))
  #lines(x_val, r_median[,1], type="o", lty=1, col="orange", pch=2)

  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*4+2], rev(r_high[,(j-1)*4+2])), border=NA, col=alpha("green", 0.2))
  #lines(x_val, r_median[,2], type="o", lty=1, col="green", pch=3)

  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*4+3], rev(r_high[,(j-1)*4+3])), border=NA, col=alpha("red", 0.2))
  #lines(x_val, r_median[,3], type="o", lty=1, col="red", pch=4)

  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*4+4], rev(r_high[,(j-1)*4+4])), border=NA, col=alpha("deepskyblue", 0.2))
  #lines(x_val, r_median[,4], type="o", lty=1, col="deepskyblue", pch=5)

  lines(x_val, r_median[,(j-1)*4+1], type="l", lty=2, col="orange")
  lines(x_val, r_median[,(j-1)*4+2], type="l", lty=3, col="green")
  lines(x_val, r_median[,(j-1)*4+3], type="l", lty=4, col="red")
  lines(x_val, r_median[,(j-1)*4+4], type="l", lty=5, col="deepskyblue")

  legend("topleft", paste(legend_names[j], "_R", sep=""), bty="n", cex=0.8)

  plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*4+1], rev(f_high[,(j-1)*4+1])), border=NA, col=alpha("orange", 0.2))
  #lines(x_val, f_median[,1], type="o", lty=1, col="orange", pch=2)

  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*4+2], rev(f_high[,(j-1)*4+2])), border=NA, col=alpha("green", 0.2))
  #lines(x_val, f_median[,2], type="o", lty=1, col="green", pch=3)

  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*4+3], rev(f_high[,(j-1)*4+3])), border=NA, col=alpha("red", 0.2))
  #lines(x_val, f_median[,3], type="o", lty=1, col="red", pch=4)

  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*4+4], rev(f_high[,(j-1)*4+4])), border=NA, col=alpha("deepskyblue", 0.2))
  #lines(x_val, f_median[,4], type="o", lty=1, col="deepskyblue", pch=5)

  lines(x_val, f_median[,(j-1)*4+1], type="l", lty=2, col="orange")
  lines(x_val, f_median[,(j-1)*4+2], type="l", lty=3, col="green")
  lines(x_val, f_median[,(j-1)*4+3], type="l", lty=4, col="red")
  lines(x_val, f_median[,(j-1)*4+4], type="l", lty=5, col="deepskyblue")

  legend("topleft", paste(legend_names[j], "_F", sep=""), bty="n", cex=0.8)

  plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*4+1], rev(ssbratio_high[,(j-1)*4+1])), border=NA, col=alpha("orange", 0.2))
  #lines(x_val, ssbratio_median[,1], type="o", lty=1, col="orange", pch=2)

  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*4+2], rev(ssbratio_high[,(j-1)*4+2])), border=NA, col=alpha("green", 0.2))
  #lines(x_val, ssbratio_median[,2], type="o", lty=1, col="green", pch=3)

  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*4+3], rev(ssbratio_high[,(j-1)*4+3])), border=NA, col=alpha("red", 0.2))
  #lines(x_val, ssbratio_median[,3], type="o", lty=1, col="red", pch=4)

  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*4+4], rev(ssbratio_high[,(j-1)*4+4])), border=NA, col=alpha("deepskyblue", 0.2))
  #lines(x_val, ssbratio_median[,4], type="o", lty=1, col="deepskyblue", pch=5)

  lines(x_val, ssbratio_median[,(j-1)*4+1], type="l", lty=2, col="orange")
  lines(x_val, ssbratio_median[,(j-1)*4+2], type="l", lty=3, col="green")
  lines(x_val, ssbratio_median[,(j-1)*4+3], type="l", lty=4, col="red")
  lines(x_val, ssbratio_median[,(j-1)*4+4], type="l", lty=5, col="deepskyblue")

  legend("topleft", paste(legend_names[j], "_SSBratio", sep=""), bty="n", cex=0.8)

  plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*4+1], rev(fratio_high[,(j-1)*4+1])), border=NA, col=alpha("orange", 0.2))
  #lines(x_val, fratio_median[,1], type="o", lty=1, col="orange", pch=2)

  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*4+2], rev(fratio_high[,(j-1)*4+2])), border=NA, col=alpha("green", 0.2))
  #lines(x_val, fratio_median[,2], type="o", lty=1, col="green", pch=3)

  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*4+3], rev(fratio_high[,(j-1)*4+3])), border=NA, col=alpha("red", 0.2))
  #lines(x_val, fratio_median[,3], type="o", lty=1, col="red", pch=4)

  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*4+4], rev(fratio_high[,(j-1)*4+4])), border=NA, col=alpha("deepskyblue", 0.2))
  #lines(x_val, fratio_median[,4], type="o", lty=1, col="deepskyblue", pch=5)

  lines(x_val, fratio_median[,(j-1)*4+1], type="l", lty=2, col="orange")
  lines(x_val, fratio_median[,(j-1)*4+2], type="l", lty=3, col="green")
  lines(x_val, fratio_median[,(j-1)*4+3], type="l", lty=4, col="red")
  lines(x_val, fratio_median[,(j-1)*4+4], type="l", lty=5, col="deepskyblue")

  legend("topleft", paste(legend_names[j], "_Fratio", sep=""), bty="n", cex=0.8)
}

mtext(text="Year",side=1,line=2.5,outer=TRUE)
mtext(text="RE",side=2,line=2.5,outer=TRUE)

legend(x=30.5,y=0.3,
       legend=c("AMAK", "ASAP", "BAM", "SS"),
       lty=c(2, 3, 4, 5),
       col=c("orange", "green", "red", "deepskyblue3"),
       cex=1,
       bty="n",
       xpd=NA)
dev.off()
#### Quantile (C0 and C12) ####
maindir_list  <- c("C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case0",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case12")

em_num <- 4
legend_names <- c("C0", "C12")

load(file.path(maindir_list[1], "output", "performance_measures_output.RData"))
ssb_median <- ssb_low <- ssb_high <-
  r_median <- r_low <- r_high <-
  f_median <- f_low <- f_high <-
  ssbratio_median <- ssbratio_low <- ssbratio_high <-
  fratio_median <- fratio_low <- fratio_high <-
  matrix(NA, nrow=nrow(pm_list[[1]]$ssb$re), ncol=em_num*length(maindir_list))

for(j in 1:length(maindir_list)){
  load(file.path(maindir_list[j], "output", "performance_measures_output.RData"))
  if (j==1){
    ssb_re <- r_re <- f_re <-
      ssbratio_re <- fratio_re <-
      matrix(NA, nrow=length(pm_list), ncol=nrow(pm_list[[1]]$ssb$re))
  }
  for (k in 1:em_num){
    for (i in 1:length(pm_list)){
      ssb_re[i,] <- as.matrix(pm_list[[i]]$ssb$re[,k])
      r_re[i,] <- as.matrix(pm_list[[i]]$recruit$re[,k])
      f_re[i,] <- as.matrix(pm_list[[i]]$Ftot$re[,k])
      ssbratio_re[i,] <- as.matrix(pm_list[[i]]$ssbratio$re[,k])
      fratio_re[i,] <- as.matrix(pm_list[[i]]$fratio$re[,k])
    }
    ssb_median[,((j-1)*4+k)] <- sapply(1:ncol(ssb_re), function(x) quantile(ssb_re[,x], c(0.5)))
    ssb_low[,((j-1)*4+k)] <- sapply(1:ncol(ssb_re), function(x) quantile(ssb_re[,x], c(0.1)))
    ssb_high[,((j-1)*4+k)] <- sapply(1:ncol(ssb_re), function(x) quantile(ssb_re[,x], c(0.9)))

    r_median[,((j-1)*4+k)] <- sapply(1:ncol(r_re), function(x) quantile(r_re[,x], c(0.5)))
    r_low[,((j-1)*4+k)] <- sapply(1:ncol(r_re), function(x) quantile(r_re[,x], c(0.1)))
    r_high[,((j-1)*4+k)] <- sapply(1:ncol(r_re), function(x) quantile(r_re[,x], c(0.9)))

    f_median[,((j-1)*4+k)] <- sapply(1:ncol(f_re), function(x) quantile(f_re[,x], c(0.5)))
    f_low[,((j-1)*4+k)] <- sapply(1:ncol(f_re), function(x) quantile(f_re[,x], c(0.1)))
    f_high[,((j-1)*4+k)] <- sapply(1:ncol(f_re), function(x) quantile(f_re[,x], c(0.9)))

    ssbratio_median[,((j-1)*4+k)] <- sapply(1:ncol(ssbratio_re), function(x) quantile(ssbratio_re[,x], c(0.5)))
    ssbratio_low[,((j-1)*4+k)] <- sapply(1:ncol(ssbratio_re), function(x) quantile(ssbratio_re[,x], c(0.1)))
    ssbratio_high[,((j-1)*4+k)] <- sapply(1:ncol(ssbratio_re), function(x) quantile(ssbratio_re[,x], c(0.9)))

    fratio_median[,((j-1)*4+k)] <- sapply(1:ncol(fratio_re), function(x) quantile(fratio_re[,x], c(0.5)))
    fratio_low[,((j-1)*4+k)] <- sapply(1:ncol(fratio_re), function(x) quantile(fratio_re[,x], c(0.1)))
    fratio_high[,((j-1)*4+k)] <- sapply(1:ncol(fratio_re), function(x) quantile(fratio_re[,x], c(0.9)))
  }
}

jpeg(file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/SSB_R_F_RE_two_surveys_quantile.jpg", width=130, height=100, units="mm", res=600)

op<-par(no.readonly=TRUE)
par(op)
par(oma=c(5,5,0,5),mar=c(0,0,0,0.2),mfrow=c(2,5), pch=16)

x_val <- 1:length(ssb_median[,1])
ylim=c(-0.35, 0.35)

for (j in 1:length(maindir_list)){
  plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  axis(2, at=c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3), labels=c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3), las=2)
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*4+1], rev(ssb_high[,(j-1)*4+1])), border=NA, col=alpha("orange", 0.2))
  #lines(x_val, ssb_median[,1], type="o", lty=1, col="orange", pch=2)

  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*4+2], rev(ssb_high[,(j-1)*4+2])), border=NA, col=alpha("green", 0.2))
  #lines(x_val, ssb_median[,2], type="o", lty=1, col="green", pch=3)

  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*4+3], rev(ssb_high[,(j-1)*4+3])), border=NA, col=alpha("red", 0.2))
  #lines(x_val, ssb_median[,3], type="o", lty=1, col="red", pch=4)

  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*4+4], rev(ssb_high[,(j-1)*4+4])), border=NA, col=alpha("deepskyblue", 0.2))
  #lines(x_val, ssb_median[,4], type="o", lty=1, col="deepskyblue", pch=5)

  lines(x_val, ssb_median[,(j-1)*4+1], type="l", lty=2, col="orange")
  lines(x_val, ssb_median[,(j-1)*4+2], type="l", lty=3, col="green")
  lines(x_val, ssb_median[,(j-1)*4+3], type="l", lty=4, col="red")
  lines(x_val, ssb_median[,(j-1)*4+4], type="l", lty=5, col="deepskyblue")

  legend("topleft", paste(legend_names[j], "_SSB", sep=""), bty="n", cex=0.8)

  plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*4+1], rev(r_high[,(j-1)*4+1])), border=NA, col=alpha("orange", 0.2))
  #lines(x_val, r_median[,1], type="o", lty=1, col="orange", pch=2)

  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*4+2], rev(r_high[,(j-1)*4+2])), border=NA, col=alpha("green", 0.2))
  #lines(x_val, r_median[,2], type="o", lty=1, col="green", pch=3)

  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*4+3], rev(r_high[,(j-1)*4+3])), border=NA, col=alpha("red", 0.2))
  #lines(x_val, r_median[,3], type="o", lty=1, col="red", pch=4)

  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*4+4], rev(r_high[,(j-1)*4+4])), border=NA, col=alpha("deepskyblue", 0.2))
  #lines(x_val, r_median[,4], type="o", lty=1, col="deepskyblue", pch=5)

  lines(x_val, r_median[,(j-1)*4+1], type="l", lty=2, col="orange")
  lines(x_val, r_median[,(j-1)*4+2], type="l", lty=3, col="green")
  lines(x_val, r_median[,(j-1)*4+3], type="l", lty=4, col="red")
  lines(x_val, r_median[,(j-1)*4+4], type="l", lty=5, col="deepskyblue")

  legend("topleft", paste(legend_names[j], "_R", sep=""), bty="n", cex=0.8)

  plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*4+1], rev(f_high[,(j-1)*4+1])), border=NA, col=alpha("orange", 0.2))
  #lines(x_val, f_median[,1], type="o", lty=1, col="orange", pch=2)

  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*4+2], rev(f_high[,(j-1)*4+2])), border=NA, col=alpha("green", 0.2))
  #lines(x_val, f_median[,2], type="o", lty=1, col="green", pch=3)

  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*4+3], rev(f_high[,(j-1)*4+3])), border=NA, col=alpha("red", 0.2))
  #lines(x_val, f_median[,3], type="o", lty=1, col="red", pch=4)

  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*4+4], rev(f_high[,(j-1)*4+4])), border=NA, col=alpha("deepskyblue", 0.2))
  #lines(x_val, f_median[,4], type="o", lty=1, col="deepskyblue", pch=5)

  lines(x_val, f_median[,(j-1)*4+1], type="l", lty=2, col="orange")
  lines(x_val, f_median[,(j-1)*4+2], type="l", lty=3, col="green")
  lines(x_val, f_median[,(j-1)*4+3], type="l", lty=4, col="red")
  lines(x_val, f_median[,(j-1)*4+4], type="l", lty=5, col="deepskyblue")

  legend("topleft", paste(legend_names[j], "_F", sep=""), bty="n", cex=0.8)

  plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*4+1], rev(ssbratio_high[,(j-1)*4+1])), border=NA, col=alpha("orange", 0.2))
  #lines(x_val, ssbratio_median[,1], type="o", lty=1, col="orange", pch=2)

  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*4+2], rev(ssbratio_high[,(j-1)*4+2])), border=NA, col=alpha("green", 0.2))
  #lines(x_val, ssbratio_median[,2], type="o", lty=1, col="green", pch=3)

  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*4+3], rev(ssbratio_high[,(j-1)*4+3])), border=NA, col=alpha("red", 0.2))
  #lines(x_val, ssbratio_median[,3], type="o", lty=1, col="red", pch=4)

  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*4+4], rev(ssbratio_high[,(j-1)*4+4])), border=NA, col=alpha("deepskyblue", 0.2))
  #lines(x_val, ssbratio_median[,4], type="o", lty=1, col="deepskyblue", pch=5)

  lines(x_val, ssbratio_median[,(j-1)*4+1], type="l", lty=2, col="orange")
  lines(x_val, ssbratio_median[,(j-1)*4+2], type="l", lty=3, col="green")
  lines(x_val, ssbratio_median[,(j-1)*4+3], type="l", lty=4, col="red")
  lines(x_val, ssbratio_median[,(j-1)*4+4], type="l", lty=5, col="deepskyblue")

  legend("topleft", paste(legend_names[j], "_SSBratio", sep=""), bty="n", cex=0.8)

  plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*4+1], rev(fratio_high[,(j-1)*4+1])), border=NA, col=alpha("orange", 0.2))
  #lines(x_val, fratio_median[,1], type="o", lty=1, col="orange", pch=2)

  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*4+2], rev(fratio_high[,(j-1)*4+2])), border=NA, col=alpha("green", 0.2))
  #lines(x_val, fratio_median[,2], type="o", lty=1, col="green", pch=3)

  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*4+3], rev(fratio_high[,(j-1)*4+3])), border=NA, col=alpha("red", 0.2))
  #lines(x_val, fratio_median[,3], type="o", lty=1, col="red", pch=4)

  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*4+4], rev(fratio_high[,(j-1)*4+4])), border=NA, col=alpha("deepskyblue", 0.2))
  #lines(x_val, fratio_median[,4], type="o", lty=1, col="deepskyblue", pch=5)

  lines(x_val, fratio_median[,(j-1)*4+1], type="l", lty=2, col="orange")
  lines(x_val, fratio_median[,(j-1)*4+2], type="l", lty=3, col="green")
  lines(x_val, fratio_median[,(j-1)*4+3], type="l", lty=4, col="red")
  lines(x_val, fratio_median[,(j-1)*4+4], type="l", lty=5, col="deepskyblue")

  legend("topleft", paste(legend_names[j], "_Fratio", sep=""), bty="n", cex=0.8)
}

mtext(text="Year",side=1,line=2.5,outer=TRUE)
mtext(text="RE",side=2,line=2.5,outer=TRUE)

legend(x=30.5,y=1,
       legend=c("AMAK", "ASAP", "BAM", "SS"),
       lty=c(2, 3, 4, 5),
       col=c("orange", "green", "red", "deepskyblue3"),
       cex=1,
       bty="n",
       xpd=NA)
dev.off()
#### Quantile (C0, C0_AC, C12, and C12_AC) ####
maindir_list  <- c("C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case0",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case0_F0_0.01_1yr",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case12",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case12_two_surveys_F0_0.01")

em_num <- 4
legend_names <- c("C0", "C0_UP", "C12", "C12_UP")

load(file.path(maindir_list[1], "output", "performance_measures_output.RData"))
ssb_median <- ssb_low <- ssb_high <-
  r_median <- r_low <- r_high <-
  f_median <- f_low <- f_high <-
  ssbratio_median <- ssbratio_low <- ssbratio_high <-
  fratio_median <- fratio_low <- fratio_high <-
  matrix(NA, nrow=nrow(pm_list[[1]]$ssb$re), ncol=em_num*length(maindir_list))

for(j in 1:length(maindir_list)){
  load(file.path(maindir_list[j], "output", "performance_measures_output.RData"))
  if (j==1){
    ssb_re <- r_re <- f_re <-
      ssbratio_re <- fratio_re <-
      matrix(NA, nrow=length(pm_list), ncol=nrow(pm_list[[1]]$ssb$re))
  }
  for (k in 1:em_num){
    for (i in 1:length(pm_list)){
      ssb_re[i,] <- as.matrix(pm_list[[i]]$ssb$re[,k])
      r_re[i,] <- as.matrix(pm_list[[i]]$recruit$re[,k])
      f_re[i,] <- as.matrix(pm_list[[i]]$Ftot$re[,k])
      ssbratio_re[i,] <- as.matrix(pm_list[[i]]$ssbratio$re[,k])
      fratio_re[i,] <- as.matrix(pm_list[[i]]$fratio$re[,k])
    }
    ssb_median[,((j-1)*4+k)] <- sapply(1:ncol(ssb_re), function(x) quantile(ssb_re[,x], c(0.5)))
    ssb_low[,((j-1)*4+k)] <- sapply(1:ncol(ssb_re), function(x) quantile(ssb_re[,x], c(0.1)))
    ssb_high[,((j-1)*4+k)] <- sapply(1:ncol(ssb_re), function(x) quantile(ssb_re[,x], c(0.9)))

    r_median[,((j-1)*4+k)] <- sapply(1:ncol(r_re), function(x) quantile(r_re[,x], c(0.5)))
    r_low[,((j-1)*4+k)] <- sapply(1:ncol(r_re), function(x) quantile(r_re[,x], c(0.1)))
    r_high[,((j-1)*4+k)] <- sapply(1:ncol(r_re), function(x) quantile(r_re[,x], c(0.9)))

    f_median[,((j-1)*4+k)] <- sapply(1:ncol(f_re), function(x) quantile(f_re[,x], c(0.5)))
    f_low[,((j-1)*4+k)] <- sapply(1:ncol(f_re), function(x) quantile(f_re[,x], c(0.1)))
    f_high[,((j-1)*4+k)] <- sapply(1:ncol(f_re), function(x) quantile(f_re[,x], c(0.9)))

    ssbratio_median[,((j-1)*4+k)] <- sapply(1:ncol(ssbratio_re), function(x) quantile(ssbratio_re[,x], c(0.5)))
    ssbratio_low[,((j-1)*4+k)] <- sapply(1:ncol(ssbratio_re), function(x) quantile(ssbratio_re[,x], c(0.1)))
    ssbratio_high[,((j-1)*4+k)] <- sapply(1:ncol(ssbratio_re), function(x) quantile(ssbratio_re[,x], c(0.9)))

    fratio_median[,((j-1)*4+k)] <- sapply(1:ncol(fratio_re), function(x) quantile(fratio_re[,x], c(0.5)))
    fratio_low[,((j-1)*4+k)] <- sapply(1:ncol(fratio_re), function(x) quantile(fratio_re[,x], c(0.1)))
    fratio_high[,((j-1)*4+k)] <- sapply(1:ncol(fratio_re), function(x) quantile(fratio_re[,x], c(0.9)))
  }
}

jpeg(file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/SSB_R_F_RE_two_surveys_quantile_UP.jpg", width=130, height=100, units="mm", res=600)

op<-par(no.readonly=TRUE)
par(op)
par(oma=c(5,5,0,5),mar=c(0,0,0,0.2),mfrow=c(4,5), pch=16)

x_val <- 1:length(ssb_median[,1])
ylim=c(-0.35, 0.35)

for (j in 1:length(maindir_list)){
  #plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="white", xlab="", ylab="", axes=F)
  box()
  axis(2, at=c(-0.2, -0.1, 0, 0.1, 0.2), labels=c(-0.2, -0.1, 0, 0.1, 0.2), las=2)
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*4+1], rev(ssb_high[,(j-1)*4+1])), border=NA, col=alpha("orange", 0.2))
  #lines(x_val, ssb_median[,1], type="o", lty=1, col="orange", pch=2)

  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*4+2], rev(ssb_high[,(j-1)*4+2])), border=NA, col=alpha("green", 0.2))
  #lines(x_val, ssb_median[,2], type="o", lty=1, col="green", pch=3)

  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*4+3], rev(ssb_high[,(j-1)*4+3])), border=NA, col=alpha("red", 0.2))
  #lines(x_val, ssb_median[,3], type="o", lty=1, col="red", pch=4)

  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*4+4], rev(ssb_high[,(j-1)*4+4])), border=NA, col=alpha("deepskyblue", 0.2))
  #lines(x_val, ssb_median[,4], type="o", lty=1, col="deepskyblue", pch=5)

  lines(x_val, ssb_median[,(j-1)*4+1], type="l", lty=2, col="orange")
  lines(x_val, ssb_median[,(j-1)*4+2], type="l", lty=3, col="green")
  lines(x_val, ssb_median[,(j-1)*4+3], type="l", lty=4, col="red")
  lines(x_val, ssb_median[,(j-1)*4+4], type="l", lty=5, col="deepskyblue")

  legend("topleft", paste(legend_names[j], "_SSB", sep=""), bty="n", cex=0.6)

  #plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="white", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*4+1], rev(r_high[,(j-1)*4+1])), border=NA, col=alpha("orange", 0.2))
  #lines(x_val, r_median[,1], type="o", lty=1, col="orange", pch=2)

  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*4+2], rev(r_high[,(j-1)*4+2])), border=NA, col=alpha("green", 0.2))
  #lines(x_val, r_median[,2], type="o", lty=1, col="green", pch=3)

  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*4+3], rev(r_high[,(j-1)*4+3])), border=NA, col=alpha("red", 0.2))
  #lines(x_val, r_median[,3], type="o", lty=1, col="red", pch=4)

  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*4+4], rev(r_high[,(j-1)*4+4])), border=NA, col=alpha("deepskyblue", 0.2))
  #lines(x_val, r_median[,4], type="o", lty=1, col="deepskyblue", pch=5)

  lines(x_val, r_median[,(j-1)*4+1], type="l", lty=2, col="orange")
  lines(x_val, r_median[,(j-1)*4+2], type="l", lty=3, col="green")
  lines(x_val, r_median[,(j-1)*4+3], type="l", lty=4, col="red")
  lines(x_val, r_median[,(j-1)*4+4], type="l", lty=5, col="deepskyblue")

  legend("topleft", paste(legend_names[j], "_R", sep=""), bty="n", cex=0.6)

  #plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="white", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*4+1], rev(f_high[,(j-1)*4+1])), border=NA, col=alpha("orange", 0.2))
  #lines(x_val, f_median[,1], type="o", lty=1, col="orange", pch=2)

  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*4+2], rev(f_high[,(j-1)*4+2])), border=NA, col=alpha("green", 0.2))
  #lines(x_val, f_median[,2], type="o", lty=1, col="green", pch=3)

  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*4+3], rev(f_high[,(j-1)*4+3])), border=NA, col=alpha("red", 0.2))
  #lines(x_val, f_median[,3], type="o", lty=1, col="red", pch=4)

  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*4+4], rev(f_high[,(j-1)*4+4])), border=NA, col=alpha("deepskyblue", 0.2))
  #lines(x_val, f_median[,4], type="o", lty=1, col="deepskyblue", pch=5)

  lines(x_val, f_median[,(j-1)*4+1], type="l", lty=2, col="orange")
  lines(x_val, f_median[,(j-1)*4+2], type="l", lty=3, col="green")
  lines(x_val, f_median[,(j-1)*4+3], type="l", lty=4, col="red")
  lines(x_val, f_median[,(j-1)*4+4], type="l", lty=5, col="deepskyblue")

  legend("topleft", paste(legend_names[j], "_F", sep=""), bty="n", cex=0.6)

  #plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="white", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*4+1], rev(ssbratio_high[,(j-1)*4+1])), border=NA, col=alpha("orange", 0.2))
  #lines(x_val, ssbratio_median[,1], type="o", lty=1, col="orange", pch=2)

  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*4+2], rev(ssbratio_high[,(j-1)*4+2])), border=NA, col=alpha("green", 0.2))
  #lines(x_val, ssbratio_median[,2], type="o", lty=1, col="green", pch=3)

  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*4+3], rev(ssbratio_high[,(j-1)*4+3])), border=NA, col=alpha("red", 0.2))
  #lines(x_val, ssbratio_median[,3], type="o", lty=1, col="red", pch=4)

  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*4+4], rev(ssbratio_high[,(j-1)*4+4])), border=NA, col=alpha("deepskyblue", 0.2))
  #lines(x_val, ssbratio_median[,4], type="o", lty=1, col="deepskyblue", pch=5)

  lines(x_val, ssbratio_median[,(j-1)*4+1], type="l", lty=2, col="orange")
  lines(x_val, ssbratio_median[,(j-1)*4+2], type="l", lty=3, col="green")
  lines(x_val, ssbratio_median[,(j-1)*4+3], type="l", lty=4, col="red")
  lines(x_val, ssbratio_median[,(j-1)*4+4], type="l", lty=5, col="deepskyblue")

  legend("topleft", paste(legend_names[j], "_SSBratio", sep=""), bty="n", cex=0.6)

  plot(x_val, rep(0, times=length(x_val)), ylim=ylim, type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*4+1], rev(fratio_high[,(j-1)*4+1])), border=NA, col=alpha("orange", 0.2))
  #lines(x_val, fratio_median[,1], type="o", lty=1, col="orange", pch=2)

  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*4+2], rev(fratio_high[,(j-1)*4+2])), border=NA, col=alpha("green", 0.2))
  #lines(x_val, fratio_median[,2], type="o", lty=1, col="green", pch=3)

  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*4+3], rev(fratio_high[,(j-1)*4+3])), border=NA, col=alpha("red", 0.2))
  #lines(x_val, fratio_median[,3], type="o", lty=1, col="red", pch=4)

  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*4+4], rev(fratio_high[,(j-1)*4+4])), border=NA, col=alpha("deepskyblue", 0.2))
  #lines(x_val, fratio_median[,4], type="o", lty=1, col="deepskyblue", pch=5)

  lines(x_val, fratio_median[,(j-1)*4+1], type="l", lty=2, col="orange")
  lines(x_val, fratio_median[,(j-1)*4+2], type="l", lty=3, col="green")
  lines(x_val, fratio_median[,(j-1)*4+3], type="l", lty=4, col="red")
  lines(x_val, fratio_median[,(j-1)*4+4], type="l", lty=5, col="deepskyblue")

  legend("topleft", paste(legend_names[j], "_Fratio", sep=""), bty="n", cex=0.6)
}

mtext(text="Year",side=1,line=2.5,outer=TRUE)
mtext(text="RE",side=2,line=2.5,outer=TRUE)

legend(x=30.5,y=1,
       legend=c("AMAK", "ASAP", "BAM", "SS"),
       lty=c(2, 3, 4, 5),
       col=c("orange", "green", "red", "deepskyblue3"),
       cex=1,
       bty="n",
       xpd=NA)
dev.off()
#### MSY (C0, C0_AC, C12, and C12_AC) ####
maindir_list  <- c("C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case0",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case0_F0_0.01_1yr",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case12",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case12_two_surveys_F0_0.01")
#### MSY based estimates plot ####
em_num <- 4
legend_names <- c("C0", "C0_UP", "C12", "C12_UP")

for(j in 1:length(maindir_list)){
  load(file.path(maindir_list[j], "output", "performance_measures_output.RData"))
  if (j==1){
    msy_re <- matrix(NA, nrow=length(pm_list), ncol=em_num*length(maindir_list))
    fmsy_re <- matrix(NA, nrow=length(pm_list), ncol=em_num*length(maindir_list))
    ssbmsy_re <- matrix(NA, nrow=length(pm_list), ncol=em_num*length(maindir_list))
  }
  for (i in 1:length(pm_list)){
    msy_re[i,((j-1)*4+1):((j-1)*4+4)] <- as.matrix(pm_list[[i]]$msy$re)
    fmsy_re[i,((j-1)*4+1):((j-1)*4+4)] <- as.matrix(pm_list[[i]]$fmsy$re)
    ssbmsy_re[i,((j-1)*4+1):((j-1)*4+4)] <- as.matrix(pm_list[[i]]$ssbmsy$re)
  }
}

msy_stat <- matrix(NA, nrow=3, ncol=em_num*length(maindir_list))
fmsy_stat <- matrix(NA, nrow=3, ncol=em_num*length(maindir_list))
ssbmsy_stat <- matrix(NA, nrow=3, ncol=em_num*length(maindir_list))

sapply(1:ncol(msy_stat), function(x) {
  msy_stat[,x] <<- boxplot.stats(msy_re[,x])$`stats`[c(1,3,5)]
  fmsy_stat[,x] <<- boxplot.stats(fmsy_re[,x])$`stats`[c(1,3,5)]
  ssbmsy_stat[,x] <<- boxplot.stats(ssbmsy_re[,x])$`stats`[c(1,3,5)]
})

# round(msy_stat[2,c(37:38, 41:42)], digits = 2)
# round(fmsy_stat[2,c(37:38, 41:42)], digits = 2)
# round(ssbmsy_stat[2,c(37:38, 41:42)], digits = 2)

jpeg(file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/MSY_RE_two_surveys.jpg", width=150, height=150, units="mm", res=300)
op<-par(no.readonly=TRUE)
par(op)
par(oma=c(2,2,0,4),mar=c(4,0,0,0.5),mfrow=c(1,3),pch=16)
color=c("orange", "green", "red", "deepskyblue3")
xlim=c(-0.5, 0.5)

plot(msy_stat[2,], rev(1:ncol(msy_stat)),
     pch=rep(((1:em_num)+1), times=length(maindir_list)),
     col=color,
     xlim=xlim,
     xlab="RE in MSY", ylab="",
     axes=F)
segments(msy_stat[1,], rev(1:ncol(msy_stat)), msy_stat[3,], rev(1:ncol(msy_stat)), col=rep(color, times=length(maindir_list)))
abline(v=0, col="gray", lty=2)
abline(h=seq(em_num+0.5, ncol(msy_stat), by=em_num), col="gray30", lty=2)
box()
axis(1)
text(x=rep(min(xlim)*0.8, times=length(maindir_list)),
     y=seq(em_num, ncol(msy_stat), by=4),
     rev(legend_names))

plot(fmsy_stat[2,], rev(1:ncol(fmsy_stat)),
     pch=rep(((1:em_num)+1), times=length(maindir_list)),
     col=color,
     xlim=xlim,
     xlab="RE in FMSY", ylab="",
     axes=F)
segments(fmsy_stat[1,], rev(1:ncol(fmsy_stat)), fmsy_stat[3,], rev(1:ncol(fmsy_stat)), col=rep(color, times=length(maindir_list)))
abline(v=0, col="gray", lty=2)
abline(h=seq(em_num+0.5, ncol(fmsy_stat), by=em_num), col="gray30", lty=2)
box()
axis(1)
text(x=rep(min(xlim)*0.8, times=length(maindir_list)),
     y=seq(em_num, ncol(fmsy_stat), by=4),
     rev(legend_names))

plot(ssbmsy_stat[2,], rev(1:ncol(ssbmsy_stat)),
     pch=rep(((1:em_num)+1), times=length(maindir_list)),
     col=color,
     xlim=xlim,
     xlab="RE in SSBMSY", ylab="",
     axes=F)
segments(ssbmsy_stat[1,], rev(1:ncol(ssbmsy_stat)), ssbmsy_stat[3,], rev(1:ncol(ssbmsy_stat)), col=rep(color, times=length(maindir_list)))
abline(v=0, col="gray", lty=2)
abline(h=seq(em_num+0.5, ncol(ssbmsy_stat), by=em_num), col="gray30", lty=2)
box()
axis(1)
text(x=rep(min(xlim)*0.8, times=length(maindir_list)),
     y=seq(em_num, ncol(ssbmsy_stat), by=4),
     rev(legend_names))

legend(x=max(xlim)*1.05,y=ncol(msy_stat)/2*1.2,
       legend=c("AMAK", "ASAP", "BAM", "SS"),
       pch=c(2, 3, 4, 5),
       lty=c(1, 1, 1, 1),
       col=c("orange", "green", "red", "deepskyblue3"),
       cex=0.9,
       bty="n",
       xpd=NA)
dev.off()

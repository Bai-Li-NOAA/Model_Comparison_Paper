#### SSB, R, and F ####
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
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case0_F0_0.01_1yr")

em_num <- 4
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
      ssbratio_re <- fratio_re <- matrix(NA, nrow=length(pm_list), ncol=nrow(pm_list[[1]]$ssb$re))
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

#### RE Figure ####
jpeg(file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/SSB_R_F_RE.jpg", width=130, height=150, units="mm", res=1200)

op<-par(no.readonly=TRUE)
par(op)
par(oma=c(5,5,0,5),mar=c(0,0,0,0.2),mfrow=c(13,5), pch=16)

x_val <- 1:length(ssb_median[,1])
ylim=c(-0.7, 0.7)

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

  if (j == length(maindir_list)) {
    legend("topleft", paste("C12_SSB", sep=""), bty="n", cex=0.8)
  } else {
    legend("topleft", paste("C", j-1, "_SSB", sep=""), bty="n", cex=0.8)
  }

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

  if (j == length(maindir_list)) {
    legend("topleft", paste("C12_R", sep=""), bty="n", cex=0.8)
  } else {
    legend("topleft", paste("C", j-1, "_R", sep=""), bty="n", cex=0.8)
  }
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
  if (j == length(maindir_list)) {
    legend("topleft", paste("C12_F", sep=""), bty="n", cex=0.8)
  } else {
    legend("topleft", paste("C", j-1, "_F", sep=""), bty="n", cex=0.8)
  }
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

  if (j == length(maindir_list)) {
    legend("topleft", paste("C12_SSBratio", sep=""), bty="n", cex=0.8)
  } else {
    legend("topleft", paste("C", j-1, "_SSBratio", sep=""), bty="n", cex=0.8)
  }

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

  if (j == length(maindir_list)) {
    legend("topleft", paste("C12_Fratio", sep=""), bty="n", cex=0.8)
  } else {
    legend("topleft", paste("C", j-1, "_Fratio", sep=""), bty="n", cex=0.8)
  }
}

mtext(text="Year",side=1,line=2.5,outer=TRUE)
mtext(text="RE",side=2,line=2.5,outer=TRUE)

legend(x=30.5,y=6,
       legend=c("AMAK", "ASAP", "BAM", "SS"),
       lty=c(2, 3, 4, 5),
       col=c("orange", "green", "red", "deepskyblue3"),
       cex=1,
       bty="n",
       xpd=NA)
dev.off()


#### SSB, R, and F vioplot ####
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
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case0_F0_0.01_1yr")

em_num <- 4

ssb_re_list <- r_re_list <- f_re_list <- ssbratio_re_list <- fratio_re_list <- list()

for(j in 1:length(maindir_list)){
  load(file.path(maindir_list[j], "output", "performance_measures_output.RData"))
  if (j==1){
    ssb_re <- r_re <- f_re <-
      ssbratio_re <- fratio_re <- matrix(NA, nrow=length(pm_list), ncol=nrow(pm_list[[1]]$ssb$re))
    ssb <- r <- f <- ssbratio <- fratio <- matrix(NA, nrow=length(pm_list)*nrow(pm_list[[1]]$ssb$re), ncol=em_num)
  }
  for (k in 1:em_num){
    for (i in 1:length(pm_list)){
      ssb_re[i,] <- as.matrix(pm_list[[i]]$ssb$re[,k])
      r_re[i,] <- as.matrix(pm_list[[i]]$recruit$re[,k])
      f_re[i,] <- as.matrix(pm_list[[i]]$Ftot$re[,k])
      ssbratio_re[i,] <- as.matrix(pm_list[[i]]$ssbratio$re[,k])
      fratio_re[i,] <- as.matrix(pm_list[[i]]$fratio$re[,k])
    }
    ssb[,k] <- as.vector(ssb_re)
    r[,k] <- as.vector(r_re)
    f[,k] <- as.vector(f_re)
    ssbratio[,k] <- as.vector(ssbratio_re)
    fratio[,k] <- as.vector(fratio_re)
  }
  ssb_re_list[[j]] <- ssb
  r_re_list[[j]] <- r
  f_re_list[[j]] <- f
  ssbratio_re_list[[j]] <- ssbratio
  fratio_re_list[[j]] <- fratio
}


library(vioplot)

jpeg(file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/violinSSB_R_F_RE.jpg", width=130, height=150, units="mm", res=1200)
op<-par(no.readonly=TRUE)
par(op)
par(oma=c(5,5,0,5),mar=c(0,0,0,0.2),mfrow=c(13,5), pch=16)

ylim=c(-0.8, 0.8)

for (i in 1:length(maindir_list)){
  plot(NA,
       type='n',
       xlim=range(0:(em_num+1)),
       ylim=ylim,
       axes=FALSE, ann=FALSE
  )
  axis(2, at=c(-0.4, 0, 0.4), labels=c(-0.4, 0, 0.4), las=2)
  vioplot(ssb_re_list[[i]][,1], ssb_re_list[[i]][,2], ssb_re_list[[i]][,3], ssb_re_list[[i]][,4], add=TRUE)
  abline(h=0, col="gray30", lty=2)
  if(i == length(maindir_list)) axis(1, at=1:em_num, labels = c("AMAK", "ASAP", "BAM", "SS"), las=2, cex=0.6)

  if(i == length(maindir_list)) {
    legend("topleft", paste("C12_SSB", sep=""), bty="n", cex=0.7)
  } else {
    legend("topleft", paste("C", i-1, "_SSB", sep=""), bty="n", cex=0.7)
  }

  plot(NA,
       type='n',
       xlim=range(0:(em_num+1)),
       ylim=ylim,
       axes=FALSE, ann=FALSE
  )
  vioplot(r_re_list[[i]][,1], r_re_list[[i]][,2], r_re_list[[i]][,3], r_re_list[[i]][,4], add=TRUE)
  abline(h=0, col="gray30", lty=2)
  if(i == length(maindir_list)) axis(1, at=1:em_num, labels = c("AMAK", "ASAP", "BAM", "SS"), las=2, cex=0.6)

  if(i == length(maindir_list)) {
    legend("topleft", paste("C12_R", sep=""), bty="n", cex=0.7)
  } else {
    legend("topleft", paste("C", i-1, "_R", sep=""), bty="n", cex=0.7)
  }
  plot(NA,
       type='n',
       xlim=range(0:(em_num+1)),
       ylim=ylim,
       axes=FALSE, ann=FALSE
  )
  vioplot(f_re_list[[i]][,1], f_re_list[[i]][,2], f_re_list[[i]][,3], f_re_list[[i]][,4], add=TRUE)
  abline(h=0, col="gray30", lty=2)
  if(i == length(maindir_list)) axis(1, at=1:em_num, labels = c("AMAK", "ASAP", "BAM", "SS"), las=2, cex=0.6)
  if(i == length(maindir_list)) {
    legend("topleft", paste("C12_F", sep=""), bty="n", cex=0.7)
  } else {
    legend("topleft", paste("C", i-1, "_F", sep=""), bty="n", cex=0.7)
  }
  plot(NA,
       type='n',
       xlim=range(0:(em_num+1)),
       ylim=ylim,
       axes=FALSE, ann=FALSE
  )
  vioplot(ssbratio_re_list[[i]][,1], ssbratio_re_list[[i]][,2], ssbratio_re_list[[i]][,3], ssbratio_re_list[[i]][,4], add=TRUE)
  abline(h=0, col="gray30", lty=2)
  if(i == length(maindir_list)) axis(1, at=1:em_num, labels = c("AMAK", "ASAP", "BAM", "SS"), las=2, cex=0.6)
  if(i == length(maindir_list)) {
    legend("topleft", paste("C12_SSBratio", sep=""), bty="n", cex=0.7)
  } else {
    legend("topleft", paste("C", i-1, "_SSBratio", sep=""), bty="n", cex=0.7)
  }
  plot(NA,
       type='n',
       xlim=range(0:(em_num+1)),
       ylim=ylim,
       axes=FALSE, ann=FALSE
  )
  vioplot(fratio_re_list[[i]][,1], fratio_re_list[[i]][,2], fratio_re_list[[i]][,3], fratio_re_list[[i]][,4], add=TRUE)
  abline(h=0, col="gray30", lty=2)
  if(i == length(maindir_list)) axis(1, at=1:em_num, labels = c("AMAK", "ASAP", "BAM", "SS"), las=2, cex=0.6)
  if(i == length(maindir_list)) {
    legend("topleft", paste("C12_Fratio", sep=""), bty="n", cex=0.7)
  } else {
    legend("topleft", paste("C", i-1, "_Fratio", sep=""), bty="n", cex=0.7)
  }
}

mtext(text="RE",side=2,line=2.5,outer=TRUE)

dev.off()




#### SSB, R, and F time series ####
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
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case0_F0_0.01_1yr")

em_num <- 4
load(file.path(maindir_list[1], "output", "om_em_output.RData"))
ssb_median <- ssb_low <- ssb_high <-
  r_median <- r_low <- r_high <-
  f_median <- f_low <- f_high <-
  ssbratio_median <- ssbratio_low <- ssbratio_high <-
  fratio_median <- fratio_low <- fratio_high <-
  matrix(NA, nrow=nrow(om_list$ssb), ncol=(em_num+1)*length(maindir_list))

for(j in 1:length(maindir_list)){
  load(file.path(maindir_list[j], "output", "om_em_output.RData"))
  data_list <- list(om_list, amak_list, asap_list, bam_list, ss_list)
  for (k in 1:(em_num+1)){
    ssb_median[,((j-1)*(em_num+1)+k)] <- sapply(1:nrow(data_list[[k]]$ssb), function(x) boxplot.stats(data_list[[k]]$ssb[x,])$`stats`[3])
    ssb_low[,((j-1)*(em_num+1)+k)] <- sapply(1:nrow(data_list[[k]]$ssb), function(x) boxplot.stats(data_list[[k]]$ssb[x,])$`stats`[1])
    ssb_high[,((j-1)*(em_num+1)+k)] <- sapply(1:nrow(data_list[[k]]$ssb), function(x) boxplot.stats(data_list[[k]]$ssb[x,])$`stats`[5])

    r_median[,((j-1)*(em_num+1)+k)] <- sapply(1:nrow(data_list[[k]]$recruit), function(x) boxplot.stats(data_list[[k]]$recruit[x,])$`stats`[3])
    r_low[,((j-1)*(em_num+1)+k)] <- sapply(1:nrow(data_list[[k]]$recruit), function(x) boxplot.stats(data_list[[k]]$recruit[x,])$`stats`[1])
    r_high[,((j-1)*(em_num+1)+k)] <- sapply(1:nrow(data_list[[k]]$recruit), function(x) boxplot.stats(data_list[[k]]$recruit[x,])$`stats`[5])

    f_median[,((j-1)*(em_num+1)+k)] <- sapply(1:nrow(data_list[[k]]$Ftot), function(x) boxplot.stats(data_list[[k]]$Ftot[x,])$`stats`[3])
    f_low[,((j-1)*(em_num+1)+k)] <- sapply(1:nrow(data_list[[k]]$Ftot), function(x) boxplot.stats(data_list[[k]]$Ftot[x,])$`stats`[1])
    f_high[,((j-1)*(em_num+1)+k)] <- sapply(1:nrow(data_list[[k]]$Ftot), function(x) boxplot.stats(data_list[[k]]$Ftot[x,])$`stats`[5])

    ssbratio_median[,((j-1)*(em_num+1)+k)] <- sapply(1:nrow(data_list[[k]]$ssbratio), function(x) boxplot.stats(data_list[[k]]$ssbratio[x,])$`stats`[3])
    ssbratio_low[,((j-1)*(em_num+1)+k)] <- sapply(1:nrow(data_list[[k]]$ssbratio), function(x) boxplot.stats(data_list[[k]]$ssbratio[x,])$`stats`[1])
    ssbratio_high[,((j-1)*(em_num+1)+k)] <- sapply(1:nrow(data_list[[k]]$ssbratio), function(x) boxplot.stats(data_list[[k]]$ssbratio[x,])$`stats`[5])

    fratio_median[,((j-1)*(em_num+1)+k)] <- sapply(1:nrow(data_list[[k]]$fratio), function(x) boxplot.stats(data_list[[k]]$fratio[x,])$`stats`[3])
    fratio_low[,((j-1)*(em_num+1)+k)] <- sapply(1:nrow(data_list[[k]]$fratio), function(x) boxplot.stats(data_list[[k]]$fratio[x,])$`stats`[1])
    fratio_high[,((j-1)*(em_num+1)+k)] <- sapply(1:nrow(data_list[[k]]$fratio), function(x) boxplot.stats(data_list[[k]]$fratio[x,])$`stats`[5])
  }
}

#### SSB ####
jpeg(file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/SSB.jpg", width=130, height=150, units="mm", res=1200)
op<-par(no.readonly=TRUE)
par(op)
par(oma=c(5,5,0,5),mar=c(0,0,0,0.2),mfrow=c(13,4), pch=16)

x_val <- 1:30

for (j in 1:length(maindir_list)){
  plot(x_val, rep(0, times=length(x_val)), ylim=range(ssb_low, ssb_high), type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  box()
  axis(2, at=c(3000, 7000, 11000), labels=c(3000, 7000, 11000), las=2)
  if (j == length(maindir_list)) {
    legend("topright", paste("C12: AMAK", sep=""), bty="n", cex=0.8)
  } else {
    legend("topright", paste("C", j-1, ": AMAK", sep=""), bty="n", cex=0.8)
  }
  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*(em_num+1)+1], rev(ssb_high[,(j-1)*(em_num+1)+1])), border=NA, col=alpha("gray30", 0.3))
  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*(em_num+1)+2], rev(ssb_high[,(j-1)*(em_num+1)+2])), border=NA, col=alpha("orange", 0.3))
  lines(x_val, ssb_median[,(j-1)*(em_num+1)+2], type="o", lty=1, col="orange", pch=19, cex=0.2)
  lines(x_val, ssb_median[,(j-1)*(em_num+1)+1], type="o", lty=1, col="gray30", pch=19, cex=0.2)


  plot(x_val, rep(0, times=length(x_val)), ylim=range(ssb_low, ssb_high), type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  if (j == length(maindir_list)) {
    legend("topright", paste("C12: ASAP", sep=""), bty="n", cex=0.8)
  } else {
    legend("topright", paste("C", j-1, ": ASAP", sep=""), bty="n", cex=0.8)
  }
  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*(em_num+1)+1], rev(ssb_high[,(j-1)*(em_num+1)+1])), border=NA, col=alpha("gray30", 0.4))
  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*(em_num+1)+3], rev(ssb_high[,(j-1)*(em_num+1)+3])), border=NA, col=alpha("green", 0.4))
  lines(x_val, ssb_median[,(j-1)*(em_num+1)+3], type="o", lty=1, col="green", pch=19, cex=0.2)
  lines(x_val, ssb_median[,(j-1)*(em_num+1)+1], type="o", lty=1, col="gray30", pch=19, cex=0.2)

  plot(x_val, rep(0, times=length(x_val)), ylim=range(ssb_low, ssb_high), type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  if (j == length(maindir_list)) {
    legend("topright", paste("C12: BAM", sep=""), bty="n", cex=0.8)
  } else {
    legend("topright", paste("C", j-1, ": BAM", sep=""), bty="n", cex=0.8)
  }
  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*(em_num+1)+1], rev(ssb_high[,(j-1)*(em_num+1)+1])), border=NA, col=alpha("gray30", 0.4))
  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*(em_num+1)+4], rev(ssb_high[,(j-1)*(em_num+1)+4])), border=NA, col=alpha("red", 0.4))
  lines(x_val, ssb_median[,(j-1)*(em_num+1)+4], type="o", lty=1, col="red", pch=19, cex=0.2)
  lines(x_val, ssb_median[,(j-1)*(em_num+1)+1], type="o", lty=1, col="gray30", pch=19, cex=0.2)

  plot(x_val, rep(0, times=length(x_val)), ylim=range(ssb_low, ssb_high), type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  if (j == length(maindir_list)) {
    legend("topright", paste("C12: SS", sep=""), bty="n", cex=0.8)
  } else {
    legend("topright", paste("C", j-1, ": SS", sep=""), bty="n", cex=0.8)
  }
  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*(em_num+1)+1], rev(ssb_high[,(j-1)*(em_num+1)+1])), border=NA, col=alpha("gray30", 0.4))
  polygon(c(x_val, rev(x_val)), c(ssb_low[,(j-1)*(em_num+1)+5], rev(ssb_high[,(j-1)*(em_num+1)+5])), border=NA, col=alpha("deepskyblue3", 0.4))
  lines(x_val, ssb_median[,(j-1)*(em_num+1)+5], type="o", lty=1, col="deepskyblue3", pch=19, cex=0.2)
  lines(x_val, ssb_median[,(j-1)*(em_num+1)+1], type="o", lty=1, col="gray30", pch=19, cex=0.2)
}
mtext(text="Year",side=1,line=2.5,outer=TRUE)
mtext(text="SSB (mt)",side=2,line=3.5,outer=TRUE)

legend(x=30.5,y=80000,
       legend=c("OM", "AMAK", "ASAP", "BAM", "SS"),
       # lty=c(1, 1, 1, 1),
       # pch=c(19, 19, 19, 19),
       col=c("gray30", "orange", "green", "red", "deepskyblue3"),
       fill = c("gray30", "orange", "green", "red", "deepskyblue3"),
       cex=0.9,
       bty="n",
       xpd=NA)
dev.off()

#### R ####
jpeg(file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/R.jpg", width=130, height=150, units="mm", res=1200)
op<-par(no.readonly=TRUE)
par(op)
par(oma=c(5,5,0,5),mar=c(0,0,0,0.2),mfrow=c(13,4), pch=16)

x_val <- 1:30

for (j in 1:length(maindir_list)){
  plot(x_val, rep(0, times=length(x_val)), ylim=range(r_low, r_high), type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  box()
  axis(2, at=c(800, 1600, 2400), labels=c(800, 1600, 2400), las=2)
  if (j == length(maindir_list)) {
    legend("topright", paste("C12: AMAK", sep=""), bty="n", cex=0.8)
  } else {
    legend("topright", paste("C", j-1, ": AMAK", sep=""), bty="n", cex=0.8)
  }
  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*(em_num+1)+1], rev(r_high[,(j-1)*(em_num+1)+1])), border=NA, col=alpha("gray30", 0.3))
  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*(em_num+1)+2], rev(r_high[,(j-1)*(em_num+1)+2])), border=NA, col=alpha("orange", 0.3))
  lines(x_val, r_median[,(j-1)*(em_num+1)+2], type="o", lty=1, col="orange", pch=19, cex=0.2)
  lines(x_val, r_median[,(j-1)*(em_num+1)+1], type="o", lty=1, col="gray30", pch=19, cex=0.2)


  plot(x_val, rep(0, times=length(x_val)), ylim=range(r_low, r_high), type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  if (j == length(maindir_list)) {
    legend("topright", paste("C12: ASAP", sep=""), bty="n", cex=0.8)
  } else {
    legend("topright", paste("C", j-1, ": ASAP", sep=""), bty="n", cex=0.8)
  }
  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*(em_num+1)+1], rev(r_high[,(j-1)*(em_num+1)+1])), border=NA, col=alpha("gray30", 0.4))
  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*(em_num+1)+3], rev(r_high[,(j-1)*(em_num+1)+3])), border=NA, col=alpha("green", 0.4))
  lines(x_val, r_median[,(j-1)*(em_num+1)+3], type="o", lty=1, col="green", pch=19, cex=0.2)
  lines(x_val, r_median[,(j-1)*(em_num+1)+1], type="o", lty=1, col="gray30", pch=19, cex=0.2)

  plot(x_val, rep(0, times=length(x_val)), ylim=range(r_low, r_high), type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  if (j == length(maindir_list)) {
    legend("topright", paste("C12: BAM", sep=""), bty="n", cex=0.8)
  } else {
    legend("topright", paste("C", j-1, ": BAM", sep=""), bty="n", cex=0.8)
  }
  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*(em_num+1)+1], rev(r_high[,(j-1)*(em_num+1)+1])), border=NA, col=alpha("gray30", 0.4))
  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*(em_num+1)+4], rev(r_high[,(j-1)*(em_num+1)+4])), border=NA, col=alpha("red", 0.4))
  lines(x_val, r_median[,(j-1)*(em_num+1)+4], type="o", lty=1, col="red", pch=19, cex=0.2)
  lines(x_val, r_median[,(j-1)*(em_num+1)+1], type="o", lty=1, col="gray30", pch=19, cex=0.2)

  plot(x_val, rep(0, times=length(x_val)), ylim=range(r_low, r_high), type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  if (j == length(maindir_list)) {
    legend("topright", paste("C12: SS", sep=""), bty="n", cex=0.8)
  } else {
    legend("topright", paste("C", j-1, ": SS", sep=""), bty="n", cex=0.8)
  }
  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*(em_num+1)+1], rev(r_high[,(j-1)*(em_num+1)+1])), border=NA, col=alpha("gray30", 0.4))
  polygon(c(x_val, rev(x_val)), c(r_low[,(j-1)*(em_num+1)+5], rev(r_high[,(j-1)*(em_num+1)+5])), border=NA, col=alpha("deepskyblue3", 0.4))
  lines(x_val, r_median[,(j-1)*(em_num+1)+5], type="o", lty=1, col="deepskyblue3", pch=19, cex=0.2)
  lines(x_val, r_median[,(j-1)*(em_num+1)+1], type="o", lty=1, col="gray30", pch=19, cex=0.2)
}
mtext(text="Year",side=1,line=2.5,outer=TRUE)
mtext(text="R (×1000 fish)",side=2,line=3.5,outer=TRUE)

legend(x=30.5,y=20000,
       legend=c("OM", "AMAK", "ASAP", "BAM", "SS"),
       # lty=c(1, 1, 1, 1),
       # pch=c(19, 19, 19, 19),
       col=c("gray30", "orange", "green", "red", "deepskyblue3"),
       fill = c("gray30", "orange", "green", "red", "deepskyblue3"),
       cex=0.9,
       bty="n",
       xpd=NA)
dev.off()

#### F ####
jpeg(file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/F.jpg", width=130, height=150, units="mm", res=1200)
op<-par(no.readonly=TRUE)
par(op)
par(oma=c(5,5,0,5),mar=c(0,0,0,0.2),mfrow=c(13,4), pch=16)

x_val <- 1:30

for (j in 1:length(maindir_list)){
  plot(x_val, rep(0, times=length(x_val)), ylim=range(f_low, f_high), type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  box()
  axis(2, at=c(0.2, 0.4, 0.6), labels=c(0.2, 0.4, 0.6), las=2)
  if (j == length(maindir_list)) {
    legend("topright", paste("C12: AMAK", sep=""), bty="n", cex=0.8)
  } else {
    legend("topright", paste("C", j-1, ": AMAK", sep=""), bty="n", cex=0.8)
  }
  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*(em_num+1)+1], rev(f_high[,(j-1)*(em_num+1)+1])), border=NA, col=alpha("gray30", 0.3))
  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*(em_num+1)+2], rev(f_high[,(j-1)*(em_num+1)+2])), border=NA, col=alpha("orange", 0.3))
  lines(x_val, f_median[,(j-1)*(em_num+1)+2], type="o", lty=1, col="orange", pch=19, cex=0.2)
  lines(x_val, f_median[,(j-1)*(em_num+1)+1], type="o", lty=1, col="gray30", pch=19, cex=0.2)


  plot(x_val, rep(0, times=length(x_val)), ylim=range(f_low, f_high), type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  if (j == length(maindir_list)) {
    legend("topright", paste("C12: ASAP", sep=""), bty="n", cex=0.8)
  } else {
    legend("topright", paste("C", j-1, ": ASAP", sep=""), bty="n", cex=0.8)
  }
  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*(em_num+1)+1], rev(f_high[,(j-1)*(em_num+1)+1])), border=NA, col=alpha("gray30", 0.4))
  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*(em_num+1)+3], rev(f_high[,(j-1)*(em_num+1)+3])), border=NA, col=alpha("green", 0.4))
  lines(x_val, f_median[,(j-1)*(em_num+1)+3], type="o", lty=1, col="green", pch=19, cex=0.2)
  lines(x_val, f_median[,(j-1)*(em_num+1)+1], type="o", lty=1, col="gray30", pch=19, cex=0.2)

  plot(x_val, rep(0, times=length(x_val)), ylim=range(f_low, f_high), type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  if (j == length(maindir_list)) {
    legend("topright", paste("C12: BAM", sep=""), bty="n", cex=0.8)
  } else {
    legend("topright", paste("C", j-1, ": BAM", sep=""), bty="n", cex=0.8)
  }
  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*(em_num+1)+1], rev(f_high[,(j-1)*(em_num+1)+1])), border=NA, col=alpha("gray30", 0.4))
  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*(em_num+1)+4], rev(f_high[,(j-1)*(em_num+1)+4])), border=NA, col=alpha("red", 0.4))
  lines(x_val, f_median[,(j-1)*(em_num+1)+4], type="o", lty=1, col="red", pch=19, cex=0.2)
  lines(x_val, f_median[,(j-1)*(em_num+1)+1], type="o", lty=1, col="gray30", pch=19, cex=0.2)

  plot(x_val, rep(0, times=length(x_val)), ylim=range(f_low, f_high), type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  if (j == length(maindir_list)) {
    legend("topright", paste("C12: SS", sep=""), bty="n", cex=0.8)
  } else {
    legend("topright", paste("C", j-1, ": SS", sep=""), bty="n", cex=0.8)
  }
  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*(em_num+1)+1], rev(f_high[,(j-1)*(em_num+1)+1])), border=NA, col=alpha("gray30", 0.4))
  polygon(c(x_val, rev(x_val)), c(f_low[,(j-1)*(em_num+1)+5], rev(f_high[,(j-1)*(em_num+1)+5])), border=NA, col=alpha("deepskyblue3", 0.4))
  lines(x_val, f_median[,(j-1)*(em_num+1)+5], type="o", lty=1, col="deepskyblue3", pch=19, cex=0.2)
  lines(x_val, f_median[,(j-1)*(em_num+1)+1], type="o", lty=1, col="gray30", pch=19, cex=0.2)
}
mtext(text="Year",side=1,line=2.5,outer=TRUE)
mtext(text="F",side=2,line=3.5,outer=TRUE)

legend(x=30.5,y=5,
       legend=c("OM", "AMAK", "ASAP", "BAM", "SS"),
       # lty=c(1, 1, 1, 1),
       # pch=c(19, 19, 19, 19),
       col=c("gray30", "orange", "green", "red", "deepskyblue3"),
       fill = c("gray30", "orange", "green", "red", "deepskyblue3"),
       cex=0.9,
       bty="n",
       xpd=NA)
dev.off()

#### SSBratio ####
jpeg(file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/SSBratio.jpg", width=130, height=150, units="mm", res=1200)
op<-par(no.readonly=TRUE)
par(op)
par(oma=c(5,5,0,5),mar=c(0,0,0,0.2),mfrow=c(13,4), pch=16)

x_val <- 1:30

for (j in 1:length(maindir_list)){
  plot(x_val, rep(0, times=length(x_val)), ylim=range(ssbratio_low, ssbratio_high), type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  box()
  axis(2, at=c(1.5, 3, 4.5), labels=c(1.5, 3, 4.5), las=2)
  if (j == length(maindir_list)) {
    legend("topright", paste("C12: AMAK", sep=""), bty="n", cex=0.8)
  } else {
    legend("topright", paste("C", j-1, ": AMAK", sep=""), bty="n", cex=0.8)
  }
  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*(em_num+1)+1], rev(ssbratio_high[,(j-1)*(em_num+1)+1])), border=NA, col=alpha("gray30", 0.3))
  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*(em_num+1)+2], rev(ssbratio_high[,(j-1)*(em_num+1)+2])), border=NA, col=alpha("orange", 0.3))
  lines(x_val, ssbratio_median[,(j-1)*(em_num+1)+2], type="o", lty=1, col="orange", pch=19, cex=0.2)
  lines(x_val, ssbratio_median[,(j-1)*(em_num+1)+1], type="o", lty=1, col="gray30", pch=19, cex=0.2)
  abline(h=1, col="gray", lty=2)


  plot(x_val, rep(0, times=length(x_val)), ylim=range(ssbratio_low, ssbratio_high), type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  if (j == length(maindir_list)) {
    legend("topright", paste("C12: ASAP", sep=""), bty="n", cex=0.8)
  } else {
    legend("topright", paste("C", j-1, ": ASAP", sep=""), bty="n", cex=0.8)
  }
  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*(em_num+1)+1], rev(ssbratio_high[,(j-1)*(em_num+1)+1])), border=NA, col=alpha("gray30", 0.4))
  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*(em_num+1)+3], rev(ssbratio_high[,(j-1)*(em_num+1)+3])), border=NA, col=alpha("green", 0.4))
  lines(x_val, ssbratio_median[,(j-1)*(em_num+1)+3], type="o", lty=1, col="green", pch=19, cex=0.2)
  lines(x_val, ssbratio_median[,(j-1)*(em_num+1)+1], type="o", lty=1, col="gray30", pch=19, cex=0.2)
  abline(h=1, col="gray", lty=2)

  plot(x_val, rep(0, times=length(x_val)), ylim=range(ssbratio_low, ssbratio_high), type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  if (j == length(maindir_list)) {
    legend("topright", paste("C12: BAM", sep=""), bty="n", cex=0.8)
  } else {
    legend("topright", paste("C", j-1, ": BAM", sep=""), bty="n", cex=0.8)
  }
  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*(em_num+1)+1], rev(ssbratio_high[,(j-1)*(em_num+1)+1])), border=NA, col=alpha("gray30", 0.4))
  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*(em_num+1)+4], rev(ssbratio_high[,(j-1)*(em_num+1)+4])), border=NA, col=alpha("red", 0.4))
  lines(x_val, ssbratio_median[,(j-1)*(em_num+1)+4], type="o", lty=1, col="red", pch=19, cex=0.2)
  lines(x_val, ssbratio_median[,(j-1)*(em_num+1)+1], type="o", lty=1, col="gray30", pch=19, cex=0.2)
  abline(h=1, col="gray", lty=2)

  plot(x_val, rep(0, times=length(x_val)), ylim=range(ssbratio_low, ssbratio_high), type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  if (j == length(maindir_list)) {
    legend("topright", paste("C12: SS", sep=""), bty="n", cex=0.8)
  } else {
    legend("topright", paste("C", j-1, ": SS", sep=""), bty="n", cex=0.8)
  }
  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*(em_num+1)+1], rev(ssbratio_high[,(j-1)*(em_num+1)+1])), border=NA, col=alpha("gray30", 0.4))
  polygon(c(x_val, rev(x_val)), c(ssbratio_low[,(j-1)*(em_num+1)+5], rev(ssbratio_high[,(j-1)*(em_num+1)+5])), border=NA, col=alpha("deepskyblue3", 0.4))
  lines(x_val, ssbratio_median[,(j-1)*(em_num+1)+5], type="o", lty=1, col="deepskyblue3", pch=19, cex=0.2)
  lines(x_val, ssbratio_median[,(j-1)*(em_num+1)+1], type="o", lty=1, col="gray30", pch=19, cex=0.2)
  abline(h=1, col="gray", lty=2)
}
mtext(text="Year",side=1,line=2.5,outer=TRUE)
mtext(text="SSB/SSBMSY",side=2,line=3.5,outer=TRUE)

legend(x=30.5,y=35,
       legend=c("OM", "AMAK", "ASAP", "BAM", "SS"),
       # lty=c(1, 1, 1, 1),
       # pch=c(19, 19, 19, 19),
       col=c("gray30", "orange", "green", "red", "deepskyblue3"),
       fill = c("gray30", "orange", "green", "red", "deepskyblue3"),
       cex=0.9,
       bty="n",
       xpd=NA)
dev.off()
#### Fratio ####
jpeg(file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/Fratio.jpg", width=130, height=150, units="mm", res=1200)
op<-par(no.readonly=TRUE)
par(op)
par(oma=c(5,5,0,5),mar=c(0,0,0,0.2),mfrow=c(13,4), pch=16)

x_val <- 1:30

for (j in 1:length(maindir_list)){
  plot(x_val, rep(0, times=length(x_val)), ylim=range(fratio_low, fratio_high), type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  box()
  axis(2, at=c(0.5, 2, 3.5), labels=c(0.5, 2, 3.5), las=2)
  if (j == length(maindir_list)) {
    legend("topright", paste("C12: AMAK", sep=""), bty="n", cex=0.8)
  } else {
    legend("topright", paste("C", j-1, ": AMAK", sep=""), bty="n", cex=0.8)
  }
  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*(em_num+1)+1], rev(fratio_high[,(j-1)*(em_num+1)+1])), border=NA, col=alpha("gray30", 0.3))
  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*(em_num+1)+2], rev(fratio_high[,(j-1)*(em_num+1)+2])), border=NA, col=alpha("orange", 0.3))
  lines(x_val, fratio_median[,(j-1)*(em_num+1)+2], type="o", lty=1, col="orange", pch=19, cex=0.2)
  lines(x_val, fratio_median[,(j-1)*(em_num+1)+1], type="o", lty=1, col="gray30", pch=19, cex=0.2)
  abline(h=1, col="gray", lty=2)


  plot(x_val, rep(0, times=length(x_val)), ylim=range(fratio_low, fratio_high), type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  if (j == length(maindir_list)) {
    legend("topright", paste("C12: ASAP", sep=""), bty="n", cex=0.8)
  } else {
    legend("topright", paste("C", j-1, ": ASAP", sep=""), bty="n", cex=0.8)
  }
  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*(em_num+1)+1], rev(fratio_high[,(j-1)*(em_num+1)+1])), border=NA, col=alpha("gray30", 0.4))
  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*(em_num+1)+3], rev(fratio_high[,(j-1)*(em_num+1)+3])), border=NA, col=alpha("green", 0.4))
  lines(x_val, fratio_median[,(j-1)*(em_num+1)+3], type="o", lty=1, col="green", pch=19, cex=0.2)
  lines(x_val, fratio_median[,(j-1)*(em_num+1)+1], type="o", lty=1, col="gray30", pch=19, cex=0.2)
  abline(h=1, col="gray", lty=2)

  plot(x_val, rep(0, times=length(x_val)), ylim=range(fratio_low, fratio_high), type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  if (j == length(maindir_list)) {
    legend("topright", paste("C12: BAM", sep=""), bty="n", cex=0.8)
  } else {
    legend("topright", paste("C", j-1, ": BAM", sep=""), bty="n", cex=0.8)
  }
  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*(em_num+1)+1], rev(fratio_high[,(j-1)*(em_num+1)+1])), border=NA, col=alpha("gray30", 0.4))
  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*(em_num+1)+4], rev(fratio_high[,(j-1)*(em_num+1)+4])), border=NA, col=alpha("red", 0.4))
  lines(x_val, fratio_median[,(j-1)*(em_num+1)+4], type="o", lty=1, col="red", pch=19, cex=0.2)
  lines(x_val, fratio_median[,(j-1)*(em_num+1)+1], type="o", lty=1, col="gray30", pch=19, cex=0.2)
  abline(h=1, col="gray", lty=2)

  plot(x_val, rep(0, times=length(x_val)), ylim=range(fratio_low, fratio_high), type="l", lty=1, col="gray30", xlab="", ylab="", axes=F)
  box()
  if(j==length(maindir_list)) axis(1, at=c(5, 15, 25), labels = c(5, 15, 25))
  if (j == length(maindir_list)) {
    legend("topright", paste("C12: SS", sep=""), bty="n", cex=0.8)
  } else {
    legend("topright", paste("C", j-1, ": SS", sep=""), bty="n", cex=0.8)
  }
  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*(em_num+1)+1], rev(fratio_high[,(j-1)*(em_num+1)+1])), border=NA, col=alpha("gray30", 0.4))
  polygon(c(x_val, rev(x_val)), c(fratio_low[,(j-1)*(em_num+1)+5], rev(fratio_high[,(j-1)*(em_num+1)+5])), border=NA, col=alpha("deepskyblue3", 0.4))
  lines(x_val, fratio_median[,(j-1)*(em_num+1)+5], type="o", lty=1, col="deepskyblue3", pch=19, cex=0.2)
  lines(x_val, fratio_median[,(j-1)*(em_num+1)+1], type="o", lty=1, col="gray30", pch=19, cex=0.2)
  abline(h=1, col="gray", lty=2)
}
mtext(text="Year",side=1,line=2.5,outer=TRUE)
mtext(text="F/FMSY",side=2,line=3.5,outer=TRUE)

legend(x=30.5,y=30,
       legend=c("OM", "AMAK", "ASAP", "BAM", "SS"),
       # lty=c(1, 1, 1, 1),
       # pch=c(19, 19, 19, 19),
       col=c("gray30", "orange", "green", "red", "deepskyblue3"),
       fill = c("gray30", "orange", "green", "red", "deepskyblue3"),
       cex=0.9,
       bty="n",
       xpd=NA)
dev.off()

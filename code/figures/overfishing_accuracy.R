library(RColorBrewer)
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


yr=1:30
em_num <- 4

overfishing_accuracy <- overfished_accuracy <- matrix (NA, nrow=em_num*length(maindir_list), ncol=length(yr))


for (j in 1:length(maindir_list)){
  overfishing_accuracy[((j-1)*4+1):((j-1)*4+4), ] <- t(as.matrix(read.csv(file.path(maindir_list[j], "output", "overfishing_performance.csv"))[, 2:(1+em_num)]))
  overfished_accuracy[((j-1)*4+1):((j-1)*4+4), ] <- t(as.matrix(read.csv(file.path(maindir_list[j], "output", "overfished_performance.csv"))[, 2:(1+em_num)]))
}


#### Line plot ####
jpeg(file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/overfishing_accuracy_line.jpeg", width=150, height=120, units="mm", res=1200)

op<-par(no.readonly=TRUE)
par(op)
par(oma=c(4,4,0,5),mar=c(0.5,0,0.1,0.5),mfrow=c(5,3),pch=16)

cases = c("case0", "case1", "case2", "case3", "case4", "case5", "case6", "case7", "case8", "case9", "case10_BAM", "case11_SS", "case0_F0_0.01_1yr")

for (i in 1:length(cases)){
  overfishing = read.csv(paste("C:/Users/bai.li/Desktop/mcp_results_r/final_cases/", cases[i], "/output/overfishing_performance.csv", sep=""))
  plot(overfishing$X, overfishing[,2], ylim=c(0, 100), type="o", xlab="", ylab="", axes=F, lty=2, pch=2, cex=0.5)
  box()
  if (i %in% c(length(cases):(length(cases)-2))) axis(1)
  if(i %in% seq(1, length(cases), by=3)) axis(2, las=2)
  lines(overfishing$X, overfishing[,3], lty=3, pch=3, type="o", cex=0.5)
  lines(overfishing$X, overfishing[,4], lty=4, pch=4, type="o", cex=0.5)
  lines(overfishing$X, overfishing[,5], lty=5, pch=5, type="o", cex=0.5)

  if (i==length(cases)) {
    legend("bottomleft", paste("C12"), bty="n")
  } else {
    legend("bottomleft", paste("C", i-1, sep=""), bty="n")
  }
}
mtext(text="Year",side=1,line=2,outer=TRUE)
mtext(text="F/Flimit Accuracy (%)",side=2,line=2,outer=TRUE)

plot.new()
plot.new()
legend("bottom",
       c("AMAK", "ASAP", "BAM", "SS"),
       lty=c(2,3,4,5),
       pch=c(2,3,4,5),
       bty="n",
       cex=1,
       xpd=NA)
dev.off()



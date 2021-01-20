
# Yield over F --------------------------------------------------------------------------------

maindir  <- c("C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case0")
load(file.path(maindir, "output", "OM", paste("OM", 1, ".RData", sep="")))
yield0 <- om_output$msy$L_eq
f0 <- om_output$msy$f_seq
msy0 <- om_output$msy$msy
fmsy0 <- om_output$msy$Fmsy

maindir  <- c("C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case10")
load(file.path(maindir, "output", "OM", paste("OM", 1, ".RData", sep="")))
yield10 <- om_output$msy$L_eq
f10 <- om_output$msy$f_seq
msy10 <- om_output$msy$msy
fmsy10 <- om_output$msy$Fmsy

maindir  <- c("C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case11")
load(file.path(maindir, "output", "OM", paste("OM", 1, ".RData", sep="")))
yield11 <- om_output$msy$L_eq
f11 <- om_output$msy$f_seq
msy11 <- om_output$msy$msy
fmsy11 <- om_output$msy$Fmsy

jpeg(file=file.path("C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures", paste("yield_over_f_bc.jpg", sep="")), width=120, height=110, units="mm", res=1200)
plot(f0, yield0,
     xlim=c(0,1), ylim=c(0, 1500),
     xlab="F", ylab="Yield",
     pch=19, cex=0.5, axes=F)
axis(1)
axis(2, las=2)
box()
lines(f10, yield10, col="red")
lines(f11, yield11, col="deepskyblue3")
abline(h=msy0, col="gray", lty=2)
abline(h=msy10, col="gray", lty=2)
abline(h=msy11, col="gray", lty=2)
abline(v=fmsy0, col="gray", lty=2)
abline(v=fmsy10, col="gray", lty=2)
abline(v=fmsy11, col="gray", lty=2)
legend("right",
       c("C0: No Bias Adjustment",
         "C10: BAM Method",
         "C11: SS Method"),
       pch=c(19, NA, NA),
       lty=c(NA, 1, 1),
       col=c("black", "red", "deepskyblue3"),
       bty="n",
       cex=0.5)
dev.off()


# RE in MSY -----------------------------------------------------------------------------------

maindir_list  <- c("C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case10/output/original",
                   "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case10/output/adhoc")
#### MSY based estimates plot ####
em_num <- 4

for(j in 1:2){
  load(file.path(maindir_list[j], "performance_measures_output.RData"))
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

jpeg(file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/MSY_RE_bc.jpg", width=150, height=100, units="mm", res=1200)
op<-par(no.readonly=TRUE)
par(op)
par(oma=c(2,2,0,4),mar=c(4,0,0,0.5),mfrow=c(1,3),pch=16)
# color=c("orange", "green", "red", "deepskyblue3")
color <- rep("black", times=em_num)
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
     rev(c("A)", "B)")))

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
     rev(c("A)", "B)")))

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
     rev(c("A)", "B)")))

legend(x=max(xlim)*1.05,y=ncol(msy_stat)/2*1.2,
       legend=c("AMAK", "ASAP", "BAM", "SS"),
       pch=c(2, 3, 4, 5),
       lty=c(1, 1, 1, 1),
       # col=c("orange", "green", "red", "deepskyblue3"),
       cex=0.9,
       bty="n",
       xpd=NA)
dev.off()


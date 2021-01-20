
maindir  <- c("C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case0")
load(file.path(maindir, "output", "OM", paste("OM", 1, ".RData", sep="")))

ratio = 0.8
msy_ratio = ratio*max(om_output$msy$L_eq)
msy_points = sort(abs(om_output$msy$L_eq-msy_ratio))[1:5]
f_points = om_output$msy$f_seq[which(abs(om_output$msy$L_eq-msy_ratio) %in% msy_points)]
flow = min(f_points)
fhigh = max(f_points)

jpeg(file=file.path("C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures", paste("yield_over_f.jpg", sep="")), width=105, height=75, units="mm", res=1200)
par(mar=c(4,5,1,1))
plot(om_output$msy$f_seq, om_output$msy$L_eq,
     xlab="", ylab="",
     type="l",
     xlim=c(0, 1), ylim=c(0, max(om_output$msy$L_eq)*1.1),
     xaxs="i", yaxs="i",
     axes=F)
box()
abline(h=max(om_output$msy$L_eq))
abline(h=msy_ratio, lty=2, col="gray30")
abline(v=om_output$msy$Fmsy)
lines(x=c(flow, flow), y=c(0, msy_ratio), lty=2, col="gray30")
lines(x=c(fhigh, fhigh), y=c(0, msy_ratio), lty=2, col="gray30")
axis(1, at=c(flow, om_output$msy$Fmsy, fhigh), labels = c("Flow", "FMSY", "Fhigh"), las=2)
axis(2, at=c(om_output$msy$msy, msy_ratio), labels = c("MSY", "0.8MSY"), las=1)
dev.off()

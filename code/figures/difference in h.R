convertSRparms <- function(R0, h, phi, sigmaR, mean2med){
  BC <- ifelse(mean2med == TRUE, exp(-0.5 * sigmaR^2), exp(0.5 * sigmaR^2))
  S0BC <- (BC * 0.8 * R0 * h * phi - 0.2 * phi * R0 * (1 - h)) / (h - 0.2)
  R0BC <-  S0BC / phi
  Rnew <- BC * 0.8 * R0 * h * 0.2 * S0BC / (0.2 * phi * R0 * (1 - h) + (h - 0.2) * 0.2 * S0BC)
  hBC <- Rnew / R0BC
  return(list(S0BC = S0BC, R0BC = R0BC, hBC = hBC))
}

R0=1000000
phi=0.01025625

sigmaR=seq(0.01, 2, by=0.01)
h=seq(0.21, 1, by=0.01)

R0diff <- hdiff <- matrix(NA, nrow=length(sigmaR), ncol=length(h))

for(i in 1:length(sigmaR)){
  for(j in 1:length(h)){
    R0diff[i,j] <- (convertSRparms(R0=R0, h=h[j], phi=phi, sigmaR=sigmaR[i], mean2med=FALSE)$R0BC-R0)/R0*100
    hdiff[i,j] <- (convertSRparms(R0=R0, h=h[j], phi=phi, sigmaR=sigmaR[i], mean2med=FALSE)$hBC-h[j])
  }
}

#### With Color ####
jpeg(file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/diff_in_R0.jpg", width=170, height=100, units="mm", res=1200)
filled.contour(sigmaR, h, R0diff,
               color = terrain.colors,
               xlab="SigmaR", ylab="Median h",
               plot.axes = {
                 axis(1); axis(2);
                 contour(sigmaR, h, R0diff,
                         col = "black",
                         add = TRUE,
                         levels = c(10, 50, 100, 500, 1000, 5000),
                         method="edge")
                 },
               key.title = {
                 par(cex.main=0.8);
                 title(main="Relative Diff\nin R0 (%)")
                 })
legend("topleft", "A)", bty="n", box.col = "")
dev.off()

jpeg(file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/diff_in_h.jpg", width=150, height=100, units="mm", res=1200)
filled.contour(sigmaR, h, hdiff,
               #color = terrain.colors,
               color = white,
               xlab="SigmaR", ylab="Median h",
               plot.axes = {
                 axis(1); axis(2);
                 contour(sigmaR, h, hdiff,
                         col = "black",
                         add = TRUE,
                         levels = c(0.02, 0.1, 0.2, 0.3, 0.4),
                         method="edge")
                 },
               key.title = {
                 par(cex.main=0.8);
                 title(main="Diff in h")
               })
legend("topleft", "B)", bty="n", box.col = "")
dev.off()

#### Without Color ####
jpeg(file="C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/diff_in_R0andh.jpg", width=200, height=100, units="mm", res=1200)
par(mfrow=c(1,2))
plot(NA,xlim=range(sigmaR),
     ylim=range(h),xlab="SigmaR",ylab="Median-unbiased h",
     frame=FALSE,axes=F,xaxs="i",yaxs="i", main="Relative Diff in R0 (%)")
contour(sigmaR, h, R0diff,
        col = "black",
        add = TRUE,
        levels = c(10, 50, 100, 500, 1000, 5000),
        method="edge")
axis(1); axis(2, las=2);
legend("topleft", "A)", bty="n")
box()

plot(NA,xlim=range(sigmaR),
     ylim=range(h),xlab="SigmaR",ylab="Median-unbiased h",
     frame=FALSE,axes=F,xaxs="i",yaxs="i", main="Diff in h")
contour(sigmaR, h, hdiff,
        col = "black",
        add = TRUE,
        levels = c(0.02, 0.1, 0.2, 0.3, 0.4),
        method="edge")
axis(1); axis(2, las=2);
legend("topleft", "B)", bty="n")
box()

dev.off()

#install.packages("jpeg")
library(jpeg) # for reading in PNGs

# example image
par(mar=rep(0,4)) # no margins

# layout the plots into a matrix w/ 12 columns, by row
layout(matrix(1:2, ncol=1, byrow=TRUE))

# do the plotting
img <- readJPEG("C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/diff_in_R0.jpg", native = FALSE)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(img,0,0,1.05,1.1)

img <- readJPEG("C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/diff_in_h.jpg", native = FALSE)
plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n")
rasterImage(img,0,0,1.05,1.1)

# write to PDF
dev.print(pdf, "C:/Users/bai.li/Desktop/mcp_results_r/cases/manuscript_figures/diff_in_R0_h.pdf")

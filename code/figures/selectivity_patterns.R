#### Plot base case ####
maindir = "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case0/output"

setwd(maindir)
load(paste("./OM/OM", 1, ".RData", sep=""))
load("./keep_sim_id.RData")
if(F){
  amak_selectivity_fish <- matrix(NA, nrow=om_sim_num, ncol=om_input$nages)
  amak_selectivity_survey <- matrix(NA, nrow=om_sim_num, ncol=om_input$nages)
  asap_selectivity_fish <- matrix(NA, nrow=om_sim_num, ncol=om_input$nages)
  asap_selectivity_survey <- matrix(NA, nrow=om_sim_num, ncol=om_input$nages)
  bam_selectivity_fish <- matrix(NA, nrow=om_sim_num, ncol=om_input$nages)
  bam_selectivity_survey <- matrix(NA, nrow=om_sim_num, ncol=om_input$nages)
  ss_selectivity_fish <- matrix(NA, nrow=om_sim_num, ncol=om_input$nages)
  ss_selectivity_survey <- matrix(NA, nrow=om_sim_num, ncol=om_input$nages)

  sapply(1:om_sim_num, function(om_sim) {
    setwd(file.path(maindir, "AMAK", paste("s", keep_sim_id[om_sim], sep="")))
    amak_output <- readRep("For_R", suffix = ".rep")
    amak_selectivity_fish[om_sim, ] <<- amak_output$sel_fsh_1[1, 3:ncol(amak_output$sel_fsh_1)]
    amak_selectivity_survey[om_sim, ] <<- amak_output$sel_ind_1[1, 3:ncol(amak_output$sel_fsh_1)]

    asap_output <- dget(file.path(maindir, "ASAP", paste("s", keep_sim_id[om_sim], sep=""), "asap3.rdat"))
    asap_selectivity_fish[om_sim, ] <<- asap_output$fleet.sel.mats$sel.m.fleet1[1,]
    asap_selectivity_survey[om_sim, ] <<- asap_output$index.sel

    bam_output <- dget(file.path(maindir, "BAM", paste("s", keep_sim_id[om_sim], sep=""), "bam-sim.rdat"))
    bam_selectivity_fish[om_sim, ] <<- bam_output$sel.age$sel.m.fleet1[1, ]
    bam_selectivity_survey[om_sim, ] <<- bam_output$sel.age$sel.m.survey1[1, ]

    ss_output <- SS_output(dir=file.path(maindir, "SS", paste("s", keep_sim_id[om_sim], sep="")), ncols = 300, verbose=F, printstats=F)
    ss_selectivity_fish[om_sim, ] <<- as.matrix(ss_output$ageselex[which(ss_output$ageselex$Label=="1_1Asel"), 9:ncol(ss_output$ageselex)])
    ss_selectivity_survey[om_sim, ] <<- as.matrix(ss_output$ageselex[which(ss_output$ageselex$Label=="1_2Asel"), 9:ncol(ss_output$ageselex)])
  })

  save(amak_selectivity_fish,
       amak_selectivity_survey,
       asap_selectivity_fish,
       asap_selectivity_survey,
       bam_selectivity_fish,
       bam_selectivity_survey,
       ss_selectivity_fish,
       ss_selectivity_survey,
       file="C:\\Users\\bai.li\\Desktop\\mcp_results_r\\cases\\manuscript_figures\\case0_selectivity.RData")
} #Use "TRUE" for the first time run of the code, otherwise use FALSE

load("C:\\Users\\bai.li\\Desktop\\mcp_results_r\\cases\\manuscript_figures\\case0_selectivity.RData")

selectivity_quantile_polygon <- function(ages, selectivity_data, quantile_value, color, lty, pch){
  low_bound <- c()
  sapply(1:length(ages), function(x) low_bound[x] <<- quantile(selectivity_data[,x], quantile_value))
  high_bound <- c()
  sapply(1:length(ages), function(x) high_bound[x] <<- quantile(selectivity_data[,x], 1-quantile_value))
  lines(ages, apply(selectivity_data, 2, median), col=alpha(color, 0.8), lty=lty, pch=pch, type="o", cex=0.6)
  #lines(ages, low_bound, col=color, lty=2)
  #lines(ages, high_bound, col=color, lty=2)
}

jpeg(file="C:\\Users\\bai.li\\Desktop\\mcp_results_r\\cases\\manuscript_figures\\selectivity_patterns.jpg", width=150, height=100, units="mm", res=1200)
op<-par(no.readonly=TRUE)
par(op)
par(oma=c(2,2,0,4),mar=c(3,3,0,0),mfrow=c(2,2),pch=16)

plot(om_input$ages, om_input$selex_fleet$fleet1, type="o", col=alpha("black", 0.5),
     xlab="", ylab="", pch=19, panel.first=grid(lty=1), cex=0.6, ylim=c(0,1), axes=F)
axis(1)
axis(2, las=2)
box()
selectivity_quantile_polygon(ages = om_input$ages,
                             selectivity_data = amak_selectivity_fish,
                             quantile_value = 0.1,
                             color = "orange",
                             lty=2,
                             pch=2)
selectivity_quantile_polygon(ages = om_input$ages,
                             selectivity_data = asap_selectivity_fish,
                             quantile_value = 0.1,
                             color = "green",
                             lty=3,
                             pch=3)
selectivity_quantile_polygon(ages = om_input$ages,
                             selectivity_data = bam_selectivity_fish,
                             quantile_value = 0.1,
                             color = "red",
                             lty=4,
                             pch=4)
selectivity_quantile_polygon(ages = om_input$ages,
                             selectivity_data = ss_selectivity_fish,
                             quantile_value = 0.1,
                             color = "deepskyblue3",
                             lty=5,
                             pch=5)
legend("bottomright",
       legend="C0: Landing Selectivity",
       cex=0.8,
       bty="n")
legend("topleft", "A)", cex=0.8, bty="n")


plot(om_input$ages, om_input$selex_survey$survey1, type="o", col=alpha("black", 0.5),
     xlab="", ylab="", pch=19, panel.first=grid(lty=1), cex=0.6, ylim=c(0,1), axes=F)
axis(1)
axis(2, las=2)
box()
selectivity_quantile_polygon(ages = om_input$ages,
                             selectivity_data = amak_selectivity_survey,
                             quantile_value = 0.1,
                             color = "orange",
                             lty=2,
                             pch=2)
selectivity_quantile_polygon(ages = om_input$ages,
                             selectivity_data = asap_selectivity_survey,
                             quantile_value = 0.1,
                             color = "green",
                             lty=3,
                             pch=3)
selectivity_quantile_polygon(ages = om_input$ages,
                             selectivity_data = bam_selectivity_survey,
                             quantile_value = 0.1,
                             color = "red",
                             lty=4,
                             pch=4)
selectivity_quantile_polygon(ages = om_input$ages,
                             selectivity_data = ss_selectivity_survey,
                             quantile_value = 0.1,
                             color = "deepskyblue3",
                             lty=5,
                             pch=5)
legend("bottomright",
       legend="C0: Survey Selectivity",
       cex=0.8,
       bty="n")
legend("topleft", "B)", cex=0.8, bty="n")

#### Plot case 5_2 ####
maindir = "C:/Users/bai.li/Desktop/mcp_results_r/final_cases/case8/output"
setwd(maindir)

load(paste("./OM/OM", 1, ".RData", sep=""))
load("./keep_sim_id.RData")
if(F){
  amak_selectivity_fish <- matrix(NA, nrow=om_sim_num, ncol=om_input$nages)
  amak_selectivity_survey <- matrix(NA, nrow=om_sim_num, ncol=om_input$nages)
  asap_selectivity_fish <- matrix(NA, nrow=om_sim_num, ncol=om_input$nages)
  asap_selectivity_survey <- matrix(NA, nrow=om_sim_num, ncol=om_input$nages)
  bam_selectivity_fish <- matrix(NA, nrow=om_sim_num, ncol=om_input$nages)
  bam_selectivity_survey <- matrix(NA, nrow=om_sim_num, ncol=om_input$nages)
  ss_selectivity_fish <- matrix(NA, nrow=om_sim_num, ncol=om_input$nages)
  ss_selectivity_survey <- matrix(NA, nrow=om_sim_num, ncol=om_input$nages)


  sapply(1:om_sim_num, function(om_sim) {
    setwd(file.path(maindir, "AMAK", paste("s", keep_sim_id[om_sim], sep="")))
    amak_output <- readRep("For_R", suffix = ".rep")
    amak_selectivity_fish[om_sim, ] <<- amak_output$sel_fsh_1[1, 3:ncol(amak_output$sel_fsh_1)]
    amak_selectivity_survey[om_sim, ] <<- amak_output$sel_ind_1[1, 3:ncol(amak_output$sel_fsh_1)]

    asap_output <- dget(file.path(maindir, "ASAP", paste("s", keep_sim_id[om_sim], sep=""), "asap3.rdat"))
    asap_selectivity_fish[om_sim, ] <<- asap_output$fleet.sel.mats$sel.m.fleet1[1,]
    asap_selectivity_survey[om_sim, ] <<- asap_output$index.sel

    bam_output <- dget(file.path(maindir, "BAM", paste("s", keep_sim_id[om_sim], sep=""), "bam-sim.rdat"))
    bam_selectivity_fish[om_sim, ] <<- bam_output$sel.age$sel.m.fleet1[1, ]
    bam_selectivity_survey[om_sim, ] <<- bam_output$sel.age$sel.m.survey1[1, ]

    ss_output <- SS_output(dir=file.path(maindir, "SS", paste("s", keep_sim_id[om_sim], sep="")), ncols = 300, verbose=F, printstats=F)
    ss_selectivity_fish[om_sim, ] <<- as.matrix(ss_output$ageselex[which(ss_output$ageselex$Label=="1_1Asel"), 9:ncol(ss_output$ageselex)])
    ss_selectivity_survey[om_sim, ] <<- as.matrix(ss_output$ageselex[which(ss_output$ageselex$Label=="1_2Asel"), 9:ncol(ss_output$ageselex)])
  })

  save(amak_selectivity_fish,
       amak_selectivity_survey,
       asap_selectivity_fish,
       asap_selectivity_survey,
       bam_selectivity_fish,
       bam_selectivity_survey,
       ss_selectivity_fish,
       ss_selectivity_survey,
       file="C:\\Users\\bai.li\\Desktop\\mcp_results_r\\cases\\manuscript_figures\\case8_selectivity.RData")
} #Use "TRUE" for the first time run of the code, otherwise use FALSE

load("C:\\Users\\bai.li\\Desktop\\mcp_results_r\\cases\\manuscript_figures\\case8_selectivity.RData")
plot(om_input$ages, om_input$selex_fleet$fleet1, type="o", col=alpha("black", 0.5),
     xlab="", ylab="", pch=19, panel.first=grid(lty=1), cex=0.6, ylim=c(0,1), axes=F)
axis(1)
axis(2, las=2)
box()
selectivity_quantile_polygon(ages = om_input$ages,
                             selectivity_data = amak_selectivity_fish,
                             quantile_value = 0.1,
                             color = "orange",
                             lty=2,
                             pch=2)
selectivity_quantile_polygon(ages = om_input$ages,
                             selectivity_data = asap_selectivity_fish,
                             quantile_value = 0.1,
                             color = "green",
                             lty=3,
                             pch=3)
selectivity_quantile_polygon(ages = om_input$ages,
                             selectivity_data = bam_selectivity_fish,
                             quantile_value = 0.1,
                             color = "red",
                             lty=4,
                             pch=4)
selectivity_quantile_polygon(ages = om_input$ages,
                             selectivity_data = ss_selectivity_fish,
                             quantile_value = 0.1,
                             color = "deepskyblue3",
                             lty=5,
                             pch=5)
legend("bottomright",
       legend="C8: Landing Selectivity",
       cex=0.8,
       bty="n")
legend("topleft", "C)", cex=0.8, bty="n")

plot(om_input$ages, om_input$selex_survey$survey1, type="o", col=alpha("black", 0.5),
     xlab="", ylab="", pch=19, panel.first=grid(lty=1), cex=0.6, ylim=c(0,1), axes=F)
axis(1)
axis(2, las=2)
box()
selectivity_quantile_polygon(ages = om_input$ages,
                             selectivity_data = amak_selectivity_survey,
                             quantile_value = 0.1,
                             color = "orange",
                             lty=2,
                             pch=2)
selectivity_quantile_polygon(ages = om_input$ages,
                             selectivity_data = asap_selectivity_survey,
                             quantile_value = 0.1,
                             color = "green",
                             lty=3,
                             pch=3)
selectivity_quantile_polygon(ages = om_input$ages,
                             selectivity_data = bam_selectivity_survey,
                             quantile_value = 0.1,
                             color = "red",
                             lty=4,
                             pch=4)
selectivity_quantile_polygon(ages = om_input$ages,
                             selectivity_data = ss_selectivity_survey,
                             quantile_value = 0.1,
                             color = "deepskyblue3",
                             lty=5,
                             pch=5)
legend("bottomright",
       legend="C8: Survey Selectivity",
       cex=0.8,
       bty="n")
legend("topleft", "D)", cex=0.8, bty="n")

mtext(text="Age",side=1,line=0,outer=TRUE)
mtext(text="Selectivity",side=2,line=0,outer=TRUE)

legend(x=12.5,y=1.5,
       legend=c("OM", "AMAK", "ASAP", "BAM", "SS"),
       pch=c(19, 2, 3, 4, 5),
       lty=c(1, 2, 3, 4, 5),
       col=c("black", "orange", "green", "red", "deepskyblue3"),
       cex=0.8,
       bty="n",
       xpd=NA)
dev.off()

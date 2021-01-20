OM_Run <- function(maindir=maindir, subdir="OM", om_sim_num=160, r_dev_matrix=r_dev_matrix, f_dev_matrix=f_dev_matrix, f_matrix=f_matrix){
  
  setwd(paste(maindir))
  
  for(om_sim in 1:om_sim_num){
    logR.resid <<- r_dev_matrix[om_sim,]
    logf.resid <<- f_dev_matrix[om_sim,]
    f <<- f_matrix[om_sim,]
    
    source(file= file.path(maindir, "OM_Settings.R")) 
    
    em_input=ObsModel(L=om_output$L.mt, 
                      survey=om_output$survey_index,
                      L.age=om_output$L.age, 
                      survey.age=om_output$survey_age_comp,
                      cv.L=cv.L, 
                      cv.survey=cv.survey, 
                      n.L=n.L, 
                      n.survey=n.survey)
  
    em_input$cv.L <- input.cv.L
    em_input$cv.survey <- input.cv.survey
    
    save(om_input, om_output, em_input, file=file.path(maindir, "output", subdir, paste("OM", om_sim, ".RData", sep="")))
  }
}







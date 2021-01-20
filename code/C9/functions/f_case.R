f_case <- function(f_pattern=1, f_min=0.01, f_max=0.39, f_end=NULL, f_mean=NULL, start_val=NULL, start_year=NULL, nyr=30, om_sim_num=160, f_dev_matrix=f_dev_matrix){
  
  ## f_pattern: 1. increase
  ##            2. decrease
  ##            3. increase first, then decrease
  ##            4. decrease first, then increase
  ##            5. fluctuate around a value 
  ##            6. ramp + fluctuate around a value
  ## f_min: min value for the main shape of fishing mortality
  ## f_max: max value for the main shape of fishing mortality
  ## f_end: ending F for patterns 3 and 4
  ## f_mean: F value for patterns 5 and 6
  ## start_val: starting F for pattern 6
  ## start_year: starting year for f_mean
  ## nyr: number of years
  ## om_sim_num: number of iterations per case
  
  f_matrix <- matrix(NA, nrow=om_sim_num, ncol=nyr)
  
  for(om_sim in 1:om_sim_num){
    if(f_pattern == 1) f_matrix[om_sim,]=seq(f_min, f_max, length=nyr)*exp(f_dev_matrix[om_sim,])
    if(f_pattern == 2) f_matrix[om_sim,]=seq(f_max, f_min, length=nyr)*exp(f_dev_matrix[om_sim,])
    if(f_pattern == 3) f_matrix[om_sim,]=c(seq(f_min, f_max, length=round(nyr/5*4)), seq(f_max, f_end, length=nyr-round(nyr/5*4)))*exp(f_dev_matrix[om_sim,])
    if(f_pattern == 4) f_matrix[om_sim,]=c(seq(f_max, f_min, length=round(nyr/5*4)), seq(f_min, f_end, length=nyr-round(nyr/5*4)))*exp(f_dev_matrix[om_sim,])
    if(f_pattern == 5) f_matrix[om_sim,]=seq(f_mean, f_mean, length=nyr)*exp(f_dev_matrix[om_sim,])
    if(f_pattern == 6) f_matrix[om_sim,]=c(rep(start_val, length=start_year-1), rep(f_mean, length=nyr-start_year+1))*exp(f_dev_matrix[om_sim,])
  }
  
  return(f_matrix)
}

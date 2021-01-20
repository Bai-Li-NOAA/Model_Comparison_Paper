r_dev_case <- function(r_dev_change=TRUE, nyr=30, logR_sd=0.1, om_sim_num=160){
  
  ## rec_dev_change:1. TRUE: generate different recruitment deviations per iteration
  ##                2. FALSE: generate same recruitment deviations per iteration
  ## nyr: number of years
  ## logR.sd: recruitment SD in log space
  ## om_sim_num: number of iterations per case
  
  r_dev_matrix <- matrix(NA, nrow=om_sim_num, ncol=nyr)
  
  for(om_sim in 1:om_sim_num){
    ## Simulate recruitment deviations
    if(r_dev_change==TRUE){
      r_dev_matrix[om_sim,] <- rnorm(nyr, mean=0, sd=logR_sd) 
    } else {
      if(om_sim==1){
        r_dev_matrix[om_sim,] <- rnorm(nyr, mean=0, sd=logR_sd)
      } else {
        r_dev_matrix[om_sim,] <- r_dev_matrix[(om_sim-1),]
      }
    }
  }
  
  return(r_dev_matrix)
}

f_dev_case <- function(f_dev_change=TRUE, nyr=30, logf_sd=0.1, om_sim_num=160){
  
  ## f_dev_change:1. TRUE: generate different annual F deviations per iteration
  ##              2. FALSE: generate same annual F deviations per iteration
  ## nyr: number of years
  ## logf.sd: SD of fishing mortality in log space
  ## om_sim_num: number of iterations per case
  
  f_dev_matrix <- matrix(NA, nrow=om_sim_num, ncol=nyr)
  
  for(om_sim in 1:om_sim_num){
    ## Simulate F deviations
    if(f_dev_change==TRUE){
      f_dev_matrix[om_sim,] <- rnorm(nyr, mean=0, sd=logf_sd) 
    } else {
      if(om_sim==1){
        f_dev_matrix[om_sim,] <- rnorm(nyr, mean=0, sd=logf_sd)
      } else {
        f_dev_matrix[om_sim,] <- f_dev_matrix[(om_sim-1),]
      }
    }
  }
  
  return(f_dev_matrix)
}


convertSRparms <- function(R0, h, phi, sigmaR, mean2med){
  BC <- ifelse(mean2med == TRUE, exp(-0.5 * sigmaR^2), exp(0.5 * sigmaR^2))
  S0BC <- (BC * 0.8 * R0 * h * phi - 0.2 * phi * R0 * (1 - h)) / (h - 0.2)
  R0BC <-  S0BC / phi
  Rnew <- BC * 0.8 * R0 * h * 0.2 * S0BC / (0.2 * phi * R0 * (1 - h) + (h - 0.2) * 0.2 * S0BC)
  hBC <- Rnew / R0BC
  return(list(S0BC = S0BC, R0BC = R0BC, hBC = hBC))
}
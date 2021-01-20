msy_calcs<-function(steep, R0, M, wgt, prop.f=0.5, selL, selD, selZ, mat.f=NULL, mat.m=NULL, sigma=0, maxF=1.0, step=0.01, verbose=FALSE){

	## Language:          R
	## Contains:          Function "msy.calcs"
	## Programmer:        Kyle Shertzer
	## First coded:       March, 2006
	## Last revised:      March, 2006
	## Purpose:           Compute MSY benchmarks
	##
  ## Note: Benchmarks based on Bev-Holt S-R fcn with steepness parameterization.
  ##        SSB may be based on females only, males only, or both.
  ##        The measure is controlled by which maturity
  ##        vectors are supplied (mat.f, mat.m, or both).
	##
	## INPUT PARAMETERS:
	##    steep = steepness parameter
	##    R0    = virgin recruitment parameter
	##    M     = natural mortality, may be constant or vector
	##    wgt   = vector of weight at age
  ##    prop.f= proportion female at age, default = 0.5
  ##    selL  = selectivity at age for landings
  ##    selD  = selectivity at age for dead discards
  ##    selZ  = selectivity at age for dead fish
  ##    OPTIONAL INPUT
	##    mat.f = vector of proportion females mature at age
	##    mat.m = vector of proportion males mature at age
	##    sigma = lognormal bias correction -- exp(sigma^2/2)
  ##    maxF  = maximum F examined
	##    step  = accuracy in MSY calculations (default is 0.01)
	##
  ## OUTPUT:
  ##    MSY   = maximum sustainable yield
  ##    Dmsy  = dead discards at MSY
  ##    Fmsy  = fishing rate at MSY
  ##    SSBmsy= spawning stock biomass at msy
  ##    Rmsy  = equilibrium recruitment at msy
  ##    Bmsy  = total biomass (male and female) at msy
  ##    Emsy  = exploitation rate at msy (total catch / number of fish)
  ##    spr_msy= spawners per recruit at msy
  ##    SPRmsy= spawning potential ratio (spr_msy/spr_virgin)

#  if (! is.numeric(amin)) stop ("Non-numeric value for minimum age!")
#  if (! is.numeric(amax)) stop ("Non-numeric value for maximum age!")
  if (! is.numeric(steep)) stop ("Non-numeric value for steepnes!")
  if (! is.numeric(R0)) stop ("Non-numeric value for R0!")
  if (! is.numeric(M)) stop ("Non-numeric value for M!")
  if (! is.numeric(wgt)) stop ("Non-numeric value for weight at age!")
  if (! is.numeric(prop.f)) stop ("Non-numeric value for proportion female at age!")
  if (! is.numeric(selL)) stop ("Non-numeric value for selectivity at age!")
  if (is.null(mat.f) & is.null(mat.m)) stop ("Specify maturity schedule -- male, female, or both!")

  if (verbose){
  	if (sigma == 0) {cat("*** MSY NOTE: Estimates contain no bias correction.\n")}
  	else {cat("*** NOTE: Estimates contain bias correction.\n")}
  	if (is.null(mat.f)){cat("*** MSY NOTE: SSB based on males only.\n")}
  	if (is.null(mat.m)){cat("*** MSY NOTE: SSB based on females only.\n")}
  	if (!is.null(mat.f) & !is.null(mat.m)){cat("*** MSY NOTE: SSB based on both sexes.\n")}
  }
  ##INITIALIZATION
  BC=exp(sigma^2/2.0)              #multiplicative bias correction
	nages=length(wgt)
	if (length(M)>1){M_age=M}        #natural mortality at age (may be constant)
	else M_age=rep(M,nages)
	prop.m=1.0-prop.f                #proportion male

  mu.f=rep(0.0,nages)              #proportion females mature at age, initialized at zero
  mu.m=rep(0.0,nages)              #proportion males mature at age, initialized at zero
	if (is.null(mat.f)){mu.m=mat.m}
	if (is.null(mat.m)){mu.f=mat.f}
	if (!is.null(mat.f) & !is.null(mat.m)){mu.f=mat.f; mu.m=mat.m}
	reprod=wgt*(prop.f*mu.f + prop.m *mu.m) #constant vector multiplied by abundance to get SSB at age


	f=seq(0.0,maxF, by=step)
	spr=rep(0.0,length(f))   #equilibrium spr at F
	S_eq=rep(0.0,length(f))  #equilibrium SSB at F
	R_eq=rep(0.0,length(f))  #equilibrium recruitment at F
	B_eq=rep(0.0,length(f))  #equilibrium biomass at F
	L_eq=rep(0.0,length(f))  #equilibrium landings at F
	D_eq=rep(0.0,length(f))  #equilibrium dead discards at F
	E_eq=rep(0.0,length(f))  #equilibrium exploitation rate at F (landings only)

  L_age=rep(0.0,nages)     #landings at age
  D_age=rep(0.0,nages)     #dead discards at age
  F_age=rep(0.0,nages)     #F at age
  Z_age=rep(0.0,nages)     #Z at age

  ## Compute virgin spr
  N0=rep(1.0,times=nages)
  for (iage in 2:nages){N0[iage]=N0[iage-1]*exp(-1.0*M_age[iage-1])}
  N0[nages]=N0[nages-1]*exp(-1.*M_age[nages-1])/(1.-exp(-1.0*M_age[nages]))
  spr_F0=sum(N0*reprod)

  ## BEGIN ALGORITHM
  for (i in 1:length(f)){
    FL_age=f[i]*selL
    FD_age=f[i]*selD
    Z_age=M_age+f[i]*selZ

    N_age=rep(1.0,nages)     #N at age
    for (iage in 2:nages){
      N_age[iage]=N_age[iage-1]*exp(-1.0*Z_age[iage-1])
    }
    #last age is pooled
    N_age[nages]=N_age[nages-1]*exp(-1.*Z_age[nages-1])/
                  (1.-exp(-1.0*Z_age[nages]));


    spr[i]=sum(N_age*reprod)
    R_eq[i]=(R0/((5.0*steep-1.0)*spr[i]))*
            (BC*4.0*steep*spr[i]-spr_F0*(1.-steep))
    if (R_eq[i]<0.0000001) R_eq[i]=0.0000001

    N_age=R_eq[i]*N_age

    S_eq[i]=sum(N_age*reprod)
    B_eq[i]=sum(N_age*wgt)

    for(iage in 1:nages){
      L_age[iage]=N_age[iage]*
                  (FL_age[iage]/Z_age[iage])*(1.-exp(-1.0*Z_age[iage]))
      D_age[iage]=N_age[iage]*
        (FD_age[iage]/Z_age[iage])*(1.-exp(-1.0*Z_age[iage]))
    }
    L_eq[i]=sum(L_age*wgt)
    D_eq[i]=sum(D_age*wgt)
    E_eq[i]=sum(L_age)/sum(N_age)

  } #END F loop

  msy_out=max(L_eq)
  F_msy_out=f[L_eq==msy_out]
  spr_msy_out=spr[L_eq==msy_out]
  SR_msy_out=spr_msy_out/spr_F0
  D_msy_out=D_eq[L_eq==msy_out]
  R_msy_out=R_eq[L_eq==msy_out]
  S_msy_out=S_eq[L_eq==msy_out]
  B_msy_out=B_eq[L_eq==msy_out]
  E_msy_out=E_eq[L_eq==msy_out]

	if (F_msy_out==maxF){cat("*** Fmsy reached a bound.\n")}

  return(list(msy=msy_out, Fmsy=F_msy_out, Dmsy=D_msy_out, spr_msy=spr_msy_out, SPRmsy=SR_msy_out, SSBmsy=S_msy_out, Rmsy=R_msy_out, Bmsy=B_msy_out, Emsy=E_msy_out, f_seq=f, L_eq=L_eq, D_eq=D_eq, SSB_eq=S_eq, R_eq=R_eq, spr=spr))

}


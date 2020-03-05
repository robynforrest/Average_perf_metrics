######################################################################
# A conventional age-structured model to run out for nyr years under a vector of Fs
# Gets called in a loop below as the "Long model"
# Runs over a sequence of Fs and values of M1
PopDyn <- function( obj, ft ,nyrs=100 )
{
  # Extract needed model parameters.
  nages <- obj$nages
  sage  <- obj$sage
  M     <- obj$M
  fa    <- obj$fa
  CR    <- obj$CR
  ro    <- obj$ro
  sel   <- obj$sel
  wa    <- obj$wa

  Ft <- ft

  #Vectors at age
  lx <- vector(length=nages)

  #Vectors at time
  SBt <- vector(length=nyrs)     # spawning biomass
  Mt  <- vector(length=nyrs)     # total mortality

  #Matrices
  Nt <- matrix(   nrow=nyrs, ncol=nages )     # Numbers at age
  Zt <- matrix( M,nrow=nyrs, ncol=nages )
  Ct <- matrix(   nrow=nyrs, ncol=nages )     # Catch at age

  #unfished survivorship
  #surv <- exp(-M0)
  surv <- exp(-M)
  lx[1] <- 1
  for(j in 2:nages) lx[j] <- surv*lx[j-1]
  lx[nages] <- lx[nages]/(1-surv)

  #unfished eqm eggs per recruit and recruitment parameters
  phie <- sum(lx*fa)
  bhalpha = CR/(phie)
  bhbeta=( bhalpha*phie - 1.)/(ro*phie)

  #First year unfished
  Nt[1,] <- ro*lx
  #Mt[1] <- M0 + (M1-M0)*(1-Bt[1]/Bt[1])
  Mt[1] <- M #fixed M for now - add random walk or DDM later
  Zt[1,] <- sel*Ft[1] + Mt[1]
  SBt[1] <- sum(Nt[1,]*fa)


  for(i in 2:nyrs){
    if(i<=sage) Nt[i,1]    <-   (bhalpha* SBt[i-1])/(1+ bhbeta*SBt[i-1]) #assume eqm SBt for first sage
    if(i>sage) Nt[i,1]    <-   (bhalpha* SBt[i-sage])/(1+ bhbeta*SBt[i-sage])

    #This is not ADMB indices must start at 1, even though age of recruit may be older.
    #This is accounted for above and in the calculation of nages in the data section ... just remember for figs etc
    for(j in 2:nages)  Nt[i,j]   <-  Nt[i-1,j-1] * exp(-Zt[i-1,j-1])
    Nt[i,nages] <- Nt[i,nages] + Nt[i-1,nages] * exp(-Zt[i-1,nages]);  #plus group

    #update biomass vectors
    SBt[i] <- sum(Nt[i,]*fa)

    #update M and Z
    Mt[i] <- M #M0 + (M1-M0)*(1-Bt[i]/Bt[1])
    Zt[i,] <- sel*Ft[i]+Mt[i]
  } #end i

  #get catch
  for(i in 1:nyrs){
    for(j in 1:nages) {
      Ct[i,j] <- ((sel[j]*Ft[i])/(sel[j]*Ft[i] + Mt[i])) * (1-exp(-(sel[j]*Ft[i] + Mt[i]))) * Nt[i,j] * wa[j]
    }#end j
  }#end i

  #return equilibrium values
  out<- list()
  out$Ctot <- rowSums(Ct) #Total catch all years
  out$SBt <- SBt
  out$Re <- Nt[nyrs,1]	#Recruits (Recruits are sage years old)

    #Final SPR at F
  SSBPRX <- SBt[nyr]/Nt[nyrs,1]
  SSBPR0 <- SBt[1]/Nt[1,1]
  out$SPRX <- SSBPRX/SSBPR0

  return(out)
}     # END function PopDyn

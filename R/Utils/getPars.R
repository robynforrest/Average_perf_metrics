# Make an input list.
#getPars <- function( xRep )
getPars <- function( xPar, xRep, matEqSel=FALSE, nagenew, Mnew, vbknew, linfnew ) # xPar and xRep are lists containing par and rep file results.
{
  # xRep is a list of parameters read from iSCAM rep file.
  val <- list()

  #fraction of spawning assumed to take place prior to fishery - important!
  #this is in iscam control file (.ctl) - hardwired here by lazy me (RF)
  val$sp_frac <- 1.0

  yr  <- xRep$yr
  nyr <- length(yr)

  # Leading parameters (MLEs)
  val$h     <- xRep$steepness
  val$CR    <- 4.0 * val$h / (1.0-val$h)
  val$ro    <- xRep$ro
  val$M     <- Mnew  #This is the estimated base M - time varying (spline) M is a function of this

  #get mean M
  val$Mtot  <- xRep$M_tot[,1]  #Estimated annual M (constant across age so just take 1st col)
  val$M     <- Mnew #All years: used in iscam calcStockRecruitment

  val$sage  <- min(xRep$age)
  val$nage  <- nagenew
  val$ages  <- val$sage:nagenew
  val$nages <- length(val$ages)
  # Selectivity - these values come from *.std file sel_par[3,4].

  # NOTE: May not need these if selectivity ogive from *.rep.
  #val$selah <- exp(1.39824905505 )         # HG Gear 2 age at 50% selectivity mle
  #val$selsd <- exp(-0.469757200741)        # SD in selectivity ogive HG

  # Extract logistic selectivity parameters for gear 2 and calculate selectivity ogive.
  selPars <- matrix( NA, nrow=1, ncol=2 )
  #selPars[ 1, ] <- exp( xPar$'sel_par[1]:' )
  selPars[ 1, ] <- exp( xPar$'sel_par[2]:' )
  #selPars[ 3, ] <- exp( xPar$'sel_par[3]:' )
  #selPars[ 4, ] <- exp( xPar$'sel_par[4]:' )
  #selPars[ 5, ] <- exp( xPar$'sel_par[5]:' )
  val$selPars <- selPars

  selAge  <- matrix( NA, nrow=nrow(selPars), ncol=val$nages )
  for ( i in 1:nrow(selPars) )
    selAge[ i, ] <- plogis( val$ages,val$selPars[ i,1 ],val$selPars[ i,2 ] )

  val$selAge  <- selAge

  # Now extract gear 2 for equilibrium calculations.
  #val$sel <- val$selAge[ 2, ]
  val$sel <- val$selAge[ 1, ]
  val$linf  <- linfnew        # L infinity
  val$vbk   <- vbknew      # von Bert K
  val$lwa   <- 4.5e-06   # length-weight a parameter
  val$lwb   <- 3.127     # length-weight b parameter
  val$amat  <- 2.055     # age at 50% maturity
  val$matsd <- 0.05      # sd in maturity ogive
  val$ma    <-  plogis( val$ages,val$amat,val$matsd ) # Maturity ogive from *.dat
  val$lt    <- val$linf * (1 - exp(-val$vbk * (val$ages)))   # Length at age
  val$wa    <- val$lwa * val$lt^val$lwb                      # Weight at age

  # This replaces val$ma, same except for age-3, i.e., ma[2] to match assessment.
  #val$ma    <- xRep$mat  # Fixed input maturity at age from assessment.

  # Set maturity equal to selectivity?
  if ( matEqSel )
    val$sel <- val$ma

  val$lt    <- val$linf * (1 - exp(-val$vbk * (val$ages)))   # Length at age
  val$wa    <- val$lwa * val$lt^val$lwb                      # Weight at age

  # Get average weight at age and average M for eqm calcs (as per iscam)
  # This replaces val$wa
  #val$obswa    <- xRep$wt_obs[1:nyr,]   #Empirical weight at age 1951-2016 - don't take projection year
  #for(i in 1:val$nages) val$wa[i] <- mean(val$obswa[,i]) #mean weight at age for eqm calcs

  #Fecundity at age
  val$fa    <- val$wa * val$ma

  # Stock and recruitment "data".
  sbt <- xRep$sbt
  rt  <- xRep$rt

  # Age-2 is the first year in the model.  An age-2 spawner in 1951
  # can't generate recruits until 3 years later. Eggs in 1951 are age-1
  # in 1952 and age-2 in 1953.
  val$ObsYear <- yr[ (val$sage+1):nyr ]
  val$ObsSB   <- sbt[ (val$sage+1):nyr ]
  val$ObsRec  <- rt
  val$Nobs    <- length(rt)

  return( val )
}     # END function getPars

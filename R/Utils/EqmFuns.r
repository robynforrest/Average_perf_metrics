# AUTHOR: Robyn Forrest, Pacific Biological Station
# Implementation: October 21, 2016 (REF)
# Revisions:      October 27, 2016 (ARK) Added passed input list to functions.
#                 October 31, 2016 (ARK) Added function calcRefCurves.
#                 October 31, 2016 (ARK) Added F0.1.

# Equilibrium functions to support herringLRP.r

# Let:
# Fe    = eqm fishing mortality
# Me    = eqm Mortality
# Phie  = Unfished eqm spawning biomass per recruit
# Phif  = Fished eqm spawning biomass per recruit
# Phib  = Unfished eqm total biomass per recruit
# Phibf = Fished eqm total biomass per recruit

# With leading parameters Ro and CR (the recruitment compensation ratio):
# Equilibrium Recruitment Re = Ro*(CR-phie/phif)/(CR-1.)

# For a given value of Fe and Me:
# B0 = Phib*Ro   // unfished biomass
# Be = Phibf*Re  // eqm fished biomass

# Also includes a population dynamics model to run out for n years at
# steady state to test eqm calcs.

######################################################################
# Functions
######################################################################

######################################################################
# EQUILIBRIUM FUNCTIONS

## calc_eqm
# Function to calculate per recruit functions and partial derivatives
# needed for fmsy.  See Table 2 in Martell et al. 2008
# Arguments: fe=fishing rate and equilibrium mortality
# update Nov 4 - now accounts for spawn timing
calc_eqm <- function( obj, fe )
{
  nages <- obj$nages
  M     <- obj$M
  sel   <- obj$sel
  ma    <- obj$ma
  fa    <- obj$fa
  ro    <- obj$ro
  CR    <- obj$CR
  wa    <- obj$wa

  #iscam spawning adjustment parameter (set to 1 in herring assessment)
  #this is the fraction of total mortality that takes place prior to spawning
  #i.e. the herring assessment assumes all fishing and natural mortality occurs before spawning
  sp_frac <- obj$sp_frac

  #eqm unfished survivorship
  lx <- vector(length=nages)
  lxx <- vector(length=nages) #adjusted for spawn timing

  #eqm fished survivorship
  lz <- vector(length=nages)
  lzz <- vector(length=nages) #adjusted for spawn timing

  sa <- exp(-M) #unfished survival rate
  za <- M + sel*fe
  zf <- 1.-exp(-za)
  qa <- (sel/za)*zf  #elem_prod(elem_div(vtot(nyr),za),zf);
  dlz_df <- 0

  # unfished survivorship
  lx[1]=1.0;
  lxx[1]=1.0*exp(-M*sp_frac) #spawn timing adjustment

 for(i in 2:nages){
    lx[i]=lx[i-1]*sa;

    #spawn timing adjustment
    lxx[i] <- lx[i] * exp(-M*sp_frac);
  }

  lx[nages] = lx[nages]/(1.-sa);
  lxx[nages] = lxx[nages]/(1.-sa);
  phie = sum(lxx*fa);

  #initialize eqm fished survivorship
  lz[1]=1.0;
  lzz[1]=1.0*exp(-M*sp_frac) #spawn timing adjustment
  dphiq_df=0; dphif_df=0;

  for(i in 1:nages)
  {
    if (i>1) lz[i]  <- lz[i-1]*exp(-za[i-1]);
    if (i>1) dlz_df <- dlz_df*exp(-za[i-1]) - lz[i-1]*sel[i-1]*exp(-za[i-1]);

    #plus group
    if(i==nages)
    {
      lz[i]  <- lz[i]/(1.-exp(-za[i]));
      dlz_df <- dlz_df/(1.-exp(-za[i]))-lz[i-1]*exp(-za[i-1])*sel[i]*exp(-za[i])/((1.-exp(-za[i]))*(1.-exp(-za[i])));
    }

    #spawn timing adjustment
    lzz[i] <- lz[i] * exp(-M*sp_frac);

    # Partial derivatives - fished eqm spawning and vulnerable biomass wrt f
    dphif_df <- dphif_df+(fa[i])*dlz_df;
    dphiq_df <- dphiq_df+(wa[i]*qa[i]*dlz_df+((lz[i]*wa[i]*sel[i]*sel[i])/za[i])*(exp(-za[i])-zf[i]/za[i]));
  }

  phif <- sum(lzz*fa);
  phiq <- sum(lzz*wa*qa);
  re <- ro*(CR-phie/phif)/(CR-1.);
  be <- phif*re;
  bo <- phie*ro;
  #partial derivative - fished eqm recruitment wrt f
  dRe_df <- (ro/(CR-1.))*(phie/phif^2)*dphif_df

  #get SR parameters while we are here
  bhalpha = CR/(phie)
  bhbeta=( bhalpha*phie - 1.)/(ro*phie)

  out <- list()
  out$phie     <- phie         # Unfished eqm SSB per recruit.
  out$phif     <- phif         # Fished eqm SSB per recruit.
  out$phiq     <- phiq
  out$dre_df   <- dRe_df
  out$dphiq_df <- dphiq_df     # Fished eqm vulnerable biomass wrt f.
  out$dphif_df <- dphif_df     # Fished eqm spawning biomass wrt f.
  out$re       <- re           # Eqm recruitment.
  out$be       <- be           # Eqm fished biomass.
  out$bo       <- bo           # Eqm fished biomass.
  out$bhalpha  <- bhalpha      # Beverton-Holt S-R alpha.
  out$bhbeta   <- bhbeta       # Beverton-Holt S-R beta.
  out$ye       <- fe*re*phiq   # Yield at eqm
  out$ypr      <- fe*phiq   # Yield at eqm

  #out$yetest <- sum( lz*wa*sel*fe*(1. - exp(-M-sel*fe))/(M+sel*fe) )*re #Test: is the same as ye
  out$yetest2 <-  re * sum(lz*wa*sel)*fe*(1.-exp(-M-fe))/(M+fe) #Test: is the same as ye

  return(out)
}     # END calc_eqm

# calcRefCurves
# Purpose:     Calculate all equilibrium relationships to fishing mortality.
# Parameters:  obj=list of all operating model parameters; nFs=the number of
#              fishing mortality points over which to compute the functions
# Returns:     a list with vectors of fishing mortality (f) and equilibrium
#              functions (these are later fitted with splines for finding
#              references points via root-finding algorithms)
# Source:      S.P. Cox, modified by ARK for herring.
calcRefCurves <- function( obj, nFs=100 )
{
  maxF <- 10

  f <- seq( from=0.0, to=maxF, length=nFs )

  recruits <- rep( NA, length=nFs )
  ssbpr    <- rep( NA, length=nFs )
  ssb      <- rep( NA, length=nFs )
  exb      <- rep( NA, length=nFs )
  ypr      <- rep( NA, length=nFs )
  yield    <- rep( NA, length=nFs )

  for( i in 1:length(f) )
  {
    tmp        <- calc_eqm( obj=obj, fe=f[i] )
    recruits[i]<- tmp$re
    ssbpr[i]   <- tmp$phif
    ssb[i]     <- tmp$be
    exb[i]     <- tmp$dphiq_df
    ypr[i]     <- tmp$ypr
    yield[i]   <- tmp$ye
  }

  refCurves <- obj
    refCurves$F        <- f
    refCurves$ypr      <- ypr
    refCurves$ssbpr    <- ssbpr
    refCurves$yield    <- yield
    refCurves$ssb      <- ssb
    refCurves$exb      <- exb
    refCurves$recruits <- recruits
  refCurves
}     # END function calcRefCurves


#################################################################################################
# Code stolen from old versions of iscam (Steve Martell)
# Reference: Martell et al. 2008
# Uses Newton_Raphson method to determine Fmsy
calc_fmsy <- function( obj )
{
  #initial guess for fmsy
  fe <- 0.3

  # iteratively solve for Fmsy using Newton algorithm
  # find value of fe that sets derivative of yield function to zero, i.e., max yield
  for( i in 1:150 )
  {
    xx <- calc_eqm( obj,fe )
    dye_df = xx$re*xx$phiq+fe*xx$phiq*xx$dre_df+fe*xx$re*xx$dphiq_df;
    ddye_df = xx$phiq*xx$dre_df + xx$re*xx$dphiq_df;
    fe = fe - dye_df/ddye_df;
    if(abs(dye_df)<1.e-4)
      break;
  }

  fmsy=fe;

  #plug back into calc_eqm to get bmsy and rmsy
  xx <- calc_eqm(obj,fmsy)

  out <- list()
  out$fmsy=fmsy
  out$bmsy=xx$be
  out$rmsy=xx$re
  out$msy=xx$ye
  return(out)
}     # END function calc_fmsy

###################################################################################################
# Function that returns spawners per recruit at fx.
# called BY calc_fx and calc_frep functions below
calc_ssbpr <- function(fe,obj)
{
  x <- calc_eqm(obj,fe)
  return(x$phif)
}     # END function calc_ssbpr

# calc_fx (called by plot_SR()) CODE STOLEN FROM MSER
# Purpose:     fit a spline function to f vs ssbpr, then use a root finder to get FX%. Can
#              get any FX% by changing the value of "target"
# Returns:     a list with all equilibrium quantities for F40%
# Source:      S.P. Cox (modified from .getF40 by K.Holt on 11-Feb-2010)
calc_fx <- function( obj, x=40 )
{
  flist <- seq(0,150, by=0.01)
  maxF <- max( flist )

  ssbpr <- sapply(flist,calc_ssbpr, obj)

  fssbprSplineFun <- splinefun( x=flist,y=ssbpr )
  ssbprAtOrigin   <- fssbprSplineFun( x=0 )

  ssbprRatio <- function( fin )
  {
    f2 <- fssbprSplineFun( fin )
    SPR <- f2/ssbprAtOrigin
    target <- x/100
    #print(target)
    return(SPR - target)
  }     # END function ssbprRatio

  fsprx <- try(uniroot( f=ssbprRatio, interval=c(0,maxF) )$root)
  if ( class( fsprx ) !="try-error" )
  {
    fsprx <- min( fsprx, maxF )

    #plug back into calc_eqm to get bmsy and rmsy
    xx <- calc_eqm(obj, fsprx)
    x0 <- calc_eqm(obj, 0)
    ssbprx <- xx$phif
    ssbpr0 <- x0$phif
    SPR <- ssbprx/ssbpr0

    out <- list()
    out$spr <- SPR
    out$fsprx <-fsprx
    out$bsprx <-xx$be
    out$rsprx <-xx$re
  } else
  {
    out <- fsprx
  }
  return(out)
}     # END function calc_fx

#calc_frep (called by calcRefPoints)
# Purpose:     fit a spline function to f vs recruits per spawner (rps), then use a root finder to get Frep.
#             #Slope of Frep line. Mace and Sissenwine 1993 p 103:
#             "The replacement line is defined as the line with a slope equal
              #to the observed median survival ratio"
# Target rps is calculated in calcRefPoints
# Returns:     a list with all equilibrium quantities for Frep
calc_frep <- function(obj, rps)
{
  flist <- seq(0,150, by=0.01)
  maxF <- max( flist )

  ssbpr <- sapply( flist,calc_ssbpr, obj)

  fssbprSplineFun <- splinefun( x=flist,y=ssbpr )
  ssbprAtOrigin   <- fssbprSplineFun( x=0 )
  ssbprRatio <- function( fin )
  {
    f2 <- fssbprSplineFun( fin )
    #solving for F that gives recruits per spawner corresponding to Frep
    #calc_ssbpr returns spawners per recruit so take inverse for recruits per spawner
    RPS <- 1/f2
    target <- rps
    #print(target)
    return(RPS - target)
  }

  frep <- try(uniroot( f=ssbprRatio,interval=c(0,maxF) )$root)

  if ( class( frep ) !="try-error" )
  {
      frep <- min( frep, maxF )

      # Plug back into calc_eqm to get brep and rrep.
      xx <- calc_eqm(obj, frep)

      out <- list()
      out$frep <-frep
      out$rps <- 1/xx$phif
      out$brep <-xx$be
      out$rrep <-xx$re
  }
  else{
    out <-frep
  }
  return(out)
}     # END function calc_frep

#calc_fext (called by calcRefPoints)
# Purpose:     fit a spline function to f vs recruits per spawner (rps), then use a root finder to get Fext = alpha.
#
# Target rps is calculated in calcRefPoints
# Returns:     a list with all equilibrium quantities for Frep
# DOES NOT WORK FOR VERY HIGH VALUES OF M (>0.6) - UNIROOT CRASHES EVEN WITH FMAX=100000
calc_fext <- function(obj, alpha)
{
  flist <- seq(0,150, by=0.01)
  maxF  <- max( flist )

  ssbpr <- sapply( flist,calc_ssbpr, obj )

  fssbprSplineFun <- splinefun( x=flist,y=ssbpr )
  ssbprAtOrigin   <- fssbprSplineFun( x=0 )

  ssbprRatio <- function( fin )
  {
    f2 <- fssbprSplineFun( fin )
    #solving for F that gives recruits per spawner corresponding to Frep
    #calc_ssbpr returns spawners per recruit so take inverse for recruits per spawner
    RPS <- 1/f2
    target <- alpha
    #print(fin)
    return(RPS - target)
  }

  fext <- try( uniroot( f=ssbprRatio, interval=c(0,maxF) )$root )

  if ( class( fext ) !="try-error" )
  {
    fext <- min( fext, maxF )

    #plug back into calc_eqm to get bmsy and rmsy
    xx       <- calc_eqm(obj, fext)

    out      <- list()
    out$fext <- fext
    out$rps  <- 1/xx$phif
  }
  else
  {
    out <-fext
  }

  return(out)
}     # END function calc_fext

# Code to calculate Fmax
# Fmax maximises YPR.
calc_fmax <- function( obj )
{
  xx <- calcRefCurves(obj, 1000)
  f   <- xx$F
  ypr <- xx$ypr

  plot(f,ypr,type="l", xlab="F", ylab="Yield per Recruit", las=1)

  yprmax <- max(ypr)
  ismax <- match(ypr,yprmax,nomatch=0)

  #get Fmax
  for (i in 1:length(ismax)) if (ismax[i] == 1) Fmax <- f[i]

  if (Fmax == f[length(ismax)]) {
    Fmax <- NA
    bmax <- NA
    rmax <- NA
    cat("\nMSG (calc_fmax) Fmax at upper F boundary, no local maximum \n")
  } else {

    #plug back into calc_eqm to get bmsy and rmsy
    xx       <- calc_eqm(obj, Fmax)
    bmax <-xx$be
    rmax <-xx$re
  }

  out <- list()
  out$fmax<-Fmax
  out$bmax<-bmax
  out$rmax<-rmax

  return(out)
}     # END function calc_fmax

# Code to calculate F01
# F0.1 is 10% of slope of YPR~F curve at origin
calc_f01 <- function( obj )
{
  xx <- calcRefCurves(obj, 1000)
  f <- xx$F
  ypr <- xx$ypr

  #Get slope at origin
  nf <- length(f)
  nf1 <- (nf-1)
  incF <- f[2]-f[1]
  dydF <- vector(length=nf-1)
  for (i in 1:nf1) dydF[i] <- (ypr[i+1]-ypr[i])/incF
  dydF0 <- dydF[1]/10.0    # 10% of slope at origin

  findf01 <- function(obsn)  {
    dydf <- dydF[obsn]
    diff   <- dydF[obsn] - dydF0
    return(diff)
  }

  F01ind <- try( uniroot( f=findf01, interval=c(1,nf1 ))$root )

  if ( class( F01ind ) !="try-error" )
  {
    F01 <- f[F01ind]

    #plug back into calc_eqm to get bmsy and rmsy
    xx       <- calc_eqm(obj, F01)
    b01 <-xx$be
    r01 <-xx$re
    cat("\nMSG (calc_f01) Verifying F0.1 calculation. Two numbers should be very close \n")
    cat( "\nMSG (calc_f01) 10% slope = ", dydF0, " check = ", dydF[F01ind],"\n" )


  }
  else
  {
    F01 <- NA
    b01 <-NA
    r01 <-NA
  }

  out <- list()
  out$f01=F01
  out$b01<-b01
  out$r01<-r01

  return(out)
}     # END function calc_f01

#This function calls all the equilibrium functions and plots them on S-R graph
#Also produces table of eqm values
calcRefPoints <- function( inp )
{
  # Dimension output table.
  #colNames <- c( "Fext","FSPR30","FMSY","Frep","FSPR40","Fmax","F0.1","Unfished" )
  colNames <- c("FMSY","Unfished")
  rowNames <- c( "F", "SSB", "SSB/B0" )
  FTable   <- as.data.frame( matrix( ncol=length(colNames), nrow=length(rowNames)) )
  names( FTable )     <- colNames
  row.names( FTable ) <- rowNames

  refPts <- list()

  # Get parameters of S-R relationship, b0 and r0
  cat( "\nMSG (calcRefPoints) Getting S-R parameters...\n" )
  params <- calc_eqm( inp, 0 )     # Run eqm model with ft=0
  Alpha  <- params$bhalpha
  Beta   <- params$bhbeta

  #these should be the same as bo and ro when ft=0
  #but is this true when spawn timing accounted for?
  be <- params$be
  re <- params$re

  b0     <- params$be
  r0     <- params$re         # Check this should be the same as ro.

  #cat( "\nMSG (calcRefPoints) ro = ",inp$ro," Eqm re = ", re, "\n" )
  #cat( "\nMSG (calcRefPoints)  bo = ",b0,"Eqm be = ", be, "\n" )

  refPts$Alpha <- Alpha
  refPts$Beta  <- Beta
  refPts$b0    <- b0
  refPts$r0    <- r0

  # Get fext by solving for Fext, where slope  R/S = alpha
  # Get frep, brep and rrep

  #cat( "\nMSG (calcRefPoints) Getting Fext, may take a minute...\n" )
  #cat( "\nMSG (calcRefPoints) If this crashes, either increase Fmax in calc_fext or reduce M\n")

  # Solve for Fext.
  #cat( "\nMSG (calcRefPoints) Solving for Extinction reference points...\n" )
  #x   <- calc_fext( inp, Alpha )

  # Did the uniroot fail, i.e., captured by "try" in calc_fext?
  #if ( !is.character(x) )
  #{
  #  rps <- x$rps

  #  cat( "\nMSG (calcRefPoints) Verifying slope at Fext, rps and Alpha should match:\n" )
  #  cat( "\nMSG (calcRefPoints) rps = ",rps," Alpha = ",Alpha,"\n" )
  #  refPts$fext <- x$fext
  #}
  #else
  #{
  #  cat( "\nERR (calcRefPoints) ",x, "\n" )
  #  refPts$fext <- NA
  #}
  # Add to output list.
  #refPts$sext <- 0
  #refPts$rext <- 0
  #refPts$dext <- 0

  # SSB and SSB/B0 are 0 at Fext.
  # FTable[ 1,1 ] <- refPts$fext
  # FTable[ 2,1 ] <- refPts$sext
  # FTable[ 3,1 ] <- refPts$dext

  # Solve for Fmsy, Bmsy and Rmsy
  cat( "\nMSG (calcRefPoints) Solving for MSY reference points...\n" )
  x <- calc_fmsy( inp )

  # Add to output list.
  refPts$msy <- x$msy
  refPts$fmsy <- x$fmsy
  refPts$smsy <- x$bmsy
  refPts$rmsy <- x$rmsy
  refPts$dmsy <- x$bmsy / b0

  # F, SSB and depletion at MSY.
  FTable[1,1] <- refPts$fmsy
  FTable[2,1] <- refPts$smsy
  FTable[3,1] <- refPts$dmsy

  # Solve for F30, B30 and R30
  #cat( "\nMSG (calcRefPoints) Solving for SPR30 reference points...\n" )
 #  x <- calc_fx( inp, 30 )
 #
 #  # Did the uniroot fail, i.e., captured by "try" in calc_fsprx?
 #  if ( !is.character(x) )
 #  {
 #    # Add to output list.
 #    refPts$spr30   <- x$spr
 #    refPts$fsprx30 <- x$fsprx
 #    refPts$ssprx30 <- x$bsprx
 #    refPts$rsprx30 <- x$rsprx
 #    refPts$dsprx30 <- x$bsprx / b0
 #    FTable[1,2] <- refPts$fsprx30
 #    FTable[2,2] <- refPts$ssprx30
 #    FTable[3,2] <- refPts$dsprx30
 #
 #   # cat( "\nMSG (calcRefPts) Spawning Potential Ratio X=0.3 is ", refPts$spr30,"\n" )
 #  }
 #  else
 #  {
 #    #cat( "\nERR (calcRefPoints) ",x, "\n" )
 #    FTable[1,2] <- NA
 #    FTable[2,2] <- NA
 #    FTable[3,2] <- NA
 #  }
 #
 #  # Solve for F40, B40 and R40
 #  #cat( "\nMSG (calcRefPoints) Solving for SPR40 reference points...\n" )
 #  x <- calc_fx( inp, 40 )
 #
 #  # Did the uniroot fail, i.e., captured by "try" in calc_sprx?
 #  if ( !is.character(x) )
 #  {
 #    # Add to output list.
 #    refPts$spr40   <- x$spr
 #    refPts$fsprx40 <- x$fsprx
 #    refPts$ssprx40 <- x$bsprx
 #    refPts$rsprx40 <- x$rsprx
 #    refPts$dsprx40 <- x$bsprx / b0
 #
 #    cat( "\nMSG (calcRefPts) Spawning Potential Ratio X=0.4 is ",refPts$spr40,"\n" )
 #
 #    # F, SSB, and depletion at spr40.
 #    FTable[1,5] <- refPts$fsprx40
 #    FTable[2,5] <- refPts$ssprx40
 #    FTable[3,5] <- refPts$dsprx40
 #  }
 #  else
 #  {
 #    FTable[1,5] <- NA
 #    FTable[2,5] <- NA
 #    FTable[3,5] <- NA
 #  }
 #
 #  #Get Fmax
 #  #Probably won't be a local maximum for this stock
 #  #calc_fmax returns NA if no local maximum
 #  #cat( "\nMSG (calcRefPoints) Solving for Fmax reference points...\n" )
 #  x <- calc_fmax(inp)
 #
 #  # Add to output list.
 #  refPts$fmax <- x$fmax
 #  refPts$smax <- x$bmax
 #  refPts$rmax <- x$rmax
 #  refPts$dmax <- x$bmax / b0
 #
 #  FTable[1,6] <- refPts$fmax
 #  FTable[2,6] <- refPts$smax
 #  FTable[3,6] <- refPts$dmax
 #
 #  #Get F0.1
 # # cat( "\nMSG (calcRefPoints) Solving for F0.1 reference points...\n" )
 #  x <- calc_f01(inp)
 #
 #  # Add to output list.
 #  refPts$f01 <- x$f01
 #  refPts$s01 <- x$b01
 #  refPts$r01 <- x$r01
 #  refPts$d01 <- x$b01 / b0
 #
 #  FTable[1,7] <- refPts$f01
 #  FTable[2,7] <- refPts$s01
 #  FTable[3,7] <- refPts$d01
 #

  # Get Beverton-Holt S-R curve from parameters.
  # For Herring, B0 is quite a bit lower than max obsSSB so draw curve out
  # to end of observations
  S <- seq( 0,max(inp$ObsSB), by=1 ) #RF reverted to this Dec 18
  #S <- seq( 0,b0, by=1 )
  R <- (Alpha*S) / (1+Beta*S)

  # Slope of Frep line. Mace and Sissenwine 1993 p 103:
  # "The replacement line is defined as the line with a slope equal to
  #  the observed median survival ratio."

  ObsRec     <- inp$ObsRec
  ObsSB      <- inp$ObsSB
  SurvRatios <- ObsRec/ObsSB         # R[2:length(R)]/S[2:length(S)]
  FrepSlope  <- median( SurvRatios )  #Dec 18, changed back to median

  #write.csv(cbind(inp$ObsRec,inp$ObsSB),file="ObsSR.csv", row.names=F)
  #cat( "\nMSG (calcRefPoints) printing survival ratios...\n" )
  #print(SurvRatios)

  # Solve for Frep, Brep and Rrep
  #cat( "\nMSG (calcRefPoints) Solving for Replacement reference points...\n" )
  # x   <- calc_frep( inp, FrepSlope )
  #
  # # Did the uniroot fail, i.e., captured by "try" in calc_sprx?
  # if ( !is.character(x) )
  # {
  #   rps <- x$rps
  #
  #   # Add to output list.
  #   refPts$rpsrep <- x$rps
  #   refPts$frep   <- x$frep
  #   refPts$srep   <- x$brep
  #   refPts$rrep   <- x$rrep
  #   refPts$drep   <- x$brep / b0
  #
  #   cat( "\nMSG (calcRefPoints) Verifying slope at Frep, rps and FrepSlope should match:\n" )
  #   cat( "\nMSG (calcRefPoints) rps = ", rps, " FrepSlope = ", FrepSlope,"\n" )
  #
  #   # F, SSB and depletion at Frep.
  #   FTable[1,4] <- refPts$frep
  #   FTable[2,4] <- refPts$srep
  #   FTable[3,4] <- refPts$drep
  # }
  # else
  # {
  #   FTable[1,4] <- NA
  #   FTable[2,4] <- NA
  #   FTable[3,4] <- NA
  # }

  # Unfished.
  refPts$f0 <- 0
  refPts$s0 <- b0
  refPts$d0 <- 1

  # F, SSB, and depletion at unfished.
  FTable[1,2] <- refPts$f0
  FTable[2,2] <- refPts$s0
  FTable[3,2] <- refPts$d0

  refPts$ObsRec <- ObsRec
  refPts$ObsSB  <- ObsSB
  refPts$R      <- R
  refPts$S      <- S

  # Slope of unfished line (unfished juv survival rate)
  refPts$UnfSlope  <- r0 / b0
  # Slope of Fmsy line (juv survival rate at fmsy)
  refPts$fmsySlope <- refPts$rmsy / refPts$smsy
  # Slope of F40 line (juv survival rate at f40)
  #refPts$f40Slope  <- refPts$rsprx40 / refPts$ssprx40
  # Slope of F30 line (juv survival rate at f30)
  #refPts$f30Slope  <- refPts$rsprx30 / refPts$ssprx30
  # Slope of F0.1 line (juv survival rate at f0.1)
  #refPts$f01Slope  <- refPts$r01 / refPts$s01
  # Slope of Frep line (juv survival rate at frep)
  #refPts$FrepSlope <- FrepSlope

  val <- list( refPts=refPts, FTable=FTable )

  return( val )
}     # END function calcRefPts.


plotSsbF <- function( x, curves,
                      refPtCols=c("red","black","green","orange","purple", "darkgray","blue") )
{
  F   <- curves$F
  ssb <- curves$ssb

  xLim <- c( 0, max(ssb ) )
  yLim <- c( 0, max(F) )
  plot( ssb, F, type="n", axes=FALSE, xlab="", xlim=xLim, ylim=yLim )
  lines( ssb, F, col="black", lty=1, lwd=2 )

  points( x$sext,    x$fext,    bg=refPtCols[1], cex=1.2, pch=21 )
  points( x$srep,    x$frep,    bg=refPtCols[2], cex=1.2, pch=21 )
  points( x$smsy,    x$fmsy,    bg=refPtCols[3], cex=1.2, pch=21 )
  points( x$ssprx40, x$fsprx40, bg=refPtCols[4], cex=1.2, pch=21 )
  points( x$ssprx30, x$fsprx30, bg=refPtCols[5], cex=1.2, pch=21 )
  points( x$b0,      x$f0,      bg=refPtCols[6], cex=1.2, pch=21 )

  axis( side=1 )
  axis( side=2, las=2 )
  axis( side=3, labels=FALSE )
  axis( side=4, labels=FALSE )
  box()

  panLegend( 0.70,0.95,
             legTxt=c( "Fext","Frep","Fmsy","F30","F40","F0" ),
             bg="white", pt.bg=refPtCols, pt.cex=1.2, lty=NA, lwd=2, pch=21  )

  mtext( side=1, line=2, "Spawning Biomass (000s t)" )
  mtext( side=2, line=3, "Fishing mortality (F)")
}     # END function plotSsbF


plotSelGear123withMaturity <- function( obj, label="Herring", scale=TRUE )
{
  # Plot 2d selectivity curve for Gear 1,2,3, overlay with Maturity ogive
  ages  <- obj$ages
  sel   <- obj$selAge[1:3,]
  mat   <- obj$ma

  ng    <- nrow(sel)
  icols <- c( "red","green","blue","black" )

  #yy <- matrix( NA, nrow=3, ncol=length(xx) )

  # This removes the first column (gear index) and takes the first row for gear i.
  #for ( ig in 1:ng )
  #{
  #tmp       <- logSel[ logSel[,1]==ig, ]          # Subset matrix of selectivities for gear ig.
  #tmp       <- as.vector( tmp[ 1,2:ncol(tmp) ] )  # Selectivity is constant, 1st row good, remove 1st col.
  #yy[ ig, ] <- exp( tmp )                         # log selectivity to selectivity.

  # if scale, then scale so that max selectivity is one (#JSC 15Aug13)
  #  if ( scale )
  #    yy[ig,] <- yy[ig,] / max(yy[ig,])
  #}

  xLim <- range( ages )
  #yLim <- c( 0,max(yy) )
  yLim <- c( 0,max(sel) )

  plot( xLim, yLim, type="n",axes=FALSE, xlab="",xlim=xLim, ylab="", ylim=yLim )
  for ( ig in 1:ng )
    lines( ages,sel[ig,], lty=1, lwd=3, col=icols[ig] )

  # Overlay maturity ogive
  lines( ages,mat, lty=1, lwd=3, col=icols[4] )

  axis( side=1 )
  axis( side=2, las=2 )
  axis( side=3, labels=FALSE )
  axis( side=4, labels=FALSE )
  box()

  panLab( 0.8, 0.8, cex=1.1, label )

  txt = c( "Winter fishery","Roe seine","Roe gillnet","Maturity-at-Age" )
  legend( "bottomright",legend=txt, col=icols, bty='n', lty=1, lwd=3, pch=-1, ncol=1 )

}     # END function plotSelGear123withMaturity


plotSR <- function( x, label="Herring",
                    refPtCols=c("red","black","green","orange","purple","pink","blue"),
                    timeCols=FALSE )
{
  # Observed spawning biomass and estimated recruitments (Age-2)
  ObsSB  <- x$ObsSB
  ObsRec <- x$ObsRec

  # Stock-recruitment curve from parameters.
  S      <- x$S
  R      <- x$R

  # Plot limits.
  xLim <- c( 0,max(ObsSB)*1.05 )
  #xLim <- c(0,x$b0)
  yLim <- c( 0,max(ObsRec)*1.05 )

  plot( ObsSB, ObsRec, type="n", axes=FALSE,main=label, xlab="", xlim=xLim, ylab="", ylim=yLim )
  clip( xLim[1], xLim[2], yLim[1],yLim[2] )

  lines( S,R, col="black", lty=1, lwd=3 )

  # Juvenile survival near the origin.
  if ( !is.na(x$Alpha) )
    abline( 0,x$Alpha,     lty=2, lwd=2, col=refPtCols[1] )

  abline( 0,x$FrepSlope, lty=2, lwd=2, col=refPtCols[2] )
  abline( 0,x$fmsySlope, lty=2, lwd=2, col=refPtCols[3] )
  abline( 0,x$f40Slope,  lty=2, lwd=2, col=refPtCols[4] )
  abline( 0,x$f30Slope,  lty=2, lwd=2, col=refPtCols[5] )
  abline( 0,x$f01Slope,  lty=2, lwd=2, col=refPtCols[6] )
  abline( 0,x$UnfSlope,  lty=2, lwd=2, col=refPtCols[7] )

  #segments( x$srep,    0, x$srep,    x$rrep,    lty=2, lwd=2, col=refPtCols[2] )
  #segments( x$smsy,    0, x$smsy,    x$rmsy,    lty=2, lwd=2, col=refPtCols[3] )
  #segments( x$ssprx40, 0, x$ssprx40, x$rsprx40, lty=2, lwd=2, col=refPtCols[4] )
  #segments( x$ssprx30, 0, x$ssprx30, x$rsprx30, lty=2, lwd=2, col=refPtCols[5] )
  #segments( x$b0,      0, x$b0,      x$r0,      lty=2, lwd=2, col=refPtCols[6] )

  # abline( h=0, col="black", lty=3 )

  # Maybe grayscale the order?
  colVec <- "white"
  if ( timeCols )
    colVec <- rev( gray.colors( length(ObsSB) ) )

  points( ObsSB, ObsRec, bg=colVec, col="black", cex=1.2, pch=21 )

  #panLab( 0.1, 0.9, cex=1.2, label )

  panLegend( 0.70,0.95,
             legTxt=c( "Fext","Frep","Fmsy","F40","F30","F0.1","F0" ),
             bg="white", col=refPtCols, lty=2, lwd=2  )

  axis( side=1 )
  axis( side=2, las=2 )
  axis( side=3, labels=FALSE )
  axis( side=4, labels=FALSE )
  box()

  mtext( side=1, line=2, "Spawning Biomass (000s t)" )
  mtext( side=2, line=3, "Age-2 Recruits (000,000s)")

  return( invisible() )
}     # END function plotSR.


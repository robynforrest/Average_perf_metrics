# AUTHORS: Robyn Forrest (REF) and Rob Kronlund (ARK), Pacific Biological Station
# Based on Herring LRP work -- use Herring iscam outputs for parameters
# Repurposed for Groundfish MP Framework to test alternative approaches for
#   calculating performance metrics in an MSE
# Function for testing FMSY calculations and other ref points for Herring
# Original:  Nov 4-7, 2016 (REF and ARK)
# Revisions: March 4, 2020

# Source the equilibrium functions
source(here::here( "R/Utils/EqmFuns.r"))
# Source functions needed to read iSCAM files.
source( here::here( "R/Utils/read.admb.r"  ))
# Source the function to extract parameters from iscam file
source(here::here( "R/Utils/getPars.r"))
# Source the pop dynamics function
source(here::here( "R/Utils/PopDyn.r"))
# Source various R tools for plotting.
source( here::here( "R/Utils/rTools.r"  ))

#######################################################################################

#--------------------------------------------------------------------#
# Main Code: Gets data, calls eqm calcs, plots, optional long run.
#--------------------------------------------------------------------#
#Get the herring iscam output files
# Set the model case folder (AM1=1 or AM2=2)
caseFld  <- "Data"
baseName <- "WCVIHerring2016qEst"
caseName <- "AM1"

# Set maturity equal to selectivity? Herring selectivity is so far to the right of maturity
# That the stock is extremely resilient to high F
matEqSel <-  TRUE  #TRUE

fileList <- list.files( caseFld )

# mgmtAreas is patched and used to name the repList of *.rep file contents.
# just one area for now
mgmtAreas <- "WCVI"
repFiles <- paste( baseName, ".rep", sep="" )

# Read in *.rep files for each stock.
parList <- as.list( 1:length( mgmtAreas ) )    # ADMB par file.
repList <- as.list( 1:length( mgmtAreas ) )    # ADMB rep file.
for ( i in 1:length(repList) )
{
  parList[[i]] <- lisread( file.path( caseFld, fileList[i], "iscam.par" ) )
  repList[[i]] <- read.rep( file.path( caseFld, fileList[i], repFiles[i] ) )
}

mgmtAreas <- c( "WCVI" )
names( parList ) <- mgmtAreas
names( repList ) <- mgmtAreas

#Set the stocks for reference point calculations e.g., stock<-"HG".
stocks <- mgmtAreas
FTables <- vector( "list", length(stocks) )
names( FTables ) <- stocks

# Indexes
syr  <- 1
nyr  <- 100     # Number of years to run out the long model
yrs  <- syr:nyr
nyrs <- length(yrs)
nrep <- 100

#Set short term and long term
LT <- 70:100
ST <- 1:30

#Outputs
Ftrep <- matrix(nrow=nyr, ncol=nrep)
Ctrep <- matrix(nrow=nyr, ncol=nrep)
Btrep <- matrix(nrow=nyr, ncol=nrep)
BBmsy <- matrix(nrow=nyr, ncol=nrep)
FFmsy <- matrix(nrow=nyr, ncol=nrep)
Cmsy <- matrix(nrow=nyr, ncol=nrep)
USR <- vector(length=nrep)
LRP <- vector(length=nrep)
FMSY <- vector(length=nrep)
BMSY <- vector(length=nrep)
MSY <- vector(length=nrep)
aboveLRP <- matrix(nrow=nyr, ncol=nrep)
aboveUSR <- matrix(nrow=nyr, ncol=nrep)
belowFMSY <- matrix(nrow=nyr, ncol=nrep)
aboveMSY <- matrix(nrow=nyr, ncol=nrep)

for ( i in 1:length(stocks))    #1:length(stocks)
{
  #Get the parameters for Herring but set new values for age, growth and M
  inp <- getPars( parList[[ stocks[i] ]], repList[[ stocks[i] ]], matEqSel=matEqSel,
                  nagenew = 50,
                  vbknew = 0.15,
                  linfnew = 50,
                  Mnew = 0.15)
  for ( j in 1:nrep)
  {
      #Add some noise to F and M
      #names(inp)
      inp$M <- inp$M + rnorm(1,0,0.01)
      inp$M

      # Now run Population dynamics model for nyr years with n replicates
      # Parameters to loop over
      Ft <- seq( 0.1, 0.15, length.out=nyr )
      #Add some noise
      Fnoise <- rnorm(nyr,0,0.055)
      Ftrep[,j] = Ft + Fnoise
      for(k in 1:nyr) { if(Ftrep[k,j]<=0) Ftrep[k,j]=0.1}

      # Run the equilibrum models to get ref points.
      result <- calcRefPoints( inp )

      MSY[j] <- result$refPts$msy
      FMSY[j] <- result$refPts$fmsy
      BMSY[j] <- result$refPt$smsy
      LRP[j] <- 0.4*BMSY[j]
      USR[j] <- 0.8*BMSY[j]

      cat("\nRunning popDyn model ... \n")

      # Run PopDynPopDyn  function
        Ft <- Ftrep[,j]

        #  Call the model and get catch under each F
        x <- PopDyn( inp, Ft, nyrs )

        # get the catch and biomass
        Ctrep[,j]  <- x$Ctot
        Btrep[,j] <- x$SBt

        BBmsy[,j] <- Btrep[,j]/BMSY[j]
        FFmsy[,j] <- Ftrep[,j]/FMSY[j]
        Cmsy[,j] <- Ctrep[,j]/MSY[j]

        aboveLRP[,j] <- ifelse(Btrep[,j]>LRP[j],1,0)
        aboveUSR[,j] <- ifelse(Btrep[,j]>USR[j],1,0)
        belowFMSY[,j] <- ifelse(FFmsy[,j]<1,1,0)
        aboveMSY[,j] <- ifelse(Cmsy[,j]>1,1,0)

      }
}

#Three ways of getting performance metrics
#1. probability is calculated for each and every year
#2. probability is calculated within replicates then averaged
#3. probability is averaged across all replicates and years

ntot <- nyr * nrep
ntotLT <- length(LT)*nrep
pLTLRP <- vector(length=3)
meanLTBBmsy <- vector(length=2)
meanCt <- vector(length=2)

#Method 1
# ICES risk3 is the maximum probability that B < Blim, where maximum is taken across years
# risk3 should be less than 5%, therefore p(B > Blim) should be greater than 95%
#http://ices.dk/sites/pub/Publication%20Reports/Guidelines%20and%20Policies/12.04.10_Criteria_for_defining_multi-annual_plans_as_precautionary.pdf
# We will take the complement
pLTLRP_yr <- apply(aboveLRP[LT,],1,mean)
risk3_complement <- min(pLTLRP_yr) #we will take the minimum prob that ssb is above blim
pLTLRP[1] <- risk3_complement

#Method 2 - mseR way. Calculate probs for each replicate then average
pLTLRP_rep <- apply(aboveLRP[LT,],2,mean) #get mean Bt > Blim in each replicate
pLTLRP[2] <- mean(pLTLRP_rep) #take the average across replicates

meanLTBBmsy_rep <- apply(BBmsy[LT,],2,mean)
meanLTBBmsy[1] <- mean(meanLTBBmsy_rep)

meanCt_rep <- apply(Ctrep,2,mean)
meanCt[1] <- mean(meanCt_rep)

#Method 3. MPF way. Calculate probs across reps and years
pLTLRP[3] <- mean(aboveLRP[LT,])
meanLTBBmsy[2] <- mean(BBmsy[LT,])
meanCt[2] <- mean(Ctrep)


matplot(BBmsy)
print(pLTLRP)
print(meanLTBBmsy)
print(meanCt)

#Yes as Rob suspected the means are the same, but what are the distributions
 #of the BBmsy and catches?
#B/Bmsy
dmeanBBmsy1 <- density(meanLTBBmsy_rep)
dmeanBBmsy2 <- density(BBmsy)
#Catches
dmeanCt1 <- density(meanCt_rep)
dmeanCt2 <- density(Ctrep)

plot(dmeanBBmsy1)
plot(dmeanBBmsy2)
plot(dmeanCt1)
plot(dmeanCt2)

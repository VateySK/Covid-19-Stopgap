graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
library(ggplot2)
library(ggpubr)
library(ks)
library(rjags)
library(runjags)
library(nimble)
library("PerformanceAnalytics")
library(psych)
library(GGally)
library(summarytools)
library(knitr)
library(dplyr)
library(data.table)
library(splitstackshape)
library(infotheo)



setwd("D:/RMIT Master of Analytics/semester 3/MATH2269 Applied Bayesian Statistics/project")
#setwd("C:/Users/jayma/Desktop/Millie/RMIT Master of Analytics/semester 3")
source("DBDA2E-utilities.R")      


#===============PRELIMINARY FUNCTIONS FOR POSTERIOR INFERENCES====================
genMCMC = function( data , zName="z" , NName="N" , sName="s" , cName="c" ,
                    numSavedSteps=5000 , saveName=NULL , thinSteps=1 ,
                    runjagsMethod=runjagsMethodDefault ,   useRunjags = TRUE,
                    nChains=nChainsDefault,   burnInSteps = 500 , adaptSteps = 500 , forInits = NULL   ) { 
  require(rjags)
  require(runjags)
  #-----------------------------------------------------------------------------
  # THE DATA.
  # N.B.: This function expects the data to be a data frame, 
  # with one component z being a vector of integer # successes,
  # one component N being a vector of integer # attempts,
  # one component s being a factor of subject identifiers,
  # and one component c being a factor of category identifiers, with
  # subjects nested in categories.
  print(cName)
  z = data[[zName]]
  N = data[[NName]]
  s = data[[sName]]
  c = data[[cName]]
  Nsubj = length(unique(s))
  Ncat =  length(unique(c))
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    z = z ,
    # N = N ,
    c = as.numeric(c) , # c in JAGS is numeric, in R is possibly factor
    Nsubj = Nsubj ,
    Ncat = Ncat
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  model {
    for ( sIdx in 1:Nsubj ) {
      z[sIdx] ~ dbern( theta[sIdx])
      theta[sIdx] ~ dbeta( omega[c[sIdx]]*(kappa[c[sIdx]]-2)+1 , 
                           (1-omega[c[sIdx]])*(kappa[c[sIdx]]-2)+1 ) 
    }
    

    for ( cIdx in 1:Ncat ) {
      omega[cIdx] ~ dbeta( omegaO*(kappaO-2)+1 , 
                           (1-omegaO)*(kappaO-2)+1 )
      kappa[cIdx] <- kappaMinusTwo[cIdx] + 2
      kappaMinusTwo[cIdx] ~ dgamma( 0.01 , 0.01 ) # mean=1 , sd=10 (generic vague)
    }
    #omegaO ~ dbeta( 1.0 , 1.0 ) 
    #omegaO ~ dbeta( 1 , 99 ) # mode=0 , concentration=100
    #omegaO ~ dbeta( 10.8 , 89.2 ) # mode=0.1 , concentration=100
    #omegaO ~ dbeta( 50 , 50 ) # mode=0.5 , concentration=100
    omegaO ~ dbeta( 89.2 , 10.8 ) # mode=0.9 , concentration=100
    
    kappaO <- kappaMinusTwoO + 2
    #kappaMinusTwoO ~ dgamma( 0.01 , 0.01 )  # mean=1 , sd=10 (generic vague)
    kappaMinusTwoO ~ dgamma( 1.01005 , 0.01005012 )  # mode=1 , sd=100
    # kappaMinusTwoO ~ dgamma( 2.6 , 26.96 )  # mode=10 , sd=2
    #kappaMinusTwoO ~ dgamma( 1.105125 , 0.1051249 )  # mode=1 , sd=10
    #kappaMinusTwoO ~ dgamma( 1.105125 , 0.01051249 )  # mode=10 , sd=100
    #kappaMinusTwoO ~ dgamma( 5.05 , 101.99 )  # mode=20 , sd=2
    #kappaMinusTwoO ~ dgamma( 0.01 , 1.22)  # mode=20 , sd=100

    #kappaMinusTwoO ~ dgamma( 10.02 , 402 )  # mode=40 , sd=2
    #kappaMinusTwoO ~ dgamma( 0.01 , 1.49 )  # mode=40 , sd=100    
  }
  " # close quote for modelString
  writeLines( modelString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  # initsList = function() {
  #   thetaInit = rep(NA,Nsubj)
  #   for ( sIdx in 1:Nsubj ) { # for each subject
  #     resampledZ = rbinom(1, size=N[sIdx] , prob=z[sIdx]/N[sIdx] )
  #     thetaInit[sIdx] = resampledZ/N[sIdx]
  #   }
  #   thetaInit = 0.001+0.998*thetaInit # keep away from 0,1    
  #   kappaInit = 100 # lazy, start high and let burn-in find better value
  #   return( list( theta=thetaInit , 
  #                 omega=aggregate(thetaInit,by=list(c),FUN=mean)$x ,
  #                 omegaO=mean(thetaInit) ,
  #                 kappaMinusTwo=rep(kappaInit-2,Ncat) ,
  #                 kappaMinusTwoO=kappaInit-2 ) )
  # }
  initsList = list( theta=rep(0.5,891) , 
                    omega=rep(0.5,3),
                    omegaO=0.5 ,
                    kappaMinusTwoO=3,
                    kappaMinusTwo=rep(3,3))
  
  # initsList = list( theta=forInits[9:899,3], 
  #                   omega=forInits[1:3,3],
  #                   omegaO=forInits[4,3] ,
  #                   kappaMinusTwoO=forInits[8,3],
  #                   kappaMinusTwo=forInits[5:7,3])
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "theta","omega","kappa","omegaO","kappaO") 
  
  
  if ( useRunjags ) {
    runJagsOut <- run.jags( method=runjagsMethod ,
                            model="TEMPmodel.txt" , 
                            monitor=parameters , 
                            data=dataList ,  
                            # inits=initsList , 
                            n.chains=nChains ,
                            adapt=adaptSteps ,
                            burnin=burnInSteps , 
                            sample=ceiling(numSavedSteps/nChains) ,
                            thin=thinSteps ,
                            summarise=FALSE ,
                            plots=FALSE )
    codaSamples = as.mcmc.list( runJagsOut )
  } else {
    # Create, initialize, and adapt the model:
    jagsModel = jags.model( "TEMPmodel.txt" , data=dataList , inits=initsList , 
                            n.chains=nChains , n.adapt=adaptSteps )
    # Burn-in:
    cat( "Burning in the MCMC chain...\n" )
    update( jagsModel , n.iter=burnInSteps )
    # The saved MCMC chain:
    cat( "Sampling final MCMC chain...\n" )
    codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                                n.iter=ceiling(numSavedSteps*thinSteps/nChains), 
                                thin=thinSteps )
  }  
  
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
} # end function

#======================================================================================

smryMCMC = function(  codaSamples , compVal=0.5 , rope=NULL , 
                      diffSVec=NULL , diffCVec=NULL , 
                      compValDiff=0.0 , ropeDiff=NULL , 
                      saveName=NULL ) {
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  summaryInfo = NULL
  rowIdx = 0
  # omega:
  for ( parName in grep("omega",colnames(mcmcMat),value=TRUE) ) {
    summaryInfo = rbind( summaryInfo , 
                         summarizePost( mcmcMat[,parName] ,
                                        compVal=compVal , ROPE=rope ) )
    rowIdx = rowIdx+1
    rownames(summaryInfo)[rowIdx] = parName
  }
  # kappa:
  for ( parName in grep("kappa",colnames(mcmcMat),value=TRUE) ) {
    summaryInfo = rbind( summaryInfo , 
                         summarizePost( mcmcMat[,parName] ,
                                        compVal=NULL , ROPE=NULL ) )
    rowIdx = rowIdx+1
    rownames(summaryInfo)[rowIdx] = parName
  }
  # theta:
  for ( parName in grep("theta",colnames(mcmcMat),value=TRUE) ) {
    summaryInfo = rbind( summaryInfo , 
                         summarizePost( mcmcMat[,parName] ,
                                        compVal=compVal , ROPE=rope ) )
    rowIdx = rowIdx+1
    rownames(summaryInfo)[rowIdx] = parName
  }
  # differences of theta's:
  if ( !is.null(diffSVec) ) {
    Nidx = length(diffSVec)
    for ( t1Idx in 1:(Nidx-1) ) {
      for ( t2Idx in (t1Idx+1):Nidx ) {
        parName1 = paste0("theta[",diffSVec[t1Idx],"]")
        parName2 = paste0("theta[",diffSVec[t2Idx],"]")
        summaryInfo = rbind( summaryInfo , 
                             summarizePost( mcmcMat[,parName1]-mcmcMat[,parName2] ,
                                            compVal=compValDiff , ROPE=ropeDiff ) )
        rowIdx = rowIdx+1
        rownames(summaryInfo)[rowIdx] = paste0(parName1,"-",parName2)
      }
    }
  }
  # differences of omega's:
  if ( !is.null(diffCVec) ) {
    Nidx = length(diffCVec)
    for ( t1Idx in 1:(Nidx-1) ) {
      for ( t2Idx in (t1Idx+1):Nidx ) {
        parName1 = paste0("omega[",diffCVec[t1Idx],"]")
        parName2 = paste0("omega[",diffCVec[t2Idx],"]")
        summaryInfo = rbind( summaryInfo , 
                             summarizePost( mcmcMat[,parName1]-mcmcMat[,parName2] ,
                                            compVal=compValDiff , ROPE=ropeDiff ) )
        rowIdx = rowIdx+1
        rownames(summaryInfo)[rowIdx] = paste0(parName1,"-",parName2)
      }
    }
  }
  # save:
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  show( summaryInfo )
  return( summaryInfo )
}

#===============================================================================

plotMCMC = function( codaSamples , 
                     data , zName="z" , NName="N" , sName="s" , cName="c" ,                      compVal=0.5 , rope=NULL , 
                     diffSList=NULL , diffCList=NULL , 
                     compValDiff=0.0 , ropeDiff=NULL , 
                     saveName=NULL , saveType="jpg" ) {
  #-----------------------------------------------------------------------------
  # N.B.: This function expects the data to be a data frame, 
  # with one component z being a vector of integer # successes,
  # one component N being a vector of integer # attempts,
  # one component s being a factor of subject identifiers,
  # and one component c being a factor of category identifiers, with
  # subjects nested in categories.
  z = data[[zName]]
  N = data[[NName]]
  s = data[[sName]]
  c = data[[cName]]
  Nsubj = length(unique(s))
  Ncat =  length(unique(c))
  # Now plot the posterior:
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  
  # kappa:
  parNames = sort(grep("kappa",colnames(mcmcMat),value=TRUE))
  nPanels = length(parNames)
  nCols = 3
  nRows = ceiling(nPanels/nCols)
  openGraph(width=2.5*nCols,height=2.0*nRows)
  par( mfcol=c(nRows,nCols) )
  par( mar=c(3.5,4,3.5,4) , mgp=c(2.0,0.7,0) )
  #xLim = range( mcmcMat[,parNames] )
  xLim=quantile(mcmcMat[,parNames],probs=c(0.000,0.995))
  mainLab = c(paste( levels(factor(data[[cName]]))),"Overall")
  print(paste("levels=", levels(factor(data[[cName]]))))
  
  mainIdx = 0
  for ( parName in parNames ) {
    mainIdx = mainIdx+1
    print(paste("Kappa Title=", mainLab[mainIdx]))
    
    postInfo = plotPost( mcmcMat[,parName] , compVal=compVal , ROPE=rope ,
                         xlab=bquote(.(parName)) , cex.lab=1.25 , 
                         main=mainLab[mainIdx] , cex.main=1.5 ,
                         xlim=xLim , border="skyblue" )
  }  
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"Kappa",sep=""), type=saveType)
  }
  
  # omega:
  parNames = sort(grep("omega",colnames(mcmcMat),value=TRUE))
  nPanels = length(parNames)
  nCols = 3
  nRows = ceiling(nPanels/nCols)
  openGraph(width=2.5*nCols,height=2.0*nRows)
  par( mfcol=c(nRows,nCols) )
  par( mar=c(3.5,4,3.5,4) , mgp=c(2.0,0.7,0) )
  #xLim = range( mcmcMat[,parNames] )
  xLim=quantile(mcmcMat[,parNames],probs=c(0.001,0.999))
  mainLab = c(paste(levels(factor(data[[cName]]))),"Overall")
  mainIdx = 0
  for ( parName in parNames ) {
    
    mainIdx = mainIdx+1
    print(paste("Omega Title=", mainLab[mainIdx]))
    
    postInfo = plotPost( mcmcMat[,parName] , compVal=compVal , ROPE=rope ,
                         xlab=bquote(.(parName)) , cex.lab=1.25 , 
                         main=mainLab[mainIdx] , cex.main=1.5 ,
                         xlim=xLim , border="skyblue" )
  }  
  if ( !is.null(saveName) ) {
    saveGraph( file=paste(saveName,"Omega",sep=""), type=saveType)
  }
  
  # Plot individual omega's and differences:
  if ( !is.null(diffCList) ) {
    for ( compIdx in 1:length(diffCList)) {
      diffCVec = diffCList[[compIdx]]
      Nidx = length(diffCVec)
      temp=NULL
      mainLab = c(paste(levels(factor(data[[cName]]))),"Overall")
      
      for ( i in 1:Nidx ) {
        temp = c( temp , which(levels(factor((data[[cName]])))==diffCVec[i]) )
      }
      #diffCVec = temp
      print(paste("Nidx=", Nidx, "diffCList", diffCList, "diffCVec", diffCVec))
       
      openGraph(width=2.5*Nidx,height=2.0*Nidx)
      par( mfrow=c(Nidx,Nidx) )
      xLim = range(c( compVal, rope,
                      mcmcMat[,paste0("omega[",diffCVec,"]")] ))
      for ( t1Idx in 1:Nidx ) {
        for ( t2Idx in 1:Nidx ) {
          parName1 = paste0("omega[",diffCVec[t1Idx],"]")
          parName2 = paste0("omega[",diffCVec[t2Idx],"]")
          if ( t1Idx > t2Idx) {  
            # plot.new() # empty plot, advance to next
            par( mar=c(3,1,3,1) , mgp=c(2.0,0.7,0) , pty="s" )
            nToPlot = 700
            ptIdx = round(seq(1,chainLength,length=nToPlot))
            
            print(paste("scatter plot x=", diffCVec[t2Idx]))
            print(paste("scatter plot y=", diffCVec[t1Idx]))
            
            print(paste("scatter plot x=", levels(factor(data[[cName]]))[diffCVec[t2Idx]]))
            print(paste("scatter plot y=", levels(factor(data[[cName]]))[diffCVec[t1Idx]]))

            plot ( mcmcMat[ptIdx,parName2] , mcmcMat[ptIdx,parName1] , 
                   cex.main=1.25 , cex.lab=1.25 , 
                     xlab=paste(levels(factor(data[[cName]]))[diffCVec[t2Idx]]) , 
                     ylab=paste(levels(factor(data[[cName]]))[diffCVec[t1Idx]]) , 
                   col="skyblue" )
            abline(0,1,lty="dotted")
          } else if ( t1Idx == t2Idx ) {
            par( mar=c(3,3.5,3,1.5) , mgp=c(2.0,0.7,0) , pty="m" )
            postInfo = plotPost( mcmcMat[,parName1] , 
                                 compVal=compVal , ROPE=rope , 
                                 cex.main=1.25 , cex.lab=1.25 , 
                                 xlab=bquote(.(parName1)) ,
                                 main=paste(levels(factor(data[[cName]]))[diffCVec[t1Idx]]) ,  
                                 xlim=xLim )
          } else if ( t1Idx < t2Idx ) {
            par( mar=c(3,1.5,3,1.5) , mgp=c(2.0,0.7,0) , pty="m" )
            postInfo = plotPost( mcmcMat[,parName1]-mcmcMat[,parName2] , 
                                 compVal=compValDiff , ROPE=ropeDiff , 
                                 cex.main=1.0 , cex.lab=1.25 , 
                                 xlab=bquote("Difference of "*omega*"'s"), 
                                  main=paste0( 
                                     levels(factor(data[[cName]]))[diffCVec[t1Idx]] ,
                                     "\n - ",
                                     levels(factor(data[[cName]]))[diffCVec[t2Idx]]))
          }
        }
      }
      if ( !is.null(saveName) ) {
        saveGraph( file=paste0(saveName,"OmegaDiff",compIdx), type=saveType)
      }
    }
  }
  
  # Plot individual theta's and differences:
  if ( !is.null(diffSList) ) {
    for ( compIdx in 1:length(diffSList) ) {
      diffSVec = diffSList[[compIdx]]
      Nidx = length(diffSVec)
      temp=NULL
      for ( i in 1:Nidx ) {

            temp = c( temp , data[[sName]]==diffSVec[i])
      }
      #diffSVec = temp
      openGraph(width=2.5*Nidx,height=2.0*Nidx)
      par( mfrow=c(Nidx,Nidx) )
      xLim = range(c(compVal,rope,mcmcMat[,paste0("theta[",diffSVec,"]")],
                     z[diffSVec]/N[diffSVec]))
      for ( t1Idx in 1:Nidx ) {
        for ( t2Idx in 1:Nidx ) {
          parName1 = paste0("theta[",diffSVec[t1Idx],"]")
          parName2 = paste0("theta[",diffSVec[t2Idx],"]")
          if ( t1Idx > t2Idx) {  
            # plot.new() # empty plot, advance to next
            par( mar=c(3,3,3,1) , mgp=c(2.0,0.7,0) , pty="s" )
            nToPlot = 700
            ptIdx = round(seq(1,chainLength,length=nToPlot))
            plot ( mcmcMat[ptIdx,parName2] , mcmcMat[ptIdx,parName1] , cex.lab=1.25 ,
                   xlab=s[diffSVec[t2Idx]] , 
                   ylab=s[diffSVec[t1Idx]] , 
                   col="skyblue" )
            abline(0,1,lty="dotted")
          } else if ( t1Idx == t2Idx ) {
            par( mar=c(3,3.5,3,1.5) , mgp=c(2.0,0.7,0) , pty="m" )
            postInfo = plotPost( mcmcMat[,parName1] , 
                                 compVal=compVal , ROPE=rope , 
                                 cex.main=1.0 , cex.lab=1.25 , 
                                 xlab=bquote(.(parName1)) ,
                                 main=paste0( s[diffSVec[t1Idx]], 
                                              " (",c[diffSVec[t1Idx]],")") ,  
                                 xlim=xLim )
            # points( z[diffSVec]/N[diffSVec] , 0 , 
            #         pch="+" , col="red" , cex=3 )
            # text( z[diffSVec[t1Idx]]/N[diffSVec[t1Idx]] , 0 , 
            #       bquote(list( z==.(z[diffSVec[t1Idx]]) ,
            #                    N==.(N[diffSVec[t1Idx]]) )) , 
            #       adj=c( (z[diffSVec[t1Idx]]/N[diffSVec[t1Idx]]-xLim[1])/
            #                (xLim[2]-xLim[1]),-3.25) , col="red" )
          } else if ( t1Idx < t2Idx ) {
            par( mar=c(3,1.5,3,1.5) , mgp=c(2.0,0.7,0) , pty="m" )
            postInfo = plotPost( mcmcMat[,parName1]-mcmcMat[,parName2] , 
                                 compVal=compValDiff , ROPE=ropeDiff , 
                                 cex.main=1.0 , cex.lab=1.25 , 
                                 xlab=bquote("Difference of "*theta*"'s") , 
                                 main=paste( 
                                   s[diffSVec[t1Idx]] , 
                                   " (",c[diffSVec[t1Idx]],")" ,
                                   "\n -",
                                   s[diffSVec[t2Idx]] , 
                                   " (",c[diffSVec[t2Idx]],")" ) )
            # points( z[diffSVec[t1Idx]]/N[diffSVec[t1Idx]]
            #         - z[diffSVec[t2Idx]]/N[diffSVec[t2Idx]] , 0 , 
            #         pch="+" , col="red" , cex=3 )
          }
        }
      }
      if ( !is.null(saveName) ) {
        saveGraph( file=paste0(saveName,"ThetaDiff",compIdx), type=saveType)
      }
    }
  }
}

#===============PRELIMINARY FUNCTIONS FOR POSTERIOR INFERENCES====================

myData <- read.csv("covid_brazil_final.csv")
describe(myData)
#x = as.matrix(myData[,c("Patient.age.quantile","Hemoglobin","Hematocrit","Mean.platelet.volume", "Lymphocytes", "Mean.corpuscular.hemoglobin.concentration",  "Neutrophils", "Monocytes", "Eosinophils", "Basophils")])

# Some more descriptives
# cat("\nCORRELATION MATRIX OF PREDICTORS:\n ")
# show( round(cor(x),3) )
# cat("\n")
# 
# 
# #Pairplots for all variables
# chart.Correlation(x, histogram=TRUE, pch=19)
# 
# 
# 
# 
# #Take away Hematocrit and Neutrophils due to high correlation
# myData2 = myData[,c("patient.ID", "Patient.age.quantile","Hemoglobin","Mean.platelet.volume", "Lymphocytes", "Mean.corpuscular.hemoglobin.concentration",  "Monocytes", "Eosinophils", "Basophils", "SARS.Cov.2.exam.result")]
names(myData)[3] <- "covid.results"
myData$covid.results <- as.numeric(as.factor(myData$covid.results)) - 1 # To get 0/1 instead of 1/2; positive = 1; negative = 0

#check for missing values
sum(is.na(myData))

# #Descriptive look
# p1 <- ggplot(myData, aes(x=Patient.age.quantile, y = covid.results)) +
#   geom_point()
# 
# p2 <- ggplot(myData, aes(x=Hemoglobin, y = covid.results)) +
#   geom_point()
# 
# p3 <- ggplot(myData, aes(x=Mean.platelet.volume, y = covid.results)) +
#   geom_point()
# 
# p4 <- ggplot(myData, aes(x=Lymphocytes, y = covid.results)) +
#   geom_point()
# 
# p5 <- ggplot(myData, aes(x=Mean.corpuscular.hemoglobin.concentration, y = covid.results)) +
#   geom_point()
# 
# p6 <- ggplot(myData, aes(x=Monocytes, y = covid.results)) +
#   geom_point()
# 
# p7 <- ggplot(myData, aes(x=Eosinophils, y = covid.results)) +
#   geom_point()
# 
# p8 <- ggplot(myData, aes(x=Basophils, y = covid.results)) +
#   geom_point()
# 
# figure <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 4, ncol = 2)
# figure <- annotate_figure(figure,
#                           top = text_grob("Covid test results vs independent variables", face = "bold", size = 14))
# 
# figure


# THE DATA.
# set.seed(999)
# trainData <- stratified(myData2, "covid.results", .7)
# testData <- setdiff(myData2, trainData)


#x = as.matrix(myData2[,c("Patient.age.quantile","Hemoglobin","Mean.platelet.volume", "Lymphocytes", "Mean.corpuscular.hemoglobin.concentration",  "Monocytes", "Eosinophils", "Basophils")])
#y = unlist(myData2[, "covid.results"])
zName = "covid.results" # column name for 0,1 values
sName = "Patient.ID" # column name for subject ID
cName = "Age.group"

myData$Age.group= 0
myData$Age.group[which(myData$Patient.age.quantile < 5)] <- 1 
myData$Age.group[which(myData$Patient.age.quantile > 4)] <- 2 
myData$Age.group[which(myData$Patient.age.quantile > 9)] <- 3 
myData$Age.group[which(myData$Patient.age.quantile > 14)] <- 4

myData$Age.group <- factor(myData$Age.group,
                    levels = c(1,2,3, 4),
                    labels = c("Age Group 1", "Age Group 2", "Age Group 3", "Age Group 4"))

write.csv(myData, "myData.csv" )

tab <- table(myData$Age.group, myData$covid.results)
tab <- cbind(tab, Total = rowSums(tab))
tab

numAdaptInSteps = 8000 ; numBburnInSteps=8000; numSavedSteps=12000 ; thinSteps=3000; nChains = 3
#numAdaptInSteps = 8000 ; numBburnInSteps=8000; numSavedSteps=8000 ; thinSteps=800; nChains = 3


startTime = proc.time()


mcmcCoda = genMCMC( data=myData , 
                    zName=zName, sName=sName, cName=cName,
                    numSavedSteps=numSavedSteps ,   useRunjags = TRUE,
                    thinSteps=thinSteps , burnInSteps = numBburnInSteps , adaptSteps = numAdaptInSteps ,nChains = nChains)


stopTime = proc.time()
duration = stopTime - startTime
show(duration)


save.image(file="rEnvironment_2000Thin_12thousands_Age_trial9.RData")
#load(file="rEnvironment_300Thin_30thousands_1_withAge.RData")

load(file="rEnvironment_800Thin_8thousands_Age_trial7.RData")
#save.image(file="rEnvironment_800Thin_8thousands_Age_trial3.RData")
#load(file="rEnvironment_800Thin_8thousands_Age_trial3.RData")

# save.image(file="rEnvironment_11Thin_Fullobs.RData")
# load(file="rEnvironment_11Thin_Fullobs.RData")

#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
#
------------------------- 
parameterNames = varnames(mcmcCoda) # get all parameter names for reference
for ( parName in c("omega[1]","omega[2]","omega[3]","omega[4]","omegaO","kappa[1]","kappa[2]","kappa[3]","kappa[4]","kappaO","theta[1]","theta[2]","theta[3]") ) { 
  diagMCMC( codaObject=mcmcCoda , parName=parName) 
}

# Get summary statistics of chain:

summaryInfo = smryMCMC( mcmcCoda , compVal=NULL , 
                        diffSVec=c(9,19, 29,39) , 
                        diffCVec=c(1,2,3,4) , 
                        compValDiff=0.0)

# Display posterior information:
plotMCMC( mcmcCoda , data=myData , 
          zName=zName, sName=sName, cName=cName, 
          compVal=NULL ,
          diffCList=list( c(1,2) ,
                          c(1,3) ,
                          c(1,4),
                          c(2,3),
                          c(2,4),
                          c(3,4)) ,
          diffSList=list( c(8,7),
                          c(8,9) ,
                          c(9,11)) ,
                          compValDiff=0.0) #ropeDiff = c(-0.05,0.05))


write.csv(summaryInfo, "summaryInfoHierarchicalAge_trial16.csv" )




# ============ Predictive check ============
confusionMatrix <- function(resp, pred){
  #if only one class is predicted, there is no accuracy at all
  if(dim(table(pred))==1)
  {  
    return(list(accuracy=0, AUC = 0))
  }  
  else
  {  
    
    classRes <- data.frame(response = resp , predicted = pred)
    conf = xtabs(~ predicted + response, data = classRes)
    
    accuracy = sum(diag(conf))/sum(conf)
    accuracy
    precision = conf[1,1]/(conf[1,1]+conf[1,2])
    precision
    recall = conf[1,1]/(conf[1,1]+conf[2,1])
    recall
    Fscore = 2*((precision*recall)/(precision+recall))
    Fscore
    tpr = conf[1,1]/(conf[1,1]+conf[2,1])
    tnr = conf[2,2]/(conf[2,2]+conf[1,2])
    AUC = (tpr + tnr) /2
    return(list(accuracy = accuracy, precision = precision, recall = recall, Fscore = Fscore, AUC=AUC,conf = conf ))
    
  }
}

# Searching threshold between 0.1 to 0.95 
thresholds <- seq(0.1, 0.95, 0.05)
cf <- array(NA,dim =c(length(thresholds),3))
for (i in 1:(length(thresholds))){
  predProbs <- summaryInfo[11:608,3]#summaryInfo[9:208,3]
  predcovidResults <- array(1, 598)#array(1, 200)
  preds <- data.frame(name = myData$Patient.ID, probPositive = predProbs, covidResultPred = predcovidResults, covidResultActual = myData$covid.results)
  preds$covidResultPred[which(preds$probPositive < thresholds[i])] <- 0 # Apply each threshold to estimate those not survived
  #write.csv(preds, paste("pred_", i, ".csv" ))
  print(paste("i=", i, "thresholds", thresholds[i], "dimension", dim(table(preds$covidResultPred))))
  confusionMatrix(resp = preds$covidResultActual, pred = preds$covidResultPred)
  cf[i,3] <- confusionMatrix(resp = preds$covidResultActual, pred = preds$covidResultPred)$accuracy
  cf[i,2] <- confusionMatrix(resp = preds$covidResultActual, pred = preds$covidResultPred)$AUC
  cf[i,1] <- thresholds[i]
}
cf
colnames(cf) <- c("Threshold", "AUC", "Accuracy")
cf

# Searching threshold between 0.01 to 0.2 
# Searching threshold between 0.01 to 0.15 - for Omega prior 

thresholds <- seq(0.01, 0.2, 0.01)
thresholds <- seq(0.01, 0.15, 0.005)

cf <- array(NA,dim =c(length(thresholds),3))
for (i in 1:(length(thresholds))){
  predProbs <- summaryInfo[11:608,3]#summaryInfo[9:208,3]
  predcovidResults <- array(1, 598)#array(1, 200)
  preds <- data.frame(name = myData$Patient.ID, probPositive = predProbs, covidResultPred = predcovidResults, covidResultActual = myData$covid.results)
  preds$covidResultPred[which(preds$probPositive < thresholds[i])] <- 0 # Apply each threshold to estimate those not survived
  #write.csv(preds, paste("pred_", i, ".csv" ))
  print(paste("i=", i, "thresholds", thresholds[i], "dimension", dim(table(preds$covidResultPred))))
  confusionMatrix(resp = preds$covidResultActual, pred = preds$covidResultPred)
  cf[i,3] <- confusionMatrix(resp = preds$covidResultActual, pred = preds$covidResultPred)$accuracy
  cf[i,2] <- confusionMatrix(resp = preds$covidResultActual, pred = preds$covidResultPred)$AUC
  cf[i,1] <- thresholds[i]
}
cf
colnames(cf) <- c("Threshold", "AUC", "Accuracy")
cf




predProbs <- summaryInfo[11:608,3]#summaryInfo[9:208,3]
predcovidResults <- array(1, 598)#array(1, 200)
preds <- data.frame(name = myData$Patient.ID, Age.group=myData$Age.group,   probPositive = predProbs, covidResultPred = predcovidResults, covidResultActual = myData$covid.results)
preds$covidResultPred[which(preds$probPositive < 0.18)] <- 0 # Apply threshold of 0.135 to estimate those not survived
confusionMatrix(resp = preds$covidResultActual, pred = preds$covidResultPred)

preds$label= 0
preds$label[which(preds$Age.group=="Age Group 1" & preds$covidResultActual == 0)] <- 1 
preds$label[which(preds$Age.group=="Age Group 1" & preds$covidResultActual == 1)] <- 2 
preds$label[which(preds$Age.group=="Age Group 2" & preds$covidResultActual == 0)] <- 3 
preds$label[which(preds$Age.group=="Age Group 2" & preds$covidResultActual == 1)] <- 4 
preds$label[which(preds$Age.group== "Age Group 3" & preds$covidResultActual == 0)] <- 5 
preds$label[which(preds$Age.group=="Age Group 3" & preds$covidResultActual == 1)] <- 6 
preds$label[which(preds$Age.group=="Age Group 4" & preds$covidResultActual == 0)] <- 7 
preds$label[which(preds$Age.group=="Age Group 4" & preds$covidResultActual == 1)] <- 8 


for(i in 1:8)
{
  minNum=min(preds[preds$label==i,]$probPositive)
  maxNum=max(preds[preds$label==i,]$probPositive)
  print(paste(minNum, maxNum))

}  

preds$label <- factor(preds$label,
                       levels = c(1,2,3,4,5,6,7,8),
                       labels = c("Age Group 1 -", "Age Group 1 +",
                                  "Age Group 2 -", "Age Group 2 +",
                                  "Age Group 3 -", "Age Group 3 +",
                                  "Age Group 4 -", "Age Group 4 +"))

write.csv(preds, paste("pred_final_trial16", ".csv" ))

a <- ggplot(preds, aes(x = probPositive))
a + geom_histogram(aes(color = label, fill = label),bins = 100,
                   alpha = 0.4, position = "identity") +
  geom_density(aes(y=..count../500,color= label), size = 1) +
  scale_fill_manual(values = c("#98c068", "#287028", "#f4d1d7", "#802939",
                               "#63bce5", "#114da9", "#fff897", "#f8e000")) +
  scale_color_manual(values = c("#98c068", "#287028", "#f4d1d7", "#802939",
                                "#63bce5", "#114da9", "#fff897", "#f8e000")) +
  ggtitle(label="Distribution of patients (mode of +ve probability) for different age groups", subtitle="(modeW=0.1, concentrationW=100, modeK=10, s.d.K=2)") +
  xlab("Positive Probability") + ylab("Count") +
  theme(plot.title = element_text(hjust = -0.3, size=16),plot.subtitle = element_text(hjust = 0.5, size=14))


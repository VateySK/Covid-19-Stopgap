
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

setwd("~/Master of Analytics/Semester 1 Year 2/MATH2269 Applied Bayesian Statistics")
source("DBDA2E-utilities.R") 

#================================================================================#
#                  PREPARING DATA FOR ANALYSIS
#================================================================================#

#================================DATA==================================

myData <- read.csv("covid_brazil_final.csv")
describe(myData)
x = as.matrix(myData[,c("Patient.age.quantile","Hematocrit", "Hemoglobin","Platelets", 
                        "Mean.platelet.volume", "Red.blood.Cells", "Lymphocytes",
                        "Mean.corpuscular.hemoglobin.concentrationA.MCHC.", "Leukocytes",
                        "Basophils", "Mean.corpuscular.hemoglobin.MCH.", "Eosinophils",
                        "Mean.corpuscular.volume.MCV.", "Monocytes", "Red.blood.cell.distribution.width.RDW.")])

#==========================Summary Statistics=====================================
# Some more descriptives
cat("\nCORRELATION MATRIX OF PREDICTORS:\n ")
show( round(cor(x),3) )
cat("\n")

x_correlation <- round(cor(x),3)
write.csv(x_correlation, "x_correlation.csv" )

#Pairplots for all variables
chart.Correlation(x, histogram=TRUE, pch=19)

#Remove Hematocrit, Red.blood.Cells and Mean.corpuscular.volume.MCV. due to high correlation

myData2 = myData[,c("Patient.ID", "Patient.age.quantile", "SARS.Cov2.exam.result", "Hemoglobin","Platelets", 
                    "Mean.platelet.volume", "Lymphocytes",
                    "Mean.corpuscular.hemoglobin.concentrationA.MCHC.", "Leukocytes",
                    "Basophils", "Mean.corpuscular.hemoglobin.MCH.", "Eosinophils",
                    "Monocytes", "Red.blood.cell.distribution.width.RDW.")]

names(myData2)[3] <- "covid.results"
describe(myData2)

myData2$covid.results <- as.numeric(as.factor(myData2$covid.results)) - 1 # To get 0/1 instead of 1/2; positive = 1; negative = 0

#check for missing values
sum(is.na(myData2))

#Descriptive look
p1 <- ggplot(myData2, aes(x=Patient.age.quantile, y = covid.results)) +
  geom_point()

p2 <- ggplot(myData2, aes(x=Hemoglobin, y = covid.results)) +
  geom_point()

p3 <- ggplot(myData2, aes(x=Platelets, y = covid.results)) +
  geom_point()

p4 <- ggplot(myData2, aes(x=Mean.platelet.volume, y = covid.results)) +
  geom_point()

p5 <- ggplot(myData2, aes(x=Lymphocytes, y = covid.results)) +
  geom_point()

p6 <- ggplot(myData2, aes(x=Mean.corpuscular.hemoglobin.concentrationA.MCHC., y = covid.results)) +
  geom_point()

p7 <- ggplot(myData2, aes(x=Leukocytes, y = covid.results)) +
  geom_point()

p8 <- ggplot(myData2, aes(x=Basophils, y = covid.results)) +
  geom_point()

p9 <- ggplot(myData2, aes(x=Mean.corpuscular.hemoglobin.MCH., y = covid.results)) +
  geom_point()

p10 <- ggplot(myData2, aes(x=Eosinophils, y = covid.results)) +
  geom_point()

p11 <- ggplot(myData2, aes(x=Monocytes, y = covid.results)) +
  geom_point()

p12 <- ggplot(myData2, aes(x=Red.blood.cell.distribution.width.RDW., y = covid.results)) +
  geom_point()


figure <- ggarrange(p1, p2, p3, p4, p5, p10, p7, p8, nrow = 4, ncol = 2)
figure <- annotate_figure(figure,
                          top = text_grob("Covid test results vs independent variables", face = "bold", size = 14))

figure

figure2 <- ggarrange(p9, p6, p11, p12, nrow = 4, ncol = 1)
figure2 <- annotate_figure(figure2,
                           top = text_grob("Covid test results vs independent variables", face = "bold", size = 14))

figure2

myData3 = myData[,c("Hemoglobin","Platelets", 
                    "Mean.platelet.volume", "Lymphocytes",
                    "Mean.corpuscular.hemoglobin.concentrationA.MCHC.", "Leukocytes",
                    "Basophils", "Mean.corpuscular.hemoglobin.MCH.", "Eosinophils",
                    "Monocytes", "Red.blood.cell.distribution.width.RDW.")] 
names(myData3)[5] <- "MCHC"
names(myData3)[8] <- "MCH"
names(myData3)[11] <- "Red Blood cell D.W"

boxplot(myData3, xaxt = "n")
text(x = 1:length(myData3), y = par("usr")[3] - 0.45, labels = names(myData3), xpd = NA, srt = 35, cex = 1.2)


#================================DATA==================================
set.seed(999)
trainData <- stratified(myData2, "covid.results", .7)
testData <- setdiff(myData2, trainData)


x = as.matrix(trainData[,c("Patient.age.quantile", "Hemoglobin","Platelets", 
                           "Mean.platelet.volume", "Lymphocytes",
                           "Mean.corpuscular.hemoglobin.concentrationA.MCHC.", "Leukocytes",
                           "Basophils", "Mean.corpuscular.hemoglobin.MCH.", "Eosinophils",
                           "Monocytes", "Red.blood.cell.distribution.width.RDW.")])


y = unlist(trainData[, "covid.results"])    

#================================================================================#
#                           LIGISTIC REGRESSION
#================================================================================#

#===============PRELIMINARY FUNCTIONS FOR POSTERIOR INFERENCES====================
genMCMC = function( x, y, numAdaptSteps=500, numBburnInSteps = 500,
                    numSavedSteps=500 ,  thinSteps=1 , saveName=NULL ,
                    runjagsMethod=runjagsMethodDefault , 
                    nChains=nChainsDefault, beta1Sens="4", beta2Sens="4", beta3Sens="4", beta4Sens="4", beta5Sens="4",
                    beta6Sens="4", beta7Sens="4", beta8Sens="4", beta9Sens="4", beta10Sens="4",
                    beta11Sens="4", beta12Sens="4"
) { 
  require(runjags)
  #-----------------------------------------------------------------------------
  # THE DATA.
  # Do some checking that data make sense:
  if ( any( !is.finite(y) ) ) { stop("All y values must be finite.") }
  if ( any( !is.finite(x) ) ) { stop("All x values must be finite.") }
  cat("\nCORRELATION MATRIX OF PREDICTORS:\n ")
  show( round(cor(x),3) )
  cat("\n")
  flush.console()
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    x = x ,
    y = y ,
    Nx = dim(x)[2] ,
    Ntotal = dim(x)[1]
  )
  #-----------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
#  Specify the model for standardized data:
  model {
    for ( i in 1:Ntotal ) {
      # In JAGS, ilogit is logistic:
      y[i] ~ dbern( mu[i] )
      mu[i] <- ( guess*(1/2) 
                 + (1.0-guess)*ilogit(beta0+sum(beta[1:Nx]*x[i,1:Nx])) )
    }    
    "
  priorString0=paste("beta0 ~ dnorm( 0 , 1/2^2 )", "\n")
  priorString1= paste("beta[1] ~ dnorm(0 , 1/", beta1Sens, ")", "\n")
  priorString2= paste("beta[2] ~ dnorm(0 , 1/", beta2Sens, ")", "\n")
  priorString3= paste("beta[3] ~ dnorm(0 , 1/", beta3Sens, ")", "\n")
  priorString4= paste("beta[4] ~ dnorm(0 , 1/", beta4Sens, ")", "\n")
  priorString5= paste("beta[5] ~ dnorm(0 , 1/", beta5Sens, ")", "\n")
  priorString6= paste("beta[6] ~ dnorm(0 , 1/", beta6Sens, ")", "\n")
  priorString7= paste("beta[7] ~ dnorm(0 , 1/", beta7Sens, ")", "\n")
  priorString8= paste("beta[8] ~ dnorm(0 , 1/", beta8Sens, ")", "\n")
  priorString9= paste("beta[9] ~ dnorm(0 , 1/", beta9Sens, ")", "\n")
  priorString10= paste("beta[10] ~ dnorm(0 , 1/", beta10Sens, ")", "\n")
  priorString11= paste("beta[11] ~ dnorm(0 , 1/", beta11Sens, ")", "\n")
  priorString12= paste("beta[12] ~ dnorm(0 , 1/", beta12Sens, ")", "\n")
  
  priorString = paste(priorString0, priorString1,priorString2, priorString3, priorString4, priorString5, priorString6,
                      priorString7, priorString8,priorString9,priorString10,priorString11, priorString12) 
  
  # Priors vague on standardized scale:
  guessString="
  
    guess ~ dbeta(1,9)
    # Transform to original scale:
    # beta[1:Nx] <- zbeta[1:Nx] / xsd[1:Nx] 
    # beta0 <- zbeta0 - sum( zbeta[1:Nx] * xm[1:Nx] / xsd[1:Nx] )
  }
  " # close quote for modelString
  # Write out modelString to a text file
  finalString = paste(modelString, priorString, guessString)
  writeLines( finalString , con="TEMPmodel.txt" )
  #-----------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Let JAGS do it...
  
  #-----------------------------------------------------------------------------
  # RUN THE CHAINS
  parameters = c( "beta0" ,  "beta" ,  
                  "guess" )
  runJagsOut <- run.jags( method=runjagsMethod ,
                          model="TEMPmodel.txt" , 
                          monitor=parameters , 
                          data=dataList ,  
                          #inits=initsList , 
                          n.chains=nChains ,
                          adapt=numAdaptSteps ,
                          burnin=numBburnInSteps , 
                          sample=ceiling(numSavedSteps/nChains) ,
                          thin=thinSteps ,
                          summarise=FALSE ,
                          plots=FALSE )
  codaSamples = as.mcmc.list( runJagsOut )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  if ( !is.null(saveName) ) {
    save( codaSamples , file=paste(saveName,"Mcmc.Rdata",sep="") )
  }
  return( codaSamples )
} # end function

#======================================================================================

smryMCMC_HD = function(  codaSamples , compVal = NULL,  saveName=NULL) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  paramName = colnames(mcmcMat)
  for ( pName in paramName ) {
    if (pName %in% colnames(compVal)){
      if (!is.na(compVal[pName])) {
        summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] , 
                                                          compVal = as.numeric(compVal[pName]) ))
      }
      else {
        summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] ) )
      }
    } else {
      summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] ) )
    }
  }
  rownames(summaryInfo) = paramName
  
  # summaryInfo = rbind( summaryInfo , 
  #                      "tau" = summarizePost( mcmcMat[,"tau"] ) )
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================

plotMCMC_HD = function( codaSamples , data , xName="x" , yName="y", preds = FALSE ,
                        showCurve=FALSE ,  pairsPlot=FALSE , compVal = NULL, 
                        saveName=NULL , saveType="jpg" ) {
  # showCurve is TRUE or FALSE and indicates whether the posterior should
  #   be displayed as a histogram (by default) or by an approximate curve.
  # pairsPlot is TRUE or FALSE and indicates whether scatterplots of pairs
  #   of parameters should be displayed.
  #-----------------------------------------------------------------------------
  y = data[,yName]
  x = as.matrix(data[,xName])
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  # zbeta0 = mcmcMat[,"zbeta0"]
  # zbeta  = mcmcMat[,grep("^zbeta$|^zbeta\\[",colnames(mcmcMat))]
  # if ( ncol(x)==1 ) { zbeta = matrix( zbeta , ncol=1 ) }
  beta0 = mcmcMat[,"beta0"]
  beta  = mcmcMat[,grep("^beta$|^beta\\[",colnames(mcmcMat))]
  if ( ncol(x)==1 ) { beta = matrix( beta , ncol=1 ) }
  if (preds){
    pred = mcmcMat[,grep("^pred$|^pred\\[",colnames(mcmcMat))]
  } # Added by Demirhan
  guess = mcmcMat[,"guess"]
  #-----------------------------------------------------------------------------
  # Compute R^2 for credible parameters:
  YcorX = cor( y , x ) # correlation of y with each x predictor
  Rsq = beta %*% matrix( YcorX , ncol=1 )
  #-----------------------------------------------------------------------------
  if ( pairsPlot ) {
    # Plot the parameters pairwise, to see correlations:
    openGraph()
    nPtToPlot = 1000
    plotIdx = floor(seq(1,chainLength,by=chainLength/nPtToPlot))
    panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
      usr = par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r = (cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex=1.25 ) # was cex=cex.cor*r
    }
    pairs( cbind( beta0 , beta )[plotIdx,] ,
           labels=c( "beta[0]" , 
                     paste0("beta[",1:ncol(beta),"]\n",xName) , 
                     expression(tau) ) , 
           lower.panel=panel.cor , col="skyblue" )
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"PostPairs",sep=""), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
  # Marginal histograms:
  
  decideOpenGraph = function( panelCount , saveName , finished=FALSE , 
                              nRow=2 , nCol=3 ) {
    # If finishing a set:
    if ( finished==TRUE ) {
      if ( !is.null(saveName) ) {
        saveGraph( file=paste0(saveName,ceiling((panelCount-1)/(nRow*nCol))), 
                   type=saveType)
      }
      panelCount = 1 # re-set panelCount
      return(panelCount)
    } else {
      # If this is first panel of a graph:
      if ( ( panelCount %% (nRow*nCol) ) == 1 ) {
        # If previous graph was open, save previous one:
        if ( panelCount>1 & !is.null(saveName) ) {
          saveGraph( file=paste0(saveName,(panelCount%/%(nRow*nCol))), 
                     type=saveType)
        }
        # Open new graph
        openGraph(width=nCol*7.0/3,height=nRow*2.0)
        layout( matrix( 1:(nRow*nCol) , nrow=nRow, byrow=TRUE ) )
        par( mar=c(4,4,2.5,0.5) , mgp=c(2.5,0.7,0) )
      }
      # Increment and return panel count:
      panelCount = panelCount+1
      return(panelCount)
    }
  }
  
  # Original scale:
  panelCount = 1
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( beta0 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(beta[0]) , main="Intercept", compVal = as.numeric(compVal["beta0"] ))
  for ( bIdx in 1:ncol(beta) ) {
    panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
    if (!is.na(compVal[paste0("beta[",bIdx,"]")])){
      histInfo = plotPost( beta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                           xlab=bquote(beta[.(bIdx)]) , main=xName[bIdx],
                           compVal = as.numeric(compVal[paste0("beta[",bIdx,"]")]))
    } else{
      histInfo = plotPost( beta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                           xlab=bquote(beta[.(bIdx)]) , main=xName[bIdx])
    }
  }
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( Rsq , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(R^2) , main=paste("Prop Var Accntd") )
  
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( guess , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab="Guess parameter" , main=paste("Prop Var Accntd"))
  
  panelCount = 1
  if ( preds){
    
    for ( pIdx in 1:ncol(pred) ) {
      panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
      histInfo = plotPost( pred[,pIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                           xlab=bquote(pred[.(pIdx)]) , main=paste0("Prediction ",pIdx) )
    }
  }
  # Standardized scale:
  panelCount = 1
  # panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
  # histInfo = plotPost( zbeta0 , cex.lab = 1.75 , showCurve=showCurve ,
  #                      xlab=bquote(z*beta[0]) , main="Intercept" )
  # for ( bIdx in 1:ncol(beta) ) {
  #   panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
  #   histInfo = plotPost( zbeta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
  #                        xlab=bquote(z*beta[.(bIdx)]) , main=xName[bIdx] )
  # }
  # panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
  # histInfo = plotPost( Rsq , cex.lab = 1.75 , showCurve=showCurve ,
  #                      xlab=bquote(R^2) , main=paste("Prop Var Accntd") )
  # panelCount = decideOpenGraph( panelCount , finished=TRUE , saveName=paste0(saveName,"PostMargZ") )
  
  #-----------------------------------------------------------------------------
}

#===============PRELIMINARY FUNCTIONS FOR POSTERIOR INFERENCES====================

numAdaptSteps = 1000 ; numBburnInSteps=1000; numSavedSteps=3000 ; thinSteps=40; nChains = 3

startTime = proc.time()
mcmcCoda = genMCMC( x, y, numAdaptSteps=numAdaptSteps, numBburnInSteps= numBburnInSteps,
                    numSavedSteps=numSavedSteps , thinSteps=thinSteps, beta12Sens=0.00000001,
                    nChains = nChains )



stopTime = proc.time()
duration = stopTime - startTime
show(duration)

#save.image(file="rEnvironment_40Thin_3000_logistic.RData")
# load(file="rEnvironment_40thin.RData")

#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = c("beta0"  ,  "beta[1]",  "beta[2]" ,  "beta[3]" ,  "beta[4]" ,  "beta[5]"  , "beta[6]"  ,
                   "beta[7]"  , "beta[8]", "beta[9]", "beta[10]", "beta[11]", "beta[12]", "guess") #varnames(mcmcCoda) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName) 
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:

compVal <- data.frame("beta0" = 0, "beta[1]" = 0, "beta[2]" = 0, "beta[3]" = 0, "beta[4]" =  0,  "beta[5]" =  0, 
                      "beta[6]" =  0, "beta[7]" = 0, "beta[8]" = 0, "beta[9]" = 0,"beta[10]" = 0,
                      "beta[11]" = 0, "beta[12]" = 0,check.names=FALSE)

summaryInfo <- smryMCMC_HD( codaSamples = mcmcCoda , compVal = compVal )
print(summaryInfo)

plotMCMC_HD( codaSamples = mcmcCoda , data = myData2, xName=c("Patient.age.quantile", "Hemoglobin","Platelets", 
                                                              "Mean.platelet.volume", "Lymphocytes",
                                                              "Mean.corpuscular.hemoglobin.concentrationA.MCHC.", "Leukocytes",
                                                              "Basophils", "Mean.corpuscular.hemoglobin.MCH.", "Eosinophils",
                                                              "Monocytes", "Red.blood.cell.distribution.width.RDW."),
             yName="covid.results", compVal = compVal, preds= FALSE, pairsPlot=TRUE )

write.csv(summaryInfo, "summaryInfoLogistic_b12_try2.csv" )

# ============ Predictive check ============

coeffs <- as.vector(summaryInfo[2:14,3])
X <- as.matrix(cbind(rep(1,nrow(testData)),testData[, c(2,4:14)]))
#X <- as.matrix(cbind(rep(1,nrow(trainData)),trainData[, c(1,3:13)]))
predProbs <- 1/(1+exp(-(X%*%coeffs)))


trueClass <- unlist(testData[,3])
#trueClass <- unlist(trainData[,3])

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

thresholds <- seq(0.05, 0.95, 0.05)
cf <- array(NA,dim =c(length(thresholds),3))
for (i in 1:(length(thresholds))){
  predClass <- as.numeric(predProbs>thresholds[i])
  cf[i,3] <- confusionMatrix(resp = trueClass, pred = predClass)$accuracy
  cf[i,2] <- confusionMatrix(resp = trueClass, pred = predClass)$AUC
  cf[i,1] <- thresholds[i]
}

colnames(cf) <- c("Threshold", "AUC", "Accuracy")
cf
write.csv(cf, "cf_b12_try2.csv" )

# Best performance using AUC
threshold <- 0.2
predClass <- as.numeric(predProbs>threshold)
confusionMatrix(resp = trueClass, pred = predClass)
preds <- data.frame(name = testData$Patient.ID, probPositive = predProbs, resp = trueClass, pred = predClass)
write.csv(preds, "preds_logistic_AUC_b12_try2.csv" )

# Best performance using Accuracy
threshold <- 0.5
predClass <- as.numeric(predProbs>threshold)
confusionMatrix(resp = trueClass, pred = predClass)
preds <- data.frame(name = testData$Patient.ID, probPositive = predProbs, resp = trueClass, pred = predClass)
write.csv(preds, "preds_logistic_Accuracy_b12_try2.csv" )


a <- ggplot(preds, aes(x = probPositive))
a + geom_histogram(aes(color = as.factor(pred), fill = as.factor(pred)),bins = 100,
                   alpha = 0.4, position = "identity")


#================================================================================#
#                       Model Comparison of 4 logistic regression models
#================================================================================#

#===============PRELIMINARY FUNCTIONS FOR POSTERIOR INFERENCES====================
fileNameRoot="Mcom-" # for output filenames
smryMCMC_HD = function(  codaSamples , compVal = NULL,  saveName=NULL) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  paramName = colnames(mcmcMat)
  for ( pName in paramName ) {
    if (pName %in% colnames(compVal)){
      if (!is.na(compVal[pName])) {
        summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] , 
                                                          compVal = as.numeric(compVal[pName]) ))
      }
      else {
        summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] ) )
      }
    } else {
      summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] ) )
    }
  }
  rownames(summaryInfo) = paramName
  
  # summaryInfo = rbind( summaryInfo , 
  #                      "tau" = summarizePost( mcmcMat[,"tau"] ) )
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================
plotMCMC_HD = function( codaSamples , data , xName="x" , yName="y" ,
                        showCurve=FALSE ,  pairsPlot=FALSE , compVal = NULL,
                        saveName=NULL , saveType="jpg" ) {
  # showCurve is TRUE or FALSE and indicates whether the posterior should
  #   be displayed as a histogram (by default) or by an approximate curve.
  # pairsPlot is TRUE or FALSE and indicates whether scatterplots of pairs
  #   of parameters should be displayed.
  #-----------------------------------------------------------------------------
  y = data[,yName]
  x = as.matrix(data[,xName])
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  beta0 = mcmcMat[,"beta0"]
  beta  = mcmcMat[,grep("^beta$|^beta\\[",colnames(mcmcMat))]
  if ( ncol(x)==1 ) { beta = matrix( beta , ncol=1 ) }
  #-----------------------------------------------------------------------------
  # Compute R^2 for credible parameters:
  YcorX = cor( y , x ) # correlation of y with each x predictor
  Rsq = beta %*% matrix( YcorX , ncol=1 )
  #-----------------------------------------------------------------------------
  # Marginal histograms:
  
  decideOpenGraph = function( panelCount , saveName , finished=FALSE , 
                              nRow=2 , nCol=3 ) {
    # If finishing a set:
    if ( finished==TRUE ) {
      if ( !is.null(saveName) ) {
        saveGraph( file=paste0(saveName,ceiling((panelCount-1)/(nRow*nCol))), 
                   type=saveType)
      }
      panelCount = 1 # re-set panelCount
      return(panelCount)
    } else {
      # If this is first panel of a graph:
      if ( ( panelCount %% (nRow*nCol) ) == 1 ) {
        # If previous graph was open, save previous one:
        if ( panelCount>1 & !is.null(saveName) ) {
          saveGraph( file=paste0(saveName,(panelCount%/%(nRow*nCol))), 
                     type=saveType)
        }
        # Open new graph
        openGraph(width=nCol*7.0/3,height=nRow*2.0)
        layout( matrix( 1:(nRow*nCol) , nrow=nRow, byrow=TRUE ) )
        par( mar=c(4,4,2.5,0.5) , mgp=c(2.5,0.7,0) )
      }
      # Increment and return panel count:
      panelCount = panelCount+1
      return(panelCount)
    }
  }
  
  # Original scale:
  panelCount = 1
  if (!is.na(compVal["beta0"])){
    panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
    histInfo = plotPost( beta0 , cex.lab = 1.75 , showCurve=showCurve ,
                         xlab=bquote(beta[0]) , main="Intercept", compVal = as.numeric(compVal["beta0"] ))
  } else {  
    histInfo = plotPost( beta0 , cex.lab = 1.75 , showCurve=showCurve ,
                         xlab=bquote(beta[0]) , main="Intercept")
  }
  for ( bIdx in 1:ncol(beta) ) {
    panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
    if (!is.na(compVal[paste0("beta[",bIdx,"]")])) {
      histInfo = plotPost( beta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                           xlab=bquote(beta[.(bIdx)]) , main=xName[bIdx],
                           compVal = as.numeric(compVal[paste0("beta[",bIdx,"]")]))
    } else{
      histInfo = plotPost( beta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                           xlab=bquote(beta[.(bIdx)]) , main=xName[bIdx])
    }
  }
  panelCount = decideOpenGraph( panelCount , finished=TRUE , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( Rsq , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(R^2) , main=paste("Prop Var Accntd") )
  
  
}

#===============================================================================

# Specify the data in a list, for later shipment to JAGS:
dataList <- list(
  x = x ,
  y = y ,
  Nx = dim(x)[2] ,
  Ntotal = dim(x)[1]
)
#------------------LOGISTIC REGRESSION 4 MODELS -------------------------
modelString = "
# data {
#   for ( j in 1:Nx ) {
#     xm[j]  <- mean(x[,j])
#     xsd[j] <-   sd(x[,j])
#     for ( i in 1:Ntotal ) {
#       zx[i,j] <- ( x[i,j] - xm[j] ) / xsd[j]
#     }
#   }
# }
model {
  for ( i in 1:Ntotal ) {
    # In JAGS, ilogit is logistic:
    y[i] ~ dbern( mu[i] )
      mu[i] <- ifelse(m == 1, (guess*(1/2) + (1.0-guess)*ilogit(beta0+beta[1]*x[i,1]+beta[2]*x[i,2]+beta[3]*x[i,3]+beta[4]*x[i,4]+beta[5]*x[i,5]+beta[6]*x[i,6]+beta[7]*x[i,7]+beta[8]*x[i,8]+beta[9]*x[i,9]+beta[10]*x[i,10]+beta[11]*x[i,11]+beta[12]*x[i,12])),
                              ifelse( m == 2, (guess*(1/2) + (1.0-guess)*ilogit(beta0+beta[1]*x[i,1]+beta[2]*x[i,2]+beta[3]*x[i,3]+beta[5]*x[i,5]+beta[6]*x[i,6]+beta[7]*x[i,7]+beta[8]*x[i,8]+beta[9]*x[i,9]+beta[10]*x[i,10]+beta[12]*x[i,12])),
                              ifelse( m == 3, (guess*(1/2) + (1.0-guess)*ilogit(beta0+beta[1]*x[i,1]+beta[2]*x[i,2]+beta[3]*x[i,3]+beta[7]*x[i,7]+beta[8]*x[i,8]+beta[9]*x[i,9]+beta[10]*x[i,10]+beta[12]*x[i,12])),
                                              (guess*(1/2) + (1.0-guess)*ilogit(beta0+beta[1]*x[i,1]+beta[2]*x[i,2]+beta[7]*x[i,7]+beta[10]*x[i,10])))))
  }
  # Priors vague on standardized scale:
  beta0 ~ dnorm( 0 , 1/2^2 )
  # non-informative run
  for ( j in 1:Nx ) {
    beta[j] ~ dnorm( 0 , 1/2^2 )
  }
  guess ~ dbeta(1,9)
  # Prior for model
  m ~ dcat( mPriorProb[] )
  mPriorProb[1] <- .25
  mPriorProb[2] <- .25
  mPriorProb[3] <- .25
  mPriorProb[4] <- .25
  # Transform to original scale:
  # beta[1:Nx] <- zbeta[1:Nx] / xsd[1:Nx]
  # beta0 <- zbeta0 - sum( zbeta[1:Nx] * xm[1:Nx] / xsd[1:Nx] )
  # Predictions
  # for ( k in 1:Npred){
  #   pred[k] <- ilogit(beta0 + sum(beta[1:Nx] * xPred[k,1:Nx]))
  # }
}
"
writeLines( modelString , con="TEMPmodel.txt" )

parameters = c( "beta0" ,  "beta", "m") # Here beta is a vector!

adaptSteps = 1000  # Number of steps to "tune" the samplers
burnInSteps = 1000
nChains = 3 
thinSteps = 42 #40 
numSavedSteps = 3000
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

# Parallel run
startTime = proc.time()
runJagsOut <- run.jags( method="parallel" ,
                        model="TEMPmodel.txt" ,
                        monitor=parameters  ,
                        data=dataList ,
                        #inits=initsList ,
                        n.chains=nChains ,
                        adapt=adaptSteps ,
                        burnin=burnInSteps ,
                        sample=numSavedSteps ,
                        thin=thinSteps , summarise=FALSE , plots=FALSE )
codaSamples = as.mcmc.list( runJagsOut )
stopTime = proc.time()
elapsedTime = stopTime - startTime
show(elapsedTime)

# save.image("REvironmentLojReg4Models42thins.RData")
# load("REvironmentLojReg4Models42thins.RData")

parameterNames = varnames(codaSamples) # get all parameter names
for ( parName in parameterNames  ) {
  diagMCMC( codaSamples , parName=parName ,
            saveName=fileNameRoot , saveType="jpg" )
}


summaryInfo <- smryMCMC_HD( codaSamples = codaSamples  )
print(summaryInfo)


# Convert coda-object codaSamples to matrix object for easier handling.
mcmcMat = as.matrix( codaSamples , chains=TRUE )
m = mcmcMat[,"m"]
beta0 = mcmcMat[,"beta0"]
beta1 = mcmcMat[,"beta[1]"]
beta2 = mcmcMat[,"beta[2]"]
beta3 = mcmcMat[,"beta[3]"]
beta4 = mcmcMat[,"beta[4]"]
beta5 = mcmcMat[,"beta[5]"]
beta6 = mcmcMat[,"beta[6]"]
beta7 = mcmcMat[,"beta[7]"]
beta8 = mcmcMat[,"beta[8]"]
beta9 = mcmcMat[,"beta[9]"]
beta10 = mcmcMat[,"beta[10]"]
beta11 = mcmcMat[,"beta[11]"]
beta12 = mcmcMat[,"beta[12]"]

# Compute the proportion of m at each index value:
pM1 = sum( m == 1 ) / length( m )
pM2 = sum( m == 2 ) / length( m )
pM3 = sum( m == 3 ) / length( m )
pM4 = 1 - (pM1 + pM2 + pM3)

priorM1 = 0.25
priorM2 = 0.25
priorM3 = 0.25
priorM4 = 0.25

# BF = (prior odds) / (posterior odds)
BF_1vs2 = (pM1/pM2) / (priorM1/priorM2) # H0: beta = beta_{m1}; H1: beta = beta_{m2}
BF_1vs3 = (pM1/pM3) / (priorM1/priorM3) # H0: beta = beta_{m1}; H1: beta = beta_{m3}
BF_1vs4 = (pM1/pM4) / (priorM1/priorM4) # H0: beta = beta_{m1}; H1: beta = beta_{m4}
BF_2vs3 = (pM2/pM3) / (priorM2/priorM3) # H0: beta = beta_{m2}; H1: beta = beta_{m3}
BF_2vs4 = (pM2/pM4) / (priorM2/priorM4) # H0: beta = beta_{m2}; H1: beta = beta_{m4}
BF_3vs4 = (pM3/pM4) / (priorM3/priorM4) # H0: beta = beta_{m3}; H1: beta = beta_{m4}

# Extract theta values for each model index:

beta1M3 = beta1[ m == 3 ]
beta2M3 = beta2[ m == 3 ]
beta3M3 = beta3[ m == 3 ]
beta4M3 = beta4[ m == 3 ]
beta5M3 = beta5[ m == 3 ]
beta6M3 = beta6[ m == 3 ]
beta7M3 = beta7[ m == 3 ]
beta8M3 = beta8[ m == 3 ]
beta9M3 = beta9[ m == 3 ]
beta10M3 = beta10[ m == 3 ]
beta11M3 = beta11[ m == 3 ]
beta12M3 = beta12[ m == 3 ]

beta1M4 = beta1[ m == 4 ]
beta2M4 = beta2[ m == 4 ]
beta3M4 = beta3[ m == 4 ]
beta4M4 = beta4[ m == 4 ]
beta5M4 = beta5[ m == 4 ]
beta6M4 = beta6[ m == 4 ]
beta7M4 = beta7[ m == 4 ]
beta8M4 = beta8[ m == 4 ]
beta9M4 = beta9[ m == 4 ]
beta10M4 = beta10[ m == 4 ]
beta11M4 = beta11[ m == 4 ]
beta12M4 = beta12[ m == 4 ]


openGraph(width=5,height=5)
plotPost( m , breaks=seq(0.9,43.1,0.2) , cenTend="mean" , xlab="m" , main="Model Index" )

openGraph(width=7,height=5)
par( mar=0.5+c(3,1,2,1) , mgp=c(2.0,0.7,0) )
layout( matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=FALSE) , widths=c(2,2) )
plotPost( beta1M3 , 
          main=bquote( beta1*" when m=3" * " ; p(m=3|D)" == .(signif(pM3,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[1])  )
plotPost( beta2M3 , 
          main=bquote( beta2*" when m=3" * " ; p(m=3|D)" == .(signif(pM3,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[2])  )
plotPost( beta3M3 , 
          main=bquote( beta3*" when m=3" * " ; p(m=3|D)" == .(signif(pM3,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[3])  )
plotPost( beta4M3 , 
          main=bquote( beta4*" when m=3" * " ; p(m=3|D)" == .(signif(pM3,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[4])  )


openGraph(width=7,height=5)
par( mar=0.5+c(3,1,2,1) , mgp=c(2.0,0.7,0) )
layout( matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=FALSE) , widths=c(2,2) )
plotPost( beta5M3 , 
          main=bquote( beta5*" when m=3" * " ; p(m=3|D)" == .(signif(pM3,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[5])  )
plotPost( beta6M3 , 
          main=bquote( beta6*" when m=3" * " ; p(m=3|D)" == .(signif(pM3,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[6])  )
plotPost( beta7M3 , 
          main=bquote( beta7*" when m=3" * " ; p(m=3|D)" == .(signif(pM3,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[7])  )
plotPost( beta8M3 , 
          main=bquote( beta8*" when m=3" * " ; p(m=3|D)" == .(signif(pM3,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[8])  )


openGraph(width=7,height=5)
par( mar=0.5+c(3,1,2,1) , mgp=c(2.0,0.7,0) )
layout( matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=FALSE) , widths=c(2,2) )
plotPost( beta9M3 , 
          main=bquote( beta9*" when m=3" * " ; p(m=3|D)" == .(signif(pM3,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[9])  )
plotPost( beta10M3 , 
          main=bquote( beta10*" when m=3" * " ; p(m=3|D)" == .(signif(pM3,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[10])  )
plotPost( beta11M3 , 
          main=bquote( beta11*" when m=3" * " ; p(m=3|D)" == .(signif(pM3,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[11])  )
plotPost( beta12M3 , 
          main=bquote( beta12*" when m=3" * " ; p(m=3|D)" == .(signif(pM3,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[12])  )


openGraph(width=5,height=5)
plotPost( m , breaks=seq(0.9,43.1,0.2) , cenTend="mean" , xlab="m" , main="Model Index" )

openGraph(width=7,height=5)
par( mar=0.5+c(3,1,2,1) , mgp=c(2.0,0.7,0) )
layout( matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=FALSE) , widths=c(2,2) )
plotPost( beta1M4 , 
          main=bquote( beta1*" when m=4" * " ; p(m=4|D)" == .(signif(pM4,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[1])  )
plotPost( beta2M4 , 
          main=bquote( beta2*" when m=4" * " ; p(m=4|D)" == .(signif(pM4,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[2])  )
plotPost( beta3M4 , 
          main=bquote( beta3*" when m=4" * " ; p(m=4|D)" == .(signif(pM4,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[3])  )
plotPost( beta4M4 , 
          main=bquote( beta4*" when m=4" * " ; p(m=4|D)" == .(signif(pM4,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[4])  )

openGraph(width=7,height=5)
par( mar=0.5+c(3,1,2,1) , mgp=c(2.0,0.7,0) )
layout( matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=FALSE) , widths=c(2,2) )
plotPost( beta5M4 , 
          main=bquote( beta5*" when m=4" * " ; p(m=4|D)" == .(signif(pM4,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[5])  )
plotPost( beta6M4 , 
          main=bquote( beta6*" when m=4" * " ; p(m=4|D)" == .(signif(pM4,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[6])  )
plotPost( beta7M4 , 
          main=bquote( beta7*" when m=4" * " ; p(m=4|D)" == .(signif(pM4,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[7])  )
plotPost( beta8M4 , 
          main=bquote( beta8*" when m=4" * " ; p(m=4|D)" == .(signif(pM4,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[8])  )

openGraph(width=7,height=5)
par( mar=0.5+c(3,1,2,1) , mgp=c(2.0,0.7,0) )
layout( matrix(c(1,2,3,4),nrow=2,ncol=2,byrow=FALSE) , widths=c(2,2) )
plotPost( beta9M4 , 
          main=bquote( beta9*" when m=4" * " ; p(m=4|D)" == .(signif(pM4,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[9])  )
plotPost( beta10M4 , 
          main=bquote( beta10*" when m=4" * " ; p(m=4|D)" == .(signif(pM4,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[10])  )
plotPost( beta11M4 , 
          main=bquote( beta11*" when m=4" * " ; p(m=4|D)" == .(signif(pM4,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[11])  )
plotPost( beta12M4 , 
          main=bquote( beta12*" when m=4" * " ; p(m=4|D)" == .(signif(pM4,3)) ) , 
          cex.main=1.75 , xlab=bquote(beta[12])  )

# ============ Predictive check ============
coeffs <- as.vector(summaryInfo[2:14,3])

X <- as.matrix(cbind(rep(1,nrow(testData)),testData[, c(1,3:13)]))

predProbs <- 1/(1+exp(-(X%*%coeffs)))


trueClass <- unlist(testData[,2])
#trueClass <- unlist(trainData[,3])

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

thresholds <- seq(0.1, 0.95, 0.05)
cf <- array(NA,dim =c(length(thresholds),3))
for (i in 1:(length(thresholds))){
  predClass <- as.numeric(predProbs>thresholds[i])
  cf[i,3] <- confusionMatrix(resp = trueClass, pred = predClass)$accuracy
  cf[i,2] <- confusionMatrix(resp = trueClass, pred = predClass)$AUC
  cf[i,1] <- thresholds[i]
}

colnames(cf) <- c("Threshold", "AUC", "Accuracy")
cf
write.csv(cf, "mcf.csv" )



# Best performance using AUC
threshold <- 0.45
predClass <- as.numeric(predProbs>=threshold)
confusionMatrix(resp = trueClass, pred = predClass)


# Best performance using Accuracy
threshold <- 0.65
predClass <- as.numeric(predProbs>=threshold)
confusionMatrix(resp = trueClass, pred = predClass)
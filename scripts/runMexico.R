#######################
### Mexico Run File ###
##### 12/31/2015 ######

## Read in libraries
library(ggplot2)
library(reshape2)
library(gridExtra)
library(pander)
library(tidyverse)

## path start -- where to save outputs
pathstart <- "~/Dropbox/Tracey's EcoA Work/2018/Mexico model/outputs/"

## Read in data input file
dataInput = read.csv("inputs/dataInputMexico_1028.csv", header=TRUE,stringsAsFactors=FALSE)


# ## Read in parameter text descriptions
# parameterText = read.csv("parameterTextMexico.csv", header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)

## Read in functions file
source("scripts/functionsMexico0114.R")



################################################
#### Set up model run
## Enter management scenarios to loop over - select 1 through 6
scenarios = c(1,2,3,4)
## Enter number of legal harvest thetas to loop over - must be greater than 1
thetas = 1
## Enter whether or not to loop over catch share cost and price scalars - put "yes" or "no"
catchShareLoop = "yes"
## Enter whether or not to loop over eliminating illegal fishing - put "yes" or "no"
illegalLoop = "yes"
################################################

################################################
### Implementation vector
delayVec <- c(2:20)

summaryData = matrix(NA,nrow=nrow(dataInput),ncol=35)
colnames(summaryData)=c("Species",
                        "B_current",
                        "B_year20_SQ",
                        "B_year20_CSL1",
                        "B_year20_CSL2",
                        "B_year20_BPM",
                        "BDelta_absolute_SQ_BPM",
                        "BDelta_percent_loss_SQ_BPM",
                        "BDelta_absolute_CSL2_BPM",
                        "BDelta_percent_loss_CSL2_BPM",
                        "H_current",
                        "H_year20_SQ",
                        "H_year20_CSL1",
                        "H_year20_CSL2",
                        "H_year20_BPM",
                        "HDelta_absolute_SQ_BPM",
                        "HDelta_percent_loss_SQ_BPM",
                        "HDelta_absolute_CSL2_BPM",
                        "HDelta_percent_loss_CSL2_BPM",
                        "P_current",
                        "P_year20_SQ",
                        "P_year20_CSL1",
                        "P_year20_CSL2",
                        "P_year20_BPM",
                        "PDelta_absolute_SQ_BPM",
                        "PDelta_percent_loss_SQ_BPM",
                        "PDelta_absolute_CSL2_BPM",
                        "PDelta_percent_loss_CSL2_BPM",
                        "Gain_SQ_BPM",
                        "Gain_CSL1_BPM",
                        "recTime_SQ",
                        "recTime_CSL1",
                        "recTime_CSL2",
                        "recTime_BPM",
                        "r")


## Loop over all fisheries
for (i in 1:2)
# for (i in 1:nrow(dataInput))
{
    outputs = projectionModel(dataInput[i,],scenarios,thetas,catchShareLoop,illegalLoop)
    summaryData[i,] = c(dataInput[i,]$Species,
                        outputs$BProjections[1,1,1,1,1,1,1],
                        outputs$BProjections[1,1,1,1,1,1,20],
                        outputs$BProjections[1,1,1,2,1,1,20],
                        outputs$BProjections[2,1,1,2,1,1,20],
                        outputs$BProjections[2,1,1,2,2,1,20],
                        outputs$BProjections[2,1,1,2,2,1,20] - outputs$BProjections[1,1,1,1,1,1,20],
                        (outputs$BProjections[2,1,1,2,2,1,20] - outputs$BProjections[1,1,1,1,1,1,20])/outputs$BProjections[2,1,1,2,2,1,20]*100,
                        outputs$BProjections[2,1,1,2,2,1,20] - outputs$BProjections[2,1,1,2,1,1,20],
                        (outputs$BProjections[2,1,1,2,2,1,20] - outputs$BProjections[2,1,1,2,1,1,20])/outputs$BProjections[2,1,1,2,2,1,20]*100,
                        outputs$HIntProjections[1,1,1,1,1,1,1],
                        outputs$HIntProjections[1,1,1,1,1,1,20],
                        outputs$HIntProjections[1,1,1,2,1,1,20],
                        outputs$HIntProjections[2,1,1,2,1,1,20],
                        outputs$HIntProjections[2,1,1,2,2,1,20],
                        outputs$HIntProjections[2,1,1,2,2,1,20] - outputs$HIntProjections[1,1,1,1,1,1,20],
                        (outputs$HIntProjections[2,1,1,2,2,1,20] - outputs$HIntProjections[1,1,1,1,1,1,20])/outputs$HIntProjections[2,1,1,2,2,1,20]*100,
                        outputs$HIntProjections[2,1,1,2,2,1,20] - outputs$HIntProjections[2,1,1,2,1,1,20],
                        (outputs$HIntProjections[2,1,1,2,2,1,20] - outputs$HIntProjections[2,1,1,2,1,1,20])/outputs$HIntProjections[2,1,1,2,2,1,20]*100,
                        outputs$profitProjections[1,1,1,1,1,1,1],
                        outputs$profitProjections[1,1,1,1,1,1,20],
                        outputs$profitProjections[1,1,1,2,1,1,20],
                        outputs$profitProjections[2,1,1,2,1,1,20],
                        outputs$profitProjections[2,1,1,2,2,1,20],
                        outputs$profitProjections[2,1,1,2,2,1,20] - outputs$profitProjections[1,1,1,1,1,1,20],
                        (outputs$profitProjections[2,1,1,2,2,1,20] - outputs$profitProjections[1,1,1,1,1,1,20])/outputs$profitProjections[2,1,1,2,2,1,20]*100,
                        outputs$profitProjections[2,1,1,2,2,1,20] - outputs$profitProjections[2,1,1,2,1,1,20],
                        (outputs$profitProjections[2,1,1,2,2,1,20] - outputs$profitProjections[2,1,1,2,1,1,20])/outputs$profitProjections[2,1,1,2,2,1,20]*100,
                        (outputs$profitProjections[2,1,1,2,2,1,20] - outputs$profitProjections[1,1,1,1,1,1,20])/(outputs$HIntProjections[2,1,1,2,2,1,20] - outputs$HIntProjections[1,1,1,1,1,1,20]),
                        (outputs$profitProjections[2,1,1,2,2,1,20] - outputs$profitProjections[2,1,1,2,1,1,20])/(outputs$HIntProjections[2,1,1,2,2,1,20] - outputs$HIntProjections[2,1,1,2,1,1,20]),
                        outputs$timeToRecovery[1,1,1,1,1,1],
                        outputs$timeToRecovery[1,1,1,2,1,1],
                        outputs$timeToRecovery[2,1,1,2,1,1],
                        outputs$timeToRecovery[2,1,1,2,2,1],
                        dataInput[i,]$g_expected*(dataInput[i,]$phi_expected+1)/dataInput[i,]$phi_expected)
    
    
    masterOutputi = cbind(rep(dataInput[i,]$Species,nrow(melt.array(outputs$BProjections))),
                          melt.array(outputs$BProjections),
                          melt.array(outputs$HIntProjections)$value,
                          melt.array(outputs$profitProjections)$value,
                          melt.array(outputs$bProjections)$value,
                          melt.array(outputs$fTotalProjections)$value)[,-c(3,4)]
    
    colnames(masterOutputi) = c("species","management","catchShare","illegalFishing","implementYear", "time","biomass",
                                "harvest","profit","BvBMSY","FvFMSY")
    
    if (i == 1) {
        masterOutput = masterOutputi
    } else {
        masterOutput = rbind(masterOutput,masterOutputi)
    }
    
    masterOutput$management[masterOutput$management == 1] = "SQ"
    masterOutput$management[masterOutput$management == 2] = "FMSY"
    masterOutput$management[masterOutput$management == 3] = "minRec"
    masterOutput$management[masterOutput$management == 4] = "econOpt"
    masterOutput$management[masterOutput$management == 5] = "close"
    masterOutput$management[masterOutput$management == 6] = "openA"
    masterOutput$catchShare[masterOutput$catchShare == 1] = "no_CS"
    masterOutput$catchShare[masterOutput$catchShare == 2] = "CS"
    masterOutput$illegalFishing[masterOutput$illegalFishing == 1] = "illegal_fishing"
    masterOutput$illegalFishing[masterOutput$illegalFishing == 2] = "no_illegal_fishing"

    recoveryOutputi = cbind(rep(dataInput[i,]$Species,nrow(melt.array(outputs$timeToRecovery))),
                            melt.array(outputs$timeToRecovery))[,-c(3,4)]
    
    colnames(recoveryOutputi) = c("species","management","catchShare","illegalFishing","implementYear", "recTime")
    
    if (i == 1) {
      recoveryOutput = recoveryOutputi
    } else {
      recoveryOutput = rbind(recoveryOutput,recoveryOutputi)
    }
    
    masterOutput$management[masterOutput$management == 1] = "SQ"
    masterOutput$management[masterOutput$management == 2] = "FMSY"
    masterOutput$management[masterOutput$management == 3] = "minRec"
    masterOutput$management[masterOutput$management == 4] = "econOpt"
    masterOutput$management[masterOutput$management == 5] = "close"
    masterOutput$management[masterOutput$management == 6] = "openA"
    masterOutput$catchShare[masterOutput$catchShare == 1] = "no_CS"
    masterOutput$catchShare[masterOutput$catchShare == 2] = "CS"
    recoveryOutput$illegalFishing[recoveryOutput$illegalFishing == 1] = "illegal_fishing"
    recoveryOutput$illegalFishing[recoveryOutput$illegalFishing == 2] = "no_illegal_fishing"
    
        
}

masterOutput$implementYear = (masterOutput$implementYear + 1)
write.csv(masterOutput, file = paste0(pathstart, "master_output.csv"))
write.csv(summaryData, file = paste0(pathstart, "summary_data.csv"))

recoveryOutput$implementYear = (recoveryOutput$implementYear + 1)
write.csv(recoveryOutput,file = paste0(pathstart, "recoveryOutput.csv"))

## opt at the end of csv names indicates results from model with policy is optimized once

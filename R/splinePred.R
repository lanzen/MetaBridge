 ###################################################################################
 # Quantile spline regression spline using polynomial models                       # 
 # Follows the approach as closely as possible of Anderson (2008) "Animal-sediment #
 # relationships re-visited: Characterising species' distributions along  an       #
 # environmental gradient using canonical analysis and quantile regression splines"#
 #                                                                                 # 
 # Anders Lanzen Oct 2019            
 # updated 2020 to take into account AIC as in Keeley et al 2018
 ###################################################################################

require(quantreg)
require(splines)
require(irr)
require(vegan)
require(TITAN2)
 
 
 ## Return TITAN2 results from features and a specific variable (impact)
 ## dvTrain - dependent variable
 ## trainingData - a set of features like taxa or OTUs 
 titanInd = function(dvTrain, trainingData, minPrevalence=.2, ncpus=16){
   require(TITAN2)
   
   nFeatures = dim(trainingData)[1]
   
   ## Prevalance filter
   trainingData.pa = decostand(trainingData, method="pa")
   print(summary(colSums(trainingData.pa)))
   trainingData.prevFiltered = trainingData[,colSums(trainingData.pa)>nFeatures*minPrevalence]
   print(dim(trainingData.prevFiltered))
   
   ot = titan(dvTrain,trainingData.prevFiltered ,nBoot=100,numPerm=100,ncpus=ncpus)
   #goodInd = (ot$sppmax[,"obsiv.prob"]<=.1 & ot$sppmax[,"purity"]>=.9 & ot$sppmax[,"reliability"]>=.8)
   return(as.data.frame(ot$sppmax))
 }
 
 ## Predict splines from indicator list based on specific index

splinePredict = function(goodTITANInds, t, otus.train, md.train, imageOutDir = NA) {
  
  groups = data.frame(row.names=names(otus.train),group=rep(0,dim(otus.train)[2]))
  
  for (ti in c(1:dim(goodTITANInds)[1])){
    
    sp = row.names(goodTITANInds)[ti]
    titanGrp = goodTITANInds$maxgrp[ti]
    if (t=="Redox" | t=="Redox5yAvg" | t=="Redox10yAvg"| t=="NSI"| t=="ISI") titanGrp = 3-titanGrp
    
    ## Model distribution using 4th degree polynomial splines based on 95% quantiles 
    ## (Andersson et al 2008 and Nigel, except not sure how to check AIC and compare to 3 and 5 df)
    mt = data.frame(ab=otus.train[,sp],t=md.train[,t])
    
    bestAIC = 1E6
    bestModel = NA
    bestDF = NA
    for (df in c(2:5)){
      try = tryCatch({
        bsp <- rq(ab ~ bs(t,degree=df),data=mt,tau=.9)
        #print(AIC(bsp))
        aic = AIC(bsp, k=df)[1]
        if (aic<bestAIC-2){ 
          bestAIC = aic
          bestModel=bsp
          bestDF = df
        }
      },
      error= function(what_condition_my_condition_is_in){
        print(paste("Warning: degree",df,"spline failed for",sp))
        try=NULL
      })
    }
    
    print(paste("Best prediction for",df,"d.f. AIC =",bestAIC))
    
    if (bestAIC < 1E6 ){
      st = seq(min(md.train[,t]),max(md.train[,t]),length.out=1000)
      pred<-predict(bestModel,data.frame(t=st))
      maxT = st[pred==max(pred)]
      if(length(maxT)>1) maxT = maxT[1] # It can happen that two maxima are returned
      ecoGroup = getBIGroupFromValue(maxT,bi=t)
      
      ## Plot distribution and check indicator status manually if needed 
      
      if (!is.na(imageOutDir) & !is.na(ecoGroup)){
        
        change=NA
        if (goodTITANInds$maxgrp[ti]==2) {
          change <- "+"
        } else {
          if (goodTITANInds$maxgrp[ti]==1) change <- "-"
        }
        
        png(paste(imageOutDir,"/",sp,"_",t,".png",sep=""),width=600,height=600)
        
        plot(otus.train[,sp]~md.train[,t],
             main=paste(sp,change,"@",round(goodTITANInds$ienv.cp[ti],1)),
             sub=paste("Ecogroup",ecoGroup,"max @",round(maxT,1)),
             xlab=t,ylab="Abundance")
        abline(v=goodTITANInds$ienv.cp[ti], col="blue")
        abline(v=maxT, col="red", lty=2)
        lines(st,pred,col="orange")
        dev.off()
      }
      #print(paste0("titanGrp: ",titanGrp," ecoGrp: ",ecoGroup)) #DEBUG
      if ((titanGrp==1 & ecoGroup>3) | (titanGrp==2 & ecoGroup<3)){
        print(paste("Warning: conflict for",sp,"TITAN =",titanGrp,"but predicted EC =",ecoGroup))
        print("Setting EC to zero")
        ecoGroup = 0
      }
    }
    else {
        print(paste("Warning: cannot predict spline for",sp))
        ecoGroup = 0
    } 
    
    groups[sp,"group"] = ecoGroup
    
  }
  return (groups)
}


## Get TITAN QRS based BI based on features and specific variable (dvTrain) of specific type (dvType, determine using ec_and_plot)
## Returns dataframe with indicator group (indGrp), indicator value
titanQRS = function(metadata, dvName, trainingData, imageOutDir=NA, titanMinPrevalence=.2, 
                    titan_threshold=.95, ncpus=16){
  
  titanResults = titanInd(metadata[,dvName], trainingData, minPrevalence=titanMinPrevalence , ncpus=ncpus)
  tind = titanResults[titanResults$reliability>=titan_threshold & titanResults$purity>=titan_threshold,]
  tind.egs = splinePredict(goodTITANInds = tind, t = dvName, otus.train = trainingData, md.train = metadata, imageOutDir = imageOutDir)
  denovoBI = tind.egs[tind.egs$group>0,,drop=F]
  denovoBI$titanIndVal = tind[row.names(denovoBI),"IndVal"]
  #print(denovoBI) #DEBUG
  denovoBI$weight = NA
  for(i in c(1:dim(denovoBI)[1])){
    denovoBI$weight[i] = getWeight(denovoBI[i,"group"], bi = dvName)
  }
  return(denovoBI)
}

predictStatus = function(bioticIndex, disTable) {
  statusPred = rep(NA,dim(disTable)[1])
  
  # the features (taxa) in distable that match the biotic index
  bioticIndex_matching = row.names(bioticIndex)[row.names(bioticIndex) %in% names(disTable)]
  if(length(bioticIndex_matching)==0){
    write("Warning: No features in validation table correspond to index. \n",stderr())
    print(names(disTable)[1:10])
    print(row.names(bioticIndex))
    return(statusPred)
  }
  features_w_Ind = disTable[,bioticIndex_matching]
  
  # number of features (taxa) present in index per sample
  nFeatures = rowSums(decostand(features_w_Ind, method="pa"))
  
  # set predicted status as 0 initially for all samples that have at least one feature
  # in biotic index. Leave the rest as NA
  statusPred[nFeatures>0] = 0
  
  # Iterate through present featurs and add weight multiplied by abundance
  for(f in names(features_w_Ind)){
    statusPred = statusPred + bioticIndex[f,"weight"]*features_w_Ind[,f]
  }
  
  # Normalise by total abundance of present features matching index
  statusPred.norm = statusPred/rowSums(features_w_Ind)
  return(statusPred.norm)
}

## Cross-validation function for TITAN/QRS
titanQRSXVal = function(disTable, metadata, depvar, xvalParameter, titanMinPrevalence=.2, 
                        titan_threshold=.95, ncpus=16) {
  
  # Check that row names in disTable == metadata and dimensions are same
  if(dim(disTable)[1] != dim(metadata)[1] | sum(row.names(disTable) != row.names(metadata))>0){
    write("Distribution table and metadata table have different number or names of rows. Breaking.\n",stderr())
    return()
  }
  
  # Check that xval_parameter and depvar are present and have no NAs
  if((!depvar %in% names(metadata)) || sum(is.na(metadata[,depvar]))>0){
    write(paste0("The dependent variable ",depvar," is missing or has NAs. Breaking.\n"),stderr())
    return()
  }
  if((!xvalParameter %in% names(metadata)) || sum(is.na(metadata[,xvalParameter]))>0){
    write(paste0("The cross-valdiation thingie ",xvalParameter," is missing or has NAs. Breaking.\n"),stderr())
    return()
  }
  
  # Make depvar var for convenience
  dv = metadata[,depvar]
  print(summary(dv)) #DEBUG
  print("") #DEBUG
  
  # Make xval var and transform into a factor with no zero count values
  xval = droplevels(as.factor(metadata[,xvalParameter]))
  print(summary(xval)) #DEBUG
  print("") #DEBUG
  
  # Make empty list of predictions
  preds = rep(NA,dim(metadata)[1])
  
  # Iterating over all unique values of xval, train an RF model without data
  # matching that value and validate on the data that does
  for (leaveOut in unique(xval)){
    trainWith = (xval != leaveOut)
    trainingData = disTable[trainWith,]
    names(trainingData) = make.names(names(trainingData))
    valData = disTable[!trainWith,]
    names(valData) = make.names(names(valData))
    
    print(paste0("Validating on ",leaveOut,". Training data ",depvar," range:"))
    print(summary(dv[trainWith]))
    print(paste0("Validation data ",leaveOut," range:"))
    print(summary(dv[!trainWith]))
    print("") #DEBUG
    
    # Build partial de novo BI with TITAN and QRS
    xvalBI = titanQRS(metadata = metadata[trainWith,], dvName = depvar, trainingData = trainingData,
                       titanMinPrevalence = titanMinPrevalence, 
                       titan_threshold = titan_threshold, ncpus = ncpus)
    
    
    # Use on validiation data
    xvalPreds = predictStatus(xvalBI, valData)
    
    # Fill in predictions for validation (leave out) data
    preds[!trainWith] = xvalPreds
    
  }
  return(preds)
}

## Make biotic index on traininig dataset and validate on other
titanQRSTrainAndVal = function(disTableTrain, disTableVal, mdTrain, depvar, 
                               titanMinPrevalence=.2, titan_threshold=.95, ncpus=16) {
  
  # Check that row names in disTableTrain == mdTrain and dimensions are same
  if(dim(disTableTrain)[1] != dim(mdTrain)[1] | sum(row.names(disTableTrain) != row.names(mdTrain))>0){
    write("Distribution table and metadata table have different number or names of rows. Breaking.\n",stderr())
    return()
  }
  
  # Check that depvar are present and have no NAs
  if((!depvar %in% names(mdTrain)) || sum(is.na(mdTrain[,depvar]))>0){
    write(paste0("The dependent variable ",depvar," is missing or has NAs. Breaking.\n"),stderr())
    return()
  }
  
  names(disTableTrain) = make.names(names(disTableTrain))
  names(disTableVal) = make.names(names(disTableVal))
  
  print(paste0("Training data ",depvar," range:"))
  print(summary(mdTrain))
  
  # Make de novo BI based on disTableTrain
  trainBI = titanQRS(metadata=mdTrain, dvName = depvar, trainingData=disTableTrain,
                     titanMinPrevalence=titanMinPrevalence, titan_threshold=titan_threshold, 
                     ncpus=ncpus)
  
  # Predict based on disTableVal
  preds = predictStatus(trainBI, disTableVal)
  
  return(preds)  
}

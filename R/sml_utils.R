######################################################################################
## Supervised machine learning utilities (inspired by Tristan Cordiers sml_compo.R) ##
##               Anders Lanzen, Tassie, Feb 2024                                    ##
######################################################################################

# Load necessary packages
require(ranger)
require(irr)


## Helper function for training an RF model
trainRF = function(dvTrain, trainingData, mtry=NA, numTrees=300, algo="RF"){
  
  
  # Set some RF options according to algo and mtry
  notus = dim(trainingData)[2]
  
  # Set mtry
  if(is.na(mtry)){
    mtry = floor(notus/3)
  }
  
  if(algo=="RFProb"){
    splitrule = "extratrees"
    split.select.weights = c(rep(1/notus,notus-1),1)
  }
  else{
    splitrule = NULL
    split.select.weights = NULL
  }
  
  return(ranger(dvTrain ~ ., data=trainingData,
                mtry=mtry, num.trees = numTrees, 
                importance= "impurity", write.forest = T,
                splitrule = splitrule, split.select.weights = split.select.weights))
  #min.node.size=5 in Tristans script seems to not improve and used inconsistently
}
  
# Helper function for TITAN. Returns TRUE/FALSE vector for inclusion
titanFilter = function(dvTrain, trainingData){
  require(TITAN2)
  trainingData.pa = decostand(trainingData, method="pa")
  trainingData.prevFiltered = trainingData[,colSums(trainingData.pa)>50]
  ot = titan(dvTrain,trainingData.prevFiltered ,nBoot=50,numPerm=50,ncpus=12)
  #table(names(trainingData.prevFiltered) == row.names(ot$sppmax))
  goodInd = (ot$sppmax[,"obsiv.prob"]<=.1 & ot$sppmax[,"purity"]>=.9 & ot$sppmax[,"reliability"]>=.8)
  return(goodInd)
}

## rfXval returns predictions based on remove-one insteance cross-validatio.
## If not yet, metada is trasnformed to a factor and unused values removed. Only RF 
## supported so far. The issue with ranger requiring ugly hacking of global env
## should be far gone so this is a reimplementation using the normal way of 
## doing things. Tristans script caused inexplicable and scary behaviour probably
## bc this, earlier. Option optim_overfit removed due to being slightly like cheating
## (TODO implement method that checks this independently of crossval itself)
##
##
## disTable - the distributions to use (OTUs or taxa, normally) across samples to use for classification 
## metadata - metadata where var name of xval_parameter and index_parameter are present and w same row.names
## depvar - name of the dependent (target) variable to be predicted, 
##          e.g. a pressure or biotic index (must be numerical)
## xvalParameter - name of the parameter to do remove-one-class cross validation on (string)
## numTrees - the number of trees to use (default 300, >1000 should never be needed)
## algo - Random Forest algorithm: RF (default) or RFProb (experimental and slow using extreme randomised trees) 
## mtry - number of variables to possibly split per node. (default is #variables/3, from Tristan)
## filterByTITAN - Use TITAN2 (50 iterations, 50 bootstraps, min prevalence=50 to check if features have
##                 significant threshold with depvar and limit RF model to those that do (very slow and experimental)

rfXval = function(disTable, metadata, depvar, xvalParameter, 
                  numTrees=300, algo="RF", mtry=NA, filterByTITAN=FALSE,
                  filterByRF = FALSE, selectedFeatures=100) { 
  
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
  #print(summary(dv)) #DEBUG
  #print("") #DEBUG
  
  # Make xval var and transform into a factor with no zero count values
  #xval = droplevels(as.factor(metadata[,xvalParameter]))
  xval = as.factor(metadata[,xvalParameter])
  print(summary(xval)) #DEBUG
  #print("") #DEBUG
  
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
    # print(paste0("Validation data ",leaveOut," range:"))
    # print(summary(dv[!trainWith]))
    print("")

    # Filters training data by TITAN if specified (not recommended) and then val data likewise
    if(filterByTITAN){
      goodInd = titanFilter(dvTrain = dv[trainWith], trainingData = trainingData)
      if(sum(goodInd)<2) {
        write(paste0("TITAN2 failed and returned only ",sum(goodInd), "filtered features \n"),stderr())
      }
      else{
        trainingData = trainingData[,goodInd]
        valData = valData[,goodInd]
      }
    }
    if(filterByRF){
      ## FEATURE SELECTION USING RF AND RETRAINING WITH 100 MOST IMPORTANT
      # Train an RF model with all data as would be used for prediction (not allowing
      # algo RFProb)
      rfAll = trainRF(dvTrain=dv[trainWith], trainingData=trainingData, 
                      mtry=mtry, numTrees = numTrees, 
                      algo = "RF")
      #print(rfAll) #DEBUG
      # Select the 100 taxa with highest variable importance
      imp = tail(sort(rfAll$variable.importance), selectedFeatures)
      #print(imp) #DEBUG
      trainingData = trainingData[,names(imp)]
    }
    
    #  why is this the default? It does seem to be a sweetspot from testing so keeping
    if(is.na(mtry)){
      mtry = floor(dim(trainingData)[2]/3)
    }
    # Train RF model
    rfModel = trainRF(dvTrain=dv[trainWith], trainingData=trainingData, 
                      mtry=mtry, numTrees = numTrees, 
                      algo = algo)
    

    # Use on validiation data
    rfPred = predict(rfModel, valData)$predictions
    
    # Fill in predictions for validation (leave out) data
    preds[!trainWith] = rfPred
    
  }
  return(preds)
}

## rfTrainAndVal:
## Train on disTableTrain, predict on disTableVal
## 

rfTrainAndVal = function(disTableTrain, disTableVal, mdTrain,  
                         depvar, numTrees=300, algo="RF", mtry=NA,
                         filterByRF = FALSE, selectedFeatures=100) { 
  
  # Check that row names in disTableTrain == mdTrain and dimensions are same
  if(dim(disTableTrain)[1] != dim(mdTrain)[1] | sum(row.names(disTableTrain) != row.names(mdTrain))>0){
    write("Distribution table and metadata table have different number or names of rows. Breaking.\n",stderr())
    return()
  }
  
  # Check that row names in disTableTrain and distTableVal are same
  if(dim(disTableTrain)[2] != dim(disTableVal)[2] | sum(names(disTableTrain) != names(disTableVal))>0){
    write("Training and validation datasets have different names of rows. Breaking.\n",stderr())
    print(dim(disTableTrain))
    print(dim(disTableVal))
    print(sum(names(disTableTrain) != names(disTableVal)))
    return()
  }
  
  # Check that depvar are present and have no NAs
  if((!depvar %in% names(mdTrain)) || sum(is.na(mdTrain[,depvar]))>0){
    write(paste0("The dependent variable ",depvar," is missing or has NAs. Breaking.\n"),stderr())
    return()
  }
  
  dv = mdTrain[,depvar]
  
  if(filterByRF){
    ## FEATURE SELECTION USING RF AND RETRAINING WITH 100 MOST IMPORTANT
    # Train an RF model with all data as would be used for prediction (not allowing
    # algo RFProb)
    rfAll = trainRF(dvTrain=dv, trainingData=disTableTrain, 
                    mtry=mtry, numTrees = numTrees, 
                    algo = "RF")
    #print(rfAll) #DEBUG
    # Select the 100 taxa with highest variable importance
    imp = tail(sort(rfAll$variable.importance), selectedFeatures)
    #print(imp) #DEBUG
    disTableTrain = disTableTrain[,names(imp)]
  }
  
  if(is.na(mtry)){
    mtry = floor(dim(disTableTrain)[2]/3)
  }
  
  names(disTableTrain) = make.names(names(disTableTrain))
  names(disTableVal) = make.names(names(disTableVal))
    
  print(paste0("Training data ",depvar," range:"))
  print(summary(mdTrain))
  
  
  # Make RF model based on disTableTrain
  rfModel = trainRF(dvTrain=dv, trainingData=disTableTrain, 
                      mtry=mtry, numTrees = numTrees, 
                      algo = algo)
  
  # Predict based on disTableVal
  preds = predict(rfModel, disTableVal)$predictions
  
  return(preds)  
}
  
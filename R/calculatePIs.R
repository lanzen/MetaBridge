###################################################################################
# Calculates abiotic pressure indices based on Aylagas et al (2017), but without  #
# including redox potential. Metabrdige version.                                  #
# 2025-04-30: Using thresholds established in Metamon                             # 
##################################o####-##########y################################

require(stringr)

# Helper functions to calculate PIs
floorZeroSix = function(x){
  zeroFloor = max(0,x)
  return(zeroFloor)
}

roofSix = function(x){
  roofSix = min(6,x) # Inappropriate since v high single metal contamination counts
  return(roofSix)
}

## PI value for a single variable
partialPI = function(md, i, limitDF, verb=FALSE){
  partialMD = md[i,row.names(limitDF)]
  if(sum(is.na(partialMD))<dim(limitDF)[1]) {
    if(verb) print("Values:")
    if(verb) print(partialMD)
    if(verb) print("Relative limit:")
    if(verb) print(partialMD/t(limitDF))
    piValues = log2(2*partialMD/t(limitDF)) +1
    if(verb) print(piValues)
    piValues = sapply(piValues, floorZeroSix)
    if(verb) print(piValues)
    ppi = mean(piValues,na.rm=TRUE)
    ppi = min(6,ppi)
    return(ppi)
  }
  else return(NA)
}


# Calculate and add PIs to metadata, returning new dataframe
calculatePIs = function(md){
  
  mdx = md # Make copy of md metadata to add variables that we later dont need
  #mdx$modRedox = 600 - mdx$Redox
  
  # Add undefined PIs or replace those that are there
  md$piMetals = NA # Heavy metals
  md$piHC = NA     # Hydrocarbons only
  #md$piEutro = NA     # Organinc material, TOC and NH3/NH4 as proxy of nutrient enrichment
  #md$piBI = NA     # Biotic index or expert assessment ranging from 0 - 6 where 0 is best
  md$pi   = NA     # Average of all the above
  #md$pi_physchem   = NA     # Average of all the above
  
  # Hardcoded (still bit hacky) limits

  # Max in Ref. stations x2 
  metalLimits = data.frame(row.names=c("Hg","Cd", "Cu", "Ba", "Pb"),
                            limits = c(.052,.15,70,250, 32))

  # Pb: 32, Cr: 24
  # metalLimits = data.frame(row.names=c("Cu", "Ba"),
  #                          limits = c(84,464))
  
  # hcLimits = data.frame(row.names=c("TotalPCB","TotalDDT","TotalHCH", "ADE", "TotalPAH","TotalHC"), 
  #                       limits = c(22.7,3.89,0.32, 6, 1607, 538))
  
  hcLimits = data.frame(row.names=c("TotalPAH","THC"), #Same as established in Metamon
                        limits = c(2,72))
  
  # omLimits = data.frame(row.names=c("TOC"), 
  #                        limits = c(1)) 
  
  # omLimits = data.frame(row.names=c("OM","TOC","NH4_pw","TotalN","TotalP", "modRedox"), 
  #                       limits = c(2,1,217,550,600,300)) 

  #includedBIs = c("AMBI","OtherAsAMBI") #,"MacrophyteAsAMBI", "AusPI") awaiting confirmation
  
  # Ensure that all variables are there else put them as NA
  allLimits = rbind(metalLimits, hcLimits)#, omLimits)
  missingVars = row.names(allLimits)[!row.names(allLimits) %in% names(mdx)]
  mdx[,missingVars] = NA
  
  # # Same for BIs
  # missingBIs = includedBIs[!includedBIs %in% names(md)]
  # mdx[,missingBIs] = NA
  
  # Normalise hcLimits and hcLimits except TotalHC
  mdx$TOC.norm = mdx$TOC
  mdx$TOC.norm[is.na(mdx$TOC.norm)] = 1
  mdx$TOC.norm[mdx$TOC.norm>10] = 10
  mdx$TOC.norm[mdx$TOC.norm<.5] = .5
  #hcMinusTHC = row.names(hcLimits)[c(1:5)] # Parameters for normalising to 1% TOC
  hcMinusTHC = "TotalPAH"
  #mdx[,row.names(hcLimits)] = mdx[,row.names(hcLimits)] / mdx$TOC.norm
  mdx[hcMinusTHC] = mdx[hcMinusTHC] / mdx$TOC.norm
  
  # Iterate over all rows
  for(i in c(1:dim(md)[1])) {
   
    md[i,"piMetals"] = partialPI(mdx,i,metalLimits)
    md[i,"piHC"] = partialPI(mdx,i,hcLimits)
    #md[i,"piEutro"] = partialPI(mdx,i,omLimits)
    
    #If at least one PIs has value, set overall pi to average over PIs, balancing up
    if(sum(is.na(md[i,c("piMetals","piHC")]))<2){
      md[i,"pi"] = mean(c(md$piMetals[i], md$piHC[i], md$piEutro[i], md$piBI[i]), 
                      na.rm=TRUE)
    }
  }
  return(md)
}

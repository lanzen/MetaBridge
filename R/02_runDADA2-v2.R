library(dada2)
packageVersion("dada2")
library(stringr)
library(dplyr)
packageVersion("dplyr")
library(ggplot2)


runDADA2 = function(fwd, filt, rev, filt.rev, truncLen=0, minOverlap = 12,
                    maxOverlapMismatch = 0, lengthMax=FALSE, lengthMin=FALSE,
                    lib_name = "SVs",cores=TRUE){
  
  allSampleNames = names(filt)
  
  
  ## Filter and trim
  print("## Filtering and Trimming reads ##")
  filtering <- filterAndTrim(fwd=fwd, filt = filt,
                       rev=rev, filt.rev = filt.rev,
                       truncLen=truncLen, 
                       matchIDs=TRUE, maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE, 
                       compress=TRUE, multithread=cores, verbose=TRUE) 

  #print(head(filtering)) # run to DEBUG
  
  ## Update sample lists to remove empty samples
  exists <- file.exists(filt) & file.exists(filt.rev)
  #print(exists) # run to DEBUG
  filt <- filt[exists]
  filt.rev <- filt.rev[exists]
  
  ## Learn errors from filtered data
  print("## Learning Errors ##")
  errF <- learnErrors(filt, nbases=5e8, multithread=cores, randomize=TRUE, 
                      errorEstimationFunction = loessErrfun_mod, verbose = TRUE,MAX_CONSIST=10) #MAX_CONSIST=20 if not converging and want to test
  
  errR <- learnErrors(filt.rev, nbases=5e8, multithread=cores, randomize=TRUE, 
                      errorEstimationFunction = loessErrfun_mod, verbose = TRUE,MAX_CONSIST=10) #MAX_CONSIST=20 if not converging
  ## Learn errors from filtered data, and time the procedure:
  #system.time(errF <- learnErrors(filt, nbases=5e8, multithread=cores, randomize=TRUE, verbose = TRUE, errorEstimationFunction = loessErrfun_mod))
  #system.time(errR <- learnErrors(filt.rev, nbases=5e8, multithread=cores, randomize=TRUE, verbose = TRUE, errorEstimationFunction = loessErrfun_mod))

  
  ## save the error model for both fwd and rev
  saveRDS(errF, paste0("dada2/",lib_name,"_mod_errF.rds"))
  saveRDS(errR, paste0("dada2/",lib_name,"_mod_errR.rds"))
                  
  ## Plot error models
  plot.errF <- plotErrors(errF, nominalQ = T)
  ggsave(plot = plot.errF, path = "dada2/img", filename= paste0(lib_name,"_mod_plot_errF.png"), device="png", width=20, height=20, units="in")
  plot.errR <- plotErrors(errR, nominalQ = T)
  ggsave(plot = plot.errR, path = "dada2/img", filename = paste0(lib_name,"_mod_plot_errR.png"), device="png", width=20, height=20, units="in")
  
  # dereplicate reads now
  print("## Dereplicating identical sequences ## ")
  drpF <- derepFastq(filt)
  drpR <- derepFastq(filt.rev)
  
  ## Core sample inference, ie the ML algo to denoise, before overlapping
  ## use pool="pseudo" for more accurate inferal of rare variants
  print("## SV inference ## ")
  dadaFs = dada(drpF, err=errF, pool=FALSE, multithread = cores,
                errorEstimationFunction = loessErrfun_mod, verbose=TRUE, selfConsist=FALSE)
  dadaRs = dada(drpR, err=errR, pool=FALSE, multithread = cores,
                errorEstimationFunction = loessErrfun_mod, verbose=TRUE, selfConsist=FALSE)

                 
  ## Merge read pairs
  print("## Merging Sequences ## ")
  mergers <- mergePairs(dadaFs, drpF, dadaRs, drpR, 
                        maxMismatch = maxOverlapMismatch, 
                        minOverlap = minOverlap, verbose=TRUE)
  
  
  ## Construct sequence table
  print("## Making sequence table ## ")
  seqtab <- makeSequenceTable(mergers)
  #print(dim(seqtab)) # run to DEBUG
  
  print("- Inspect ASV length distribution : ")
  print(table(nchar(getSequences(seqtab))))
  
  ## Filter by length
  print(paste0("- Removing SVs outside length range [", lengthMin, ":", lengthMax, "]"))
  seqtab.length <- seqtab[,nchar(colnames(seqtab)) %in% seq(lengthMin,lengthMax)]
  
  ## Remove chimeras
  print("## Removing chi(bi)meras ## ")
  seqtab.nochim <- removeBimeraDenovo(seqtab.length, method="consensus", 
                                      multithread=cores, verbose=TRUE)
  
  print("- Percentage of non chimeric reads:")
  print(sum(seqtab.nochim)/sum(seqtab)*100)
  print("- Percentage of non chimeric ASVs :")
  print(ncol(seqtab.nochim)/ncol(seqtab)*100)
  

  print(" ## Producing final outputs ## ")
  #print(dim(seqtab.nochim)) #run to DEBUG
  
  track <- cbind(filtering, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab.length), rowSums(seqtab.nochim),rowSums(seqtab>0),rowSums(seqtab.length>0),rowSums(seqtab.nochim>0))
  colnames(track) <- c("input", "filtered", "denoised", "merged", "length.filtered", "nonchim","uniqueSVs","uniqueSVsAfterLengthFilter","uniqueSVsAfterBimera")
  rownames(track) <- sample.names
  write.csv(as.data.frame(track), file=paste0("dada2/",lib_name,"_read_counts.csv"))
  return(seqtab.nochim) # command producing the seqtab.nochim.Cat object
}

# Error model inference with binned NovaSeq basecalling
# see:
# https://github.com/benjjneb/dada2/issues/1307
# best solution is to alter loess arguments (weights, degree, span OR just weights and span) & enforce monotonicity
# but it is recommended to check a couple of solutions and compare!
loessErrfun_mod <- function (trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow = 0, ncol = length(qq))
  for (nti in c("A", "C", "G", "T")) {
    for (ntj in c("A", "C", "G", "T")) {
      if (nti != ntj) {
        errs <- trans[paste0(nti, "2", ntj), ]
        tot <- colSums(trans[paste0(nti, "2", c("A","C", "G", "T")), ])
        rlogp <- log10((errs + 1)/tot)
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q = qq, errs = errs, tot = tot,
                         rlogp = rlogp)
        # original
        #mod.lo <- loess(rlogp ~ q, df, weights=tot)
        # Guillem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        ## mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2) #uncomment for Salazar
        # hhollandmoritz' 4th solution based on JonaLim (alter weights, span, and degree).
        mod.lo <- loess(rlogp ~ q, df, weights=log10(tot), degree=1, span = 0.95) #comment for salazar
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred) > maxrli] <- pred[[maxrli]]
        pred[seq_along(pred) < minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      }
    }
  }
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-07
  est[est > MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est < MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # # enforce monotonicity, @hhollandmoritz
  # # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  ##--
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1 - colSums(est[1:3, ]), est[1:3, ],
               est[4,], 1 - colSums(est[4:6, ]), est[5:6, ],
               est[7:8, ], 1 - colSums(est[7:9, ]), est[9, ],
               est[10:12, ], 1 - colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A", "C", "G", "T"), each = 4),
                          "2", c("A", "C", "G", "T"))
  colnames(err) <- colnames(trans)
  return(err)
}

sampleName = function(filename){
  return(str_sub(filename,start = 1,end = nchar(filename)-21)) # this removes 21 characters from the end of filename
  #return(sub("_R1_trimmed\\.fastq\\.gz$", "", filename))
}

getN <- function(x){
  return(sum(getUniques(x)))
}

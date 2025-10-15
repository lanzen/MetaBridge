#!/bin/bash

# 01_demultiplex.sh lib-metadata.csv Lib1 [Lib2...] (space-separated list of library names to demultiplex) > logfile 
# example usage: (win users need to remove \r at end of lines before using lib-metadata files: sed -i 's/\r//g' lib-metadata.csv
# bash 01_demultiplex.sh lib-metadata-metabridge.csv LIB_1 LIB_8 > logfile_demult_LIB1_LIB8.log
# The 2 scripts used need to be symlinked to the WD! The script looks for all input data in the WD!
#
# Anders Lanzen for Metabridge, 2023-11-21
# Final comments/edits by Miriam Brandt, 2024-01-10
#
# Runs cutadapt for Demultiplexing paired-end reads sequentially, because it does not work (too many file handles)
# with combinatorial dual indexes (see https://cutadapt.readthedocs.io/en/stable/guide.html#combinatorial-demultiplexing).
#

# The number of cores to use - increase for faster execution on a server
cores=16

libFile=$1 # Library info file

libsToDemultiplex=${@:2}
echo "Attempting to demultiplex: $libsToDemultiplex" #DEBUG

# Activate conda environment with cutadapt installed (v3+)
# Comment if using system-wide cutadapt is preferred
eval "$(conda shell.bash hook)"
module load Miniconda3
# shellcheck source=/dev/null
source "$EBROOTMINICONDA3"/etc/profile.d/conda.sh
conda init bash
#conda activate cutadapt
conda activate /scratch/ssd/fastwork/metabridge/common/conda/cutadapt-v4.5
## --


# Read library metadata file - line by line
skip_headers=1
while IFS=, read -r col1 col2 col3 col4 col5 col6 col7
do
    if ((skip_headers))
    then
        ((skip_headers--))
    else
	#DEBUG: echo $col1 $col2 $col5 $col6 $col7
    libName=$col1 # how you want the lib folder to be named in the demult/ output folder
	fqr1=$col2 # (path to) Forward FASTQ multiplexed file
	fqr2=$col3 # (path to) Reverse FASTQ multiplexed file
	fPrimer=$col4 # Forward primer sequence (no tag)
	rPrimer=$col5 # Reverse primer sequence (no tag)
	t2sFile=$col6 # (path to) The tag2sample metadata file (3 columns: sample names, F tag name , and R tag name). Will only use the first 3 columns.
	tagSeq=$col7 # The sequences (in FASTA format, i.e. incl names) of F and R tags used to index this library
	
	echo "$libName"
		
	## Proceed with demultiplexing only if the lib of this line in the metadata corresponds to a user supplied lib name
	for lib in $libsToDemultiplex
	do
	    if [ "$libName" == "$lib" ]; then
	    
	    ## Make directory named after library under demult if it does not already exist.
		## Sample-specific re-oriented files will be placed here.
		echo "Demultiplexing $lib"
		mkdir -p demult/"$libName"

		## Make metadata file called [lib]_cutadapt.csv using 01_prepCutadaptData.py
		## to link tag sequences from FASTA to sample names from t2s
		libMeta=${libName}_cutadapt.csv

		python3 01_prepCutadaptData.py "$t2sFile" "$tagSeq" > "$libMeta"

		## Iterate over samples
		
		while IFS=, read -r j1 j2 j3
		do
		    sample=$j1
		    fTag=$j2
		    rTag=$j3
		   
		    ## Sample matching: forward tag j found in R1 (and reverse tag in R2), always in the first 8 bases (^)
		    ## Output *R1_F.fastq is R1 and R2_R is R2. No errors allowed over the 8 nt tag (-e).
			## cutadapt by default removes what it is looking for.
		    cutadapt -e 0 -O 8 --no-indels -j $cores --discard-untrimmed --max-n=0 -g ^$fTag -G ^$rTag \
			     -o ${sample}_R1_F.fastq.gz -p ${sample}_R2_R.fastq.gz ${fqr1} ${fqr2}

		    ## Sample matching: forward tag j in R2, and *R1_R.fastq output is R2
		    cutadapt -e 0 -O 8 --no-indels -j $cores --discard-untrimmed --max-n=0 -g ^$fTag -G ^$rTag \
			     -o ${sample}_R2_F.fastq.gz -p ${sample}_R1_R.fastq.gz ${fqr2} ${fqr1}
			     
		    ## Remove target primer sequences in concatenated sample. Default error rate (10%)
		    fPrimLength=${#fPrimer}
		    echo "DEBUG: forward primer is $fPrimLength nt"
		    rPrimLength=${#rPrimer}
		    echo "DEBUG: reverse primer is $rPrimLength nt"
			## set overlap in cutadapt to be the length of the shortest primer
		    if [[ $fPrimLength -gt $rPrimLength ]]; then
		    	overlap=$rPrimLength
		    else
		    	overlap=$fPrimLength
		    fi
		    echo "DEBUG: Using cutadapt overlap (-O) of $overlap nt to remove primers"
			
		    #can enforce a minimum length using cutadapt with the --minimum-length N flag to ensure no zero length (empty) sequences
		    cutadapt -j $cores -e 0.1 --no-indels --max-n=0 --discard-untrimmed -O $overlap -g $fPrimer -G $rPrimer \
			     -o  ${sample}_R1_F_clean.fastq -p ${sample}_R2_R_clean.fastq \
			     ${sample}_R1_F.fastq.gz ${sample}_R2_R.fastq.gz
			     
		    cutadapt -j $cores -e 0.1 --no-indels --max-n=0 --discard-untrimmed -O $overlap -g $fPrimer -G $rPrimer \
			     -o  ${sample}_R2_F_clean.fastq -p ${sample}_R1_R_clean.fastq \
			     ${sample}_R2_F.fastq.gz ${sample}_R1_R.fastq.gz
		
		    		     
		    ## Concatenate directions, then compress all FASTQ files and move to output directory
		    cat ${sample}_R1_F_clean.fastq ${sample}_R2_F_clean.fastq > ${sample}_cat_F_clean.fastq
		    cat ${sample}_R2_R_clean.fastq ${sample}_R1_R_clean.fastq > ${sample}_cat_R_clean.fastq
		    gzip ${sample}_cat_F_clean.fastq ${sample}_cat_R_clean.fastq
		    mv ${sample}_cat*.fastq.gz demult/$libName
		    
		    ##Comment to keep R1 and R2 derived files separately
		    rm ${sample}_R1_F_clean.fastq ${sample}_R2_F_clean.fastq ${sample}_R2_R_clean.fastq ${sample}_R1_R_clean.fastq
		    rm ${sample}_R1_F.fastq.gz ${sample}_R2_F.fastq.gz ${sample}_R2_R.fastq.gz ${sample}_R1_R.fastq.gz
		   
		done < "$libMeta"				
		#exit
	    fi
	done
    fi
done < "$libFile"
conda deactivate

# prepCutadaptData.py lib-name tag-metadata.csv tags.fasta
#
# Anders Lanzen for Metabridge 2023-10-29
#
# Helper script for demultiplex.sh. Reads tag sequences and prepares metadata file with
# tag1-seq tag2-seq
import sys
import re
import os

#libName = sys.argv[1]
metadata = sys.argv[1]
tagFile = sys.argv[2]

tags = open(tagFile,"r")
tagSeq = {}
# Read FASTA format tag sequences

for line in tags:
    if line[0]==">":            
        name = line[1:].rstrip()
    else:
        seq=line.rstrip()
        tagSeq[name] = seq
           

tags.close()

# Print data for cutadapt of library and with tag sequences
sampleData = open(metadata,"r")
header = True
for line in sampleData:
    if header:
        header=False
    else:        
        items = re.split(r'[\t,;]+', line.replace(" ","").rstrip())
        sampleName = items[0]        
        fTag = items[1]
        rTag = items[2]
        print("%s,%s,%s"%(sampleName,tagSeq[fTag],tagSeq[rTag]))
sampleData.close()

            

        

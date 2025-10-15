#!/usr/bin/env python

# Written By Anders Lanzen and Miriam Brandt, 13.02.2024

# Usage:
#I have modified 05_get_SWARM_table_from_SV_and_clustering_annotated.py and put it in the 0_scripts folder. The new file:
#•	is executable (you don't have to type "python" before it, which is better because you can use tab to select input files and you can put it in the Linux $PATH)
#•	uses an options parser with -i being representative SV sequences from SWARM or DADA2, -t being the SV table and -c being the SWARM clustering output (".swarms" file)
#•	if you want to see how it works simply use 05_get_SWARM_table_from_SV_and_clustering_annotated.py -h which produces the output below
#•	writes output FASTA and a SWARM table to the name stub given by option -o instead of writing a table to standard out (if you do not provide a name it will write it to SWARM.tsv and SWARM.fasta)
#•	uses the name SWARM_i counting from 1 and up instead of the SV names
#There is an example output in 0_scripts_resources, see DADA2SWARM_example.zip

# 05_get_SWARM_table_from_SV_and_clustering_annotated.py -i LIB_8_DADA2_SWARM.fasta -t LIB_8_DADA2.tsv -c LIB_8_DADA2_SWARM.swarms -o LIB_8_DADA2_SWARM_parsed
# Assumes tab-separated sequence variant (SV) table giving per sample read abundance, with SVs in rows and samples in columns
# Assumes clustering file format like:
#ASV_1;size=206990
#ASV_2;size=152401 ASV_545;size=1100 ASV_930;size=500 ASV_1355;size=292 ASV_3140;size=83 ASV_6693;size=13

import sys
from optparse import OptionParser

parser = OptionParser()
    

parser.add_option("-i", "--fasta-in",
                      dest="fasta",
                      type="string",
                      default=None,
                      help="Representative SV sequence input in FASTA format")

parser.add_option("-t", "--sv-table",
                      dest="svtable",
                      type="string",
                      default=None,
                      help="SV table from DADA2")

parser.add_option("-c", "--swarms",
                  dest="swarmfile",
                  type="string",
                  default=None,
                  help="SWARM cluster file")

parser.add_option("-o", "--output-stub",
                      dest="output",
                      type="string",
                      default="SWARMs",
                      help="File name (without suffix) for Name of output SWARM OTU table and FASTA files ")

(options, args) = parser.parse_args()

if not (options.svtable and options.swarmfile and options.fasta):
    sys.stderr.write("\nSV table, SWARM cluster file, and representative sequences required. Quitting.\n")
    sys.exit()


## Read all sequences and store in sv_seq
sv_seq = {}
seq=""
svn=""

svRep = open(options.fasta)
for line in svRep:
    if line.startswith(">"):

        if seq and svn:
            sv_seq[svn] = seq
        
        
        svn=line[1:line.find(";")]
        seq=""
    else:
        seq = seq + line.rstrip()

sv_seq[svn] = seq

svRep.close()

# Open files for reading and writing

repSeqOut = open(options.output + ".fasta","w")
swarmTableOut = open(options.output + ".tsv","w")

clust = open(options.swarmfile)
svTab = open(options.svtable)

# make abundance dictionary
abundances = {}
# we do not care about 1st line of the ASV table, as only headers
firstLine = True
# list all abundance values per line of ASV table:
for line in svTab:
    # Extract abundance values from tsv table, and put this in a list of strings called ab. 
    # For each line, remove end of line sign (equivalent of line.rstrip('\n')) and split lines of sv file (split by tab). Otherwise python will add extra empty lines for each end of line sign.
    ab = line[:-1].split("\t")  
    if firstLine:
        #Print first line of SV table from beginning to the end
        swarmTableOut.write(line)
        ## identify first column header entry by searching for the first tabular character, and call it tabPosition
        #tabPosition = re.search("\t",line).start()
        ##Print swarm_name and print first line of SV table from tabPosition (i.e. removing "SV_names" title) to the end, removing end of line sign (-1) ==> creates new header line 
        #print("swarm_name" + line[tabPosition:-1])
        # remove firstline i.e. headers
        firstLine = False
    else:
        # fill abundance dictionary : ab[0] is first item in the list ab, i.e. "SV_x". So for every SV_name (the key), list the abundance values
        abundances[ab[0]] = [int(i) for i in ab[1:]]
# close the sv file
svTab.close()

# SWARM output number
swarm_nb = 0
# read the clustering results to know which SVs to sum within each swarm
for line in clust:
    swarm_nb = swarm_nb+1
    # remove end of line sign, and split lines separated by spaces, as clustering.txt is a text file. Create asvs_w_size list, as names are in the form ASV_x;size=xx
    asvs_w_size = line[:-1].split(" ")
    # make new empty list that will contain asv names of each swarm OTU
    asvs = []
    # for every line (swarm OTU) of the clustering file, extract corresponding SV_names, by using ";" key (0 means before)
    for aws in asvs_w_size:
        asvs.append(aws.split(";")[0])
    # take abundance values of first asv of each swarm OTU as start and put this into aSum list
    aSum = abundances[asvs[0]]
    # for all the remaining asv names making up each swarm (a from 1 to end), extract the abundance values and put these into new list(aab)
    for a in asvs[1:]:
        aab = abundances[a]
        # in the list aSum (abundance table of first SV of each swarm), we add the abundances of the first SV in each swarm to the abundances of the other SVs making it up.
        for i in range(len(aSum)):
            aSum[i] = (aSum[i] + aab[i])
# combine swarm OTU name (the first SV, i.e. (asvs[[0]]), so the seed SV of each swarm OTU) with final swarm OTU abundance in aSum, and print in tabular format
    swarm_name = "SWARM_"+str(swarm_nb)
    swarmTableOut.write("\t".join([swarm_name]+[str(i) for i in aSum]))
    swarmTableOut.write("\n")
# combine swarm OTU new name with final OTU abundance in aSum, and print in tabular format
    #print("\t".join(["swarm_"+str(swarm_nb)]+[str(i) for i in aSum]))

    ## Write fasta rep sequence
    repSeqOut.write(">%s\n%s\n" % (swarm_name, sv_seq[asvs[0]]))

# Close file handles
clust.close()
swarmTableOut.close()
repSeqOut.close()



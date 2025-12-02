findTags_mib.py
-----------

This script helps defining the tags to be attached 5' of amplification primers. 
It assumes that the user has collected the consensus sequence of an aligned, representative set of sequences corresponding to the n positions located 5' of annealing position of the amplification primer sequence, with n being the length of the desired tags.
The only required input is one or more than one sequence based on which to design tags. 
The outputs are fasta files containing the tags (fused to the amplification primer sequence upon user querying) and possibly a log file to keep track of all selections.

This script is based on Frank L's findTags.py script, and was adapted to python3 by Miriam Brandt on 21.02.2022. 
It should run properly on MacOS >8 or windows, using a unix terminal such as MobaXterm.
The module Levenshtein should be installed, using for e.g.:
python -m pip install Levenshtein

Please do not distribute - it is not commented nor licensed! :)
© Franck Lejzerowicz - all rights reserved.


Usage:
python3.x findTags_mib.py [-h] -s [S [S ...]] [-max MAX] [-min MIN] [-d D] [--v]
                   [-log] [-i I] [-n N]

optional arguments:
  -h, --help      show this help message and exit
  -s [S [S ...]]  Consensus anti-complementary sequence(s): the nucleotides that you do NOT want (if several possibilities at the same position, then use IUPAC code).
  Note: more than one set of tags can be designed. For e.g., you can specify the forward and then reverse sequences.
  -max MAX        Maximum number of tag [default = 30]
  -min MIN        Minimum number of tag [default = 20]
  -d D            Minimum number of differences [default = 3]
  -i I            Maximum number of dinucleotide occurrences (~homopolymers)
                  in a tag [default = 0]
  -n N            Maximum number of selections to be searched before querying
                  to stop [default = 50]
  --v             Verbose mode (prints the intermediary tags) [default = On]
  --log           Write a log file with all the searches [default = On]

Notes:
- The input sequence(s) should be the consensus sequence in regular 5'-3' orientation and may include IUPAC degenerate code.
- The interaction of the -min -max and -d parameters may cause to program to look for primers in a search space from which the desired number of primers could not be attained. 
If the script searched for more than -n selections without result, it will query you to stop (press enter) and up to you to restart with different parameters (less primers, increased tag length); 
otherwise press any key + enter and possibly wait until another query after -n more searches (this may cause indefinite number of searches!).
- By default, no tag containing homopolymers will be included in the selections but since this should not be a problem for Illumina sequencing base calling, it could be modified using the -i option. Moreover, it increases the search space for highly degenerate consensus sequences. If i is set to e.g. 4, then selections will first try to include tags with 1 dinucleotide occurrence, then tags with 2 dinucleotide occurrences, etc.

Example on mac/linux/win:
python3.8 findTags_mib.py -s AAGGCGTC CGTATGCG -min 10 -max 10 -i 2
python findTags_mib.py -s MWARYYHS GYGSBYYY -min 21 -max 23 -i 2

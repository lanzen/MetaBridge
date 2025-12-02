#!/usr/bin/env python

#### Exporting 8 nt before the primers and the amplified fragment from a target database (ecoPCR)

#### script from Tristan Cordier, adapted to python3 by Miriam Brandt
#### e.g. use: python ecoPCR_and_tags.py --db SILVA_SSURef_NR99_taxsilva_U-T_noquote.fasta -f GTGYCAGCMGCCGCGGTAA -r CCGYCAATTYMTTTRAGTTT -s 100
#### careful: the databse file needs to have no U nucleotides (otherwise primers can't "anneal") and no special characters in headers.


import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import random
# Bio.Alphabet is no longer in biopython, since sept 2020 it was deprecated
#from Bio.Alphabet import generic_nucleotide
#from Bio.Alphabet import IUPAC
from Bio import SeqUtils
import random


def main():
	"""
	In silico PCR from an unaligned reference sequence database.
	IUPAC nt are allowed and primers seqs have to be in 5' to 3' direction (reverse will be reverse complement).
	Three files are being outputed:
	AMPLIF contains the fragement amplified: i.e. in silico PCR results
	FWD contains the forward primer seq and the previous 8 nt on the template (to design tags..). Forward primers anneal to the template strand (anti-sense strand, 3'-5'), so they have the same sequence as the sense strand, which runs in the 5'-3', is the same sequence as mRNA, and is the conventionally reported strand. So fwd primers can be found as is in DB sequences.
	REV contains the reverse primer and the prior 8 nt on the template (in revcomp direction): Reverse primers anneal to the sense strand (coding strand, 5'-3') and are thus complementary to it. As DBs conventionally show the sense strand, the reverse primer is found only as rev_complement in the DB sequences.
	REVcomp contains the initial 5'-3' orientation of the reverse primer and previous 8 nt on the template
	"""
	parser = argparse.ArgumentParser()
	parser.add_argument('--db', '-db', required=True, help='Unaligned database fasta file')
	parser.add_argument('--forward_primer', '-f', required=True, help='Forward primer in 5-3 direction')
	parser.add_argument('--reverse_primer', '-r', required=True, help='The reverse primer in 5-3 direction (rev comp will be done..)')
	parser.add_argument('--sampling', '-s', required=True, type=int, help='Fraction (1-100%) of the database to sample for matching the primers')
	args = parser.parse_args()

	fwd = Seq(str(args.forward_primer))
	rev = Seq(str(args.reverse_primer))
	revComp = rev.reverse_complement()
	
	print(rev)
	print(revComp)

	ofwd = open(str(args.db.split(".fasta")[0] + "_FWD.fasta"),'w')
	orev = open(str(args.db.split(".fasta")[0] + "_REV.fasta"),'w')
	orevC = open(str(args.db.split(".fasta")[0] + "_REV_revComp.fasta"),'w')
	ampli = open(str(args.db.split(".fasta")[0] + "_AMPLIF.fasta"),'w')

	for i in SeqIO.parse(open(str(args.db)), 'fasta'):
		id = i.description
		seq = i.seq
		if random.randint(1, 100) < args.sampling:
			matchFwd = SeqUtils.nt_search(str(seq), str(fwd))
			matchRev = SeqUtils.nt_search(str(seq), str(revComp))
			if len(matchFwd) > 1 and len(matchRev) > 1:
				if len(matchFwd) > 2:
					print("more than a match FWD for: " + i.id)
				if len(matchRev) > 2:
					print("more than a match REV for: " + i.id)
				idxFwd = matchFwd[1]
				idxRev = matchRev[1]
				if idxFwd >3 and idxRev >3:
					FWD_out = str(seq)[idxFwd-8:idxFwd+len(fwd)]
					REV_out = str(seq)[idxRev:idxRev+len(revComp)+8]
					REV_outC = seq[idxRev:idxRev+len(revComp)+8].reverse_complement()
					if len(FWD_out) > 0 and len(REV_out) > 0:
						ofwd.write('>%s\n%s\n' % (id, FWD_out))
						orev.write('>%s\n%s\n' % (id, REV_out))
						orevC.write('>%s\n%s\n' % (id, REV_outC))
						ampli.write('>%s\n%s\n' % (id, str(seq)[idxFwd+len(fwd):idxRev]))


	ofwd.close()
	orev.close()
	orevC.close()
	ampli.close()




if __name__ == "__main__":
	main()

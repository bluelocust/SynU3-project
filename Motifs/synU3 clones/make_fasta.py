import os
import sys
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main(script, 
         rootdir='C:\Documents and Settings\qsun\My Documents\Year 4-1\Amgen\SynU3 Project\synU3 clones'):
 
	seq_as_nts = []
	seq_names = []
	for subdir, dirs, files in os.walk(rootdir):
		for file in files:
			# only process .seq files [!CLEAN]
			if '.seq' in file:
				f=open(file, 'r')
				lines=f.readlines()
				
				seq_names.append(file)
				
				seq = ''.join([x.rstrip('\n') for x in lines])
				seq_as_nts.append(seq)
				
				f.close()
	
	# make this code a module with parts of all_motifs.py for cleanliness
	# and conservation of import statements
	# 1) make Seq_Record
	# 2) export Seq_Record as FASTA file
	seq_NT = seq_generator(seq_as_nts, IUPAC.ambiguous_dna, seq_names) 
	output = open('exp_seqs.fasta', 'w')
	
	# export as FASTA
	SeqIO.write(seq_NT, output, 'fasta')
	output.close()
		
def seq_generator(data, datatype, dataID):
	'''process a list of sequences into SeqRecord'''
	for index in xrange(len(data)):
		yield SeqRecord(Seq(data[index], datatype), 
		                id = dataID[index] + '---' + 
						     str(len(data[index])),
		                description = '')

if __name__ == '__main__':
	main(*sys.argv)
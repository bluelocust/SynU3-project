import os
import sys
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main(script, ugly_MEME = 'prelim_seqs.txt'):
	pretty_print(ugly_MEME)

def pretty_print(file):	
	'''turns box-shapes into REAL linebreaks'''
	fp = open(file, 'r')
	lines = fp.readlines()
	fp.close()
	
	fp = open(file, 'w')	
	
	# specific sequence : list of motifs
	seq_M = {}
	for line in lines[1:]:
		newline = doWhatYouWant(line)
		fp.write(newline)
		
		vars = newline.split('	')

		# what's a p-value?
		if float(vars[5]) < .001:
			motif = vars[0]
			left = min(vars[2], vars[3])
			right = max(vars[3], vars[2])
			if left == vars[3]:
				motif += "'"
			motif_nts = vars[-1]
			
			motif_info = motif + '\t' + left + '\t' + right + '\t' + motif_nts
			
			if vars[1] in seq_M:
				seq_M[vars[1]].append(motif_info)
			else:
				seq_M[vars[1]] = [motif_info]
	fp.close()
	
	fp = open('MOTIF'+file, 'w')
	
	for key in seq_M.keys():
		fp.write('>' + key + '\n')
		for element in seq_M[key]:
			fp.write(element)
		fp.write('\n')
	fp.close()
	
def doWhatYouWant(line):
	return line[0:-1] + '\n'
	
if __name__ == '__main__':
	main(*sys.argv)
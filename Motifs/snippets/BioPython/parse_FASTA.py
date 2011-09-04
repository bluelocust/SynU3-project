from Bio import SeqIO
import sys

def main(script, seq_file = 'ls_orchid.fasta'):
	for seq_record in SeqIO.parse(seq_file, "fasta"):
		print type(seq_record)
	
if __name__ == '__main__':
	main(*sys.argv)
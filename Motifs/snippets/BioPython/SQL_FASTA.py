from Bio import SeqIO
import sys

def main(script, seq_file = 'ls_orchid.fasta'):
	orchid = SeqIO.index_db("orchid.idx", seq_file, "fasta")
	print "%i sequences indexed" % len(orchid)

# def main(script, seq_file = 'gbvrl1.seq'):
	# virus = SeqIO.index_db('gbvrl1.idx', seq_file, "genbank")
	# print "%i sequences indexed" % len(virus)

if __name__ == '__main__':
	main(*sys.argv)
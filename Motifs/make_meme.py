# Sam Sun
# -----------------------------------------------
# Processes a text file of LTR enhancer parts into a motif file
# that can be read by the open-source MEME Suite 
# -----------------------------------------------
# TO DO:
# *Pull background letter frequencies from Gillespie simulations
# -----------------------------------------------

from __future__ import division
import sys

# globals 
__motifs = {}
__IUPAC = {}

def main(script, LTR_motifs = 'LTR.txt', 
				 code = 'IUPAC.txt', 
				 file0 = 'min_MEME.txt'):
	# build motif dictionary
	num_motifs = build_motifs(open(LTR_motifs))
	# build IUPAC dictionary
	build_IUPAC(code)
		
	print '# of motifs:', num_motifs
	print 'motif dictionary:', __motifs
	
	# description
	print 'How many motifs would you like MEME Suite to find? Which ones?'
	
	# choose a motif, make a MEME
	i1 = raw_input("Number of motifs to look for (1, 2, ..., or ALL): ")
	i2 = []
	if i1 == 'ALL':
		i2 = __motifs.keys()
	else:
		for j in range(int(i1)):
			i2.append(raw_input("Enter a motif name: "))
	
	# load a starter file
	fp0 = open(file0, 'r')
	contents0 = fp0.read()
	fp0.close()

	# create a new file with starter contents
	fp = open('my_MEME' + '.txt', 'w')
	fp.write(contents0)
	
	# replace with frequencies from gillespie
	nt = 'ACGT'
	freq = [.25, .25, .25, .25]
	zip = ''.join([nt[n] + ' ' + str(freq[n]) + ' ' for n in xrange(4)])
	fp.write(zip + '\n')
	
	# create an nt-probability matrix for each desired motif
	for motif in i2:
		build_MEME(fp, motif)

	# cleanup
	fp.close()

def build_motifs(fp):
	'''creates dictionary of motif names and sequences'''
	global __motifs

	num = 0
	for line in fp:
		num += 1
		name, seq = process_line(line)
		__motifs[name] = seq
	return num
		
def build_IUPAC(code):
	'''creates IUPAC code for ambiguous nucleotides using a .txt'''
	global __IUPAC
	
	fp = open(code)

	for line in fp:
		nuc, nuc_options = process_line(line)
		__IUPAC[nuc] = list(nuc_options)
	
def process_line(line):
	'''keeps track of # of lines; also, processes each line
	to isolate name from sequence'''
	new_line = line.split()
	name, seq = new_line[0], ''.join(new_line[1:])
	
	return name, seq

def build_MEME(fp, motif):
	'''creates MEME text file for the desired motif'''
	
	fp.write(''.join(['\n', 'MOTIF ', motif]))
	fp.write(''.join(['\n', 'letter-probability matrix: ']))
	fp.write("alength= 4 w= %d nsites= 100 E= 0 \n" %len(__motifs[motif]))
	
	# use IUPAC to create matrix
	build_matrix(fp, __motifs[motif])

def build_matrix(filename, seq):	
	'''turn a nucleotide sequence into numbers'''
	ref = nt_probability()
	for n in seq:
		filename.write(ref[n])
		filename.write('\n')
	
def nt_probability():	
	''' create dictionary of 'NT' and probability #s for 'AGCT'
	    e.g. KEY = 'N'
	         VALUE = ' 0.250000  0.250000  0.250000  0.250000 ' '''   
	nt_matrix = {}
	for n in __IUPAC.keys():
		nt_matrix[n] = []
		for m in 'ACGT':
			if m in __IUPAC[n]: 
				nt_matrix[n].append(num_format(1/len(__IUPAC[n])))
			else:
				nt_matrix[n].append(num_format(0))
		# concatenates into one string of probabilities
		nt_matrix[n] = ''.join(nt_matrix[n])
	return nt_matrix
	
def num_format(num, decimals = 6):
	'''takes a float number and returns it in string format with a
	fixed # of decimals'''
	return ' ' + "%.6f" % num + ' '

if __name__ == '__main__':
	main(*sys.argv)
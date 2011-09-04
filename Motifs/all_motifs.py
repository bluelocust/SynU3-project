# Sam Sun
# Promoter Library
# -----------------------------------------------
# DONE:
# *Enumerates all possible motifs
# *Implements a basic gillespie algorithm and checked result with graphs
# *Gillespie function reflects diff. sizes (PEG)
# *Densitometry on DNA gels to get distribution of fragment sizes [MATLAB]
# *Exports sequences as nts and motifs in FASTA format
# -----------------------------------------------
# TO DO:
# 1) Separate file into modules...this is getting unwieldly
# 2) Create classes and methods instead of random functions
# 3.1) Convert arbitrary DNA units into molecules, etc.
# 3.2) Compare experimental vs. computational size distribution
# 4) Train HMMs on each motif separately, and apply to fragments
# 5) Make it water-tight
# 5.1) TDD
# 5.2) Error-handling
# -----------------------------------------------

from random import choice, random
import bisect
import matplotlib.pyplot as plt
import sys
import math
from Bio.Alphabet import IUPAC, Alphabet
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from my_hist import hist_and_fit
from random import choice

# globals 
__IUPAC = {}
__rxns = {}

def main(script, motifs = 'LTR.txt', code = 'IUPAC.txt'):
	global __num_motifs
	
	fp = open(motifs)
	
	# build IUPAC dictionary
	build_IUPAC(code)
	
	# build motif dictionary
	d = dict()
	for line in fp:
		motif, seq = process_line(line)
		d[motif] = seq
	
	# enumerate all possibilities for each motif
	all = {}
	for key in d:
		seq = d[key]
		
		# ID for specific, fully-determined sequence
		# entry 1: sequence position
		# entry 2: # of nucleotide choices
		seq_ID = []
		
		# find ambiguous nucleotides
		for pos in range(len(seq)):
			if seq[pos] in 'AGCT':
				pass
			else:
				seq_ID.append((pos, len(__IUPAC[seq[pos]])))

		all[key] = generate(seq, seq_ID, [], len(seq_ID))
	
	# build rxn code for Gillespie
	rxn_lookup()
	
	# Gillespie algorithm
	time = 0					
	parts = [{0:50}, {1:400}, {}]	# 100 HP, 200 parts, 0 HP-stuff-HP
	
	# WHEN YOU GET TO THE END, I ONLY CARE ABOUT DICT3
	##full_parts = []
	full_time = []					# why is this tracked?
	for n in range(2000):
		##full_parts.append(tuple(parts)) 
		full_time.append(time)
		
		dt, parts = gillespie(parts)
		
		# end condition; COMBINE; oh, I'm just passing up right now
		if dt == -1:
			break
		time += dt

	# aggregate HP-mer-HP fragments for motif-finding & visualization
	export_seq(parts[2], all)
	
	# aggregate DNA fragments for gel display
	DNA_gel = merge_dict(parts)
		
	plt.plot(DNA_gel.keys(), DNA_gel.values(), 'bo-')
	plt.show()
	
	#print '# of motifs:', len(d.keys())
	#print '# of sequences:', len(full)
	#print 'sequences unique?:', test_uniqueness(full)
	#print 'mean, S.D.:'
	
# Rather than using an IUPAC dictionary, BioPython built-in's could be
# used, but let's only opt for BioPython when creating FASTA files since
# FASTA is a bioinformatics standard	
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
	
def generate(seq, seq_ID, ACGT_sequences, depth):
	'''processes an ambiguous sequence into a list of all 
	possible sequences'''

	# extra case, if motif has not ambiguous nucleotides, e.g. Hairpin
	if depth == 0:
		return [seq]

	unknown_pos, num_choices = extract(seq_ID, depth-1)
	
	# base case
	# 1) fill in final unknown nucleotide in the sequence
	# 2) add 2-4 determined seq's to ACGT_sequences, and return 
	# ACGT_sequences
	if depth == 1:
		for m in range(num_choices):
			full_seq = fill_in(seq, unknown_pos, m)
			ACGT_sequences.append(''.join(full_seq))
		return ACGT_sequences

	# recursive case
	else:
		for n in range(num_choices):
			new_seq = fill_in(seq, unknown_pos, n)
			generate(new_seq, seq_ID, ACGT_sequences, depth-1)
		return ACGT_sequences
		
def extract(list_of_tuples, which_tuple):
	'''extracts a given tuple's entries from a list of tuples'''
	pair = list_of_tuples[which_tuple]
	a, b = pair[0], pair[1]
	return a, b

def fill_in(seq, unknown_pos, which_nuc):
	'''takes an ambiguous sequence and fills in one unknown 
	nucleotide with ACGT'''
	# copy sequence
	new_seq = list(seq)
	
	# fill in an unknown nucleotide
	nuc_choices = __IUPAC[new_seq[unknown_pos]]
	nuc = nuc_choices[which_nuc]
	new_seq[unknown_pos] = nuc
	
	return new_seq 

def test_uniqueness(all_seq):
	'''tests if all elements of a list are unique'''
	if len(set(all_seq)) == len(all_seq):
		return True
	else:
		return False

def gillespie(parts):
	'''implements Gillespie algorithm, based on # of monomers, and #
	of hairpins'''
	# rate constant
	k = .01
	# parts
	HPmer = parts[0]
	polymer = parts[1]
	HPmerHP = parts[2]
	
	# rxns
	X = sum(HPmer.values())
	Y = sum(polymer.values())
	
	a = [k*X*(X-1), k*X*Y, k*Y*(Y-1)]
	a_0 = sum(a)
	
	# end condition
	if a_0 == 0:
		return -1, parts
	
	# find time tau
	tau = -(1/a_0) * math.log(random())
	
	# choose a rxn 
	this_rxn = random()
	ctr = 0
	for n in range(len(a)):
		if this_rxn > ctr and this_rxn < ctr + a[n]/a_0:
			rxn = n
			break
		else:
			ctr += a[n]/a_0
	
	# change time, # of molecules
	dM = list(__rxns[rxn])
	return tau, rxn_update(parts, dM)

def rxn_lookup():
	'''creates reaction, or molecule-updating code for gillespie
	algorithm'''
	global __rxns
	
	# HPmer, polymer, HPmerHP
	__rxns[0] = [-2, 0, 1]		# require two HPmers
	__rxns[1] = [-99, -1, 0]
	__rxns[2] = [0, -100, 0]	# require two polymers

def rxn_update(parts, dM):
	'''updates randomly chosen molecules based on rxn type'''
	# track reactants
	lenR = []

	# delete reactant(s)
	for m in range(len(dM)):
		while dM[m] < 0:
			# randomly select a molecule size
			all_len = parts[m].keys()
			temp = all_len[weighted_rnd(parts[m].values())]
			
			# destroy one particle
			parts[m][temp] -= 1
			# housekeeping
			if parts[m][temp] == 0:
				del parts[m][temp]
			lenR.append(temp)
			
			if dM[m] == -99:
				dM[m] += 99
			dM[m] += 1

	# create product
	for n in range(len(dM)):
		if dM[n] == 1:
			lenP = sum(lenR)			
			if lenP in parts[n]:
				parts[n][lenP] += 1
			else:
				parts[n][lenP] = 1
	
	return parts

def weighted_rnd(weights):
	'''randomly selects an element from a list, with the chances of
	each length to be selected defined by weights'''
	rnd = random() * sum(weights)
	for i, w in enumerate(weights):
		rnd -= w
		if rnd < 0:
			return i

def merge_dict(parts):
	'''merges a list of dictionaries into a single dictionary
	--> + 1 to all keys in HPmer
	--> + 2 to all keys in HPmerHP'''
	#LISTS are more efficient at this point
	full = {}
	full.update(parts[0])
	for m in range(1, len(parts)):
		for k in parts[m].keys():
			if k in full:
				full[k] += parts[m][k]
			else:
				full[k] = parts[m][k]
	print full
	return full

def export_seq(SynU3s, ACGT_motifs):
	'''processes a dictionary of SynU3 sizes (# of motifs : # of seq's) 
	into two lists, seq_as_motifs and seq_as_nt, and exports two 
	FASTA-formatted files'''
	
	motifs = ACGT_motifs.keys();
	
	# destroy key-entry '0' since HP-HP is not an LTR
	if 0 in SynU3s:
		del SynU3s[0]
	else:
		pass
	
	# create a list whose entries correspond to U3 size for all SynU3s
	l = [[key]*SynU3s[key] for key in SynU3s.keys()]
	SynU3_sizes = [item for sublist in l for item in sublist]

	seq_as_motifs = []
	seq_as_nts = []
	# --- randomly choose motif, randomly choose sequence ---
	# 1) use a generator to step through the list, returning # of mers
	# 2) for each U3, choose a motif/seq and store in separate lists
	
	for seq in xrange(len(SynU3_sizes)):
		motif_string = [choice(motifs) for m in xrange(SynU3_sizes[seq])]
		seq_as_motifs.append(''.join([x + ' ' for x in motif_string]))
		seq_as_nts.append(''.join([choice(ACGT_motifs[x]) 
		                           for x in motif_string]))
	
	# seq_NT: feed into MEME or MochiView
	# seq_M: compare MEME output to seq_M, using same file format
	seq_M = seq_generator(seq_as_motifs, Alphabet())
	seq_NT = seq_generator(seq_as_nts, IUPAC.unambiguous_dna) 
	output1 = open('sim_motifs.fasta', 'w')
	output2 = open('sim_seqs.fasta', 'w')
	
	# export as FASTA
	SeqIO.write(seq_M, output1, 'fasta')
	SeqIO.write(seq_NT, output2, 'fasta')
	output1.close()
	output2.close()

def seq_generator(data, datatype):
	'''process a list of sequences into SeqRecord'''
	for index in xrange(len(data)):
		yield SeqRecord(Seq(data[index], datatype), 
		                 id = 'SynU3---#' + str(index), 
						 description = '')

# *test with MEME and MochiView
# *optimize the p-value, or enforce non-overlapping motifs, and compare 
# against seq_M
# *give two options (non-overlapping) vs. using p-values
						 
if __name__ == '__main__':
	main(*sys.argv)
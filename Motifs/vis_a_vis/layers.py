#!/usr/bin/env python

import sys
from wx.lib.pubsub import Publisher as pub

class Model(object):
	def __init__(self, script = None, sequences = 'MOTIFprelim_seqs.txt'):
		self.make_matrices(sequences)
		
	def make_matrices(self, sequences):
		fp = open(sequences)
		self.motif_M = []	# matrix of motif names
		self.len_M = []		# matrix of motif lengths
		self.start_M = []	# matrix of motif start positions
		
		#Turn into {seq_name : seq_nt's}, or use BioPython FASTA I/O
		self.ID = []		# array of sequence IDs
		for line in fp:	
			if line[0] == '>':
				seq_motifs = []
				seq_lengths = []
				seq_starts = []
				self.ID.append(line[1:-1])
			elif line[0].isalnum():
				motif, length, start = self.process_line(line)
				# Create two matrices, motif names and lengths
				seq_motifs.append(motif)
				seq_lengths.append(int(length))
				seq_starts.append(int(start))
			else:
				self.motif_M.append(seq_motifs)
				self.len_M.append(seq_lengths)
				self.start_M.append(seq_starts)
		
		# Fill in len_M with zeroes to prevent IndexError's
		max_num_motifs = len(max(self.len_M, key = len))
		for arr in self.len_M:
			for n in xrange(len(arr), max_num_motifs):
				arr.append(0)
		
		# UGLY UGLY UGLY
		for arr in self.start_M:
			for n in xrange(len(arr), max_num_motifs):
				arr.append(0)
				
		# Encapsulate...at some point
		for arr in self.motif_M:
			for n in xrange(len(arr), max_num_motifs):
				arr.append(" ")
		
		# does this go here?
		pub.sendMessage("ID CHANGED", self.ID)
		pub.sendMessage("Seq_Lengths CHANGED", self.len_M)
		pub.sendMessage("Seq_Names CHANGED", self.motif_M)
		
	def process_line(self, line):
		'''keeps track of # of lines; also, processes each line
		to isolate name from sequence'''
		new_line = line.split()

		print new_line
		
		name = new_line[0]
		length = int(new_line[2]) - int(new_line[1])
		start = new_line[1]
		
		return name, length, start
	
if __name__ == '__main__':
	model = Model(*sys.argv)
	print 'sequence data [motifs]', model.motif_M
#!/usr/bin/env python

# clean up imports
import sys
import matplotlib
from pylab import *
from matplotlib.colors import colorConverter	

def main(script, sequences = 'my_seqs_motifs.txt'):
	fp = open(sequences)
	motif_M = []	# matrix of motif names
	len_M = []		# matrix of motif lengths
	for line in fp:	
		if line[0] == '>':
			seq_motifs = []
			seq_lengths = []
		elif line[0].isalnum():
			motif, length = process_line(line)
			# Create two matrices, motif names and lengths
			seq_motifs.append(motif)
			seq_lengths.append(int(length))
		else:
			motif_M.append(seq_motifs)
			len_M.append(seq_lengths)
		
	# Fill in len_M with zeroes to prevent IndexError's
	max_num_motifs = len(max(len_M, key = len))
	for arr in len_M:
		for n in xrange(len(arr), max_num_motifs):
			arr.append(0)
	# Encapsulate...at some point
	for arr in motif_M:
		for n in xrange(len(arr), max_num_motifs):
			arr.append(" ")

	# Extract unique motifs
	keys = {}
	for motif in [item for sublist in motif_M for item in sublist]:
		keys[motif] = 1
	motifs_fin = keys.keys()
	
	# Get some pastel shades for the colours
	colours = get_colours(len(len_M[0]))
	colours.reverse()
	
	# Zip a dictionary of motifs and colors
	color_lookup = {}
	x = 0
	for RGB in get_colours(len(motifs_fin)):
		color_lookup[motifs_fin[x]] = RGB
		x += 1
	
	cols = len(len_M[0])
	height = .8	 # the height of the bars
	ind = arange(1, len(len_M)+1) - height/2 # the y locations for the groups

	xlim(0,max([sum(len_M[x]) for x in xrange(len(len_M))]) * 1.1)
	ylim(0,len(len_M)+1)

	xoff = array([0.0] * len(len_M)) # the left values for stacked bar chart

	for col in xrange(cols):
		#barh(bottom, width, **kwargs)
		
		len_col = [len_M[j].__getitem__(col) for j in xrange(len(len_M))]
		motif_col = [motif_M[j].__getitem__(col) for j in xrange(len(motif_M))]
		
		barh(ind, len_col, height, left=xoff, 
		color = [color_lookup[m] for m in motif_col])	

		#list comprehension :)
		for row in xrange(len(len_col)):
			if len_col[row] != 0:
				#text(x, y, s, fontdict=None, **kwargs)
				matplotlib.pyplot.text(xoff[row]+.5*len_col[row],row+1,
				motif_M[row].__getitem__(col),horizontalalignment = 'center',
				verticalalignment = 'center',fontsize=10)

		xoff = xoff + len_col
		
	show()
	
def process_line(line):
	'''keeps track of # of lines; also, processes each line
	to isolate name from sequence'''
	new_line = line.split()

	name, length = new_line[0], new_line[1]
	
	return name, length
	
def pastel(colour, weight=2.4):
	''' Convert colour into a nice pastel shade '''
	rgb = asarray(colorConverter.to_rgb(colour))
	# scale colour
	maxc = max(rgb)
	if maxc < 1.0 and maxc > 0:
		# scale colour
		scale = 1.0 / maxc
		rgb = rgb * scale
	# now decrease saturation
	total = sum(rgb)
	slack = 0
	for x in rgb:
		slack += 1.0 - x

	# want to increase weight from total to weight
	# pick x s.t.  slack * x == weight - total
	# x = (weight - total) / slack
	x = (weight - total) / slack

	rgb = [c + (x * (1.0-c)) for c in rgb]

	return rgb

def get_colours(n):
	""" Return n pastel colours. """
	base = asarray([[1,0,0], [0,1,0], [0,0,1]])

	if n <= 3:
		return base[0:n]

	# how many new colours to we need to insert between
	# red and green and between green and blue?
	needed = (((n - 3) + 1) / 2, (n - 3) / 2)

	colours = []
	for start in (0, 1):
		for x in linspace(0, 1, needed[start]+2):
			colours.append((base[start] * (1.0 - x)) +
						   (base[start+1] * x))

	return [pastel(c) for c in colours[0:n]]
	
if __name__ == '__main__':
	main(*sys.argv)
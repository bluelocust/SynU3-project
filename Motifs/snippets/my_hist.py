#!/usr/bin/env python
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import sys

def main(script):
	mu, sigma = 100, 15
	x = mu + sigma*np.random.randn(10000)
	
	hist_and_fit(x)

def hist_and_fit(x):
	# the histogram of the data
	n, bins, patches = plt.hist(x, 50, normed=1, facecolor='green', alpha=0.75)

	# add a 'best fit' line
	#y = mlab.normpdf( bins, mu, sigma)
	#l = plt.plot(bins, y, 'r--', linewidth=1)

	plt.xlabel('LTR Length (nt)')
	plt.ylabel('Probability')
	plt.title('SynU3 Fragments')
	#plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
	#plt.axis([40, 160, 0, 0.03])
	plt.grid(True)

	plt.show()

if __name__ == '__main__':
	main(*sys.argv)
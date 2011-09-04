import random
import pylab
aas = "ACDEFGHIKLMNPQRSTVWXY-"
fig = pylab.figure()
xmax = 100
ymax = 50
for iX in range(xmax):
    for iY in range(ymax):
        pylab.text(iX + 0.5,
                   iY + 0.5,
                   random.choice(aas),
                   family = 'sans-serif',
                   size = 6,
                   horizontalalignment = 'center',
                   verticalalignment = 'center')
pylab.xlim(0, xmax)
pylab.ylim(0, ymax)
pylab.savefig("test.png")
pylab.show()
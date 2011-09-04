import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)

t = np.arange(0.0, 5.0, 0.01)
s = np.cos(2*np.pi*t)
line, = ax.plot(t, s, lw=2)

# ax.annotate('local max', xy=(2, 1), xytext=(2, 1.5),
            # arrowprops=dict(facecolor='black', shrink=0.05),
            # )

for j in range(19):
	# for k in range(9):
		# plt.text(.5+.5*k,-1.8+.2*j,'motif' + str(k),
		# horizontalalignment = 'center',fontsize=10)
	plt.text(.5+.5*4,-1.8+.2*j,'motif   ' * 12,
	horizontalalignment = 'center',fontsize=10)


	
ax.set_ylim(-2,2)
plt.show()
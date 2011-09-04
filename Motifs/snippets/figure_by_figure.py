
import matplotlib.pyplot as plt
import time

plt.ion()

# figure by figure
for fig_number, data in [(0, [1,2,3]), (1, [10, 20, 30]), (0, [4, 5, 6]),
                         (1, [40, 50, 60])]:
    plt.figure(fig_number)
    plt.plot(data)
    plt.draw()
    plt.draw()
    time.sleep(5)
    plt.close(fig_number)

# only clear the axes without rebuilding the whole figure
fig = plt.figure()
ax = fig.add_subplot(111)
for data in [[1,2,3], [10, 20, 30], [4, 5, 6], [40, 50, 60]]:
    ax.cla() # clear the axes
    ax.plot(data)
    plt.draw()
    time.sleep(5)

plt.ioff()
plt.show()

import matplotlib.pyplot as plt

# left click: event.button = 1
# right click: event.button = 2

# Options
# *Use xdata, ydata to find what motif what pressed (lame)
# *Build off of PickEvent, but query the mouseevent

#class matplotlib.backend_bases.PickEvent(name, canvas, mouseevent, 
#artist, guiEvent=None, **kwargs)
def on_press(event):
    print 'you pressed', event.button, event.xdata, event.ydata

fig = plt.figure()
cid = fig.canvas.mpl_connect('button_press_event', on_press)

plt.show()

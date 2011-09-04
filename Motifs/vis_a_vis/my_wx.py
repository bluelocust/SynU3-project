"""
Sam Sun
F.W. Olin College of Engineering
"""
import os
import pprint
import random
import wx
from wx.lib.pubsub import Publisher as pub
import numpy as np
from layers import *

# The recommended way to use wx with mpl is with the WXAgg backend. 
import matplotlib
matplotlib.use('WXAgg')
from matplotlib.colors import colorConverter	
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_wxagg import \
	FigureCanvasWxAgg as FigCanvas, \
	NavigationToolbar2WxAgg as NavigationToolbar


class BarsFrame(wx.Frame):
	''' The main frame of the application '''
	title = 'Demo: wxPython with matplotlib'
	
	def __init__(self, ID, len_M, motif_M, start_M):
		wx.Frame.__init__(self, None, -1, self.title)
		
		self.create_menu()
		self.create_status_bar()
		self.create_main_panel()

		self.textbox.SetValue('')
		self.draw_figure(ID, len_M, motif_M, start_M)

	def create_menu(self):
		self.menubar = wx.MenuBar()
		
		menu_file = wx.Menu()
		m_expt = menu_file.Append(-1, "&Save plot\tCtrl-S", "Save plot to file")
		self.Bind(wx.EVT_MENU, self.on_save_plot, m_expt)
		menu_file.AppendSeparator()
		m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
		self.Bind(wx.EVT_MENU, self.on_exit, m_exit)
		
		menu_help = wx.Menu()
		m_about = menu_help.Append(-1, "&About\tF1", "About the demo")
		self.Bind(wx.EVT_MENU, self.on_about, m_about)
		
		self.menubar.Append(menu_file, "&File")
		self.menubar.Append(menu_help, "&Help")
		self.SetMenuBar(self.menubar)

	def create_main_panel(self):
		""" Creates the main panel with all the controls on it:
			 * mpl canvas 
			 * mpl navigation toolbar
			 * Control panel for interaction
		"""
		self.panel = wx.Panel(self)
		
		# Create the mpl Figure and FigCanvas objects. 
		# 5x4 inches, 100 dots-per-inch
		#
		self.dpi = 100
		self.fig = Figure((5 * (1+np.sqrt(5))/2, 5.0), dpi=self.dpi)
		self.canvas = FigCanvas(self.panel, -1, self.fig)
		
		# Since we have only one plot, we can use add_axes 
		# instead of add_subplot, but then the subplot
		# configuration tool in the navigation toolbar wouldn't
		# work.
		#
		
		self.axes = self.fig.add_subplot(111)
		
		# Bind the 'pick' event for clicking on one of the bars
		#
		self.canvas.mpl_connect('pick_event', self.on_pick)
		
		self.textbox = wx.TextCtrl(
			self.panel, 
			size=(200,-1),
			style=wx.TE_PROCESS_ENTER)
		self.Bind(wx.EVT_TEXT_ENTER, self.on_text_enter, self.textbox)
		
		self.findbutton = wx.Button(self.panel, -1, "Find Motif Module!")
		self.Bind(wx.EVT_BUTTON, self.on_find_button, self.findbutton)

		self.clearbutton = wx.Button(self.panel, -1, "Clear Modules")
		self.Bind(wx.EVT_BUTTON, self.on_clear_button, self.clearbutton)
		
		# Create the navigation toolbar, tied to the canvas
		self.toolbar = NavigationToolbar(self.canvas)
		
		#
		# Layout with box sizers
		#
		
		self.vbox = wx.BoxSizer(wx.VERTICAL)
		self.vbox.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
		self.vbox.Add(self.toolbar, 0, wx.EXPAND)
		self.vbox.AddSpacer(10)
		
		self.hbox = wx.BoxSizer(wx.HORIZONTAL)
		flags = wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
		self.hbox.Add(self.textbox, 0, border=3, flag=flags)
		self.hbox.Add(self.findbutton, 0, border=3, flag=flags)
		self.hbox.Add(self.clearbutton, 0, border =3, flag=flags)
		self.hbox.AddSpacer(30)
		self.vbox.Add(self.hbox, 0, flag = wx.ALIGN_LEFT | wx.TOP)
		
		self.panel.SetSizer(self.vbox)
		self.vbox.Fit(self)
	
	def create_status_bar(self):
		self.statusbar = self.CreateStatusBar()

	#do var's need to be self.vars?
	def draw_figure(self, ID, len_M, motif_M, start_M):
		'''Redraws the figure'''
		#Look at all "NFkB-1" or "NFAT/NFkB-1"
		#str = self.textbox.GetValue()
		
		self.len_M, self.motif_M, self.ID = len_M, motif_M, ID
		
		self.start_M = start_M
		
		# clear the axes and redraw the plot anew
		self.axes.clear()		
		
		# Extract unique motifs...
		keys = {}
		for motif in [item for sublist in self.motif_M for item in sublist]:
			keys[motif] = 1
		motifs_fin = keys.keys()
		
		# ...to zip a dictionary of motifs and colors
		color_lookup = {}
		x = 0
		for RGB in self.get_colours(len(motifs_fin)):
			color_lookup[motifs_fin[x]] = RGB
			x += 1
		
		rows = len(self.len_M)
		height = .8	 # the height of the bars
		ind = np.arange(1,len(self.len_M)+1)-height/2 # y positions for seq's
		
		# make it pretty
		# self.axes.set_xlim(0,
		# max([sum(self.len_M[x]) for x in xrange(len(self.len_M))]) * 1.1)
		# self.axes.set_ylim(0,len(self.len_M)+1)
		
		self.axes.set_title('Sequence Set')
		self.axes.set_xlabel('Length (Base Pairs)')
		self.axes.set_ylabel('Sequence #')
		
		for row in xrange(rows):
			lens = len_M[row]		# motif lengths for a given sequence
			motifs = motif_M[row]	# motif names for a given sequence
			
			starts = start_M[row]
			
			# motif positional offsets for a given sequence
			#offsets = [sum(lens.__getslice__(0,z)) for z in xrange(len(lens))]
		
			#barh(bottom, width, **kwargs)
			wtf = self.axes.barh([ind[row]] * len(lens), lens, height, 
			left = starts, color=[color_lookup[m] for m in motifs], 
			picker = True, gid=self.ID[row], alpha = .3)
			
			# add text
			# for col in xrange(len(lens)):
				# if lens[col] != 0:
					# #text(x, y, s, fontdict=None, **kwargs)
					# self.axes.text(starts[col]+.5*lens[col],row+1,
					# motifs[col],horizontalalignment = 'center',
					# verticalalignment = 'center',fontsize=10)

		self.canvas.draw()
			
	def pastel(self, colour, weight=2.4):
		''' Convert colour into a nice pastel shade '''
		rgb = np.asarray(colorConverter.to_rgb(colour))
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

	def get_colours(self, n):
		""" Return n pastel colours. """
		base = np.asarray([[1,0,0], [0,1,0], [0,0,1]])

		if n <= 3:
			return base[0:n]

		# how many new colours to we need to insert between
		# red and green and between green and blue?
		needed = (((n - 3) + 1) / 2, (n - 3) / 2)

		colours = []
		for start in (0, 1):
			for x in np.linspace(0, 1, needed[start]+2):
				colours.append((base[start] * (1.0 - x)) +
							   (base[start+1] * x))

		return [self.pastel(c) for c in colours[0:n]]
	
	def on_clear_button(self, event):
		# clean! why am I passing variables to myself?
		self.draw_figure(self.ID, self.len_M, self.motif_M)
		
	# useful for re-drawing graphs:
	# 1) layer-to-layer, i.e. 10e6 to 10e2
	# 2) highlighting motifs or architectures
	def on_find_button(self, event):
		# any motifs highlighted? (gui list of lists object)
		# if no, then open a dialog asking so-and-so to click motifs
		# if yes, ask model to process stuff...err...and...?
		
		
		# re-draw figure with red outlines
		
		
		#self.draw_figure()
		print 'clicked on button'

	#use a single binding?		
	def on_text_enter(self, event):
		#self.draw_figure()
		print 'entered text'
	
	def on_pick(self, event):
		# The event received here is of the type
		# matplotlib.backend_bases.PickEvent
		#
		# It carries lots of information, of which we're using
		# only a small amount here.
		# 
		#box_points = event.artist.get_bbox().get_points()

		if event.mouseevent.button == 1:		
			seq = event.artist.get_gid()
			msg = "Sequence:\n %s " % seq
			#msg = "You clicked on the motif:\n %s \n %s" % box_points

			#option of saving as text?
			dlg = wx.MessageDialog(
				self, 
				msg, 
				"Sequence Selected",
				wx.OK | wx.ICON_INFORMATION)

			dlg.ShowModal() 
			dlg.Destroy()
		elif event.mouseevent.button == 3:
			# Positional data of Rectangle
			rect = event.artist
			x, y = rect.get_xy()
			
			# Draw red outlines
			self.red_box([(x, y, rect.get_width(), rect.get_height())])
			
			# Store 'clicked' motifs
			# *** very necessary for finding other modules
			# *** prevent highlighting across sequences
			
	'''accepts a list of tuples'''
	def red_box(self, pos):
		for x, y, w, h in pos:
			self.axes.barh(y, w, left = x, height = h, facecolor='none',
			linestyle='dashed',linewidth=3,edgecolor='r')
		self.canvas.draw()

	def on_save_plot(self, event):
		file_choices = "PNG (*.png)|*.png"
		
		dlg = wx.FileDialog(
			self, 
			message="Save plot as...",
			defaultDir=os.getcwd(),
			defaultFile="plot.png",
			wildcard=file_choices,
			style=wx.SAVE)
		
		if dlg.ShowModal() == wx.ID_OK:
			path = dlg.GetPath()
			self.canvas.print_figure(path, dpi=self.dpi)
			self.flash_status_message("Saved to %s" % path)
		
	def on_exit(self, event):
		self.Destroy()
		
	def on_about(self, event):
		msg = """ DNA visualization tool for the Schaffer Lab:
		
		 * Use the matplotlib navigation bar
		 * Save the plot to a file using the File menu
		 * Click on a bar to receive an informative message
		"""
		dlg = wx.MessageDialog(self, msg, "About", wx.OK)
		dlg.ShowModal()
		dlg.Destroy()
	
	def flash_status_message(self, msg, flash_len_ms=1500):
		self.statusbar.SetStatusText(msg)
		self.timeroff = wx.Timer(self)
		self.Bind(
			wx.EVT_TIMER, 
			self.on_flash_status_off, 
			self.timeroff)
		self.timeroff.Start(flash_len_ms, oneShot=True)
	
	def on_flash_status_off(self, event):
		self.statusbar.SetStatusText('')

class Controller:
	def __init__(self, app):
		self.mod = Model()
		
		# initialize main frame with current data
		self.view = BarsFrame(self.mod.ID,self.mod.len_M,
							  self.mod.motif_M,self.mod.start_M)
		
		# subscribe to all "DATA CHANGED" messages from the Model
		pub.subscribe(self.ID_Changed, "ID CHANGED")
		pub.subscribe(self.Lengths_Changed, "Seq_Lengths CHANGED")
		pub.subscribe(self.Names_Changed, "Seq_Names CHANGED")
		
		self.view.Show()
	
	# def DataPlease(self):
		# '''asks the model to update its data registry (using file I/O)'''
	
	def ID_Changed(self, message):
		'''pubsub will call this handler as messages are sent from the model'''
		#self.view.drag_figure(message.data)
		print 'do some stuff'
		
	def Lengths_Changed(self, message):
		print 'hello goodbye'
	
	def Names_Changed(self, message):
		print 'hi hi hi'

if __name__ == '__main__':
	app = wx.PySimpleApp()
	controller = Controller(app)
	app.MainLoop()

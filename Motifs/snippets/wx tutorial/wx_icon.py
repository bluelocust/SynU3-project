# experiment with wxPython's wx.Icon() and SetIcon()
# you can use your own icon or use the base64 encoded icon embedded in the code
# tested with Python24 and wxPython26     vegaseat     20oct2005

import wx
import base64

class MyFrame(wx.Frame):
    """make a frame, inherits wx.Frame"""
    def __init__(self):
        # create a frame/window, no parent, default to wxID_ANY
        wx.Frame.__init__(self, None, wx.ID_ANY, 'wxPython custom icon',
            pos=wx.Point(300, 150), size=wx.Size(300, 350))
        self.SetBackgroundColour('green')
        
        # this is the icon file generated from the base64 text
        # or pick an icon file (.ico) you have ...
        iconFile = "zzz77.ico"
        icon1 = wx.Icon(iconFile, wx.BITMAP_TYPE_ICO)
        self.SetIcon(icon1)

        # show the frame
        self.Show(True)
        
        
# this is a base64 encoded icon image
ico1_b64 = \
"""AAABAAEAICAQAAAAAADoAgAAFgAAACgAAAAgAAAAQAAAAAEABAAAAAAAgAIAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAACAAACAAAAAgIAAgAAAAIAAgACAgAAAgICAAMDAwAAAAP8AAP8AAAD//wD/AAAA
/wD/AP//AAD///8AqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqoACqqqqqqqqqqqqqqqqg
AAqqqqqqqqqqqqqqqqoAAKqqqqqqqqqqqqqqqqAAAAAKqqqqqqqqqqqqoAAKoAAAAKqqqqAACqqq
AAqqqqAAAAAKAKoAAACgAKqqqqqgAAAAAKqgAAAACqqqqqqqqgAAAAqqoAAAAAqqqqqqqqAAoAAK
qgAAAAAKqqqqoAAKqqqqCqqqoAAACqqgAAqqqqqqqgoAqqqgAAoACgAACqqqqqoAqqqqqgAAqqAA
AACqqqAACqqqAACqCqqgAAAACgAKqgqqAKqqqgqqoAAAAACqqqoKqqqqqqoKqqoAAAAKqqqqCgCq
qqqqCqoAAAAACqqqqgCqqqqqqgAAqqqqqgqqoAAKqqqqoAAKqqqqqqoKAAqqCqqqAAqqCqqqqqqg
AKqqqgqqAKqqqgqqqqAACgqqqqoKAKqqqqoAAAAKqqoKqqqgAKqqqqqqAAAKqqqqCqoACqqqqqqg
AAAAAAqqqgAAqqqqqqoACqqgAAAAAAAKqqqqqqoAqqqqqgAAAAAKCqqqqqqqqqqqqqoAAAAAAAqq
qqqqqqqqqqoACqoAAAAKqqqqqqqqqqAAqqqqqgAACqqqqqqqqgAKqqqqqqqgAAqqqqqqqgCqqqqq
qqqqqqqqqqqqqqoAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA==
"""
# convert back to binary to save it
ico1 = base64.b64decode(ico1_b64)
fout = open("zzz77.ico","wb")
fout.write(ico1)
fout.close()

application = wx.PySimpleApp()
# call class MyFrame
window = MyFrame()
# start the event loop
application.MainLoop()
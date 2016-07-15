#!/usr/bin/env python
import numpy as np
import Tools as misctools
import matplotlib.pyplot as plt
import Tomo_slice_manipulation_tools as slicetools



class PointBrowser:
    """
    Controls user interaction with the Mapper plot
    """
    def __init__(self):
        self.dragging = None
        self.line = None
        self.points = None
        self.linelats = None
        self.linelons = None
        self.drawinglines = None
        self.arrow = None

    def addobjs(self,mapobj,figobj):

      '''Add a map object so that we can plot on it'''

      self.mapobj = mapobj
      self.figobj = figobj

    def updatelinecoords(self):

      '''reset the box coordinates to none'''

      self.linelats = None
      self.linelons = None
      self.points = None
      self.line = None
      self.drawinglines = None

    def motion(self,event):

      '''define what happens when the user moves the mouse over the canvas'''

      lon = event.xdata
      lat = event.ydata

      if self.drawinglines:

        if self.dragging:

          print 'Dragging!'
          #print lon,lat

        if self.line:
           self.line[0].remove()
        if self.points:
          self.points[0].remove()
        if self.arrow:
          self.arrow.remove()

        self.linelats = [self.startlat,lat]
        self.linelons = [self.startlon,lon]
        xevent,yevent = self.mapobj(self.linelons,self.linelats)

        self.arrow = plt.arrow(xevent[0],yevent[0],xevent[1]-xevent[0],yevent[1]-yevent[0],fc="k", ec="k", linewidth = 2, head_width=3, head_length=3)
        self.line = self.mapobj.plot(xevent,yevent,'r-',linewidth=2,alpha=0.9)
        self.points = self.mapobj.plot(xevent,yevent,'ko')
        self.figobj.canvas.draw()

    def releasepick(self,event):

      '''define what happens when the user releases the cursor'''

      lon = event.xdata
      lat = event.ydata

      if self.dragging:

        self.dragging = None


    def returnboxcoords(self):

      '''return box coordinates to user'''

      if self.boxlats:

        return self.boxlats,self.boxlons


    def onpick(self, event):
      '''
      Define what happens when the user presses the cursor. If the user chooses to make a cross section
      the Ritsema tools are automatically called and slices though the various models are plotted. A stack of 
      tomography slices is also made and contours associated with the slab are traced 
      '''

      if self.drawinglines:

        if event.button == 3:
          print 'Create profile!'

          s = str(raw_input('Create a profile? [Y/N]: '))

          if s == 'Y':
            print 'making profile!!'
            #print self.linelats
            #print self.linelons

            midlon = self.linelons[0]
            midlat = self.linelats[0]
            lon1 = self.linelons[1]
            lat1 = self.linelats[1]

            #Determine azimuth from the start point to the selected location
            azimuth = misctools.coords_for_profile(self.linelons[0],self.linelats[0],self.linelons[1],self.linelats[1])

            #print 'Getting ready to run Ritsema codes with midlon/midlat = %g/%g and azimuth of %g' %(midlon,midlat,azimuth)
            #misctools.Ritsema_180_sections(midlon,midlat,azimuth)

            print 'Getting ready to run Becker codes with midlon/midlat = %g/%g and azimuth of %g' %(midlon,midlat,azimuth)
            misctools.Becker_slice(midlon, midlat, lon1, lat1, azimuth)

          else:
            print 'No profile selected. Continue'

          self.drawinglines = None

        if event.button == 1:
          print 'Draw desired profile now. Click mouse 3 when done with your profile'
          self.dragging = True

          lon = event.xdata
          lat = event.ydata
          self.startlon = lon
          self.startlat = lat
      
      elif event.button == 3:

        print 'proceeding to draw lines!'
        self.drawinglines = True

      elif event.button == 1:

        print 'No drawing lines!'
        print 'Use mouse 3 to start or stop drawing a profile'
        print 'Use mouse 1 to '

      else:

        print 'Not a recognized command!'





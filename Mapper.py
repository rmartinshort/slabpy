#!/usr/bin/env python 

import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from Browser import PointBrowser
from Tools import read_faults
import sys

matplotlib.rcParams.update({'font.size': 12})

Browse = PointBrowser()

class DrawableMap:

	def __init__(self,slabcat=None):
		'''Initialize the map'''

		print '\n-----------------------------------\n'
		print 'INTERACTIVE MAP START'
		print '\n-----------------------------------\n'


		self.figure = plt.figure(facecolor='white',figsize=(25,8))

		#Connect the figure canvas to user's mouse motions/key presses
		self.figure.canvas.mpl_connect('button_press_event',Browse.onpick)
		self.figure.canvas.mpl_connect('motion_notify_event', Browse.motion)
		self.figure.canvas.mpl_connect('button_release_event',Browse.releasepick)


		self.a = self.figure.add_subplot(111)

		print 'Drawing map....'

		self.map = Basemap(projection='cyl',ax=self.a,lat_0=0,lon_0=90,resolution ='l',llcrnrlon=0,llcrnrlat=-89,urcrnrlon=360,urcrnrlat=89)
		#self.map.arcgisimage(service='NatGeo_World_Map',verbose=False,xpixels=10000)
		self.map.drawparallels(np.arange(-90,90,30),labels=[1,1,0,0],linewidth=0.5,fontsize=10)
		self.map.drawmeridians(np.arange(0,360,30),labels=[0,0,0,1],linewidth=0.5,fontsize=10)
		self.map.fillcontinents()
		#self.map.shadedrelief()
		self.map.drawcountries()
		#self.map.shiftgrid()

		if isinstance(slabcat,list):
			print 'Recognized slabcat!'
		else:
			print 'A correct slab catalog was not entered!'
			sys.exit(1)

		for slab in slabcat:
			bbox = slab.Fukaoslab_details['Bounds']

			print bbox

			if len(bbox) == 4:

				minlon = bbox[0]
				minlat = bbox[1]
				maxlon = bbox[2]
				maxlat = bbox[3]

				boxlats = [minlat,maxlat,maxlat,minlat,minlat]
				boxlons = [minlon,minlon,maxlon,maxlon,minlon]

			elif len(bbox) == 8:

				crn1lon = bbox[0]
				crn1lat = bbox[1]

				crn2lon = bbox[2]
				crn2lat = bbox[3]

				crn3lon = bbox[4]
				crn3lat = bbox[5]

				crn4lon = bbox[6]
				crn4lat = bbox[7]

				boxlats = [crn1lat,crn2lat,crn3lat,crn4lat,crn1lat]
				boxlons = [crn1lon,crn2lon,crn3lon,crn4lon,crn1lon]

				minlon = crn1lon
				minlat = crn1lat



			xevent,yevent = self.map(boxlons,boxlats)

			self.map.plot(xevent,yevent,'r-',linewidth=1,alpha=0.9)
			self.a.text(minlon, minlat, '%s' %slab.name, style='italic',bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})

		#print 'Drawing plate boundaries.....'

		self.faults = read_faults('data/plates.xy')
		for i in self.faults:
			faults_lons = self.faults[i][0]
			faults_lats = self.faults[i][1]
			x,y = self.map(faults_lons, faults_lats)
			self.map.plot(x,y,'b--',linewidth=1.0)

		Browse.addobjs(self.map,self.figure)

	def showmap(self):

		plt.show()

if __name__ == '__main__':

	mymap = DrawableMap()
	mymap.showmap()






		#self.map.arcgisimage(service='NatGeo_World_Map',verbose=False,xpixels=10000)
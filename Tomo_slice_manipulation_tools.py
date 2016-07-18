#!/usr/bin/env python

#Tools for manipulation and comparision of tomography cross sections

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import glob
import os
from skimage import measure


def gmt_norm_stack(files,plot=True):
	'''Uses gmt to normalize a set of GRD files and then stack them - the aim is to see where they agree most'''

	#Normalize the files

	for grdfile in files:

		os.system('gmt grdmath %s NORM = %s' %(grdfile,grdfile[:-3]+'NORM.grd'))

	basefilename = files[0][:-5]
	normedfiles = glob.glob(basefilename+'*'+'NORM.grd')
	print normedfiles

	#Start summing the normed files

	os.system('gmt grdmath %s %s ADD = %s' %(normedfiles[0],normedfiles[1],basefilename+'SUM.nc'))

	for normedfile in normedfiles[2:]:
		print 'Adding %s to stack' %normedfile
		os.system('gmt grdmath %s %s ADD = %s' %(basefilename+'SUM.nc',normedfile,basefilename+'SUM.nc'))

	#plot a gmt figure of the stacked tomography
	if plot == True:
		os.system('./stackmap.sh %s' %basefilename+'SUM.nc')

	os.system('rm output_slices/*.grd')

	return basefilename+'SUM.nc'

def plotfigs(model=None,axobj=None,i=0,infile=None,stack=True,datatracker=None,name=None,axz=None,quakes=None):

	'''
	Plotting command for the becker slice'''

	if not axobj:
		fig = plt.figure(facecolor='white',figsize=(10,6))
		datatracker = {}
		axobj = fig


	infile = Dataset(infile, model='r')

	depths = infile.variables['y'][:]
	lengths = infile.variables['x'][:]
	data = infile.variables['z'][:][:]

	infile.close()

	datatracker[infile] = {'data' : data}
	subplotno = '11%s' %str(i)

	if i==0:
		ax = axobj.add_subplot(321)
	elif i==1:
		ax = axobj.add_subplot(322)
	elif i==2:
		ax = axobj.add_subplot(323)
	elif i==3:
		ax = axobj.add_subplot(324)
	elif i==4:
		ax = axobj.add_subplot(325)

	frame1 = plt.gca()


	image = ax.imshow(data, interpolation='nearest',aspect='auto',cmap=plt.cm.seismic_r)


	#plot on the boundaries:
	#This assumes that the bottom of the map corresponds to 2000km depth! Its possible to change this, so 
	#be careful!
	thousandkm = int(len(depths)/2.0)
	sixsixtykm = int(len(depths)*6.6/20.0)


	ax.plot([0,len(lengths)],[thousandkm,thousandkm],'k--',label='1000km')
	ax.plot([0,len(lengths)],[sixsixtykm,sixsixtykm],'k-',label='660km')

	#Get earthquake information if available and scale to the dimensions of matrix
	#This is working to some extent, but adds to the run time and may fail in certain cases
	#In the end it may not be that useful - may be better to plot earthquakes on a different plot,
	#with just the slab1.0 contours, for example

	if quakes:
		if (len(quakes[0]) > 1):

			# quakepoints = np.array(quakes[0])
			# quakedepths = -np.array(quakes[1])
			# maxlen = float(quakes[2])

			# depthinc = len(depths)/2000.0
			# xinc = len(lengths)/maxlen

			# print depthinc,xinc

			# quakedepths = quakedepths*depthinc
			# quakepoints = len(lengths) - quakepoints*xinc

			# print quakedepths
			# print quakepoints

			# ax.plot(quakepoints,quakedepths,'g.',linewidth=0.4)

			if i == 2:
				ax.text(30, 200, 'Max earthquake depth: %g km' %-min(quakes[1]), style='italic',bbox={'facecolor':'red', 'alpha':0.5, 'pad':10})

	ax.set_title('%s' %model)

	ax.set_ylim([len(depths),0])
	ax.set_xlim([0,len(lengths)])

	cbar_ax = axobj.add_axes([0.88, 0.03, 0.03, 0.25])
	colors = image.set_clim(-1.5,1.5)
	axobj.colorbar(image,cbar_ax)

	#Remove the numbers on the x and y axes - these are distrating

	frame1.axes.get_xaxis().set_visible(False)
	frame1.axes.get_yaxis().set_visible(False)

	return axobj,datatracker




def tracestackedslab(infile,depthinc=5,llinc=((6371*np.pi)/360),cval=0.5):
	'''
	Takes a netcdf file containing a stacked, normed profiles and attempts to contour what might be a slab - can use this to estimate dip, profile etc
	The values of depthinc and llinc are defaults from the Ritsema code, which creates slices over angles of 180 degrees
	and with a depth increment of 5 km

	This produces a map showing the slice and contours in stacked velocity perturbation at a chosen level (typically 0.5)
	'''

	Mantlebase = 2895

	infile = Dataset(infile, model='r')

	filevariables = infile.variables.keys()

	#Get data from the netCDF file
	depths = infile.variables['y'][:]
	lengths = infile.variables['x'][:]
	data = infile.variables['z'][:][:]

	infile.close()

	#print np.shape(data)
	#print np.shape(lengths)
	#print np.shape(depths)

	#Use image processing suite to find contours
	contours = measure.find_contours(data,cval)

	#Various plotting commands to produce the figure
	fig, ax = plt.subplots()

	thousandkm = int((Mantlebase-1000)/depthinc)
	sixsixtykm = int((Mantlebase-660)/depthinc)

	plt.set_cmap('jet_r')

	ax.imshow(data, interpolation='linear',aspect='auto')

	ax.plot([0,len(lengths)],[thousandkm,thousandkm],'k--',label='1000km')
	ax.plot([0,len(lengths)],[sixsixtykm,sixsixtykm],'k-',label='660km')

	for n, contour in enumerate(contours):
		if n == 0:
			ax.plot(contour[:, 1], contour[:, 0], 'r-', linewidth=2,label='contour at %g' %cval)
		else:
			ax.plot(contour[:, 1], contour[:, 0], 'r-', linewidth=2)


	ax.set_ylim([0,len(depths)])
	ax.set_xlim([len(lengths),0])
	ax.set_title('Stacked slab image from netCDF')
	plt.xlabel('Cross section x increment')
	plt.ylabel('Cross section depth increment')
	plt.legend(loc='best')

	#plt.gca().invert_yaxis()
	#plt.gca().invert_xaxis()

	plt.show(block=False)






if __name__ == '__main__':

	grdfiles = glob.glob('uit*.grd')
	print grdfiles

	stackedfile = gmt_norm_stack(grdfiles,plot=True)
	tracestackedslab(stackedfile)


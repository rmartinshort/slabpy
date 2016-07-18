#!/usr/bin/env python

from __future__ import print_function
import numpy as np
from scipy.io import netcdf
from slab import Slab
from catalog import Catalog
import glob
import matplotlib.pyplot as plt


def read_slabdata(slabcat=None):
	'''
	Read in the slab1.0 data as a catalog of slabs object. The argument options allows one to 
	append these slabs to an existing catalog
	'''

	slabdir='data/slab1.0/'

	#Names given to the slab 1.0 slabs
	slabbounds = {}
	slabbounds['Alaska'] = {'Filename' : 'alu_slab1.0_clip.grd'}
	slabbounds['Cascadia'] = {'Filename' :'cas_slab1.0_clip.grd'}
	slabbounds['Izu-Bonin'] = {'Filename' : 'izu_slab1.0_clip.grd'}
	slabbounds['Mexico'] = {'Filename' : 'mex_slab1.0_clip.grd'}
	slabbounds['Tonga'] = {'Filename' : 'ker_slab1.0_clip.grd'}
	slabbounds['Japan'] = {'Filename' : 'kur_slab1.0_clip.grd'}
	slabbounds['Philippines'] = {'Filename' : 'phi_slab1.0_clip.grd'}
	slabbounds['Ryukyu'] = {'Filename' : 'ryu_slab1.0_clip.grd'}
	slabbounds['Vanautu'] = {'Filename' : 'van_slab1.0_clip.grd'}
	slabbounds['Scotia'] = {'Filename' : 'sco_slab1.0_clip.grd'}
	slabbounds['Solomon'] = {'Filename' : 'sol_slab1.0_clip.grd'}
	slabbounds['South_America'] = {'Filename' : 'sam_slab1.0_clip.grd'}
	slabbounds['Sumatra'] = {'Filename' : 'sum_slab1.0_clip.grd'}

	#try:
	#	os.chdir(slabdir)
	#except:
	#	print('Data directory %s does not exist!' %slabdir)

	slab_ncs = glob.glob(slabdir+'*.grd')

	#Create a catalog object for the slab1.0 slabs
	if slabcat:
		Slabs = slabcat
	else:
		Slabs=Catalog()

	for slab in slab_ncs:

		infile = netcdf.netcdf_file(slab, 'r')

		filevariables = infile.variables.keys()

		lats = infile.variables[filevariables[0]][:]
		lons = infile.variables[filevariables[1]][:]

		#This is a lat x lon array
		depths = infile.variables[filevariables[2]][:]
		infile.close()

		minlat = min(lats)
		maxlat = max(lats)
		minlon = min(lons)
		maxlon = max(lons)

		boundingbox = [minlon,maxlon,minlat,maxlat]

		for name in slabbounds:

			if slabbounds[name]['Filename'] == slab.split('/')[-1].strip():

				#create a Slab object 
				newslab = Slab(name)
				slabbounds[name]['Lat_array'] = lats
				slabbounds[name]['Lon_array'] = lons
				slabbounds[name]['Depth_array'] = depths
				slabbounds[name]['Bounding_box'] = boundingbox

				#Add sub-dicionary to the slab object
				newslab.add_slab1_details(slabbounds[name])
				Slabs.append(newslab)
				break

	return Slabs,slabbounds




if __name__ == '__main__':

	Slabs, slabbounds = read_slabdata()
	
	slist = Slabs.slabs

	#Testing slab plotting ability

	for slab in slist:

		axisobj, mapobj = slab.map_slab1()
		plt.show()




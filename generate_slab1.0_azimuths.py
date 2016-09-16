#!/usr/bin/env python


#RMS July 2016

#Generate table of subduction parameters using information from the Syracuse paper
#This was a 'one time' code to generate a list of slab strikes by interpolating the slab1.0 database at each of the points
#from the Syracuse paper

import os
import glob
import sys
import math as mt
from scipy.io import netcdf
import Tools as tomotools
import numpy as np

#Path to models

#--------------------------------------------
#Important: The paths to these data directories must be updated on a new system
# the becker directory also includes all of T.Becker's extraction code too
#--------------------------------------------
global becker_data_directory
becker_data_directory = '/Volumes/TOSHIBA EXT/Workshops/CIDER_2016/Slab_group2/datasets/Profile_plotting'



def generate_syracuse_latlon(infile='Syracuse_latlons.dat'):

	datafile = open(infile,'r')
	lines = datafile.readlines()
	datafile.close()

	slabpoints = []

	for line in lines:
		vals = line.split()
		name = vals[0]
		lon = vals[1]
		lat = vals[2]
		slabpoints.append([name,lon,lat])

	return slabpoints


def read_slabdata():
	'''
	Read the boundaries of the slabs1.0 data 
	'''

	#Names given to the slab 1.0 slabs
	slabbounds = {}
	slabbounds['Alaska'] = {'Filename' : 'alu_slab1.0_strclip.grd'}
	slabbounds['Cascadia'] = {'Filename' :'cas_slab1.0_strclip.grd'}
	slabbounds['Izu-Bonin'] = {'Filename' : 'izu_slab1.0_strclip.grd'}
	slabbounds['Mexico'] = {'Filename' : 'mex_slab1.0_strclip.grd'}
	slabbounds['Tonga'] = {'Filename' : 'ker_slab1.0_strclip.grd'}
	slabbounds['Japan'] = {'Filename' : 'kur_slab1.0_strclip.grd'}
	slabbounds['Philippines'] = {'Filename' : 'phi_slab1.0_strclip.grd'}
	slabbounds['Ryukyu'] = {'Filename' : 'ryu_slab1.0_strclip.grd'}
	slabbounds['Vanautu'] = {'Filename' : 'van_slab1.0_strclip.grd'}
	slabbounds['Scotia'] = {'Filename' : 'sco_slab1.0_strclip.grd'}
	slabbounds['Solomon'] = {'Filename' : 'sol_slab1.0_strclip.grd'}
	slabbounds['South_America'] = {'Filename' : 'sam_slab1.0_strclip.grd'}
	slabbounds['Sumatra'] = {'Filename' : 'sum_slab1.0_strclip.grd'}


	slab_ncs = glob.glob('*strclip.grd')

	slabpoints = generate_syracuse_latlon()

	#For every slab1.0 netcdf file, get the boundaries and check which of the Syracuse points falls within them

	for slab in slab_ncs:


		infile = netcdf.netcdf_file(slab, 'r')

		filevariables = infile.variables.keys()

		lats = infile.variables[filevariables[0]][:]
		lons = infile.variables[filevariables[1]][:]

		infile.close()

		minlat = min(lats)
		maxlat = max(lats)
		minlon = min(lons)
		maxlon = max(lons)

		#Convert eveything to 0-->360 longitude 

		if minlon < 0:
		 	maxlontmp = 360+maxlon
			minlontmp = 360+minlon
			boundingbox = [minlontmp,maxlontmp,minlat,maxlat]
		else:
			boundingbox = [minlon,maxlon,minlat,maxlat]

		print '-------------------------------------\n'

		print 'Dealing with slab1.0 slab %s' %slab
		print '%s has boundingbox %s' %(slab,boundingbox)

		print '\n-------------------------------------'

		for name in slabbounds:

			if slabbounds[name]['Filename'] == slab:

				slabbounds[name]['Bounding_box'] = boundingbox
				slabbounds[name]['Syracuse_points'] = []

				for syracusepoint in slabpoints:

					lon = float(syracusepoint[1])
					lat = float(syracusepoint[2])

					#Convert eveything to 0-->360 longitude 

					if lon < 0:
						lon = 360+lon

					#Check if the coordinates are within the region defined by the slab 1.0

					if (boundingbox[0] <= lon <= boundingbox[1]) and (boundingbox[2] <= lat <= boundingbox[3]):

						slabbounds[name]['Syracuse_points'].append(syracusepoint)

				break

	return slabbounds

def generate_grdtrackpoints(slabbounds):
	'''
	Takes a slabbounds dictionary and extracts slab srike from slab1.0 database
	'''

	#Make a record of Syracuse points that are not avaialble in the slab1.0 data, for some reason
	NaNsfile = open('Undefined_points.dat','w')
	SyracuseFile = open('Syracuse_points_bathymetry.dat','w')
	SyracuseFile.write('Slice_name Start_lon Start_lat End_lon End_lat Strike Azimuth\n')

	for slab in slabbounds:

		syracusepoints = slabbounds[slab]['Syracuse_points']
		filename = slabbounds[slab]['Filename']

		if len(syracusepoints) > 0:

			ofilename='tmp.%s.dat' %slab
			ofile = open(ofilename,'w')

			#Generate input file for grdtrack program
			for element in syracusepoints:

				name = element[0]
				lon = float(element[1])
				lat = float(element[2])

				#Cascadia slab1.0 has coodinates in a different format to the other slabs. Need to account for this

				if slab != 'Cascadia':

					if (lon < 0):
						lon = 360+lon

				ofile.write('%g %g %s\n' %(lon,lat,name))

			ofile.close()

			#Run grdtrack to get the strike of the slab1.0 slab at each of the points

			command = 'gmt grdtrack %s -G%s > %s' %(ofilename,filename,ofilename[:-3]+'strikes.dat')

			os.system(command)

			ifile = open(ofilename[:-3]+'strikes.dat','r')

			lines = ifile.readlines()

			#Distace in degrees to determine the cross section over 
			tdist = 30


			for line in lines:
				vals = line.split()
				midlon = float(vals[0])
				midlat = float(vals[1])
				name = vals[2].strip()
				strike = float(vals[3])
				az1 = strike-90

				lon1 = midlon + (tdist/2)*np.sin(az1*(np.pi/180))
				lat1 = midlat + (tdist/2)*np.cos(az1*(np.pi/180))

				lon0 = midlon - (tdist/2)*np.sin(az1*(np.pi/180))
				lat0 = midlat - (tdist/2)*np.cos(az1*(np.pi/180))

				#Sometimes the point falls outside the slab1.0 bounds. If this is the case, we write it to a file so we can investigate later

				if not mt.isnan(az1):

					SyracuseFile.write('%s %g %g %g %g %g %g\n' %(name,lon0,lat0,lon1,lat1,strike,az1))

					print '--------------------------\n'
					print 'Slice: %g %g %g %s' %(midlon,midlat,az1,name)
					print '--------------------------\n'

					tomotools.Becker_slice(midlon,midlat,lon1,lat1,az1,model='all',dist=tdist,name=name)

				else:

					NaNsfile.write('%g %g %g %s\n' %(midlon,midlat,az1,name))

					SyracuseFile.write('%s %g %g %g %g %g %g\n' %(name,lon0,lat0,lon1,lat1,strike,az1))


	NaNsfile.close()
	SyracuseFile.close()



if __name__ == '__main__':

	cwd = os.getcwd()

	os.chdir('slab1.0_data')

	slabbounds = read_slabdata()

	generate_grdtrackpoints(slabbounds)

	#Clean up
	os.system('rm tmp*')

	os.chdir(cwd)


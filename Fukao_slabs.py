#!/usr/bin/env python

from __future__ import print_function
import numpy as np
from slab import Slab
from catalog import Catalog

def read_slabdata():
	'''
	Creates slab objects for the slabs identified by Fukao & Obayashi (2013)
	'''

	slabs = {}

	#This dictionary contains information about the boundaries of the slabs mentioned in the 
	#Fukao and Obayashi paper

	#The coordinates are entered as follows:
	#[minlon,minlat,maxlon,maxlat]

	slabs['Honshu'] = {'Bounds' : [109,35,159,45]}
	slabs['Bonin'] = {'Bounds' : [109,23,159,35]}
	slabs['Mariana'] = {'Bounds' : [120,8.3,152,25]}
	slabs['Java'] = {'Bounds' : [89,-10,129,19]}
	slabs['Phillippine'] = {'Bounds' : [114,5,128,18]}
	slabs['Tonga'] = {'Bounds' : [157,-26,190,-10]}
	slabs['Kermadec'] = {'Bounds' : [157,-37,190,-26]}
	slabs['Peruvian'] = {'Bounds' : [270,-15,307,6.7]}
	slabs['Chilean'] = {'Bounds' : [270,-45,307,-15]}
	slabs['Central_American'] = {'Bounds' : [249,4.5,282,33]}

	#For these two regions, we need an oblique box. The coodinates are entered as follows
	#[upper left,lower left,lower right,upper right]

	slabs['Kurile_N'] = {'Bounds' : [134.5,65.5,130.5,61.5,158.5,43,162.5,47]}
	slabs['Kurile_S'] = {'Bounds' : [127.5,58.5,123.5,54.5,151.5,36,155.5,40]}

	FukaoSlabs = Catalog()

	for slabname in slabs:

		newslab = Slab(slabname)

		#add directory to the slabs
		newslab.add_Fukao_slab_details(slabs[slabname])
		FukaoSlabs.append(newslab)

	return FukaoSlabs

if __name__ == '__main__':

	'''Testing - plot the Fukao slab boxes'''

	FukaoSlabs = read_slabdata()
	print(FukaoSlabs)

	slist = FukaoSlabs.slabs

	#for slab in slist:
	#	slab.map_Fukao_slab_box()

	FukaoSlabs.InteractiveMap(plotslabs=True)






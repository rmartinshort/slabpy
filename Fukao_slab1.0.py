#!/usr/bin/env python

#RMS July 2016 

#Script to read both the Fukao slab data and slab1.0, and plot on the interactive map

import Fukao_slabs as read_fukao
import read_slab1 as read_slab1

def main():

	Slabs1, slabbounds = read_slab1.read_slabdata() #Generate a slab satalog containing the slab1.0 data

	appendedslabs = read_fukao.read_slabdata(slabcat=Slabs1) #Append Fukaoslabs to the existing catalog\

	appendedslabs.InteractiveMap(plotslabs=True)


if __name__ == '__main__':

	main()
	
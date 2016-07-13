#!/usr/bin/env python 

from obspy.core.util.geodetics.base import gps2DistAzimuth
import os

def read_faults(input_file):
    '''
    Read in the faults file to plot the faults - used for adding tectonic plate boundaries
    '''
    count = 0
    faults = {}
    lats = []
    lons = []
    for line in open(input_file):
        # Check for empty/commented lines
        if not line or line.startswith('>'):
            # If 1st one: new block
            if count == 0:
                pass
            count += 1
        # Non empty line: add line in current(last) block
            faults[count] = [lons, lats]
            lats = []
            lons = []
            
        else:
            try:
                lon = float(line.split()[1])
                lat = float(line.split()[0])

                if lon < 0:
                    lon = 360+lon
                    
                lats.append(lat)
                lons.append(lon)
            except:
                print 'error in reading coordinate: skipped'            
            
    return faults

def Ritsema_180_sections(midlon,midlat,az1):
    '''
    Direct link to Ritsema c-shell/fortran codes for plotting a 180 degree slice through several tomography models. Makes a pdf figure.
    '''

    cwd = os.getcwd()
    cmddir = 'Ritsema_extraction_tools'
    os.chdir(cmddir)
    os.system('./mkcr_all %g %g %g' %(midlat,midlon,az1))
    os.chdir(cwd)


def coords_for_profile(lon0,lat0,lon1,lat1):
    '''
    Determine the minpoint coordinates and the two azimuths for input into program to determine the tomography slice
    '''

    #determine midpoint:

    # if lat1*lat2 < 0:
    #     midlat = (max([lat1,lat2])+min(lat1,lat2))/2
    # else:
    #     midlat = min([lat1,lat2]) + (max([lat2,lat1])-min([lat1,lat2]))/2

    # if lon2*lon1 < 0:
    #     midlon = (max([lon1,lon2])+min(lon1,lon2))/2
    # else:
    #     midlon = min([lon1,lon2]) + (max([lon1,lon2])-min(lon1,lon2))/2

    # print 'midpoint:'
    # print midlat,midlon

    #Find the azimuth from the mipoint to the map marker

    d1, az1, bz1 = gps2DistAzimuth(lat0,lon0,lat1,lon1)

    return az1





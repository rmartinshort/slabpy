#!/usr/bin/env python 

from obspy.core.util.geodetics.base import gps2DistAzimuth
import Tomo_slice_manipulation_tools as slicetools
import matplotlib.pyplot as plt
import os
import glob
import numpy as np

global becker_data_directory
becker_data_directory = '/Users/rmartinshort/Documents/Workshops/CIDER_2016/Slab_group2/datasets/Profile_plotting'



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

    wdir = os.getcwd()
    print wdir
    cmddir = 'Ritsema_extraction_tools'
    os.chdir(cmddir)
    os.system('./mkcr_all %g %g %g' %(midlat,midlon,az1))

    grdfiles = glob.glob('output_slices/uit*.grd')

    print '----------------------------------------------'
    print 'Stacking the following section files: %s' %grdfiles
    print '----------------------------------------------'

    summedfilename = slicetools.gmt_norm_stack(grdfiles,plot=True)
    print '\n\n\n'
    print summedfilename

    #Produce an image of the traced, stacked slab contours - eventually will return some information about the slab itself from these
    slicetools.tracestackedslab(summedfilename)

    os.chdir(wdir)
    dd =  os.getcwd()
    print dd


def coords_for_profile(lon0,lat0,lon1,lat1):
    '''
    Determine the minpoint coordinates and the two azimuths for input into program to determine the tomography slice
    '''

    d1, az1, bz1 = gps2DistAzimuth(lat0,lon0,lat1,lon1)

    return az1

def Becker_slice(midlon,midlat,lon1,lat1,az1,model='all'):

    '''
    Link to T. Becker's scripts to extract a cartesian .grd file from global tomography datasets. This requires a path 
    to the datesets to be set, and the datasets to be in a particular format. The start and end points define a distance (in degrees)
    which is then used to make the model slice
    '''

    halfdistance = np.sqrt((lat1-midlat)**2 + (lon1-midlon)**2)
    print halfdistance
    print 'LON LAT:'
    print midlon,midlat,lon1,lat1,az1

    #Total length of the slice, in 'degrees'
    tdist = int(halfdistance*2)

    cwd = os.getcwd()

    os.chdir(becker_data_directory)
    print os.getcwd()

    #Make a temp file containing the information to be read into the scripts
    ofile = open('subduction_params.tmp','w')
    ofile.write('#/bin/bash\n')
    ofile.write('slab=${1-1}\n')
    ofile.write('len=%g\n' %tdist)
    ofile.write('case $slab in\n1)\n\nlon=%g;lat=%g;azi=%g;name=USERSLICE\n;;\n\nesac\n' %(midlon,midlat,az1))
    ofile.write('echo $lon $lat $azi $len $name')
    ofile.close()

    os.system('chmod 755 subduction_params.tmp')

    if model == 'all':
        #Make all slices and plot
        models = ['gap_p4','mitp08','llnlg3p','tx2015','semucb-wm1']

        modelcount = 0

        for model in models:
            modeldir = os.getcwd()
            if os.path.exists(model):
                os.chdir(model)

                #run script to extract slice
                os.system('./auto_extract.sh')
                grdfile = glob.glob('tmp.prof*')

                if len(grdfile) > 1:
                    print 'Error! grdfile list should contain one entry only; contains %s' %grdfile

                else:

                    if modelcount == 0:
                        fig, datadic = slicetools.plotfigs(model=model,infile=grdfile[0],i=modelcount)
                        print fig
                    else:
                        fig, datadic = slicetools.plotfigs(model=model,axobj=fig,infile=grdfile[0],i=modelcount,datatracker=datadic)

                    os.system('mv %s %s' %(grdfile[0],grdfile[0][4:]))

                os.chdir(modeldir)
                modelcount += 1

    os.chdir(cwd)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':

    Becker_slice(179.274193548,-21.0483870968,168.024193548,-14.5161290323,299.397714924)














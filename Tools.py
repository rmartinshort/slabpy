#!/usr/bin/env python 


from obspy.fdsn import Client
from obspy.core.util.geodetics.base import gps2DistAzimuth
import Tomo_slice_manipulation_tools as slicetools
import matplotlib.pyplot as plt
import os
import glob
import numpy as np
from mpl_toolkits.basemap import Basemap
import datetime
from obspy import UTCDateTime

global becker_data_directory
becker_data_directory = '/Users/rmartinshort/Documents/Workshops/CIDER_2016/Slab_group2/datasets/Profile_plotting'
global topofile
topofile = '/Users/rmartinshort/Documents/Workshops/CIDER_2016/plumes/Bathy/gebco_08.nc'



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

def get_quakes(startcoors,endcoors,minmag=5.0):
    '''Get earthquakes within a region around the start and end coordinates. These will be plotted on the section, with x distance'''

    client = Client('USGS')

    #Only get intermediate depth quakes, or deeper
    mindepth = 60
    boxHW = 0.5

    quakefile = open('quakedata.dat','w')

    startlon = startcoors[0] - boxHW
    startlat = startcoors[1] + boxHW

    endlon = endcoors[0] + boxHW
    endlat = endcoors[1] - boxHW

    minlon = min(startlon,endlon)
    maxlon = max(startlon,endlon)

    minlat = min(startlat,endlat)
    maxlat = max(startlat,endlat)

    starttime = '1970-01-01'
    endtime = str(datetime.datetime.today()).split(' ')[0]
    
    #print startcoors,endcoors
    
    #print minlon,minlat,maxlon,maxlat

    print '---------------------\nUsing Obspy to get quakes\n---------------------'

    quakecat = client.get_events(starttime=UTCDateTime(starttime), endtime=UTCDateTime(endtime), minlongitude=minlon, maxlongitude=maxlon, minlatitude=minlat, maxlatitude=maxlat, minmagnitude=minmag,mindepth=mindepth)

    for event in quakecat:

        evlon = event.origins[0].longitude
        evlat = event.origins[0].latitude
        evdep = event.origins[0].depth

        quakefile.write('%s %s %s\n' %(evlon,evlat,evdep))

    quakefile.close()

    #Work out the distance from each quake to the profile line, and write to a file 
        #Create the section coordinates. Should make some sort of auto-decision about the spacing
    sectionname = 'tmp_toposection.dat'

    print '---------------------\nMaking section through topography\n---------------------'

    os.system('gmt project -C%g/%g -E%g/%g -G10 -Q > %s' %(minlon,minlat,maxlon,maxlat,sectionname))
    os.system('gmt grdtrack %s -G%s > gridsectiontopo.dat' %(sectionname,topofile))

    #Open the topo file and extract the longest distance. This will be used to scale the quake locations

    infile = open('gridsectiontopo.dat','r')
    lines = infile.readlines()
    infile.close()

    topoX = []
    topoY = []

    for line in lines:
        vals = line.split()
        topoX.append(float(vals[2]))
        topoY.append(float(vals[3]))

    maxdist = topoX[-1]
    topoX = np.array(topoX)
    topoY = np.array(topoY)

    print '---------------------\nGetting quake distances\n---------------------'

    #Make a file containing quakelon, quakelat, dist, and dist along profile
    os.system('gmt mapproject quakedata.dat -Lgridsectiontopo.dat/k > quake_dist.dat')

    #Reorder the columns and do another grdtrack to get distance along the profile
    os.system("awk '{print $5,$6,$1,$2,$3,$4}' quake_dist.dat > quaketmp.dat")
    os.system("rm quake_dist.dat")

    #Now, calculate distance along the profile from the start point
    os.system('gmt mapproject quaketmp.dat -G%g/%g/k > quake_points.dat' %(minlon,minlat))
    os.system('rm quaketmp.dat')

    #Now, open the newly created file and grid section file, and pull the distance data
    infile1 = open('quake_points.dat','r')
    lines1 = infile1.readlines()
    infile1.close()

    Xdistances_quakes = []
    Ydepths_quakes = []

    for line in lines1:
        vals = line.split(' ')
        try:

            evlon = float(vals[0].split('\t')[-1])
            evlat = float(vals[1])
            evdep = float(vals[2])
            evdist = float(vals[3].split('\t')[-2])
            evdistalong = float(vals[3].split('\t')[-1])

            #Only keep if the distance between the event and the profile line is less then 50km
            if evdist <= 50:
                Xdistances_quakes.append(evdistalong)
                Ydepths_quakes.append(-evdep/1000.0)

            #for some reason, some depths don't exist, so use this try; except statement
        except:
            continue

    os.system('rm quake_points.dat')

    return Xdistances_quakes, Ydepths_quakes, maxdist, topoX, topoY

def Becker_slice(midlon,midlat,lon1,lat1,az1,model='all',dist=None,name=None,showplot=False):

    '''
    Link to T. Becker's scripts to extract a cartesian .grd file from global tomography datasets. This requires a path 
    to the datesets to be set, and the datasets to be in a particular format. The start and end points define a distance (in degrees)
    which is then used to make the model slice
    '''

    if not dist:

        halfdistance = np.sqrt((lat1-midlat)**2 + (lon1-midlon)**2)
        print halfdistance
        print 'LON LAT:'
        print midlon,midlat,lon1,lat1,az1

        #Total length of the slice, in 'degrees'
        tdist = int(halfdistance*2)

    else:
        tdist = dist

    print 'Total distance in degrees: %g' %tdist

    lon0 = midlon - (tdist/2)*np.sin(az1*(np.pi/180))
    lat0 = midlat - (tdist/2)*np.cos(az1*(np.pi/180))

    #If there is some failure in the get_quakes routine, we just set the resulting variables 
    #to None type and confinue

    try:
        X,Y,maxdist,topoX,topoY = get_quakes([lon0,lat0],[lon1,lat1])

        #print topoX,topoY

        #Plot a section showing the earthquakes

        quakefig = plt.figure(facecolor='white')
        ax1 = quakefig.add_subplot(211)
        ax1.plot(topoX,topoY,'k')
        ax1.set_xlabel('Distance along profile [km]')
        ax1.set_ylabel('Surface height [m]')
        #ax1.set_xlim([min(X)-10,max(X)+10])
        ax1.set_ylim([min(topoY),max(topoY)])
        ax1.set_title('Topography')

        ax2 = quakefig.add_subplot(212)

        ax2.plot(X,Y,'k.',label='Quake [D>60km M>5.0]')
        ax2.plot([0,maxdist],[0,0],'r--',linewidth=3)
        #ax2.set_xlim([min(X)-10,max(X)+10])
        ax2.set_ylim([-700,50])
        ax2.set_xlabel('Distance along profile [km]')
        ax2.set_ylabel('Depth [km]')
        plt.tight_layout()
        plt.grid()
        plt.legend(loc='best')
        #ax2.set_title('Test earthquake profile plot')
        plt.savefig('Earthquakes_%s.jpg' %name,dpi=200)
        #plt.show()

        #x is a list of distances along profile to quake foci, y is a list of depth to hypocenters, maxdist is the length of the profile
        #according to the GMT commands
        quakeinfo = [X,Y,maxdist]

    except:
        X = None
        Y = None
        maxdist = None
        quakeinfo = None


    cwd = os.getcwd()

    os.chdir(becker_data_directory)

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
                        fig, datadic = slicetools.plotfigs(model=model,infile=grdfile[0],i=modelcount,quakes=quakeinfo)
                        print fig
                    else:
                        fig, datadic = slicetools.plotfigs(model=model,axobj=fig,infile=grdfile[0],i=modelcount,datatracker=datadic,quakes=quakeinfo)

                    os.system('mv %s %s' %(grdfile[0],grdfile[0][4:]))

                os.chdir(modeldir)
                modelcount += 1



        ax = fig.add_subplot(326)

        #if midlon > 180:
        #    midlon = midlon-360

        print name, midlon,midlat
        if midlat == 0:
            midlat = 0.1

        m = Basemap(ax=ax,llcrnrlon=(midlon-10),llcrnrlat=(midlat-10),urcrnrlon=(midlon+10),urcrnrlat=(midlat+11),resolution='l',projection='lcc',lat_0=midlat,lon_0=midlon)
        #m.fillcontinents(color='coral',lake_color='aqua')
        m.drawmeridians(np.arange((midlon-20),(midlon+16),5),labels=[0,0,0,0])
        m.drawparallels(np.arange((midlat-10),(midlat+10),5),labels=[0,0,0,0])
        m.drawmapboundary()

        #try:
        m.etopo(ax=ax)
        #except:


        x,y = m(midlon,midlat)
        m.plot(x,y,'ro',linewidth=2.0)

        try:
            xline,yline = m([midlon,lon1],[midlat,lat1])
            m.plot(xline,yline,'r--',linewidth=1.5)
            xline,yline = m([lon0,midlon],[lat0,midlat])
            m.plot(xline,yline,'r--',linewidth=1.5)
        except:
            pass

        os.chdir(cwd)
        fig.suptitle('%s with azimuth %g' %(name,az1))
        plt.tight_layout()
        plt.savefig('%s_%g.jpg' %(name,az1),dpi=300)

        if showplot == True:
            plt.show()





if __name__ == '__main__':

    Becker_slice(179.274193548,-21.0483870968,168.024193548,-14.5161290323,299.397714924,showplot=True)














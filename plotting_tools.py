#!/usr/bin/env python
#Tools to help with plotting slab data on maps 

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
                if lon < 0:
                	lon = lon + 360
                lat = float(line.split()[0])
                lats.append(lat)
                lons.append(lon)
            except:
                print 'error in reading coordinate: skipped'            
            
    return faults

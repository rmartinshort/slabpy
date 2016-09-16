#!/usr/bin/env python 

#Sanne Cottar July 2016
#plot some slab parameters

from __future__ import print_function
import numpy as np
from slab import Slab
from catalog import Catalog
import matplotlib.pyplot as plt


def read_seracuse():
    file = open('data/syracuseetal_parameters.txt','r')


    slabs=Catalog()
    for line in file.readlines():
        val = line.split()
        name = val[0]
        savei=0
        for i in range(1,4):
            try:
                float(val[i])
                break
            except:
                savei=i
                name = name + ' ' +val[i]

        slab=Slab(name)
        slab.transition_depth = float(val[savei+1])
        slab.transition_T = float(val[savei+2])
        slab.slab_T = float(val[savei+3])
        slab.moho_T = float(val[savei+4])
        slab.max_mantle_T = float(val[savei+5])
        slab.max_mantle_T_depth = float(val[savei+6])
        slab.transition_offset = float(val[savei+7])
        slab.slab_surface_T_30km = float(val[savei+8])
        slab.slab_surface_T_240km = float(val[savei+9])
        slab.min_slab_T_240km = float(val[savei+10])
        slabs.append(slab)
    return slabs


slabs = read_seracuse()
slabs.plot_two_parameters('slab_surface_T_30km','slab_surface_T_240km')
plt.show()
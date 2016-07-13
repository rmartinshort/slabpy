#!/usr/bin/env python 

from __future__ import print_function
import numpy as np

from slab import Slab
import matplotlib.pyplot as plt
from Mapper import DrawableMap

"""
Shamelessly copied implementations from Stream() and Trace() in ObsPy
"""

class Catalog(object):

    """
    Base class for all slabs.

    """

    def __init__(self,slabs=None):
        self.slabs=[]
        if slabs:
            self.slabs.extend([slabs])

    def print_slabs(self):
        '''
        For debugging - just print the slab database
        '''

        for slab in self.slabs:
            slab.print_slab1_details()


    def list(self,name):
        z=[]
        for i in range(len(self.slabs)):
            z.append(getattr(self.slabs[i],name))
        return z


    def append(self, slab):
        """
        append slab
        """
        if isinstance(slab, Slab):
            self.slabs.append(slab)
        else:
            msg = 'Append only supports a single Slab object as an argument.'
            raise TypeError(msg)
        return self

    def plot_two_parameters(self,param1,param2):
        plt.plot(self.list(param1),self.list(param2),'.k' )
        plt.xlabel(param1)
        plt.ylabel(param2)

    def InteractiveMap(self,plotslabs=False):
        '''
        Load an interactive map instance and plot the slab information, if available and desired
        '''

        if plotslabs == True:
            catalogmap = DrawableMap(slabcat=self.slabs)
        catalogmap.showmap()


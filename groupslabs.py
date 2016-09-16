#!/usr/bin/env python 

#RMS July 2016
#Base class for slab catalog object

from __future__ import print_function
import numpy as np


class catalog(object):

    """
    Base class for all slabs.

    """

    def __init__(self,slabs):



def append(self, trace):
    """
        Append a single Trace object to the current Stream object.
        
        :param trace: :class:`~obspy.core.stream.Trace` object.
        
        .. rubric:: Example
        
        >>> from obspy import read, Trace
        >>> st = read()
        >>> tr = Trace()
        >>> tr.stats.station = 'TEST'
        >>> st.append(tr)  # doctest: +ELLIPSIS
        <...Stream object at 0x...>
        >>> print(st)  # doctest: +ELLIPSIS
        4 Trace(s) in Stream:
        BW.RJOB..EHZ | 2009-08-24T00:20:03.000000Z ... | 100.0 Hz, 3000 samples
        BW.RJOB..EHN | 2009-08-24T00:20:03.000000Z ... | 100.0 Hz, 3000 samples
        BW.RJOB..EHE | 2009-08-24T00:20:03.000000Z ... | 100.0 Hz, 3000 samples
        .TEST..      | 1970-01-01T00:00:00.000000Z ... | 1.0 Hz, 0 samples
        """
            if isinstance(trace, Trace):
            self.traces.append(trace)
                else:
                    msg = 'Append only supports a single Trace object as an argument.'
                        raise TypeError(msg)
                            return self


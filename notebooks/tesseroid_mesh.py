"""
Defines a TesseroidRelief object.

A mesh discretizing a relief that ondulates around a reference level.
"""
from __future__ import division
import copy
import numpy as np
from fatiando.mesher import Tesseroid


class TesseroidRelief(object):
    """
    Implements a relief ondulating around a reference level using tesseroids.
    """
    
    def __init__(self, area, shape, relief, reference, props=None):
        assert shape[0]*shape[1] == relief.size
        assert len(area) == 4
        assert area[0] < area[1] and area[2] < area[3]
        self.area = area
        self.shape = shape
        w, e, s, n = area
        nlat, nlon = shape
        self.spacing = ((e - w)/nlon, (n - s)/nlat)
        self._relief = relief
        self.reference = reference
        self.set_top_bottom()
        if props is None:
            self.props = {}
        else:
            self.props = props
        self._i = 0
        
    def set_top_bottom(self):
        self._top = self.relief.copy()
        self._bottom = self.reference*np.ones_like(self.relief)
        isbelow = self._top <= self.reference
        self._top[isbelow] = self.reference
        self._bottom[isbelow] = self.relief[isbelow]
        
    @property
    def top(self):
        return self._top
    
    @property
    def bottom(self):
        return self._bottom
    
    @property
    def relief(self):
        return self._relief
    
    @relief.setter
    def relief(self, z):
        assert z.size == self.size
        self._relief = z
        self.set_top_bottom()
        
    @property
    def size(self):
        return self.relief.size
    
    def __len__(self):
        return self.size
    
    def __iter__(self):
        self._i = 0
        return self
    
    def next(self):
        if self._i >= self.size:
            raise StopIteration
        cell = self.__getitem__(self._i)
        self._i += 1
        return cell
    
    def __getitem__(self, index):
        nlat, nlon = self.shape
        w, e, s, n = self.area
        dlon, dlat = self.spacing
        j = index//nlon
        i = index - j*nlon
        cw = w + j*dlon
        ce = cw + dlon
        cs = s + i*dlat
        cn = cs + dlat
        top = self.top[index]
        bottom = self.bottom[index]
        props = {}
        for p in self.props:
            props[p] = self.props[p][index]
        cell = Tesseroid(cw, ce, cs, cn, top, bottom, props)
        return cell
    
    def copy(self, deep=False):
        if deep:
            other = copy.deepcopy(self)
        else:
            other = copy.copy(self)
        return other
"""
Defines classes and functions to load and manipulate the CRUST1.0 model.

Use *fetch_crust1* to load the model from the zip archive
(you can download it from http://igppweb.ucsd.edu/~gabi/crust1.html ).
The function will return a *Crust1* object.
"""
from __future__ import division
import tarfile
import hashlib
import numpy as np
from fatiando.mesher import Tesseroid
    

class Crust1(object):
    """
    The CRUST1.0 model.
    
    Allows you to get the layers, slice and get properties from the model.
    """
    
    layers = """
    water 
    ice 
    upper_sediments
    middle_sediments
    lower_sediments
    upper_crust
    middle_crust
    lower_crust
    mantle
    """.split()
    
    # The SHA256 hash of the model file downloaded from
    # http://igppweb.ucsd.edu/~gabi/crust1.html
    # Used to check if the file used is the correct one
    # or if it's been corrupted.
    sha256 = '0b41b46fc3e1a76debbbcb66ab407febaeee4dc3033db7a0a24d2bb7c7adfe3e'
    
    def __init__(self, lats, lons, topo, vp, vs, density):
        self.lons, self.lats = lons, lats
        self.topo = topo
        self.vp = vp
        self.vs = vs
        self.density = density
        self.lon, self.lat = np.meshgrid(lons, lats)
        self._make_layers()
        
    def _make_layers(self):
        for i in range(len(self.layers) - 1):
            layer = _Layer(self.lat, self.lon, self.topo[i], self.topo[i + 1],
                           vp=self.vp[i], vs=self.vs[i], density=self.density[i])
            setattr(self, self.layers[i], layer)
        layer = _Layer(self.lat, self.lon, self.topo[-1], None,
                       vp=self.vp[-1], vs=self.vs[-1], density=self.density[-1])
        setattr(self, self.layers[-1], layer)
        
    @property
    def sediment_thickness(self):
        """
        A 2D array with the total sediment thickness.
        """
        thick = (self.upper_sediments.thickness 
                 + self.middle_sediments.thickness 
                 + self.lower_sediments.thickness)
        return thick
    
    @property
    def crustal_thickness(self):
        """
        A 2D array with the total crustal thickness.
        """
        thick = (self.sediment_thickness
                 + self.upper_crust.thickness 
                 + self.middle_crust.thickness 
                 + self.lower_crust.thickness)
        return thick
    
    @property
    def moho_depth(self):
        """
        A 2D array with the Moho depth.
        """
        return -self.lower_crust.bottom
    
    def cut(self, area):
        """
        Extract a subset of the model contained in the given area.
        
        *area* should be (w, e, s, n) in degrees.
        """
        s, n, w, e = area
        imin = np.searchsorted(self.lats, s)
        imax = np.searchsorted(self.lats, n)
        jmin = np.searchsorted(self.lons, w)
        jmax = np.searchsorted(self.lons, e)
        data = self.__class__(
            lons=self.lons[jmin:jmax],
            lats=self.lats[imin:imax],
            topo=self.topo[:, imin:imax, jmin:jmax],
            vp=self.vp[:, imin:imax, jmin:jmax],
            vs=self.vs[:, imin:imax, jmin:jmax],
            density=self.density[:, imin:imax, jmin:jmax])
        return data

    
class _Layer(object):
    """
    Store a single layer of the model.
    """
    
    def __init__(self, lat, lon, top, bottom=None, **kwargs):
        self.lon = lon
        self.lat = lat
        self.top = top
        self.bottom = bottom
        self.props = kwargs
        for p in self.props:
            setattr(self, p, kwargs[p])
    
    @property
    def thickness(self):
        """
        2D array with the thickness of this layer.
        """
        if self.bottom is None:
            raise ValueError('This layer has no bottom')
        return self.top - self.bottom
    
    @property
    def tesseroids(self):
        """
        Get a tesseroid representation of this layer.
        """
        ds = (self.lon[0, 1] - self.lon[0, 0])/2
        arrays = [self.lon, self.lat, self.top, self.bottom,
                  self.vp, self.vs, self.density]
        args = zip(*[i.ravel() for i in arrays])
        gen = (Tesseroid(lon - ds, lon + ds, lat - ds, lat + ds,
                         top, bottom, dict(vp=vp, vs=vs, density=density))
               for lon, lat, top, bottom, vp, vs, density in args
               if abs(top - bottom) > 10)
        return gen


def fetch_crust1(fname):
    """
    Load the CRUST1.0 model from a file.
    
    You can download the file from http://igppweb.ucsd.edu/~gabi/crust1.html
    
    Parameters:
    
    * fname : string
        The file name (or full path) of the .zip file.
        
    Returns:
    
    * crust1 : Crust1
        The model in a *Crust1* class.
        
    """
    _check_hash_crust1(fname)
    with tarfile.open(fname, 'r:gz') as arc:
        topo = _extract_file(arc, 'bnds')
        density = _extract_file(arc, 'rho')
        vp = _extract_file(arc, 'vp')
        vs = _extract_file(arc, 'vs')
    lons = np.linspace(-179.5, 179.5, 360)
    lats = np.linspace(-89.5, 89.5, 180)
    return Crust1(lats, lons, topo, vp, vs, density)
    

def _check_hash_crust1(fname):
    "Check the hash of the file agains the one recorded in  the class."
    sha = _stream_sha(fname)
    msg = ' '.join([
        'Error reading model from file "{}".'.format(fname),
        'Possibly corrupted file',
        '(invalid SHA256 hash "{}").'.format(sha)
        ])
    assert sha == Crust1.sha256, msg
        

def _stream_sha(fname, chunksize=65536):
    "Calculate the SHA256 hash of a file in chunks"
    hasher = hashlib.sha256()
    with open(fname, 'rb') as f:
        buf = f.read(chunksize)
        while len(buf) > 0:
            hasher.update(buf)
            buf = f.read(chunksize)
    return hasher.hexdigest()

    
def _extract_file(archive, ext):
    "Extract the data from a file in the zip. Returns in numpy array."
    f = archive.extractfile('crust1.{}'.format(ext))
    shape = (9, 180, 360)
    # Convert from km, g/cm^3 and km/s to m, kg/m^3 and m/s
    # The slice is to invert the latitude axis. CRUST uses lat from +89.5 to -89.5
    # (North to South). We use South to North.
    data = 1000*np.loadtxt(f, unpack=True).reshape(shape)[:, ::-1, :]
    return data       
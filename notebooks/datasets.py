# -*- coding: utf-8 -*-
u"""
Defines classes and functions to load and manipulate the datasets we use.

All datasets are provided in the 'data' directory of this repository.

CRUST1.0 model
--------------

* ``fetch_crust1``: loads the model from the zip archive (you can download 
  it from http://igppweb.ucsd.edu/~gabi/crust1.html ). The function will 
  return a *Crust1* object which you can use to interact with the model.

Moho depth
----------

* ``fetch_assumpcao_moho_points``: loads the point data from the
  seismic crustal thickness compilation of Assumpção et al. (2012).
  
ICGEM data files
----------------

* ``load_icgem_gdf``: loads data from an ICGEM (http://icgem.gfz-potsdam.de/ICGEM/)
  ``.gdf`` file.
* ``down_sample``: down-sample the grid data.
  
"""
from __future__ import division
import tarfile
import hashlib
import numpy as np
from fatiando.mesher import Tesseroid
    
    
def down_sample(arrays, shape, every):
    """
    Down-sample data by taking every other data point.
    
    Parameters:
    
    * arrays : list of 1d-arrays
        The data arrays to down-sample
    * shape : tuple = (nlat, nlon)
        The original shape of the data.
    * every : int
        Will take every *every* data. For example, if
        ``every=3``, will take every 3 data points. 
        Equivalent of doing ``[::3, ::3]`` on a numpy
        2d-array.
     
    Returns:
    
    * array1, array2, array3, ..., shape
        The down-sampled arrays and the new shape of the
        data grid.
        
    """
    downsampled = [array.reshape(shape)[::every, ::every]
                   for array in arrays]
    newshape = downsampled[0].shape
    return [v.ravel() for v in downsampled] + [newshape]

    
def load_icgem_gdf(fname, usecols=None):
    """
    Load data from an ICGEM .gdf file.
    
    Returns:
    
    * data : dict
        A dictionary with the data from the file. 
        Reads the column data and other metadata from 
        the file. Column data are numpy arrays.
        
    """
    with open(fname) as f:
        # Read the header and extract metadata
        metadata = []
        shape = [None, None]
        size = None
        height = None
        attributes = None
        attr_line = False
        area = [None]*4
        for line in f:
            if line.strip()[:11] == 'end_of_head':
                break
            metadata.append(line)
            if not line.strip():
                attr_line = True
                continue
            if not attr_line:
                parts = line.strip().split()
                if parts[0] == 'height_over_ell':
                    height = float(parts[1])
                elif parts[0] == 'latitude_parallels':
                    shape[0] = int(parts[1])
                elif parts[0] == 'longitude_parallels':
                    shape[1] = int(parts[1])
                elif parts[0] == 'number_of_gridpoints':
                    size = int(parts[1])
                elif parts[0] == 'latlimit_south':
                    area[0] = float(parts[1])
                elif parts[0] == 'latlimit_north':
                    area[1] = float(parts[1])
                elif parts[0] == 'longlimit_west':
                    area[2] = float(parts[1])
                elif parts[0] == 'longlimit_east':
                    area[3] = float(parts[1])
            else:
                attributes = line.strip().split()
                attr_line = False
        # Read the numerical values
        rawdata = np.loadtxt(f, usecols=usecols, ndmin=2, unpack=True)
    # Sanity checks
    assert all(n is not None for n in shape), "Couldn't read shape of grid."
    assert size is not None, "Couldn't read size of grid."
    shape = tuple(shape)
    assert shape[0]*shape[1] == size, \
        "Grid shape '{}' and size '{}' mismatch.".format(shape, size)
    assert attributes is not None, "Couldn't read column names."
    if usecols is not None:
        attributes = [attributes[i] for i in usecols]
    assert len(attributes) == rawdata.shape[0], \
        "Number of attributes ({}) and data columns ({}) mismatch".format(
            len(attributes), rawdata.shape[0])
    assert all(i is not None for i in area), "Couldn't read the grid area."
    # Return the data in a dictionary with the attribute names
    # that we got from the file.
    data = dict(shape=shape, area=area, metadata=''.join(metadata))
    for attr, value in zip(attributes, rawdata):
        # Need to invert the data matrices in latitude "[::-1]"
        # because the ICGEM grid gets varies latitude from N to S
        # and the TesseroidRelief expects the opposite.
        data[attr] = value.reshape(shape)[::-1].ravel()
    if (height is not None) and ('height' not in attributes):
        data['height'] = height*np.ones(size)
    if 'latitude' in attributes and 'longitude' in attributes:
        lat, lon = data['latitude'], data['longitude']
        area = (lat.min(), lat.max(), lon.min(), lon.max())
        assert np.allclose(area, data['area']), \
            "Grid area read ({}) and calculated from attributes ({}) mismatch.".format(
                data['area'], area)
    return data        
    
    
# The SHA256 has of the Assumpção et al dataset archive file
# Used to verify that the given file is not corrupted.
ASSUMPCAO_HASH = '80bf7f3be4b4cc8899d403b95ee8d1cc874c7eb70d8e83151b2cbfa353c00179'


def fetch_assumpcao_moho_points(fname, todepth=True):
    u"""
    Extract the point seismic data of Assumpção et al. (2012) from the 
    tar.gz archive.
    
    Parameters:
    
    * fname : string
        The name (or full path) of the tar.gz archive with the data.
    * todepth : True or False
        If True, will convert the crustal thickness data to Moho depth
        (in meters)
        
    Returns:
    
    * lat, lon, height, data : 1d-arrays
        The latitude, longitude, and altitude coordinates of each data 
        point, the corresponding crustal thickness or Moho depth (in meters),
        and the corresponding uncertainty.
        
    """
    _check_hash(fname, ASSUMPCAO_HASH)
    with tarfile.open(fname, 'r:gz') as archive:
        # Need to get the data separate because there is a break in the format
        # when the Marcelo and Andres files are concatenated into the SAm file
        f1 = archive.extractfile('Moho_Map_SAm2013_data/Compilation_SAm.06NOV2012.IXYEHU.dat')
        data1 = np.genfromtxt(f1.readlines()[:754], usecols=[1, 2, 3, 4, 5])
        f2 = archive.extractfile('Moho_Map_SAm2013_data/Compilation_Andres.XYEHU.dat')
        data2 = np.genfromtxt(f2.readlines(), usecols=[0, 1, 2, 3, 4])
        data = np.vstack([data1, data2]).T
    lon, lat, height, crustal_thick, uncert = data
    # Convert from km to meters
    crustal_thick *= 1000
    uncert *= 1000
    if todepth:
        crustal_thick -= height
    return lat, lon, height, crustal_thick, uncert
    

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
    _check_hash(fname, Crust1.sha256)
    with tarfile.open(fname, 'r:gz') as arc:
        topo = _extract_file_crust1(arc, 'bnds')
        density = _extract_file_crust1(arc, 'rho')
        vp = _extract_file_crust1(arc, 'vp')
        vs = _extract_file_crust1(arc, 'vs')
    lons = np.linspace(-180, 180, 360, endpoint=False)
    lats = np.linspace(-90, 90, 180, endpoint=False)
    model = Crust1(lats, lons, topo, vp, vs, density)
    assert model.shape == (180, 360), \
        "Model shape mismatch: {}".format(model.shape)
    return model
    

def _check_hash(fname, true_hash):
    "Check the hash of the file agains the one recorded in  the class."
    sha = _stream_sha(fname)
    msg = '\n'.join([
        'Error reading model from file "{}". Possibly corrupted file.'.format(fname),
        '  - Calculated SHA256 hash: {}'.format(sha),
        '  - Known (recorded) SHA256 hash: {}'.format(true_hash)])
    assert sha == true_hash, msg
        

def _stream_sha(fname, chunksize=65536):
    "Calculate the SHA256 hash of a file in chunks"
    hasher = hashlib.sha256()
    with open(fname, 'rb') as f:
        buf = f.read(chunksize)
        while len(buf) > 0:
            hasher.update(buf)
            buf = f.read(chunksize)
    return hasher.hexdigest()

    
def _extract_file_crust1(archive, ext):
    "Extract the data from the CRUST1.0 file in the zip. Returns in numpy array."
    f = archive.extractfile('crust1.{}'.format(ext))
    shape = (9, 180, 360)
    # Convert from km, g/cm^3 and km/s to m, kg/m^3 and m/s
    # The slice is to invert the latitude axis. CRUST uses lat from +89.5 to -89.5
    # (North to South). We use South to North.
    data = 1000*np.loadtxt(f, unpack=True).reshape(shape)[:, ::-1, :]
    return data       
    
    
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
        assert np.allclose(lons[1:] - lons[0:-1], 1), "Spacing in lon is not 1"
        assert np.allclose(lats[1:] - lats[0:-1], 1), "Spacing in lat is not 1"
        # Longitude and latitude of the SW corner of the cells
        self.lons, self.lats = lons, lats
        self.lon, self.lat = np.meshgrid(lons, lats)
        # Longitude and latitude of the center of each cell
        self.clons = lons + 0.5 
        self.clats = lats + 0.5 
        self.clon = self.lon + 0.5 
        self.clat = self.lat + 0.5 
        self.topo = topo
        self.vp = vp
        self.vs = vs
        self.density = density
        self._make_layers()
        
    @property
    def shape(self):
        return (len(self.lats), len(self.lons))
        
    @property
    def area(self):
        s = self.lats[0]
        n = self.lats[-1] + 1
        w = self.lons[0]
        e = self.lons[-1] + 1
        return [s, n, w, e]
        
    def _make_layers(self):
        for i in range(len(self.layers) - 1):
            layer = _Layer(self.clat, self.clon, self.topo[i], self.topo[i + 1],
                           vp=self.vp[i], vs=self.vs[i], density=self.density[i])
            setattr(self, self.layers[i], layer)
        layer = _Layer(self.clat, self.clon, self.topo[-1], None,
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
        
        *area* should be (s, n, w, e) in degrees.
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
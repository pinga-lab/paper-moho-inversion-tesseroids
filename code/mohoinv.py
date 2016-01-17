# -*- coding: utf-8 -*-
"""
Implements the proposed Moho inversion method.

See the Jupyter notebooks 'moho-inversion-example.ipynb' and
'tesseroid-relief-example.ipynb' for instructions and example usage of these
classes.

Defines the classes:

* `MohoGravityInvSpherical`: Performs the inversion itself. Uses the classes
  from *fatiando.inversion* as a base and *fatiando.gravmag.tesseroid* for
  forward modeling.

* `TesseroidRelief`: Class to represent a relief undulating around a reference
  level and discretized into tesseroids.

"""
from __future__ import division
from fatiando.inversion import Misfit, CachedMethodPermanent, CachedMethod
from fatiando.gravmag import tesseroid
from fatiando import constants, utils, gridder
from fatiando.mesher import Tesseroid
import scipy.sparse
import numpy as np
import copy
import multiprocessing
import warnings


def split_data(data, shape, every_other):
    """
    Split the data into inversion and test datasets
    
    Parameters:
    
    * data : list of 1d-arrays
        List with the data arrays that will be divided (like lon, lat, 
        gravity, etc)
    * shape : tuple = (nlat, nlon)
        The original shape of the data grid
    * every_other : int
        Will divide by taking every other "every_other" grid point
        for the inversion dataset. The remaining points will be the 
        test data.
        
    Returns:
    
    * inversion_set, test_set, inversion_shape:
        The inversion_set and test_set are lists in the same order as "data".
        inversion_shape is the shape of the new inversion dataset
        
    """
    data = [i.reshape(shape) for i in data]
    # Take "every_other" point for the inversion
    inversion = [i[::every_other, ::every_other] for i in data]
    # mask marks the grid points I didn't take for inversion
    mask = np.ones(shape, dtype=np.bool)
    mask[::every_other, ::every_other] = False  # These are the ones I took
    test = [i[mask] for i in data]
    shape = inversion[0].shape
    # Sanity checks
    assert all(t.size + i.size == d.size for t, i, d in zip(test, inversion, data)), \
        "Number of points in inversion + test set different from original data."
    assert all(t.size == test[0].size for t in test), "Test set has differet number of points"
    assert all(i.size == inversion[0].size for i in inversion), "Inversion set has differet number of points"
    assert all(i.shape == inversion[0].shape for i in inversion), "Inversion set has differet shape"
    return [i.ravel() for i in inversion], test, shape


def make_mesh(area, shape, relief=None, reference=None):
    """
    Make a TesseroidRelief mesh for use in MohoGravityInvSpherical.

    This functions will create a mesh that has one tesseroid below each
    data point. Data is considered to be on a regular grid.

    Parameters:

    * area : list = (s, n, w, e)
        The dimensions of the data grid in degrees.
    * shape : list = (nlat, nlon)
        The number of data points in latitude and longitude, respectively.
    * relief : 1d-array
        The relief of the interface. If None, will use a constant of 1.
    * reference : float
        The reference level. If None, will use a constant of 0.

    Returns:

    * mesh : TesseroidRelief
        The generated mesh.

    """
    dlat, dlon = gridder.spacing(area, shape)
    s, n, w, e = area
    modelarea = (s - dlat/2, n + dlat/2, w - dlon/2, e + dlon/2)
    if relief is None:
        relief = np.ones(shape).ravel()
    if reference is None:
        reference = 0
    mesh = TesseroidRelief(modelarea, shape, relief, reference)
    return mesh


class MohoGravityInvSpherical(Misfit):
    """
    Gravity inversion for the relief of the Moho using tesseroids.
    """

    def __init__(self, lat, lon, height, data, mesh, njobs=1, field='gz'):
        super(MohoGravityInvSpherical, self).__init__(data=data, nparams=mesh.size,
                                                      islinear=False, cache=False)
        # Not caching the Jacobian and Hessian because
        # the Jacobian (and thus the Hessian) change depending
        # on the optimization method used. This will not be
        # much of a problem because they take little time to build.
        self.predicted = CachedMethod(self, 'predicted')
        if mesh.size != data.size:
            msg = ("The mesh size ({}) is different".format(mesh.size)
                   + " from the data size ({}).".format(data.size)
                   + " The mesh elements should be below each data point."
                   + " Make sure you know what you're doing.")
            warnings.warn(msg, RuntimeWarning)
        assert field in ['gz']
        self.lon = lon
        self.lat = lat
        self.height = height
        self.mesh = mesh.copy(deep=True)
        self.field = field
        self._kernel = getattr(tesseroid, field)
        self.njobs = njobs
        self.pool = None
        # The value that will go in the Jacobian
        self.kernel_kwargs = dict()

    def predicted(self, p):
        """
        Calculate the predicted data for a given parameter vector.
        """
        mesh = self.mesh
        mesh.relief = p
        # Need to fix the density in case the relief crossed
        # the reference level at some point (density would
        # change sign if that happens).
        self.fix_density(mesh)
        # Check if this should run in parallel
        njobs = self.njobs
        if self.pool is None:
            njobs = 1
        data = self._kernel(self.lon, self.lat, self.height, mesh,
                            njobs=njobs, pool=self.pool, **self.kernel_kwargs)
        return data

    def jacobian(self, p):
        """
        Calculate the Jacobian matrix.

        Will ignore the value of *p*.

        The Jacobian will have 2*Pi*G*rho in the diagonal
        if using Newton's method or Levemberg-Marquardt
        and 1/(2*2*Pi*G*rho) if using Steepest Descent.
        """
        method = getattr(self, 'fit_method', None)
        step = utils.si2mgal(2*np.pi*constants.G*self.mesh.props['density'])
        if method is not None and method == 'steepest':
            step = 1/(2*step)
        # Use a sparse matrix for the Jacobian in Compressed Sparse Row (CSR)
        # format. See
        # http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html
        jac = scipy.sparse.diags([np.abs(step)], [0]).tocsr()
        return jac

    def fix_density(self, mesh):
        """
        Invert the sign of any cell whose relief crossed the
        reference since the last check.
        """
        dens = mesh.props['density']
        should_b_pos = mesh.relief > mesh.reference
        dens[(should_b_pos) & (dens < 0)] *= -1
        dens[(~should_b_pos) & (dens > 0)] *= -1
        mesh.props['density'] = dens
        return mesh

    def set_reference(self, reference):
        """
        Set the reference level of the mesh.
        """
        self.mesh.reference = reference
        self.fix_density(self.mesh)
        return self

    def set_density(self, density):
        """
        Set the density constrast along the mesh.
        """
        self.mesh.props['density'] = np.ones(self.mesh.size)*density
        self.fix_density(self.mesh)
        return self

    def fit(self):
        """
        Run the inversion.
        """
        # Create the process pool if using more than one job.
        # Need to to this here so that I can reuse the pool object
        # on every call to *predicted*. Otherwise, would have to
        # create a pull everytime and that has a large overhead.
        if self.njobs > 1:
            self.pool = multiprocessing.Pool(self.njobs)
        super(MohoGravityInvSpherical, self).fit()
        if self.pool is not None:
            self.pool.close()
            self.pool = None

    def fmt_estimate(self, p):
        """
        Invert the estimated relief into the mesh and fix the density.

        Returns a TesseroidRelief with the estimated Moho.
        """
        self.mesh.relief = p
        mesh = self.fix_density(self.mesh)
        return mesh

    def config_kernel(self, **kwargs):
        """
        Configuration that will be passed on to the tesseroid
        forward modeling functions.
        """
        self.kernel_kwargs = kwargs
        return self


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
        s, n, w, e = area
        nlat, nlon = shape
        self.lons = np.linspace(w, e, nlon, endpoint=False)
        self.lats = np.linspace(s, n, nlat, endpoint=False)
        self.lon, self.lat = np.meshgrid(self.lons, self.lats)
        self.spacing = self.lats[1] - self.lats[0], self.lons[1] - self.lons[0]
        self._relief = relief
        self._reference = reference*np.ones_like(relief)
        self.set_top_bottom()
        if props is None:
            self.props = {}
        else:
            self.props = props
        self._i = 0

    @property
    def clons(self):
        dlon = self.spacing[1]
        return self.lons + dlon/2

    @property
    def clats(self):
        dlat = self.spacing[0]
        return self.lats + dlat/2

    @property
    def clon(self):
        dlon = self.spacing[1]
        return self.lon + dlon/2

    @property
    def clat(self):
        dlat = self.spacing[0]
        return self.lat + dlat/2

    def addprop(self, prop, values):
        """
        Add physical property values to the mesh.

        Different physical properties of the grid are stored in a dictionary.

        Parameters:

        * prop : str
            Name of the physical property.
        * values :  list or array
            Value of this physical property in each point of the grid

        """
        self.props[prop] = values

    def set_top_bottom(self):
        self._top = self.relief.copy()
        self._bottom = self.reference.copy()
        isbelow = self._top <= self.reference
        self._top[isbelow] = self.reference[isbelow]
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
    def reference(self):
        return self._reference

    @reference.setter
    def reference(self, reference):
        self._reference = np.ones_like(self.relief)*reference
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
        dlat, dlon = self.spacing
        i = index//nlon
        j = index - i*nlon
        w = self.lons[j]
        e = w + dlon
        s = self.lats[i]
        n = s + dlat
        top = self.top[index]
        bottom = self.bottom[index]
        props = {}
        for p in self.props:
            props[p] = self.props[p][index]
        cell = Tesseroid(w, e, s, n, top, bottom, props)
        return cell

    def copy(self, deep=False):
        if deep:
            other = copy.deepcopy(self)
        else:
            other = copy.copy(self)
        return other

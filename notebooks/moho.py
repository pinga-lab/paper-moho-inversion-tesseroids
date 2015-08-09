"""
The gravity inversion class for the relief of the Moho.

Uses the base classes from *fatiando.inversion*.
"""
from __future__ import division
from fatiando.inversion import Misfit, CachedMethodPermanent, CachedMethod
from fatiando.gravmag import tesseroid
from fatiando import constants, utils
import scipy.sparse
import numpy as np
import copy
import multiprocessing


class MohoGravityInvSpherical(Misfit):
    """
    Gravity inversion for the relief of the Moho using tesseroids.
    """
    
    def __init__(self, lon, lat, height, data, mesh, njobs=1, field='gz'):
        super(MohoGravityInvSpherical, self).__init__(data=data, nparams=mesh.size, 
                                                      islinear=False, cache=False)
        # Not caching the Jacobian and Hessian because
        # the Jacobian (and thus the Hessian) change depending
        # on the optimization method used. This will not be 
        # much of a problem because they take little time to build.
        self.predicted = CachedMethod(self, 'predicted')
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
        self.step = utils.si2mgal(2*np.pi*constants.G*self.mesh.props['density'])
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
        step = self.step
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
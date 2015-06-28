from __future__ import division
from fatiando.inversion import Misfit, CachedMethodPermanent, CachedMethod
from fatiando.gravmag import tesseroid
from fatiando import constants, utils
import scipy.sparse
import numpy as np
import copy
import multiprocessing


class MohoGravityInversion(Misfit):
    def __init__(self, lon, lat, height, data, mesh, njobs=1, field='gz'):
        super(MohoGravityInversion, self).__init__(data=data, nparams=mesh.size, islinear=False, cache=False)
        assert mesh.size == data.size
        assert field in ['gz']
        self.lon = lon
        self.lat = lat
        self.height = height
        self.mesh = mesh.copy(deep=True)
        # Not caching the Jacobian and Hessian because
        # The Jacobian (and thus the Hessian) change depending
        # on the optimization method used. This will not be 
        # much of a problem because they take little time to build.
        self.predicted = CachedMethod(self, 'predicted')
        self.field = field
        self._kernel = getattr(tesseroid, field)
        self.njobs = njobs
        self.pool = None
        # The value that will go in the Jacobian
        self.set_step()
        self.kernel_kwargs = dict()
        
    def predicted(self, p):
        mesh = self.mesh
        mesh.relief = p
        self.fix_density(mesh)
        njobs = self.njobs
        if self.pool is None:
            njobs = 1
        data = self._kernel(self.lon, self.lat, self.height, mesh, 
                            njobs=njobs, pool=self.pool, **self.kernel_kwargs)
        return data
    
    def jacobian(self, p):
        # Multiply by ones in case the step is a scalar
        step = np.ones(self.nparams)*self.step
        method = getattr(self, 'fit_method', None)
        if method is not None and method == 'steepest':
            step = 1/(2*step)
        jac = scipy.sparse.diags([np.abs(step)], [0]).tocsr()
        return jac
    
    def fix_density(self, mesh):
        dens = mesh.props['density']
        should_b_pos = mesh.relief > mesh.reference
        dens[(should_b_pos) & (dens < 0)] *= -1
        dens[(~should_b_pos) & (dens > 0)] *= -1
        mesh.props['density'] = dens
        return mesh
        
    def fit(self):        
        if self.njobs > 1:
            self.pool = multiprocessing.Pool(self.njobs)
        super(MohoGravityInversion, self).fit()
        if self.pool is not None:
            self.pool.close()
            self.pool = None
    
    def fmt_estimate(self, p):
        self.mesh.relief = p
        mesh = self.fix_density(self.mesh)
        return mesh
        
    def config_kernel(self, **kwargs):
        self.kernel_kwargs = kwargs
        return self
    
    def set_step(self, step=None):
        if step is None:
            self.step = utils.si2mgal(2*np.pi*constants.G*self.mesh.props['density'])
        else:
            self.step = step
        return self
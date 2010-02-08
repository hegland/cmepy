"""
Implementation of Solver-style interface for matrix-exponential CME solver.
"""

import numpy
import cmepy.core
import scipy.sparse
from cmepy.solver import Solver

import cmepy.core.matrix_cme



class CmeSolverMatrix(Solver):
    """
    A cmepy.solver.Solver style interface for a cmepy.core.matrix_cme based
    solver. This solver uses SciPy's VODE solver under the hood, but the
    differential equations are computed using a sparse matrix. This sparse
    matrix may be explicitly passed when initialising the CmeSolverMatrix.
    """
    
    def __init__(self, model, sparse_matrix = None):
        """
        Initialise CmeSolverMatrix using the provided model, and optional
        sparse matrix.
        """
        Solver.__init__(self, model)
        self._args = None
        self.sparse_matrix = sparse_matrix
        self.t = None
        self.p = None
        self.cme = None
        
        # load np from model, if it is present
        if ('np' in model):
            self.set_solver_params(np = model['np'])
        
    def set_solver_params(self, **args):
        if self._args is None:
            self._args = {}
        self._args.update(args)
    
    def _begin(self):
        nl = numpy.product(self._args['np'])
        if self.p0 is None:
            p0 = numpy.zeros(nl)
            p0[0] = 1.0
            self.p0 = p0
        else:
            self.p0 = numpy.reshape(self.p0, (nl, ))

        cme_args = {}
        cme_args.update(self.model)
        cme_args.update(self._args)
        
        if self.sparse_matrix is None:
            self.sparse_matrix = cmepy.core.sparse_matrix_factory(**cme_args)
        self.sparse_matrix = self.sparse_matrix.tocsr()
        self.sparse_matrix.sum_duplicates()
        
        cme_args['p0'] = self.p0
        cme_args['t0'] = self.t0
        
        # construct ode solver instance based on sparse matrix action
        self.cme = cmepy.core.matrix_cme.cme(self.sparse_matrix,
                                             **cme_args)
        
        # set initial solver state
        self.p = numpy.array(self.p0)
        self.t = self.t0
        
    def step(self, t):
        """
        Advances the current solution p(t) to the time t.
        """
        
        if self.cme is None:
            self._begin()
        assert(t >= self.cme.t)
        
        if t == self.cme.t:
            return
        
        self.cme.integrate(t)
        if not self.cme.successful():
            complaint = ('CME integration failure '+
                         '(look for messages from DVODE / vode)')
            raise RuntimeError, complaint
        self.t = t
    
    def get_p(self):
        """
        Obtain the current solution p(t) of the probability distribution.
        """
        if self.cme is None:
            self._begin()
        p_current = numpy.array(self.cme.y)
        return numpy.reshape(p_current, self._args['np'])

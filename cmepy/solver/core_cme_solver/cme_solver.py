"""
Implementation of Solver-style interface for cmepy.core cme solver.
"""

import numpy
import cmepy.core
from cmepy.solver import Solver



class CmeSolver(Solver):
    """
    A cmepy.solver.Solver style interface for the cmepy.core.cme ode-based
    solver.
    """
    
    def __init__(self, model):
        """
        Initialise CmeSolver using the provided model.
        """
        Solver.__init__(self, model)
        self._args = None
        self.cme = None
        
        # load np from model, if it is present
        if ('np' in model):
            self.set_solver_params(np = model['np'])
        
    def set_solver_params(self, **args):
        assert('np' in args), 'missing argument \'np\' : shape of state space'
        self._args = args
    
    def _begin(self):
        if self.p0 is None:
            nl = numpy.product(self._args['np'])
            p0 = numpy.zeros(nl)
            p0[0] = 1.0
            p0 = numpy.reshape(p0, self._args['np'])
            self.p0 = p0
        
        cme_args = {}
        cme_args.update(self._args)
        cme_args['propensities'] = self.model['propensities']
        cme_args['p0'] = self.p0
        cme_args['t0'] = self.t0
        self.cme = cmepy.core.cme(**cme_args)
    
    def step(self, t):
        """
        Advances the current solution p(t) to the time t.
        """
        
        if not self.cme:
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
        if not self.cme:
            self._begin()
        p_current = numpy.array(self.cme.y)
        return numpy.reshape(p_current, self._args['np'])
            
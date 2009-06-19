import unittest
from test import test_support
from cmepy.solver.core_cme_solver import cme_solver


class CmeSolverTest(unittest.TestCase):
    """
    A test class for the cme_solver module
    
    This test class is adapted from the cme_tests file from cmepy/core .
    """

    def test1DPoisson(self):
        """
        testing 1D Poisson process
        """
        from scipy.stats import poisson
        from numpy import arange, zeros
        from numpy.testing import assert_array_almost_equal
        np = (1024,)
        lam = 30.0
        pexact = poisson(lam).pmf(arange(np[0]))
        p0 = zeros(np)
        p0[0] = 1
        
        model = {'propensities' : (lambda x : lam, )}
        solver = cme_solver.CmeSolver(model)
        solver.set_solver_params(np=np)
        solver.set_initial_values(t0=0.0, p0=p0)
        solver.step(1.0)
        p = solver.get_p()
        
        assert_array_almost_equal(p[:-1], pexact[:-1])

    def test2DPoisson(self):
        """
        testing 2D Poisson process
        """
        from scipy.stats import poisson
        from numpy import arange, outer
        from numpy.testing import assert_array_almost_equal
        np = (16, 32)
        lam = (5.0, 2.3)
        pex1 = poisson(lam[0]).pmf(arange(np[0]))
        pex2 = poisson(lam[1]).pmf(arange(np[1]))
        pexact = outer(pex1, pex2)
        
        model = {'propensities' : (lambda x1, x2 : lam[0],
                                   lambda x1, x2 : lam[1]) }
        
        solver = cme_solver.CmeSolver(model)
        solver.set_solver_params(np=np)
        # use default behaviour of cmepy.core.cme when handling
        # a missing p0 value ... 
        solver.set_initial_values(t0=0.0, p0=None)
        solver.step(1.0)
        p = solver.get_p()
        
        # note that the last rows/columns are not the same as propensity is set
        # to zero there
        assert_array_almost_equal(p[:-1, :-1], pexact[:-1, :-1])

    def test1DBinomial(self):
        """
        testing 1D Binomial process
        """
        from scipy.stats import binom
        from numpy import arange, zeros
        from numpy.testing import assert_array_almost_equal
        from math import exp
        np = (32,)
        lam = 1.0
        pr = 1.0 - exp(-lam)
        pexact = binom(np[0], pr).pmf(arange(np[0]))
        p0 = zeros(np)
        p0[0] = 1
        
        model = {'propensities' : (lambda x : lam*(np[0]-x), ) }
        solver = cme_solver.CmeSolver(model)
        solver.set_solver_params(np=np)
        solver.set_initial_values(t0=0.0, p0=p0)
        solver.step(1.0)
        p = solver.get_p()
        
        assert_array_almost_equal(p[:-1], pexact[:-1])

def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(CmeSolverTest)
    return suite

def main():
    test_support.run_unittest(CmeSolverTest)

if __name__ == '__main__':
    main()

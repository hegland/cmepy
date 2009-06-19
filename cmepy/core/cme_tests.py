"""
Unit tests for core_cme.

These compare numerical solutions of certain models
against their exact solutions ...
"""

import unittest
from test import test_support
from cmepy.core import cme

def get_exact_poisson_p(lambdas, p_shape):
    """
    returns exact Poisson distribution
    """
    import scipy.stats
    import itertools
    import numpy
    poissons = []
    for (lam, dim) in itertools.izip(lambdas, p_shape):
        poissons += [scipy.stats.poisson(lam).pmf(numpy.arange(dim))]
    p_exact = reduce(numpy.multiply.outer, poissons)
    return p_exact


class CMETest(unittest.TestCase):
    """
    A test class for the CME module.
    """

    def test_1d_poisson(self):
        """
        testing 1D Poisson process
        """
        
        import numpy
        
        p_shape = (1024,)
        lambdas = (30.0,)

        p_exact = get_exact_poisson_p(lambdas, p_shape)
        
        p_0 = numpy.zeros(p_shape)
        p_0[0] = 1
        
        cme_de = cme(p_shape, (lambda x : lambdas[0],))
        cme_de.integrate(1.0)
        
        p_final = cme_de.y
        
        numpy.testing.assert_array_almost_equal(p_final[:-1],
                                                p_exact[:-1])

    def test_2d_poisson(self):
        """
        testing 2D Poisson process
        """
        import numpy
        
        p_shape = (16, 32)
        lambdas = (5.0, 2.3)
        
        p_exact = get_exact_poisson_p(lambdas, p_shape)
        
        cme_de = cme(p_shape, (lambda *x : lambdas[0],
                               lambda *x : lambdas[1]))
        cme_de.integrate(1.0)
        
        p_final = cme_de.y
        p_final = numpy.reshape(p_final, p_shape)
        
        # note that the last rows/columns are not the same as propensity is set
        # to zero there
        numpy.testing.assert_array_almost_equal(p_final[:-1, :-1],
                                                p_exact[:-1, :-1])

    def test_1d_binomial(self):
        """
        testing 1D Binomial process
        """
        import math
        import numpy
        import scipy.stats

        p_shape = (32,)
        lam = 1.0
        
        pr = 1.0 - math.exp(-lam)
        p_exact = scipy.stats.binom(p_shape[0],
                                    pr).pmf(numpy.arange(p_shape[0]))
        
        p_0 = numpy.zeros(p_shape)
        p_0[0] = 1
        
        cme_de = cme(p_shape, (lambda x : lam*(p_shape[0]-x),))
        cme_de.integrate(1.0)
        p_final = cme_de.y
        
        numpy.testing.assert_array_almost_equal(p_final[:-1],
                                                p_exact[:-1])

    def test_a_to_b_to_a(self):
        """
        testing reversible A <-> B
        with kappa1 = kappa2 = 1.0
        """
        
        from cmepy.models import A2B2A as model
        import math
        import numpy
        import scipy.stats
        
        species_counts = model['species counts']
        
        t_max = 1.0  # don't choose too large !   
        
        # exact distribution of Y0
        pi0  = (1.0 + math.exp(-2*t_max))/2.0
        
        # total number of species
        y_tot = species_counts[0](0, 0) + species_counts[1](0, 0) 
        p_exact  = scipy.stats.binom(y_tot, pi0).pmf(numpy.arange(y_tot+1))
        
        # compute numerical solution
        p_shape = model['np']
        cme_de = cme(p_shape, model['propensities'])
        cme_de.integrate(t_max)
        p_final = numpy.reshape(cme_de.y, p_shape)
        
        # determine the probability distr. for X0 - X1:
        p_0 = p_final[:, 0]
        for i in xrange(1, p_shape[0]):
            p_0[:-i] += p_final[i:, i]
        p_0 = p_0[:(y_tot+1)]
        
        # reverse order of p0, as Y0 = YTOT - X0 + X1:
        p_0 = p_0[::-1]
        
        # test against exact solution
        numpy.testing.assert_array_almost_equal(p_exact, p_0)

    def test_a_to_b_to_c(self):
        """
        testing A -> B -> C
        with kappa1 = kappa2 = 1.0
        """
        
        from cmepy.models import A2B2C as model
        import math
        import numpy
        import scipy.stats
     
        t_max = 0.01  # don't choose too large !   
        y_tot = 31
        p_shape = model['np']
        cme_de = cme(p_shape, model['propensities'])
        
        # exact distribution of Y1  --all marginal probabilities are binomial
        pi_1  = t_max*math.exp(-t_max)
        p_exact = scipy.stats.binom(y_tot, pi_1).pmf(numpy.arange(y_tot+1))
        
        # compute numerical solution
        cme_de.integrate(t_max)
        p_final = numpy.reshape(cme_de.y, p_shape)
        
        # determine the probability distr. for Y1 = X0 - X1:
        p_0 = p_final[:, 0]
        for i in xrange(1, p_shape[0]):
            p_0[:-i] += p_final[i:, i]
        p_0 = p_0[:(y_tot+1)]
        numpy.testing.assert_array_almost_equal(p_exact, p_0)

def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(CMETest)
    return suite

def main():
    test_support.run_unittest(CMETest)

if __name__ == '__main__':
    main()

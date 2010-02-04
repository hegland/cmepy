import unittest
from test import test_support

from scipy.stats import binom, poisson
from itertools import izip
import numpy
from numpy.testing import assert_almost_equal

import cmepy.recorder
import cmepy.solver
from cmepy import model

def exact_poisson(rates, shape):
    poissons = []
    assert len(rates) == len(shape)
    
    for (rate, dim) in izip(rates, shape):
        poissons.append(poisson(rate).pmf(numpy.arange(dim)))
    exact = reduce(numpy.multiply.outer, poissons)
    return exact

def exact_binomial(rate, size):
    return binom(size, 1.0 - numpy.exp(-rate)).pmf(numpy.arange(size))

def exact_monomolecular_aba(t, size):
    pi_0 = (1.0 + numpy.exp(-2.0 * t))/2.0
    return binom(size, pi_0).pmf(numpy.arange(size+1))

def exact_monomolecular_abc(t, size):
    pi_1  = t * numpy.exp(-t)
    return binom(size, pi_1).pmf(numpy.arange(size+1))

def compare_against_exact(model, t, f, p_exact):
    solver = cmepy.solver.create(model, sink = True)
    for t_prime in numpy.linspace(0.0, t, 11):
        solver.step(t_prime)
    p, p_sink = solver.y
    assert_almost_equal(f(p), p_exact)

def compare_against_poisson(rates, shape):
    # generate model & create a solver for it
    size = len(shape)
    def constant_propensity(rate):
        return lambda *x : rate
    props = tuple(constant_propensity(r) for r in rates)
    transitions = tuple((0,)*i + (1,) + (0,)*(size-1-i) for i in xrange(size))
    
    m = model.create(
        propensities = props,
        transitions = transitions,
        shape = shape,
        initial_state = (0, )*size
    )
    
    compare_against_exact(m,
                          t = 1.0,
                          f = lambda p : p.to_dense(shape),
                          p_exact = exact_poisson(rates, shape))

def compare_against_binomial(rate, shape):
    
    m = model.create(
        propensities = (lambda *x : rate*(shape[0]-x[0]), ),
        transitions = ((1, ), ),
        shape = shape,
        initial_state = (0, )
    )
        
    compare_against_exact(m,
                      t = 1.0,
                      f = lambda p : p.to_dense(shape),
                      p_exact = exact_binomial(rate, shape[0]))
    
class SolverTests(unittest.TestCase):
    def test_against_poisson_processes(self):
        """
        compare against analytic solution for Poisson process
        """
        compare_against_poisson((6.0, ), (66, ))
        compare_against_poisson((19.0, ), (512, ))
        compare_against_poisson((30.0, ), (1024, ))
        
        compare_against_poisson((5.0, 2.3, ), (16, 32, ))
        compare_against_poisson((2.1, 0.9, ), (50, 50, ))
        
        compare_against_poisson((1.1, 3.3, 2.2, ), (12, 12, 12, ))
        compare_against_poisson((7.1, 0.1, 4.2, ), (8, 3, 11, ))
    
    def test_against_binomial_processes(self):
        """
        testing 1D Binomial processes
        """

        compare_against_binomial(1.0, (32, ))
        compare_against_binomial(2.9, (47, ))
        compare_against_binomial(4.5, (6, ))
    
    def test_monomolecular_a_to_b_to_a(self):
        """
        testing reversible A <-> B
        with kappa1 = kappa2 = 1.0
        """
        
        from cmepy.models.mono_molecular import A2B2A as m
        
        # model is defined in reaction coordinates, hence t_max must be small,
        # otherwise significant quantities of probability will escape the
        # truncated domain and cause error
        t_max = 1.0   
        
        # compute net copy count
        exact_size = m.species_counts[0](0, 0) + m.species_counts[1](0, 0) 
        
        def f(p):
            recorder = cmepy.recorder.create(
                (m.species,
                 m.species_counts)
            )
            recorder.write(t_max, p)
            d = recorder['A'].distributions[-1]
            return d.to_dense((exact_size+1, ))
        
        compare_against_exact(
            m,
            t_max,
            f,
            exact_monomolecular_aba(t_max, exact_size)
        )
    
    def test_monomolecular_a_to_b_to_c(self):
        """
        testing A -> B -> C
        with kappa1 = kappa2 = 1.0
        """
        
        from cmepy.models.mono_molecular import A2B2C as m
        
        # model is defined in reaction coordinates, hence t_max must be small,
        # otherwise significant quantities of probability will escape the
        # truncated domain and cause error
        t_max = 0.01
        
        # compute net copy count
        exact_size = m.species_counts[0](0, 0) + m.species_counts[1](0, 0)
        def f(p):
            recorder = cmepy.recorder.create(
                (m.species,
                 m.species_counts)
            )
            recorder.write(t_max, p)
            d = recorder['B'].distributions[-1]
            return d.to_dense((exact_size+1, ))
        
        compare_against_exact(
            m,
            t_max,
            f,
            exact_monomolecular_abc(t_max, exact_size)
        )

def suite():
    test_suite = unittest.TestLoader().loadTestsFromTestCase(SolverTests)
    return test_suite

def main():
    test_support.run_unittest(SolverTests)

if __name__ == '__main__':
    main()

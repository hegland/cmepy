"""
Statistical utilities for working with CME solver output.
"""

import operator
import numpy
from cmepy.new_core import domain

class Distribution(dict):
    """
    Distributions are mappings that represent probability distributions
    over a discrete state space.
    """
    def __init__(self, p=None):
        """
        Distribution([p]) -> distribution
        
        Returns distribution instance initialised using mapping p.
        They keys of p define states in the state space, with the associated
        values of p defining the corresponding probabilities.
        
        If p is omitted, defaults to an empty distribution, that is, no states.
        """
        if p is None:
            p = {}
        dict.__init__(self, p)
    
    def _get_statistics(self):
        """
        return mapping of bound statistical methods keyed by method name
        """
        return {
            'expectation' : self.expectation,
            'expected_value' : self.expectation,
            'variance' : self.variance,
            'covariance' : self.covariance,
            'standard_deviation' : self.standard_deviation
        }
    statistics = property(_get_statistics)
    
    def _get_dimension(self):
        """
        infer dimension of state space of this distribution
        """
        # attempt to infer dimension of state space
        if self:
            for state in self.iterkeys():
                try:
                    # self is keyed by vector
                    # arguments, return the length of the first one
                    return len(state)
                except TypeError:
                    # self is keyed by scalar arguments,
                    # treat this as dimension 1
                    return 1
        else:
            # empty distributions have dimension zero
            return 0
    dimension = property(_get_dimension)
    
    def map(self, f, g=None):
        """
        d.map(f [, g]) -> distribution
        
        Returns a copy of the distribution d, with each key replaced by its
        image under f. Any duplicate image keys are merged, with the value of
        the merged key equal to the sum of the values.
        
        If g is supplied, it is used instead of addition to reduce the values of
        duplicate image keys.
        """
        return Distribution(map_distribution(f, self, g))
    
    def expectation(self):
        """
        d.expectation() -> mu
        
        Returns expected state of the distribution d, provided dimension > 0.
        """
        assert self.dimension > 0
        return expectation(self)
    
    def variance(self):
        """
        d.variance() -> sigma_squared
        
        Returns variance of the Distribution d, provided dimension == 1.
        """
        assert self.dimension == 1
        return variance(self)
    
    def covariance(self):
        """
        d.covariance() -> cov
        
        Returns covariance of the distribution d, provided dimension == 2.
        """
        assert self.dimension == 2
        return covariance(self)
    
    def standard_deviation(self):
        """
        d.standard_deviation() -> sigma
        
        Returns std deviation of the Distribution d, provided dimension == 1.
        """
        return numpy.sqrt(self.variance())
    

def map_distribution(f, p, g=None):
    """
    map_distribution(f, p [, g]) -> mapping
    
    Returns a copy of the mapping p, with each key replaced by its
    image under f. Any duplicate image keys are merged, with the value of
    the merged key equal to the sum of the values.
    
    If g is supplied, it is used instead of addition to reduce the values of
    duplicate image keys.
    """
    
    if g is None:
        g = operator.add
    f_p = {}
    for state, probability in p.iteritems():
        f_state = f(state)
        if f_state in f_p:
            f_p[f_state] = g(f_p[f_state], probability)
        else:
            f_p[f_state] = probability
    return f_p

def expectation(p):
    """
    expectation(p) -> mu
    
    Returns the expected value mu, treating the mapping p as a distribution
    p : states -> probabilities.
    """
    
    if type(p) is tuple:
        assert len(p) == 2
        states, probabilities = p
    else:
        states, probabilities = domain.from_mapping(p)
    weighted_states = states * probabilities[numpy.newaxis, :]
    mu = numpy.add.reduce(weighted_states, axis=1)
    
    return mu

def _metavariance(p, exponent):
    """
    _metavariance(p, exponent) -> alpha
    
    Returns alpha := E[ product_i (X_i - mu_i)**alpha ], where
    i ranges over the dimension of the keys of p.
    """
    
    states, probabilities = domain.from_mapping(p)
    mu = expectation((states, probabilities))
    diffs = (states - numpy.asarray(mu)[:, numpy.newaxis])
    if exponent != 1:
        diffs **= exponent
    product = numpy.multiply.reduce(diffs, axis = 0)
    alpha = expectation((product, probabilities))
    assert alpha.shape == (1, )
    return alpha[0]
    
    
def variance(p):
    """
    variance(p) -> sigma_squared
    
    Returns the variance sigma_squared, treating the mapping p as a distribution
    p : states -> probabilities.
    """
    return _metavariance(p, exponent=2)

def covariance(p):
    """
    covariance(p) -> cov
    
    Returns the covariance cov, treating the mapping p as a distribution
    p : states -> probabilities.
    """
    return _metavariance(p, exponent=1)

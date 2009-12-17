import numpy

import cmepy.util

def compute_run_ends(a):
    """
    maps flat array a to a boolean array e, where
    e[i] is defined by (a[i] != a[i+1])
    
    the last element of e is defined to be True
    """
    run_ends = numpy.zeros(numpy.shape(a), dtype = numpy.bool)
    run_ends[-1] = True
    run_ends[:-1] = numpy.not_equal(a[:-1], a[1:])
    return run_ends

def create_coord_sparse_marginal(dim, p, norigin=None):
    np = numpy.shape(p)
    dims = len(np)
    marginal = numpy.array(p)
    axis = 0
    for reduce_dim in xrange(dims):
        if reduce_dim == dim:
            axis += 1
            continue
        marginal = numpy.sum(marginal, axis)
    if norigin is None:
        dim_offset = 0
    else:
        dim_offset = norigin[dim]
    states = numpy.arange(dim_offset, dim_offset+np[dim])
    sparse_marginal = SparseMarginal(states,
                                     marginal)
    return sparse_marginal.compress()
    

def create_sparse_marginal(f, p, norigin=None):
    np = numpy.shape(p)
    dims = len(np)
    all_indices = (slice(None, None, None), )*dims
    states = cmepy.util.indices_ext(np, all_indices, norigin)
    states = [numpy.ravel(coord_vector) for coord_vector in states]
    
    keys = f(*states)
    values = numpy.ravel(p)
    sparse_marginal = SparseMarginal(f(*states),
                                     numpy.array(numpy.ravel(p)))
    
    return sparse_marginal.compress().sum_duplicates()

class SparseMarginal(object):
    def __init__(self, keys, values):
        assert len(keys) == len(values)
        self.keys = keys
        self.values = values
    
    def compress(self, epsilon=0.0):
        """
        removes (key, value) pairs where abs(value) <= epsilon
        
        returns self
        """
        large_value_indices = numpy.abs(self.values)>epsilon
        self.keys = self.keys[large_value_indices]
        self.values = self.values[large_value_indices]
        return self
    
    def sum_duplicates(self):
        """
        replaces keys, values with unique_keys, summed_values, where
            summed_values[i] := sum over j of values[j],
            for all j such that keys[j] == unique_keys[i]
        
        returns self
        """
        sort_indices = numpy.argsort(self.keys)
        sorted_keys = self.keys[sort_indices]
        sorted_values = self.values[sort_indices]
        run_ends = compute_run_ends(sorted_keys)
        unique_keys = sorted_keys[run_ends]
        
        # since there is some finite number of unique keys we can
        # put them into the obvious bijection with the natural numbers (incl. 0)
        # (we need to do this to interface with numpy.bincount, since
        # if bincount is to return a length n array the values of the
        # first array must be integers from 0, .., n-1
        
        sorted_bijected_keys = numpy.zeros(numpy.shape(self.keys),
                                           dtype = numpy.int)
        sorted_bijected_keys[1:] = numpy.add.accumulate(run_ends[:-1])
        
        merged_values = numpy.bincount(sorted_bijected_keys,
                                       sorted_values)
        
        self.keys = unique_keys
        self.values = merged_values
        return self
    
    def expected_value(self):
        """
        Evaluates and returns the expected value of the marginal
        """
        return numpy.add.reduce(self.keys * self.values)
    
    def standard_deviation(self, expected_value=None):
        """
        Evaluates and returns the standard deviation of the marginal
        
        Optionally, a pre-computed expected value may be supplied.
        """
        
        if expected_value is None:
            expected_value = self.expected_value()
        
        x = self.values*((self.keys - expected_value)**2)
        return numpy.sqrt(numpy.sum(x))
        
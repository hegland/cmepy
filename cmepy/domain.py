"""
Utility functions for generating arrays of states in domains
"""

import numpy
import cmepy.util

def from_rect(shape, slices=None, origin=None):
    """
    from_rect(shape [, slices [, origin]]) -> array
    
    Returns array of states in rectangular domain of shape 'shape'.
    
    Works similarly to cmepy.util.indices_ext, except all dimensions
    of the array following the first have been flattened.
    
    See cmepy.util.indices_ext for a discussion of the optional arguments.
    """
    
    indices = cmepy.util.indices_ext(shape, slices, origin)
    flat_coord_arrays_shape = (len(indices), -1)
    return numpy.reshape(indices, flat_coord_arrays_shape)

def from_iter(state_iter):
    """
    from_iter(state_iter) -> array
    
    Returns array of all states from the state iterator 'state_iter'.
    """
    
    return numpy.array(list(state_iter)).transpose()

def from_mapping(state_mapping):
    """
    from_iter(state_mapping) -> state_array, value_array
    
    'state_mapping' is assumed to be a mapping keyed by states.
    
    Returns pair of arrays ('state_array', 'value_array'),
    where the i-th state in 'state_array' is the key
    associated with the i-th element in 'value_array'.
    
    This function is intended to ease the construction
    of state arrays with corresponding value arrays from
    a state mapping. If a value array is not required,
    simply use from_iter on the state mapping instead.
    """
    states = from_iter(state_mapping.iterkeys())
    values = numpy.array(state_mapping.values())
    return states, values

def to_iter(state_array):
    """
    to_iter(state_array) -> iter
    
    Returns iterator over all states in the state array 'state_array'.
    """
    for state in state_array.transpose():
        yield tuple(state)
    return
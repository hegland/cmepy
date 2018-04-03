"""
Common utility routines for Finite State Projection.
"""

import numpy
from cmepy.cme_matrix import non_neg_states
from cmepy import lexarrayset

def grow_domain(domain_states, transitions, depth, validity_test = None):
    """
    Returns domain_states grown by depth along transitions.
    
    Resulting states are filtered by the validity_test. By default,
    only states without a negative coordinate are valid.
    """
    if numpy.size(domain_states) == 0:
        raise ValueError('there must be at least one state to expand')
    
    if validity_test is None:
        validity_test = non_neg_states
    
    expanded_states = domain_states
    for _ in xrange(depth):
        # expand support states by one along each state transition
        for transition in transitions:
            shifted_states = lexarrayset.shift(domain_states, transition)
            # filter invalid states (ie states with a negative coord)
            valid = validity_test(shifted_states)
            shifted_states = shifted_states[:, valid]
            expanded_states = lexarrayset.union(
                expanded_states,
                shifted_states
            )
        domain_states = expanded_states
    
    return expanded_states

"""
A collection of models from the literature
naming   Gou07_a, Gou07_b   where Gou07 refers to the paper
in my bibliography and a, b etc is a counter referring to
the specific model in the paper

TODO: fix doc to give a link to the actual paper ?
"""

import numpy
from cmepy import model
from cmepy.util import non_neg

def create_model_uni_dim(initial_copies = 10, rate = 0.001):
    """
    returns model for P+Q -> PQ with initial copies of P, Q and rate specified. 
    """
    m = model.create(
        name = 'Unidirectional Dimerisation',
        propensities = [lambda x : rate*(initial_copies-x)**2],
        transitions = [(1, )],
        reactions = ['P+Q -> PQ'],
        species = ['P', 'Q', 'PQ'],
        species_counts = (lambda x : initial_copies - x, )*2 + (lambda x : x, ),
        shape = (initial_copies + 1, ),
        origin = (0, )
    )
    return m

def create_model_quad_autocat(max_p=30,
                              max_q=30,
                              fixed_s=True,
                              s_0=10,
                              d_0=2,
                              vol=1.0):
    """
    Creates a species-count based model for the system of reactions:
    
        S->P
        D+P->D+P+P
        P+P->P+Q
        P+Q->Q+Q
        P->*
        Q->*
    
    The copy counts of the species S and D are assumed to be constant,
    with s_0 (default 10) copies of S and d_0 (default 2) copies of D.
    
    If fixed_s is set to False, copy count of species S will no longer be
    constant, and will be decreased by the reaction S->P.
    """
    
    model_name = 'Quadratic Autocatalator (%s S)'
    
    p = lambda *x : x[0]
    q = lambda *x : x[1]
    d = lambda *x : d_0
    
    if fixed_s:
        s = lambda *x : s_0
        model_name %= 'fixed'
        
        transitions = (
            (1, 0),
            (1, 0),
            (-1, 1),
            (-1, 1),
            (-1, 0),
            (0, -1),
        )
        
        shape = (max_p+1, max_q+1)
        origin = (0, )*2
    else:
        model_name %= 'variable'
        s = lambda *x : x[2]
        transitions = (
            (1, 0, -1),
            (1, 0, 0),
            (-1, 1, 0),
            (-1, 1, 0),
            (-1, 0, 0),
            (0, -1, 0),
        )
        
        shape = (max_p+1, max_q+1, s_0+1)
        origin = (0, 0, 3, )
    
    m = model.create(
        name = model_name,
        reactions = (
            'S->P',
            'D+P->D+2P',
            'P+P->P+Q',
            'P+Q->2Q',
            'P->*',
            'Q->*',
        ),
        propensities = (
            lambda *x : 0.002 * (vol ** -1) * s(*x),
            lambda *x : 0.001 * (vol ** -2) * d(*x) * p(*x),
            lambda *x : 0.005 * (vol ** -2) * p(*x) * non_neg(p(*x) - 1) / 2.0,
            lambda *x : 0.004 * (vol ** -2) * p(*x) * q(*x),
            lambda *x : 0.002 * (vol ** -1) * p(*x),
            lambda *x : 0.050 * (vol ** -1) * q(*x),
        ),
        transitions = transitions,
        species = ('P', 'Q', 'S', 'D',
        ),
        species_counts = ( p, q, s, d, ),
        shape = shape,
        origin = origin
    )
    return m

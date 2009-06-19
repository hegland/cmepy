"""
Some models adapted from the Discrete Stochastic Models Test Suite.

See http://code.google.com/p/dsmts/
"""

from numpy import maximum

DSMTS_001_01 = {'np' : (200, 200),
                'propensities' : (lambda x0, x1: 0.1*maximum(x0-x1+100, 0),
                                  lambda x0, x1: 0.11*maximum(x0-x1+100, 0)),
                'doc' : 'Birth-death model (001), variant 01',
                'reactions' : ['X -> 2X', 'X -> *'],
                'species' : ['X'],
                'species counts' : (lambda x0, x1: 100+x0-x1,)}

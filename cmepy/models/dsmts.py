"""
Some models adapted from the Discrete Stochastic Models Test Suite.

See http://code.google.com/p/dsmts/
"""

from cmepy.util import non_neg
from cmepy import model

DSMTS_001_01 = model.create(
    name = 'Birth-death model (001), variant 01',
    propensities = (
        lambda *x: 0.1*x[0],
        lambda *x: 0.11*x[0],
    ),
    reactions = (
        'X -> 2X',
        'X -> *',
    ),
    transitions = (
        (1, ),
        (-1, ),
    ),
    species_counts = (
        lambda *x : x[0],
    ),
    species = (
        'X',
    ),
    shape = (200, ),
    initial_state = (100, )
)
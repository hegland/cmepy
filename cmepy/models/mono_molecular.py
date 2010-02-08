"""
some mono-molecular models
"""

from cmepy.util import non_neg
from cmepy import model

A2B2C = model.create(
    name = 'simplest nontrivial mono-molecular reaction',
    propensities = (
        lambda *x: non_neg(31.0-x[0]),
        lambda *x: non_neg(x[0]-x[1]),
    ),
    transitions = (
        (1, 0),
        (0, 1),
    ),
    reactions = (
        'A->B',
        'B->C',
    ),
    species_counts = (
        lambda *x: non_neg(31-x[0]),
        lambda *x: non_neg(x[0]-x[1]),
        lambda *x: x[1],
    ),
    species = (
        'A',
        'B',
        'C',
    ),
    shape = (32, 32),
    initial_state = (0, 0)
)

A2B2A = model.create(
    name = 'simplest reversible mono-molecular reaction',
    propensities = (
        lambda *x: non_neg(31.0-x[0]+x[1]),
        lambda *x: non_neg(x[0]-x[1]),
    ),
    transitions = (
        (1, 0),
        (0, 1),
    ),
    reactions = (
        'A->B',
        'B->A',
    ),
    species_counts = (
        lambda *x: non_neg(31-x[0]+x[1]),
        lambda *x: non_neg(x[0]-x[1]),
        lambda *x: x[1],
    ),
    species = (
        'A',
        'B',
        'C',
    ),
    shape = (50, 60),
    initial_state = (0, 0)
)
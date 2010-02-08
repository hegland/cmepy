"""
Model definition constants
"""

PROPENSITIES = 'propensities'
TRANSITIONS = 'offset_vectors'
NAME = 'doc'
ORIGIN = 'norigin'
SHAPE = 'np'
SPECIES_NAMES = 'species'
SPECIES_COUNTS = 'species counts'
REACTION_NAMES = 'reactions'
ENTRIES = frozenset([
    PROPENSITIES,
    TRANSITIONS,
    NAME,
    ORIGIN,
    SHAPE,
    SPECIES_NAMES,
    SPECIES_COUNTS,
    REACTION_NAMES,
])

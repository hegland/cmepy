"""
Model definition constants
"""

import cmepy.validate

PROPENSITIES = 'propensities'
TRANSITIONS = 'transitions'
NAME = 'name'
ORIGIN = 'origin'
SHAPE = 'shape'
SPECIES_NAMES = 'species'
SPECIES_COUNTS = 'species_counts'
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

def create(**entries):
    """
    create(**entries) -> model
    """
    
    cmepy.validate.model(entries)

    return entries

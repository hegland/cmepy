"""
Model definition constants
"""

import cmepy.validate

PROPENSITIES = 'propensities'
TRANSITIONS = 'transitions'
NAME = 'name'
INITIAL_STATE = 'initial_state'
SHAPE = 'shape'
SPECIES_NAMES = 'species'
SPECIES_COUNTS = 'species_counts'
REACTION_NAMES = 'reactions'
ENTRIES = frozenset([
    PROPENSITIES,
    TRANSITIONS,
    NAME,
    INITIAL_STATE,
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

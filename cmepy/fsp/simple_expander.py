"""
A simple domain expansion routine for the FSP algorithm.
"""

import cmepy.fsp.util

class SimpleExpander(object):
    """
    Simple FSP expander that expands the entire domain.
    
    The domain is expanded along the given transitions, up to the specified
    depth.
    """
    def __init__(self, transitions, depth):
        """
        Simple FSP expander that expands the entire domain.
    
        The domain is expanded along the given transitions, up to the specified
        depth.
        """
        self.transitions = transitions
        self.depth = depth
    
    def expand(self, **kwargs):
        """
        Returns expanded domain states
        """
        return cmepy.fsp.util.grow_domain(
            kwargs['domain_states'],
            self.transitions,
            self.depth
        )

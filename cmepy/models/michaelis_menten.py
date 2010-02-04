"""
Models for Michaelis-Menten systems
"""

from cmepy.util import non_neg

def __create_mm_simple_model(np = None):
    """
    Creates a model for a simple michaelis-menten system.
    """
    # first, define functions mapping reaction counts x to species counts
    species_c = (lambda *x : non_neg(50 - x[0] - x[1]),
                 lambda *x : non_neg(10-x[0]+x[1]+x[2]),
                 lambda *x : non_neg(x[0]-x[1]-x[2]),
                 lambda *x : x[2],
                )
    # second, define reaction propensities via species counts
    props = (lambda *x : 0.01*species_c[0](*x)*species_c[1](*x),
             lambda *x : 35.0*species_c[2](*x),
             lambda *x : 30.0*species_c[2](*x),
            )
    
    # construct the model
    model = {}
    model['doc'] = 'simple Michaelis-Menten system'
    model['species'] = ('S', 'E', 'C', 'D', )
    model['species counts'] = species_c
    model['reactions'] = ('E+S->C', 'C->E+S', 'C->E+D', )
    model['propensities'] = props
    if np is None:
        model['np'] = (21, 17, 17, )
    else:
        model['np'] = np
    return model

MM_SIMPLE = __create_mm_simple_model()
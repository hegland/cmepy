"""
A collection of models from the literature
naming   Gou07_a, Gou07_b   where Gou07 refers to the paper
in my bibliography and a, b etc is a counter referring to
the specific model in the paper

TODO: fix doc to give a link to the actual paper ?
"""

GOU07_A = {'np' : (2, ),
           'propensities' : (lambda x : 0.001*(1.0-x)**2,),
           'doc' : 'unidirectional dimerisation',
           'reactions' : ['P+Q-> PQ'],
           'species' : ['P','Q', 'PQ'],
           'species counts' : (lambda x : 1.0-x,
                               lambda x : 1.0-x,
                               lambda x : x)}

GOU07_B = {'np' : (11, ),
           'propensities' : (lambda x : 0.001*(10.0-x)**2,),
           'doc' : 'unidirectional dimerisation',
           'reactions' : ['P+Q-> PQ'],
           'species' : ['P','Q', 'PQ'],
           'species counts' : (lambda x : 10.0-x,
                               lambda x : 10.0-x,
                               lambda x : x)}


def __gen_gou07_c(vol, num_d, num_s):
    """
    This appears to be the scheme used in the paper, i.e. the number of copies
    of S is assumed to be set to the constant value 10.
    """
    import numpy
    species_counts = (lambda *x : numpy.maximum(num_s-x[0], 0),
                      lambda *x : numpy.maximum(x[0]+x[1]-x[2]-x[3]-x[4], 0),
                      lambda *x : numpy.maximum(x[2]+x[3]-x[5], 0),
                      )
    propensities = (lambda *x : (0.002*(num_s)/vol),
                    lambda *x : (0.001*num_d*species_counts[1](*x)/(vol*vol)),
                    lambda *x : (0.005*species_counts[1](*x)
                                 *numpy.maximum(species_counts[1](*x)-1, 0)
                                 /(2.0*vol*vol)),
                    lambda *x : (0.004*species_counts[1](*x)
                                 *species_counts[2](*x)/(vol*vol)),
                    lambda *x : (0.002*species_counts[1](*x)/vol),
                    lambda *x : (0.05*species_counts[2](*x)/vol),
                    )
    model = {'doc' : 'quadratic autocatalator',
             'propensities' : propensities,
             'reactions' : ('S->P',
                            'D+P->D+2P',
                            'P+P->P+Q',
                            'P+Q->2Q',
                            'P->*',
                            'Q->*',
                            ),
             'species counts' : species_counts,
             'species' : ('S', 'P', 'Q',),
             'np' : (16, 11, 11, 5, 6, 11, ),
             }
    
    return model

GOU07_C = __gen_gou07_c(vol = 1.0,
                        num_d = 2,
                        num_s = 10)
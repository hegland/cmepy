"""
some mono-molecular models
"""

from numpy import maximum

A2B2C  = {'np' : (32, 32),
          'propensities' : (lambda x1, x2: maximum(31.0-x1, 0.0),
                            lambda x1, x2: maximum(x1-x2, 0.0)),
          'doc' : 'simplest nontrivial monomolecular reaction',
          'reactions' : ['A->B', 'B->C'],
          'species' : ['A', 'B', 'C'],
          'species counts' : (lambda x1, x2: maximum(31-x1,0.0),
                              lambda x1, x2: maximum(x1-x2, 0.0),
                              lambda x1, x2: x2)}
A2B2A  = {'np' : (50, 60),
          'propensities' : (lambda x1, x2: maximum(31.0-x1+x2, 0.0),
                            lambda x1, x2: maximum(x1-x2, 0.0)),
          'doc' : """simplest reversible monomolecular reaction
                     for unit testing set np = (50, 60), and 
                     propensities:
                     (lambda x1, x2: maximum(31.0-x1+x2, 0.0),
                            lambda x1, x2: maximum(x1-x2, 0.0))
                  """,
          'reactions' : ['A->B', 'B->A'],
          'species' : ['A', 'B'],
          'species counts' : (lambda x1, x2: maximum(31-x1+x2,0.0),
                              lambda x1, x2: maximum(x1-x2, 0.0))}


"""
gene toggle model

@conference{munsky2008computation,
  title={{Computation of switch time distributions in stochastic gene regulatory networks}},
  author={Munsky, B. and Khammash, M.},
  booktitle={Proc. 2008 American Control Conference.(June 2008)},
  pages={2761--2766}
}

"""

from numpy import maximum

MUNK08_A = {'np' : (55, 8),
            'propensities' : (lambda *x: 16.0,
                              lambda *x: (x[0]-x[1])),
            'doc' : ('simple precursor to Gardner\'s gene toggle model'
                     +' -- birth-death'),
            'reactions' : '*->S, S->*',
            'species' : ['S'],
            'species counts' : (lambda *x : x[0]-x[1],)}

MUNK08 = {'np' : (18, 8, 75, 30),
          'propensities' : (lambda *x: 16.0/(1.0+maximum(x[2]-x[3],0)),
                            lambda *x: maximum(x[0]-x[1],0),
                            lambda *x: 50.0/(1+maximum(x[0]-x[1],0)**2.5),
                            lambda *x: maximum(x[2]-x[3],0)),
          'doc' : 'Gardner\'s gene toggle model according to MunK08',
          'reactions' : ['*->S1', 'S1->*', '*->S2', 'S2->*'],
          'species' : ['S1', 'S2'],
          'species counts' : (lambda *x : x[0]-x[1],
                              lambda *x : x[2]-x[3])}

"""
Example implementation of the pap-pili epigenetic switch model.

This is the second example from Munsky & Khammash :

    @article{munsky2006finite,
      title={{The finite state projection algorithm for the
              solution of the chemical master equation}},
      author={Munsky, B. and Khammash, M.},
      journal={The Journal of chemical physics},
      volume={124},
      pages={044104},
      year={2006}
    }
"""

import numpy

import cmepy.solver
import cmepy.recorder
import cmepy.model
from cmepy.util import non_neg

def create_model(max_papi, papi_count, dna_count, lrp_count):
    # define mapping from state space to species counts
    c = {
        'DNA' : lambda *x : x[0],
        'DNA-LRP' : lambda *x : x[1],
        'LRP-DNA' : lambda *x : x[2],
        'LRP-DNA-LRP' : lambda *x : non_neg(dna_count - x[0] - x[1] - x[2]),
        'LRP' : lambda *x : non_neg(lrp_count - 2*non_neg(dna_count - x[0] - x[1] - x[2]) - x[1] - x[2]),
        'PapI' : lambda *x : x[3],
    }
    def f(r):
        return r / (1.0 + r)
    # define reaction propensities
    props = (
        lambda *x : 1.0*c['DNA'](*x)*c['LRP'](*x),
        lambda *x : (2.5 - 2.25*f(c['PapI'](*x)))*c['LRP-DNA'](*x),
        lambda *x : 1.0*c['DNA'](*x)*c['LRP'](*x),
        lambda *x : (1.2 - 0.20*f(c['PapI'](*x)))*c['DNA-LRP'](*x),
        lambda *x : 0.01*c['LRP-DNA'](*x)*c['LRP'](*x),
        lambda *x : (1.2 - 0.20*f(c['PapI'](*x)))*c['LRP-DNA-LRP'](*x),
        lambda *x : 0.01*c['DNA-LRP'](*x)*c['LRP'](*x),
        lambda *x : (2.5 - 2.25*f(c['PapI'](*x)))*c['LRP-DNA-LRP'](*x),
        lambda *x : 10.0*c['LRP-DNA'](*x),
        lambda *x : 1.0*c['PapI'](*x),
    )
    # define corresponding reaction state space offsets
    transitions = (
        (-1, 0, 1, 0), # LRP + DNA -> LRP-DNA
        (1, 0, -1, 0), # LRP-DNA -> LRP + DNA
        (-1, 1, 0, 0), # DNA + LRP -> DNA-LRP
        (1, -1, 0, 0), # DNA-LRP -> DNA + LRP
        (0, 0, -1, 0), # LRP-DNA + LRP -> LRP-DNA-LRP
        (0, 0, 1, 0), # LRP-DNA-LRP -> LRP-DNA + LRP
        (0, -1, 0, 0), # LRP + DNA-LRP -> LRP-DNA-LRP
        (0, 1, 0, 0), # LRP-DNA-LRP -> LRP + DNA-LRP
        (0, 0, 0, 1), # * -> PapI
        (0, 0, 0, -1), # PapI -> *
    )
    
    # define state space shape
    shape = (dna_count+1, )*3 + (max_papi+1, )
    
    species_names = []
    species_counts = []
    for key, value in c.iteritems():
        species_names.append(key)
        species_counts.append(value)
    
    initial_state = (dna_count, 0, 0, papi_count)
    return cmepy.model.create(
        name = '',
        species = species_names,
        species_counts = species_counts,
        propensities = props,
        transitions = transitions,
        shape = shape,
        initial_state = initial_state
    )
    
def main():
    m = create_model(
        max_papi = 100,
        papi_count = 5,
        dna_count = 1,
        lrp_count = 100
    )
    
    solver = cmepy.solver.create(
        m,
        sink = True
    )
    recorder = cmepy.recorder.create(
        (cmepy.model.SPECIES_NAMES,
         m[cmepy.model.SPECIES_NAMES],
         m[cmepy.model.SPECIES_COUNTS])
    )
    
    time_steps = numpy.linspace(0.0, 10.0, 101)
    for t in time_steps:
        solver.step(t)
        p, p_sink = solver.y
        print 't = %g; p_sink = %g' % (t, p_sink)
        recorder.write(t, p)
    
    cmepy.recorder.display_plots(recorder,
                                 cmepy.model.SPECIES_NAMES,
                                 title = m[cmepy.model.NAME])

if __name__ == '__main__':
    main()

"""
Munsky & Khammash's stochastic version of Gardner's gene toggle model.

@conference{munsky2008computation,
  title={{Computation of switch time distributions
          in stochastic gene regulatory networks}},
  author={Munsky, B. and Khammash, M.},
  booktitle={Proc. 2008 American Control Conference.(June 2008)},
  pages={2761--2766}
}
"""

import numpy
from cmepy import domain
import cmepy.solver
import cmepy.recorder

def create_model(max_s1_copies, max_s2_copies):
    s1_count = lambda s1, s2 : s1
    s2_count = lambda s1, s2 : s2
    
    s1_birth = lambda s1, s2 : 16.0/(1.0+s2)
    s1_death = lambda s1, s2 : 1.0*s1
    
    s2_birth = lambda s1, s2 : 50.0/(1.0+(s1**2.5))
    s2_death = lambda s1, s2 : 1.0*s2
    
    propensities = (s1_birth, s1_death, s2_birth, s2_death)
    
    offsets = ((1, 0), (-1, 0), (0, 1), (0, -1))
    
    np = (max_s1_copies+1, max_s2_copies+1)
    norigin = (0, 0)
    
    doc = 'Gardner\'s gene toggle according to Munsky & Khammash'
    model = {'doc' : doc,
             'species' : ('S1', 'S2'),
             'species counts' : (s1_count, s2_count),
             'propensities' : propensities,
             'offset_vectors' : offsets,
             'np' : np,
             'norigin' : norigin}
    return model

def main():

    max_s1_copies = 40
    max_s2_copies = 100
    model = create_model(max_s1_copies, max_s2_copies)
    
    # define domain states as union of two rectangular regions along the axes
    a_shape = (max_s1_copies, 6)
    b_shape = (10, max_s2_copies)
    domain_states_a = set(domain.to_iter(domain.from_rect(a_shape)))
    domain_states_b = set(domain.to_iter(domain.from_rect(b_shape)))
    states = domain.from_iter(domain_states_a | domain_states_b)
    
    solver = cmepy.solver.create(
        model,
        sink = True,
        domain_states = states
    )
    
    recorder = cmepy.recorder.create(
        ('species',
         model['species'],
         model['species counts'])
    )
    
    time_steps = numpy.linspace(0.0, 100.0, 101)
    for t in time_steps:
        solver.step(t)
        p, p_sink = solver.y
        
        print '\ttime = %g' % t
        print '\tsink state probability %g' % p_sink
        
        recorder.write(t, p)
    
    cmepy.recorder.display_plots(recorder, 'species', title = model['doc'])


if __name__ == '__main__':
    main()
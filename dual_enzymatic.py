"""
Example system with two enzymatic reactions.

System :
    S + E1 <--> C1 --> P + E1
    P + E2 <--> C2 --> P + E2
"""

import numpy
import cmepy.solver
from cmepy import domain

def gen_states(initial_copies):
    """
    gen_states(initial_copies) -> generator
    
    Returns generator yielding all reachable states in state space.
    """
    # state space format: (c2_copies, c1_copies, s_copies)
    
    s_0 = initial_copies['S']
    e1_0 = initial_copies['E1']
    e2_0 = initial_copies['E2']
    for s in xrange(s_0+1):
        for c1 in xrange(min(e1_0, s_0 - s)+1):
            for c2 in xrange(min(e2_0, s_0 - s - c1)+1): 
                yield (c2, c1, s)
    return

def create_model(initial_copies):
    """
    create_model(initial_copies) -> mapping
    
    Returns mapping storing model.
    """
    # state space format: (c2_copies, c1_copies, s_copies)
    
    s_copies = lambda *x : x[2]
    c1_copies = lambda *x : x[1]
    c2_copies = lambda *x : x[0]
    p_copies = lambda *x : initial_copies['S'] - x[0] - x[1] - x[2]
    e1_copies = lambda *x : initial_copies['E1'] - x[1]
    e2_copies = lambda *x : initial_copies['E2'] - x[0]
    
    species = (
        'S',
        'C1',
        'C2',
        'P',
        'E1',
        'E2',
    )
    
    species_counts = (
        s_copies,
        c1_copies,
        p_copies,
        e1_copies,
        e2_copies,
    )
    
    props = (
        lambda *x : 4.0*s_copies(*x)*e1_copies(*x),
        lambda *x : 5.0*c1_copies(*x),
        lambda *x : 1.0*c1_copies(*x),
        lambda *x : 4.0*p_copies(*x)*e2_copies(*x),
        lambda *x : 5.0*c2_copies(*x),
        lambda *x : 1.0*c2_copies(*x),
    )
    
    offsets = (
        (0, 1, -1),
        (0, -1, 1),
        (0, -1, 0),
        (1, 0, 0),
        (-1, 0, 0),
        (-1, 0, 1),
    )
    
    doc = 'dual enzymatic reactions'
    origin = (0, 0, initial_copies['S'])
    
    model = {
        'propensities' : props,
        'offset_vectors' : offsets,
        'species' : species,
        'species counts' : species_counts,
        'doc' : doc,
        'norigin' : origin,
    }
    
    return model

def save_p(p, file_name):
    f = open(file_name, 'w')
    keys = p.keys()
    keys.sort()
    for key in keys:
        f.write('%s : %g\n' % (key, p[key]))
    f.close()

def main():
    """
    Solve dual enzymatic reaction system
    """
    initial_copies = {
        'S' : 100,
        'E1' : 50,
        'E2' : 10,
    }
    model = create_model(initial_copies)
    states = domain.from_iter(gen_states(initial_copies))
    
    solver = cmepy.solver.create(
        model,
        sink = True,
        domain_states = states,
    )
    
    t_final = 10.0
    steps_per_time = 100
    time_steps = numpy.linspace(0.0, t_final, int(steps_per_time*t_final) + 1)
    for i, t in enumerate(time_steps):
        print '\tt = %f' % t
        solver.step(t)
        
        p, p_sink = solver.y
        
        if i % 10 == 0:
            # save every 10-th time step to disk
            file_name = 'dual_enzymatic_results_%f.txt' % t
            save_p(p, file_name)

if __name__ == '__main__':
    import cProfile, pstats
    PROFILE_FILE = 'dual_enzymatic.profile'
    cProfile.run('main()', PROFILE_FILE)
    stats = pstats.Stats(PROFILE_FILE)
    stats.sort_stats('cumulative').print_stats(30)

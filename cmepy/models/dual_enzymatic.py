"""
Example system with two enzymatic reactions.

System :
    S + E1 <--> C1 --> P + E1
    P + E2 <--> C2 --> P + E2
"""

from cmepy import model

def gen_states(initial_copies = None):
    """
    gen_states([initial_copies]) -> generator
    
    Returns generator yielding all reachable states in state space.
    
    NB state space format: (c2_copies, c1_copies, s_copies)
    """
    
    if initial_copies is None:
        initial_copies = {
            'S' : 100,
            'E1' : 50,
            'E2' : 10,
        }
    
    s_0 = initial_copies['S']
    e1_0 = initial_copies['E1']
    e2_0 = initial_copies['E2']
    for s in xrange(s_0+1):
        for c1 in xrange(min(e1_0, s_0 - s)+1):
            for c2 in xrange(min(e2_0, s_0 - s - c1)+1): 
                yield (c2, c1, s)
    return

def create_model(initial_copies = None):
    """
    create_model([initial_copies]) -> mapping
    
    Returns mapping storing model.
    
    NB state space format: (c2_copies, c1_copies, s_copies)
    """
    
    if initial_copies is None:
        initial_copies = {
            'S' : 100,
            'E1' : 50,
            'E2' : 10,
        }
    
    s_copies = lambda *x : x[2]
    c1_copies = lambda *x : x[1]
    c2_copies = lambda *x : x[0]
    p_copies = lambda *x : initial_copies['S'] - x[0] - x[1] - x[2]
    e1_copies = lambda *x : initial_copies['E1'] - x[1]
    e2_copies = lambda *x : initial_copies['E2'] - x[0]
    
    return model.create(
        name = 'dual enzymatic reactions',
        reactions = (
            'S+E1 -> C1',
            'C1 -> S+E1',
            'C1 -> P+E1',
            'P+E2 -> C2',
            'C2 -> P+E2',
            'C2 -> S+E2',
        ),
        propensities = (
            lambda *x : 4.0*s_copies(*x)*e1_copies(*x),
            lambda *x : 5.0*c1_copies(*x),
            lambda *x : 1.0*c1_copies(*x),
            lambda *x : 4.0*p_copies(*x)*e2_copies(*x),
            lambda *x : 5.0*c2_copies(*x),
            lambda *x : 1.0*c2_copies(*x),
        ),
        transitions = (
            (0, 1, -1),
            (0, -1, 1),
            (0, -1, 0),
            (1, 0, 0),
            (-1, 0, 0),
            (-1, 0, 1),
        ),
        species = (
            'S',
            'C1',
            'C2',
            'P',
            'E1',
            'E2',
        ),
        species_counts = (
            s_copies,
            c1_copies,
            c2_copies,
            p_copies,
            e1_copies,
            e2_copies,
        ),
        initial_state = (0, 0, initial_copies['S'])
    )

def main():
    """
    Solve dual enzymatic reaction system
    """
    
    import numpy
    import cmepy.solver
    from cmepy import domain
    
    solver = cmepy.solver.create(
        model = create_model(),
        sink = True,
        domain_states = domain.from_iter(gen_states()),
    )
    
    t_final = 10.0
    steps_per_time = 100
    time_steps = numpy.linspace(0.0, t_final, int(steps_per_time*t_final) + 1)
    
    print 'solving dual enzymatic system (this may take some time ...)'
    for t in time_steps:
        print '\t = %f' % t
        solver.step(t)

def create_catalytic_model():
    """
    creates a model for the catalytic reaction example
    """
    import cmepy.model
    
    a_initial = 15
    d_initial = 20
    
    s_1 = lambda *x : a_initial - x[0]
    s_2 = lambda *x : x[0] - x[1]
    s_3 = lambda *x : x[1]
    s_4 = lambda *x : d_initial - x[2]
    s_5 = lambda *x : x[2]

    return cmepy.model.create(
        name = 'catalytic reaction',
        initial_state = (0, )*3,
        shape = (a_initial + 1, )*2 + (d_initial + 1,),
        species = ('A', 'B', 'C', 'D', 'E'),
        species_counts = (s_1, s_2, s_3, s_4, s_5),
        reactions = ('A->B', 'B->C', 'B+D->B+E'),
        propensities = (
            lambda *x : 1.0 * s_1(*x),
            lambda *x : 1000.0 * s_2(*x),
            lambda *x : 100.0 * s_4(*x) * s_2(*x)
        ),
        transitions = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    )
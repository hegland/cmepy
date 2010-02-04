"""
model of two competing clonotypes
"""

from cmepy import model

def create_model_competing_clonotypes():
    """
    create species count state space version of the competing clonotypes model
    """
    
    shape = (50, 50)
     
    # we first define the mappings from the state space to species counts
    # this is pretty easy since we choose species count state space
    species_count_a = lambda *x : x[0]
    species_count_b = lambda *x : x[1]
    
    # we now define the reaction propensities using the species counts
    def reaction_a_birth(*x):
        """
        propensity of birth reaction for species a
        """
        s_a = species_count_a(*x)
        s_b = species_count_b(*x)
        return numpy.where(s_a + s_b > 0,
                           60.0*s_a*(numpy.divide(0.5, s_a + s_b) +
                                   numpy.divide(0.5, (s_a + 10*100))),
                           0.0)
    
    def reaction_a_decay(*x):
        return 1.0*species_count_a(*x)
    
    def reaction_b_birth(*x):
        """
        propensity of birth reaction for species b
        """
        s_a = species_count_a(*x)
        s_b = species_count_b(*x)
        return numpy.where(s_a + s_b > 0,
                           60.0*s_b*(numpy.divide(0.5, s_a + s_b) +
                                   numpy.divide(0.5, (s_b + 10*100))),
                           0.0)
    
    def reaction_b_decay(*x):
        return 1.0*species_count_b(*x)
    
    return model.create(
        name = 'T Cell clonoTypes',
        reactions = (
            '*->A',
            'A->*',
            '*->B',
            'B->*',
        ),
        propensities = (
            reaction_a_birth,
            reaction_a_decay,
            reaction_b_birth,
            reaction_b_decay,
        ),
        transitions = (
            (1, 0),
            (-1, 0),
            (0, 1),
            (0, -1),
        ),
        species = (
            'A',
            'B',
        ),
        species_counts = (
            species_count_a,
            species_count_b,
        ),
        shape = shape,
        origin = (10, 10)
    )

def create_time_deps_competing_clonotypes():
    """
    create time dependencies for the competing clonotypes model
    """
    # 0-th and 2-nd reactions are scaled by the following
    # time dependent factor
    return {frozenset([0, 2]) : lambda t : math.exp(-0.1*t)}

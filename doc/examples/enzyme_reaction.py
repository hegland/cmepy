from cmepy import model
from cmepy.util import non_neg

s_0 = 50
e_0 = 10

# species count function definitions
s = lambda *x : x[0]
e = lambda *x : non_neg(e_0 - x[1])
c = lambda *x : x[1]
p = lambda *x : non_neg(s_0 - x[0] - x[1])

# model definition, in terms of species count functions
m = model.create(
    name = 'enzyme kinetics',
    species = ('S', 'E', 'C', 'P', ),
    species_counts = (s, e, c, p, ),
    reactions = ('E+S->C', 'C->E+S', 'C->E+P', ),
    propensities = (
        lambda *x : 0.01*s(*x)*e(*x),
        lambda *x : 35.0*c(*x),
        lambda *x : 30.0*c(*x),
    ),
    transitions = (
        (-1, 1),
        (1, -1),
        (0, -1)
    ),
    shape = (s_0 + 1, min(s_0, e_0) + 1),
    initial_state = (s_0, 0)
)

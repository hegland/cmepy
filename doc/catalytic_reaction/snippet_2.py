reactions = ('A->B', 'B->C', 'B+D->B+E')

propensities = (
    lambda *x : 1.0 * s_1(*x),
    lambda *x : 1000.0 * s_2(*x),
    lambda *x : 100.0 * s_4(*x) * s_2(*x)
)

transitions = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
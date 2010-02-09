s_1 = lambda *x : a_initial - x[0]
s_2 = lambda *x : x[0] - x[1]
s_3 = lambda *x : x[1]
s_4 = lambda *x : d_initial - x[2]
s_5 = lambda *x : x[2]

species = ('A', 'B', 'C', 'D', 'E')
species_counts = (s_1, s_2, s_3, s_4, s_5)
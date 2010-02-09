import numpy
import cmepy.recorder
import cmepy.solver

def create_catalytic_model():
    """
    creates a model for the catalytic reaction example
    """
    pass # TODO CREATE MODEL HERE

m = create_catalytic_model()
solver = cmepy.solver.create(
    model = m,
    sink = False
)

recorder = cmepy.recorder.create(
    (m.species, m.species_counts)
)

time_steps = numpy.linspace(0.0, 0.5, 101)
for t in time_steps:
    solver.step(t)
    recorder.write(t, solver.y)

e_ev = recorder['E'].expected_value[-1]
print 'expected copy count of species E at final time: %g' % e_ev

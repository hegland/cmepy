import numpy

import cmepy.solver
import cmepy.recorder
from cmepy.models.gou07 import create_model_quad_autocat

def main():
    
    model = create_model_quad_autocat(fixed_s = False)
    
    solver = cmepy.solver.create(
        model,
        sink = True
    )
    recorder = cmepy.recorder.create(
        ('species',
         model['species'],
         model['species_counts'])
    )
    time_steps = numpy.linspace(0.0, 1000.0, 101)
    for t in time_steps:
        solver.step(t)
        p, p_sink = solver.y
        recorder.write(t, p)
        print 't = %g; p_sink = %g' % (t, p_sink)
    
    cmepy.recorder.display_plots(recorder, 'species', title = model['name'])

if __name__ == '__main__':
    main()
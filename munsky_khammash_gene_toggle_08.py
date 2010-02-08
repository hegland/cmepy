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
import cmepy.new_core.cme_solver as cme_solver
import cmepy.new_core.recorder as cme_recorder
import pylab

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

def display_plots(recorder, title):
    pylab.figure()
    for measurement in recorder.measurements('species'):
        pylab.plot(measurement.times,
                   measurement.expected_value,
                   label = measurement.name)
    pylab.legend()
    pylab.title(title+': species count expected value')
    
    pylab.figure()
    for measurement in recorder.measurements('species'):
        pylab.plot(measurement.times,
                   measurement.standard_deviation,
                   label = measurement.name)
    pylab.legend()
    pylab.title(title+': species count standard deviation')

def main():

    max_s1_copies = 100
    max_s2_copies = 100
    
    model = create_model(max_s1_copies, max_s2_copies)
    
    solver = cme_solver.create_cme_solver(model)
    recorder = cme_recorder.CmeRecorder(model)
    recorder.add_target('species',
                        ['expected value', 'standard deviation'],
                        model['species'],
                        model['species counts'])
    time_steps = numpy.linspace(0.0, 100.0, 101)
    for t in time_steps:
        print ('\tt = %g' % t)
        solver.step(t)
        recorder.write(t, solver.y)
        bdry_mass = numpy.add.reduce(numpy.abs(numpy.ravel(solver.y[-1, :]))) + numpy.add.reduce(numpy.abs(numpy.ravel(solver.y[:, -1])))
        print '\tbdry_mass %g' % bdry_mass
    
    display_plots(recorder, model['doc'])
    pylab.show()


if __name__ == '__main__':
    main()
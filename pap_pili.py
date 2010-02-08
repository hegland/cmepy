"""
Example implementation of the pap-pili epigenetic switch model.

This is the second example from Munsky & Khammash :

    @article{munsky2006finite,
      title={{The finite state projection algorithm for the
              solution of the chemical master equation}},
      author={Munsky, B. and Khammash, M.},
      journal={The Journal of chemical physics},
      volume={124},
      pages={044104},
      year={2006}
    }
"""

import numpy

import cmepy.new_core.cme_solver as cme_solver
import cmepy.new_core.recorder as cme_recorder

import pylab

def create_pap_pili_model(max_papi, papi_count, dna_count, lrp_count):
    # define mapping from state space to species counts
    c = {'DNA' : lambda *x : x[0],
         'DNA-LRP' : lambda *x : x[1],
         'LRP-DNA' : lambda *x : x[2],
         'LRP-DNA-LRP' : lambda *x : dna_count - x[0] - x[1] - x[2],
         'LRP' : lambda *x : lrp_count - 2*(dna_count - x[0]) + x[1] + x[2],
         'PapI' : lambda *x : x[3], }
    def f(r):
        return r / (1.0 + r)
    # define reaction propensities
    props = (lambda *x : 1.0*c['DNA'](*x)*c['LRP'](*x),
             lambda *x : (2.5 - 2.25*f(c['PapI'](*x)))*c['LRP-DNA'](*x),
             lambda *x : 1.0*c['DNA'](*x)*c['LRP'](*x),
             lambda *x : (1.2 - 0.20*f(c['PapI'](*x)))*c['DNA-LRP'](*x),
             lambda *x : 0.01*c['LRP-DNA'](*x),
             lambda *x : (1.2 - 0.20*f(c['PapI'](*x)))*c['LRP-DNA-LRP'](*x),
             lambda *x : 0.01*c['DNA-LRP'](*x),
             lambda *x : (2.5 - 2.25*f(c['PapI'](*x)))*c['LRP-DNA-LRP'](*x),
             lambda *x : 10.0*c['LRP-DNA'](*x),
             lambda *x : 1.0*c['PapI'](*x), )
    # define corresponding reaction state space offsets
    offsets = ((-1, 0, 1, 0),
               (1, 0, -1, 0),
               (-1, 1, 0, 0),
               (-1, -1, 0, 0),
               (0, 0, -1, 0),
               (0, 0, 1, 0),
               (0, -1, 0, 0),
               (0, 1, 0, 0),
               (0, 0, 0, 1),
               (0, 0, 0, -1), )
    
    # define state space shape
    np = (dna_count+1, )*3 + (max_papi+1, )
    
    species_names = []
    species_counts = []
    for key, value in c.iteritems():
        species_names.append(key)
        species_counts.append(value)
    
    model = {'doc' : 'pap-pili epigenetic switch',
             'species' : species_names,
             'species counts' : species_counts,
             'propensities' : props,
             'offset_vectors' : offsets,
             'np' : np, }
    p_0 = numpy.zeros(np)
    p_0[-1, 0, 0, papi_count] = 1.0
    return model, p_0

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

def main_pap_pili():
    model, p_0 = create_pap_pili_model(max_papi = 50,
                                       papi_count = 5,
                                       dna_count = 1,
                                       lrp_count = 100)
    
    solver = cme_solver.create_cme_solver(model, p_0)
    recorder = cme_recorder.CmeRecorder(model)
    recorder.add_target('species',
                        ['expected value', 'standard deviation'],
                        model['species'],
                        model['species counts'])
    
    time_steps = numpy.linspace(0.0, 10.0, 101)
    for t in time_steps:
        print ('stepping to t = %f' % t)
        solver.step(t)
        recorder.write(t, solver.y)
    
    display_plots(recorder, model['doc'])
    pylab.show()

if __name__ == '__main__':
    main_pap_pili()

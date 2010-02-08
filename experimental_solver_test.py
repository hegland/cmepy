"""
*** version of munsky_khammash_gene_toggle using experimental solver ***

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
import cmepy.new_core.domain as domain
import cmepy.new_core.state_enum as state_enum
import cmepy.new_core.cme_solver_experimental as cme_solver
import cmepy.new_core.recorder as cme_recorder
import pylab

from munsky_khammash_gene_toggle_08 import create_model

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
    model['norigin'] = (0, 0)
    
    
    # explicitly create our own domain_enum so we have access to the
    # enumeration order of the states in the domain
    
    # dense version
    #domain_enum = state_enum.create(domain.from_rect(model['np']))
    
    # cheapo sparse version
    
    part_a_shape = (model['np'][0], 10)
    part_b_shape = (10, model['np'][1])
    state_space_part_a = set(domain.to_iter(domain.from_rect(part_a_shape)))
    state_space_part_b = set(domain.to_iter(domain.from_rect(part_b_shape)))
    state_space = domain.from_iter(state_space_part_a | state_space_part_b)
    
    domain_enum = state_enum.create(state_space)
    solver = cme_solver.create_cme_solver(model,
                                          sink = True,
                                          domain_enum = domain_enum)
    
    print 'domain_enum.size : %d' % domain_enum.size
    
    recorder = cme_recorder.CmeRecorder(model)
    recorder.add_target('species',
                        ['expected value', 'standard deviation'],
                        model['species'],
                        model['species counts'])
    time_steps = numpy.linspace(0.0, 100.0, 101)
    for t in time_steps:
        print ('\tt = %g' % t)
        solver.step(t)
        p_flat, p_sink = solver.y
        print 'p net probability %g' % numpy.add.reduce(p_flat)
        print 'sink probability %g' % p_sink
        
        # translate solution from sparse mapping format to
        # rectangular dense array, so recorder can understand it
        
        p_sparse = domain_enum.unpack_distribution(p_flat)
        p_dense = numpy.zeros(model['np'])
        for state, probability in p_sparse.iteritems():
            p_dense[state] = probability
        
        recorder.write(t, p_dense)
        
    display_plots(recorder, model['doc'])
    pylab.show()


if __name__ == '__main__':
    main()
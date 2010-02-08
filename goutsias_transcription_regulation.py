"""
the transcription regulation example, as formulated by:

@article{goutsias2005quasiequilibrium,
  title={{Quasiequilibrium approximation of fast reaction
          kinetics in stochastic biochemical systems}},
  author={Goutsias, J.},
  journal={The Journal of chemical physics},
  volume={122},
  pages={184102},
  year={2005}
}

a time-independent approximation of this example was also considered by

@conference{burrage2006krylov,
  title={{A Krylov-based finite state projection algorithm for
          solving the chemical master equation arising in the
          discrete modelling of biological systems}},
  author={Burrage, K. and Hegland, M. and Macnamara, S. and Sidje, R.},
  booktitle={Proc. of The AA Markov 150th Anniversary Meeting},
  pages={21--37},
  year={2006}
}
"""

import itertools

import numpy

import cmepy.new_core.ode_solver as ode_solver
import cmepy.new_core.cme_solver as cme_solver
import cmepy.new_core.recorder as cme_recorder

import cmepy.util

import pylab

def create_optimised_flux_data(np, flux_data, epsilon, initial_state_indices):
    
    flat_state_indices = numpy.arange(numpy.product(np))
    state_indices = numpy.reshape(flat_state_indices, np)
    
    new_flux_data = []
    print 'removing insignificant coefficients'
    for coeffs, source, dest in flux_data:
        coeffs = numpy.ravel(numpy.array(coeffs))
        
        source = numpy.ravel(numpy.array(state_indices[source]))
        dest = numpy.ravel(numpy.array(state_indices[dest]))
        
        print '%s %s %s' % (str(numpy.shape(coeffs)),
                            str(numpy.shape(source)),
                            str(numpy.shape(dest)))
        assert numpy.shape(coeffs) == numpy.shape(source) == numpy.shape(dest)
        
        large_coeffs = coeffs > epsilon
        num_total = numpy.size(large_coeffs)
        num_removed = num_total - numpy.add.reduce(numpy.ravel(large_coeffs))
        print '\tremoved %d of %d coefficients' % (num_removed, num_total)
        new_flux_data.append((numpy.array(coeffs[large_coeffs]),
                             numpy.array(source[large_coeffs]),
                             numpy.array(dest[large_coeffs])))
    print 'reachability analysis'
    
    reachable_states = set()
    open = set(initial_state_indices)
    def reachable(s):
        return s in reachable_states
    v_reachable = numpy.vectorize(reachable)
    
    edges = {}
    for c, s, d in new_flux_data:
        for ss, dd in itertools.izip(list(s), list(d)):
            if ss in edges:
                edges[ss].add(dd)
            else:
                edges[ss] = set([dd])
    
    while len(open)>0:
        if len(reachable_states)%1000 == 0:
            print '\treachability : reachable %d open %d' % (len(reachable_states), len(open))
        state = open.pop()
        reachable_states.add(state)
        if state in edges:
            open.update(edges[state].difference(reachable_states))
    
    print 'applying compressed indexing scheme'
    # convert from set to array
    v_reachable_states = numpy.array(list(reachable_states))
    
    compressed_flux_data = []
    for c, s, d in new_flux_data:
        mask = v_reachable(s)
        c_prime = c[mask]
        print '\t%d of %d flux coeffs are to reachable states, compressing indexing' % (numpy.size(c_prime), numpy.size(c))
        s_prime = numpy.searchsorted(v_reachable_states, s[mask])
        d_prime = numpy.searchsorted(v_reachable_states, d[mask])
        compressed_flux_data.append((c_prime, s_prime, d_prime))
    
    print 'defining pack / unpack functions'
    full_mask = v_reachable(flat_state_indices)
    full_transform = numpy.searchsorted(v_reachable_states, flat_state_indices[full_mask])
    full_transform_inverse = numpy.argsort(full_transform)
    
    def compressed_pack(p):
        return (numpy.ravel(p)[full_mask])[full_transform]
    
    def compressed_unpack(y):
        p = numpy.zeros(numpy.product(np))
        p[full_mask] = y[full_transform_inverse]
        return numpy.reshape(p, np)
        
    return compressed_flux_data, compressed_pack, compressed_unpack    
        

def non_neg(x):
    return numpy.where(x>0.0, x, 0.0)

def create_transcription_regulation_time_dependencies():
    v_0 = 1.0e-15 # litres
    T = 35*60.0 # seconds
    avogadro_number = 6.0221415e23
    
    def volume(t):
        return v_0*numpy.exp(numpy.log(2.0)*t/T)
    
    def make_time_dep(k, divisor):
        if divisor:  
            def time_dep(t):
                return k/(avogadro_number*volume(t))
        else:
            def time_dep(t):
                return k
        return time_dep
    
    reaction_coefficients = [(4.3e-2, False),
                             (7.0e-4, False),
                             (7.8e-2, False),
                             (3.9e-3, False),
                             (1.2e7, True),
                             (4.791e-1, False),
                             (1.2e5, True),
                             (8.765e-12, False),
                             (1.0e8, True),
                             (0.5, False)]
    
    return [make_time_dep(k, divisor) for k, divisor in reaction_coefficients]

def create_transcription_regulation_model(dna_count):
    c = {'DNA' : lambda *x : x[0],
         'DNA-D' : lambda *x : x[1],
         'DNA-2D' : lambda *x : dna_count - x[0] - x[1],
         'RNA' : lambda *x : x[2],
         'M' : lambda *x : x[3],
         'D' : lambda *x : x[4], }
    
    props = (lambda *x : c['RNA'](*x),
             lambda *x : c['M'](*x),
             lambda *x : c['DNA-D'](*x),
             lambda *x : c['RNA'](*x), 
             lambda *x : c['D'](*x)*c['D'](*x), 
             lambda *x : c['DNA-D'](*x), 
             lambda *x : c['DNA-D'](*x)*c['D'](*x), 
             lambda *x : c['DNA-2D'](*x), 
             lambda *x : 0.5*c['M'](*x)*non_neg(c['M'](*x)-1), 
             lambda *x : c['D'](*x), )
    
    offsets = ((0, 0, 0, 1, 0),     # RNA -> RNA + M
               (0, 0, 0, -1, 0),    # M -> *
               (0, 0, 1, 0, 0),     # DNA-D -> RNA + DNA-D
               (0, 0, 1, 0, 0),     # RNA -> *
               (-1, 1, 0, 0, -1),   # DNA + D -> DNA-D
               (1, -1, 0, 0, 1),    # DNA-D -> DNA + D
               (0, -1, 0, 0, -1),   # DNA-D + D -> DNA-2D
               (0, 1, 0, 0, 1),     # DNA-2D -> DNA-D + D
               (0, 0, 0, -2, 1),    # M + M -> D
               (0, 0, 0, 2, -1), )  # D -> M + M
    
    species_names = []
    species_counts = []
    for key, value in c.iteritems():
        species_names.append(key)
        species_counts.append(value)
    
    model = {'doc' : 'Goutsias transcription regulation',
             'species' : species_names,
             'species counts' : species_counts,
             'propensities' : props,
             'offset_vectors' : offsets, }
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

def main_simple(graph=False):
    dna_count = 2
    rna_max = 15
    m_max = 15
    d_max = 15
    model = create_transcription_regulation_model(dna_count)
    np = (dna_count+1, )*2 + (rna_max+1, m_max+1, d_max+1)
    model['np']= np
    
    p_0 = numpy.zeros(np)
    p_0[-1, 0, 0, 2, 6] = 1.0
    
    time_dependencies = create_transcription_regulation_time_dependencies()
    
    flux_evaluator = cme_solver.create_time_dependent_flux_evaluator(time_dependencies)
    flux_data = cme_solver.create_flux_data(model)
    
    optimise = True
    
    if optimise:
        # attempt to optimise flux data & indexing to use less space ...
        # (cull unreachable states from state space, etc ...)
        # (introduces explicit indirect state indexing behind the scenes ...)
        # this multiplies the running time of the integration step
        # by a factor of 11/8...
        
        initial_state_indices = numpy.nonzero(numpy.ravel(p_0))[0]
        opt = create_optimised_flux_data(np,
                                         flux_data,
                                         epsilon = 0.0,
                                         initial_state_indices = initial_state_indices)
        compressed_flux_data, compressed_pack, compressed_unpack = opt
        dy_dt = cme_solver.create_diff_eqs(compressed_flux_data, flux_evaluator)
        solver = ode_solver.Solver(dy_dt, p_0)
        solver.set_packing(compressed_pack,
                           compressed_unpack,
                           transform_dy_dt=False)
    else:
        dp_dt = cme_solver.create_diff_eqs(flux_data, flux_evaluator)
        solver = ode_solver.Solver(dp_dt, p_0)
        pack, unpack = cme_solver.create_packing_functions(model)
        solver.set_packing(pack, unpack)
    
    recorder = cme_recorder.CmeRecorder(model)
    recorder.add_target('species',
                        ['expected value', 'standard deviation'],
                        model['species'],
                        model['species counts'])
    
    t_final = 1*60.0 
    time_steps = numpy.linspace(0.0, t_final, 1*60+1)
    
    for t in time_steps:
        print ('stepping to t = %f' % t)
        solver.step(t)
        recorder.write(t, solver.y)
        print 'boundary mass:'
        p = solver.y
        mass = numpy.add.reduce(numpy.ravel(p[:, :, -1, :, :]))
        print '\tRNA bdry mass : %g' % mass
        mass = numpy.add.reduce(numpy.ravel(p[:, :, :, -1, :]))
        print '\tM bdry mass : %g' % mass
        mass = numpy.add.reduce(numpy.ravel(p[:, :, :, :, -1]))
        print '\tD bdry mass : %g' % mass
        # check out how evenly distributed the mass is ...
        sorted_cumulative_mass = numpy.add.accumulate(numpy.sort(numpy.ravel(p)))
        spread = numpy.array([10.0**(-x) for x in xrange(10)])
        spread_indices = numpy.searchsorted(sorted_cumulative_mass, spread)
        pylab.figure()
        pylab.plot(numpy.size(p) - spread_indices)
        pylab.savefig('./dump/scm_%g.png' % t)
        pylab.close()
    if graph:
        display_plots(recorder, model['doc'])
        pylab.show()

if __name__ == '__main__':
    import cProfile, pstats
    PROFILE_FILE = 'goutsias_transcription_regulation.profile'
    cProfile.run('main_simple()', PROFILE_FILE)
    stats = pstats.Stats(PROFILE_FILE)
    stats.sort_stats('cumulative').print_stats(30)
import numpy

import scipy.sparse
import scipy.linalg

import cmepy.models
import cmepy.new_core.cme_solver as cme_solver
import cmepy.new_core.ode_solver as ode_solver
import cmepy.new_core.recorder as cme_recorder

import cmepy.core.matrix_cme as matrix_cme
import perturb

def approx_cme_solver(full_model, slow_reactions, fast_reactions, epsilon, k):
    reaction_subsets = (slow_reactions, fast_reactions)
    (slow_model, fast_model) = perturb.decompose_model(full_model,
                                                       reaction_subsets)
    
    slow_flux_data = cme_solver.create_flux_data(slow_model)
    slow_matrix = matrix_cme.gen_sparse_matrix(full_model['np'],
                                               slow_flux_data)

    fast_flux_data = cme_solver.create_flux_data(fast_model)
    fast_matrix = matrix_cme.gen_sparse_matrix(full_model['np'],
                                               fast_flux_data)
    fast_matrix *= epsilon

    
    # m := limit of exp(fast_matrix *t) as t --> +ive infty
    # approximate limit by using a large T
    T_INFTY = 100000.0
    dense_fast_matrix = fast_matrix.todense()
    m_matrix = scipy.linalg.expm(dense_fast_matrix*T_INFTY)
    
    print 'shape of m matrix : %s' % str(numpy.shape(m_matrix))
    
    # compute truncated svd of k largest singular values
    u_bar, s_bar, v_bar = perturb.truncated_svd(m_matrix, k)
    
    # use svd to compute aggregation and disaggregation matrices
    f = u_bar
    f_sparse = scipy.sparse.csr_matrix(f)
    e = numpy.dot(numpy.diag(s_bar), v_bar)
    e_sparse = scipy.sparse.csr_matrix(e)
    
    # define initial distribution as system with
    # maximum copies of first species and 0 copies of second species
    p_0 = numpy.zeros(full_model['np'])
    p_0[-1, 0] = 1.0
    
    a_hat = numpy.dot(e, slow_matrix*f)
    a_hat_sparse = scipy.sparse.csr_matrix(a_hat)
    
    #print 'f : '
    #print str(f)
    #print 'e : '
    #print str(e)
    #print 'a hat:'
    #print str(a_hat)
    
    pack, unpack = cme_solver.create_packing_functions(fast_model)
    
    def pack_aggregate(p):
        return e_sparse*pack(p)
    
    def deaggregate_unpack(y):
        return unpack(f_sparse*y)
    
    def dy_dt(t, y):
        return a_hat_sparse*y
    
    new_solver = ode_solver.Solver(dy_dt, p_0)
    new_solver.set_packing(pack_aggregate,
                           deaggregate_unpack,
                           transform_dy_dt=False)
    return new_solver

def test():
    copy_count = 30
    k = copy_count+1
    slow_reactions = (0, )
    fast_reactions = (1, )
    epsilon = 0.01
    
    
    model = perturb.get_simple_model(copy_count)
    
    solver = approx_cme_solver(model,
                               slow_reactions,
                               fast_reactions,
                               epsilon,
                               k)
    
    recorder = cme_recorder.CmeRecorder(model)
    recorder.add_target('species',
                        ['expected value', 'standard deviation'],
                        model['species'],
                        model['species counts'])
    
    time_steps = numpy.linspace(0.0, 10.0, 101)
    for t in time_steps:
        solver.step(t)
        recorder.write(t, solver.y)
    
    """
    
    # uncomment for graphing code ...
    
    import pylab
    pylab.figure()
    for measurement in recorder.measurements('species'):
        pylab.plot(measurement.times,
                   measurement.expected_value,
                   label = measurement.name)
    pylab.legend()
    pylab.title('species count expected value')
    pylab.figure()
    for measurement in recorder.measurements('species'):
        pylab.plot(measurement.times,
                   measurement.standard_deviation,
                   label = measurement.name)
    pylab.legend()
    pylab.title('species count standard deviation')
    pylab.show()
    """
    
if __name__ == '__main__':
    import cProfile, pstats
    PROFILE_FILE = 'rank_approx.profile'
    cProfile.run('test()', PROFILE_FILE)
    stats = pstats.Stats(PROFILE_FILE)
    stats.sort_stats('cumulative').print_stats(30)
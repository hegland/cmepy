import numpy

import scipy.sparse
import scipy.linalg

import cmepy.models
import cmepy.new_core.cme_solver as cme_solver
import cmepy.new_core.ode_solver as ode_solver
import cmepy.new_core.recorder as cme_recorder

import cmepy.core.matrix_cme as matrix_cme
import perturb

def cached_svd_trunc(matrix, k, cache={}):
    if len(cache)==0:
        (k_max, k_max) = numpy.shape(matrix)
        u, s, vh = perturb.truncated_svd(matrix, k_max)
        cache['u'] = u
        cache['s'] = s
        cache['vh'] = vh
    
    u_bar = cache['u'][:, :k]
    s_bar = cache['s'][:k]
    vh_bar = cache['vh'][:k, :]
    
    return u_bar, s_bar, vh_bar
    

def approx_cme_solver(full_model,
                      slow_reactions,
                      fast_reactions,
                      epsilon,
                      k,
                      p_0 = None):
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
    
    print 'shape of fast_matrix : %s' % str(fast_matrix.shape)
    
    # m := limit of exp(fast_matrix *t) as t --> +ive infty
    # approximate limit by using a large T
    T_INFTY = 100000.0
    dense_fast_matrix = fast_matrix.todense()
    m_matrix = scipy.linalg.expm(dense_fast_matrix*T_INFTY)
    
    # compute truncated svd of k largest singular values
    u_bar, s_bar, v_bar = cached_svd_trunc(m_matrix, k)
    
    print '  dumping s_bar info:'
    s_bar_acc = numpy.add.accumulate(s_bar)
    tol = 1.0e-15
    tail_indices = s_bar_acc + tol >= s_bar_acc[-1] 
    last_head_index = (numpy.nonzero(tail_indices)[0])[0]
    print '\tlast_head_index %d' % last_head_index
    s_view_radius = 5
    s_view_lo = max(0, last_head_index - s_view_radius)
    s_view_hi = min(numpy.size(tail_indices), last_head_index + s_view_radius)
    for i in xrange(s_view_lo, s_view_hi):
        print '\ts[%d] = %f' % (i, s_bar[i])
    
     
    # use svd to compute aggregation and disaggregation matrices
    f = u_bar
    f_sparse = scipy.sparse.csr_matrix(f)
    e = numpy.dot(numpy.diag(s_bar), v_bar)
    e_sparse = scipy.sparse.csr_matrix(e)
    
    # define initial distribution as system with
    # maximum copies of first species and 0 copies of second species
    if p_0 is None:
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

def get_catalyser_model(initial_s_count, initial_e_count, epsilon):
    species = ('S', 'E', 'C', 'D')
    species_counts = (lambda *x : x[0],
                      lambda *x : x[1],
                      lambda *x : x[2],
                      lambda *x : initial_s_count - x[0] - x[2])
    propensities = (lambda *x : 1.0*x[0]*x[1],
                    lambda *x : 1.0*x[2],
                    lambda *x : epsilon*x[2])
    offset_vectors = ((-1, -1, 1),
                      (1, 1, -1),
                      (0, 1, -1))
    np = (initial_s_count+1, initial_e_count+1, initial_s_count+1)
    model = {'species' : species,
             'species counts' : species_counts,
             'propensities' : propensities,
             'offset_vectors' : offset_vectors,
             'np' : np}
    slow_reactions = (0, 1)
    fast_reactions = (2, )
    return model, slow_reactions, fast_reactions

def test(graph = False):
    """
    copy_count = 30
    k = copy_count+1
    slow_reactions = (0, )
    fast_reactions = (1, )
    epsilon = 0.01
    model = perturb.get_simple_model(copy_count)
    """
    
    epsilon = 0.1
    s_0 = 20
    e_0 = 3
    model, slow_reactions, fast_reactions = get_catalyser_model(s_0,
                                                                e_0,
                                                                epsilon)
    
    p_0 = numpy.zeros(model['np'])
    p_0[-1, -1, 0] = 1.0
    
    #k_range_coarse = [1764, 504, 475, 450, 425, 400, 375, 350, 325, 300, 275, 250]
    
    k_range_fine = numpy.linspace(350, 500, (500-350)/5 + 1)
    
    for k in k_range_fine:
        solver = approx_cme_solver(model,
                                   slow_reactions,
                                   fast_reactions,
                                   epsilon,
                                   k,
                                   p_0)
        
        recorder = cme_recorder.CmeRecorder(model)
        recorder.add_target('species',
                            ['expected value', 'standard deviation'],
                            model['species'],
                            model['species counts'])
        
        time_steps = numpy.linspace(0.0, 5.0, 101)
        for t in time_steps:
            solver.step(t)
            recorder.write(t, solver.y)
        
        if graph:
        #  graphing code ...
            title = 'using rank %d approx for M;\n' % k
            import pylab
            pylab.figure()
            for measurement in recorder.measurements('species'):
                pylab.plot(measurement.times,
                           measurement.expected_value,
                           label = measurement.name)
            pylab.legend()
            pylab.title(title+'species count expected value')
            pylab.savefig('ev_rank_%d_approx.png' % k)
            pylab.figure()
            for measurement in recorder.measurements('species'):
                pylab.plot(measurement.times,
                           measurement.standard_deviation,
                           label = measurement.name)
            pylab.legend()
            pylab.title(title+'species count standard deviation')
            pylab.savefig('sd_rank_%d_approx.png' % k)

def profile():
    import cProfile, pstats
    PROFILE_FILE = 'rank_approx.profile'
    cProfile.run('test()', PROFILE_FILE)
    stats = pstats.Stats(PROFILE_FILE)
    stats.sort_stats('cumulative').print_stats(30)

if __name__ == '__main__':
    #test(graph=True)
    profile()
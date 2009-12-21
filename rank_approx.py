import numpy

import scipy.sparse
import scipy.linalg

import cmepy.models
import cmepy.new_core.cme_solver as cme_solver
import cmepy.new_core.ode_solver as ode_solver
import cmepy.new_core.recorder as cme_recorder

import cmepy.core.matrix_cme as matrix_cme
import perturb
import block_diagonal

def create_change_of_basis_matrices(permutation):
    n = numpy.size(permutation)
    shape = (n, )*2
    forward_row = numpy.array(permutation)
    forward_col = numpy.arange(n)
    forward_val = numpy.ones((n, ))
    forward = scipy.sparse.coo_matrix((forward_val,
                                      (forward_row,
                                       forward_col)),
                                       shape).tocsr()
    inverse_row = numpy.arange(n)
    inverse_col = numpy.array(permutation)
    inverse_val = numpy.ones((n, ))
    inverse = scipy.sparse.coo_matrix((inverse_val,
                                      (inverse_row,
                                       inverse_col)),
                                       shape).tocsr()
    return (forward, inverse)

def create_ordered_basis(model, f):
    # change basis so that states are enumerated in order of increasing
    # value under the transform f
    indices = numpy.arange((numpy.product(model['np'],)))
    states = [numpy.ravel(i) for i in numpy.indices(model['np'])]
    f_states = f(*states)
    new_order = numpy.argsort(f_states)
    basis_permutation = indices[new_order]
    return basis_permutation    

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
    slow_matrix = slow_matrix.tocsr()

    fast_flux_data = cme_solver.create_flux_data(fast_model)
    fast_matrix = matrix_cme.gen_sparse_matrix(full_model['np'],
                                               fast_flux_data)
    fast_matrix *= epsilon
    fast_matrix = fast_matrix.tocsr()
    
    """
    # compute change of basis so that fast_matrix will be block diagonal
    # to do this we index states wrt their x[0], x[2], x[1] ordering
    
    np = full_model['np']
    state_ordering = lambda *x : (numpy.product(np[1:])*x[0]
                                  + numpy.product(np[1])*x[2]
                                  + x[1])
    
    block_diag_basis = create_ordered_basis(fast_model, state_ordering)
    beta, beta_inverse = create_change_of_basis_matrices(block_diag_basis)
    
    import pylab
    pylab.figure()
    pylab.spy(beta.todense())
    pylab.figure()
    pylab.spy(beta_inverse.todense())
    pylab.show()
    """
    
    #print 'beta : %s' % repr(beta)
    print 'fast_matrix : %s' % repr(fast_matrix)
    
    bd_fast = block_diagonal.from_sparse_matrix(fast_matrix)
    #bd_fast = block_diagonal.from_sparse_matrix(beta*fast_matrix)
    
    
    # m := limit of exp(fast_matrix *t) as t --> +ive infty
    # approximate limit by using a large T
    T_INFTY = 1000.0
    assert len(bd_fast.blocks)>0
    bd_m = block_diagonal.expm(bd_fast, T_INFTY)
    assert len(bd_m.blocks)>0
    for (start, size, block) in bd_m.blocks:
        if not numpy.logical_and.reduce(numpy.ravel(numpy.isfinite(block.data))):
            print 'sub block contains nonfinite elements :('
            print 'block start : %d' % start
            print 'block size : %d' % size
            print str(block.data)
    # compute rank k approx
    bd_svd_m = block_diagonal.block_svd(bd_m)
    assert len(bd_svd_m.blocks)>0
    (e_sparse, f_sparse) = block_diagonal.to_sparse_rank_k_approx(bd_svd_m, k)
    e_sparse = e_sparse.tocsr()
    f_sparse = f_sparse.tocsr()
    print 'e_sparse : %s'% repr(e_sparse)
    print 'f_sparse : %s'% repr(f_sparse)
    # define initial distribution as system with
    # maximum copies of first species and 0 copies of second species
    if p_0 is None:
        p_0 = numpy.zeros(full_model['np'])
        p_0[-1, 0] = 1.0
    
    a_hat = e_sparse*slow_matrix*f_sparse
    #a_hat = e_sparse*beta*slow_matrix*f_sparse
    print 'a_hat : %s' % repr(a_hat)
    
    pack, unpack = cme_solver.create_packing_functions(fast_model)
    
    def pack_aggregate(p):
        return e_sparse*pack(p)
        #return e_sparse*beta*pack(p)
    
    def deaggregate_unpack(y):
        return unpack(f_sparse*y)
        #return unpack(beta_inverse*f_sparse*y)
    
    def dy_dt(t, y):
        return a_hat*y
    
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
    
    k_range_coarse = [1764, 504, 475, 450, 425, 400, 375, 350, 325, 300, 275, 250]
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
            pylab.close()
            pylab.figure()
            for measurement in recorder.measurements('species'):
                pylab.plot(measurement.times,
                           measurement.standard_deviation,
                           label = measurement.name)
            pylab.legend()
            pylab.title(title+'species count standard deviation')
            pylab.savefig('sd_rank_%d_approx.png' % k)
            pylab.close()
            
def profile():
    import cProfile, pstats
    PROFILE_FILE = 'rank_approx.profile'
    cProfile.run('test()', PROFILE_FILE)
    stats = pstats.Stats(PROFILE_FILE)
    stats.sort_stats('cumulative').print_stats(30)

if __name__ == '__main__':
    test(graph=True)
    #profile()
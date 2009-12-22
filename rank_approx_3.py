import numpy
import scipy.sparse

import cmepy.core.matrix_cme as matrix_cme
import cmepy.new_core.ode_solver as ode_solver
import cmepy.new_core.cme_solver as cme_solver
import cmepy.new_core.recorder as cme_recorder

import block_diagonal

def get_catalyser_models(initial_s_count, initial_e_count, epsilon):
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
    
    slow_model = {'propensities' : propensities[:-1],
                  'offset_vectors' : offset_vectors[:-1],
                  'np' : np}
    
    fast_model = {'propensities' : (lambda *x : epsilon*x[1], ),
                  'offset_vectors' : ((1, -1), ),
                  'np' : (initial_e_count+1, initial_s_count+1)}
    return model, slow_model, fast_model

def create_indexing(np, f):
    # xxx todo strictly there should be something about norigin in here
    # use cmepy.util(s).indices_ext ....
    states = [numpy.ravel(i) for i in numpy.indices(np)]
    f_states = f(*states)
    permutation = numpy.argsort(f_states)
    inverse_permutation = numpy.argsort(permutation)
    return permutation, inverse_permutation

def create_cme_solver(initial_s_count,
                      initial_e_count,
                      epsilon,
                      approximation_rank_ratio):
    # ==========================================================================
    # define model & decomposition into slow and fast models
    # ==========================================================================
    print 'creating models'
    full_model, slow_model, fast_model = get_catalyser_models(initial_s_count,
                                                              initial_e_count,
                                                              epsilon)
    num_fast_states = numpy.product(fast_model['np'])
    duplication = slow_model['np'][0]
    num_slow_states = num_fast_states*duplication
    assert num_slow_states == numpy.product(slow_model['np'])
    
    # ==========================================================================
    # define state space indexing scheme as follows:
    #
    # for any two states x & y,
    #   1) if keys contain S values, order by this first; if tied, continue
    #   2) order by comparing x.E + x.C with y.E + y.C ; if tied, continue
    #   3) order by comparing x.E - x.C with y.E - y.C 
    #
    # the reason this scheme is chosen is to ensure that the matrix for the
    # fast_model ends up having a really nice block diagonal structure.
    # since the indexing scheme for the slow_model doesn't really matter
    # we just define it so transforming between the slow_model scheme
    # and the fast_model scheme is simple
    # ==========================================================================
    print 'creating indexing schemes'
    # xxx todo define these bijections in a nice way without
    # the magic numbers
    f_fast = lambda *x : (x[0]+x[1])*1000 + (x[0]-x[1])
    f_slow = lambda *x : x[0]*1000000 + (x[1]+x[2])*1000 + x[1]-x[2]
    
    fast_indexing, fast_indexing_inverse = create_indexing(fast_model['np'],
                                                           f_fast)
    slow_indexing, slow_indexing_inverse = create_indexing(slow_model['np'],
                                                           f_slow)
    # important! pass the **inverse** indexing permutation to the
    # sparse matrix construction routine ....
    print 'creating matrices'
    fast_flux_data = cme_solver.create_flux_data(fast_model)
    fast_matrix = matrix_cme.gen_sparse_matrix(fast_model['np'],
                                               fast_flux_data,
                                               fast_indexing_inverse)
    fast_matrix *= epsilon
    fast_matrix = fast_matrix.tocsr()
    
    slow_flux_data = cme_solver.create_flux_data(slow_model)
    slow_matrix = matrix_cme.gen_sparse_matrix(slow_model['np'],
                                               slow_flux_data,
                                               slow_indexing_inverse)
    slow_matrix = slow_matrix.tocsr()
    
    
    
    # ==========================================================================
    # compute matrix exponential of fast_matrix
    # -- use accelerated block diagonal routines ...
    # ==========================================================================
    print 'taking exponential of reduced fast matrix'
    
    # xxx todo this is a bit of a hack to compute the limit as t goes to infty
    # if too large a time is used the sparse block diagonal expm routine
    # will a bit upset by infinities popping up in sub blocks and
    # start returning nan entries ...
    
    # on the other hand if this is too small then its not a good
    # approximation of the limit.
    # t_infty = 100000 seems to work ok for initial_s 50, initial_e 6
    t_infty = 100000.0
    bd_fast = block_diagonal.from_sparse_matrix(fast_matrix)
    bd_m = block_diagonal.expm(bd_fast, t_infty)
    
    print 'duplicating exponential for full exponential fast matrix'
    shape = (duplication*num_fast_states, )*2
    bd_m_duplicates = block_diagonal.BlockDiagonalMatrix(shape)
    for i in xrange(duplication):
        start = i*num_fast_states
        size = num_fast_states
        bd_m_duplicates.add_block_diagonal(start,
                                           size,
                                           bd_m)
    print '***'
    print ''
    print '\tbd_m : %s' % repr(bd_m)
    print '\tbd_m_duplicates : %s' % repr(bd_m_duplicates)
    print ''
    print '***'
    
    # ==========================================================================
    # compute rank k approximation via svd
    # ==========================================================================
    k = int(duplication*num_fast_states*approximation_rank_ratio)
    print 'computing rank k approximation via svd, for k = %d' % k
    bd_svd_m = block_diagonal.block_svd(bd_m_duplicates)
    assert len(bd_svd_m.blocks)>0
    (e, f) = block_diagonal.to_sparse_rank_k_approx(bd_svd_m, k)
    e = e.tocsr()
    f = f.tocsr()
    
    print '***'
    print ''
    print '\te : %s'% repr(e)
    print '\tf : %s'% repr(f)
    print ''
    print '***'
    
    
    # ==========================================================================
    # define initial conditions, differential equations, and packing functions
    # ==========================================================================
    print 'defining initial conditions'
    # initialise system with maximal copies of S & E, zero copies of C & D
    p_0 = numpy.zeros(slow_model['np'])
    p_0[-1, -1, 0] = 1.0
    
    print 'defining dy_dt function'
    
    a_hat = e*slow_matrix*f
    def dy_dt(t, y):
        return a_hat*y
    
    print 'defining pack / unpack routines'
    # flatten p
    # permute ordering to match slow indexing scheme
    def pack(p):
        return e*(numpy.ravel(p)[slow_indexing])
    
    # inverse permute to match default indexing scheme
    # inflate
    def unpack(y):
        return numpy.reshape((f*y)[slow_indexing_inverse], slow_model['np'])
    
    if False:
        import pylab
        pylab.figure()
        pylab.spy(slow_matrix.todense())
        pylab.title('slow matrix')
        pylab.figure()
        pylab.spy(e.todense())
        pylab.title('e matrix')
        pylab.figure()
        pylab.spy(f.todense())
        pylab.title('f matrix')
        pylab.figure()
        pylab.spy(a_hat.todense())
        pylab.title('a_hat matrix')
        pylab.show()
        raise RuntimeError('debug : halting after plotting')
    
    # ==========================================================================
    # construct the cme solver
    # ==========================================================================
    print 'constructing solver'
    new_solver = ode_solver.Solver(dy_dt, p_0)
    # nb -- it is important to set the transform_dy_dt flag to False here
    # otherwise the solver will assume that the provided dy_dt is actually
    # of the form dp_dt and needs to be conjugated by the pack / unpack
    # functions ...
    new_solver.set_packing(pack,
                           unpack,
                           transform_dy_dt=False)
    return new_solver, full_model

def main():
    # ==========================================================================
    # set model definition parameters
    # ==========================================================================
    
    # xxx todo if these counts are set higher than 999 bad things will happen.
    # cf dirty hacks used to define state ordering for index schemes
    initial_s_count = 100
    initial_e_count = 10
    epsilon = 0.01
    approximation_rank_ratio = 1.0 / 5.0
    
    # ==========================================================================
    # compute the approximate cme solver for this model
    # ==========================================================================
    
    # nb the (full) model is used to interface with the cme recorder ...
    solver, model = create_cme_solver(initial_s_count,
                                      initial_e_count,
                                      epsilon,
                                      approximation_rank_ratio)
    
    # ==========================================================================
    # advance the solution and plot some results
    # ==========================================================================
    recorder = cme_recorder.CmeRecorder(model)
    recorder.add_target('species',
                        ['expected value', 'standard deviation'],
                        model['species'],
                        model['species counts'])
    
    time_steps = numpy.linspace(0.0, 1.0, 21)
    for t in time_steps:
        print 'stepping to t = %g' % t
        solver.step(t)
        recorder.write(t, solver.y)
    
    graph = True
    
    #  graphing code ...
    if graph:
        title = 's_count %d e_count %d;\n' % (initial_s_count, initial_e_count)
        import pylab
        pylab.figure()
        for measurement in recorder.measurements('species'):
            pylab.plot(measurement.times,
                       measurement.expected_value,
                       label = measurement.name)
        pylab.legend()
        pylab.title(title+'species count expected value')
        pylab.savefig('ev_%d_%d_approx_rank_%4f.png' % (initial_s_count,
                                                        initial_e_count,
                                                        approximation_rank_ratio))
        pylab.close()
        pylab.figure()
        for measurement in recorder.measurements('species'):
            pylab.plot(measurement.times,
                       measurement.standard_deviation,
                       label = measurement.name)
        pylab.legend()
        pylab.title(title+'species count standard deviation')
        pylab.savefig('sd_%d_%d_approx_rank_%4f.png' % (initial_s_count,
                                                        initial_e_count,
                                                        approximation_rank_ratio))
        pylab.close()
    
def profile():
    import cProfile, pstats
    PROFILE_FILE = 'rank_approx_3.profile'
    cProfile.run('main()', PROFILE_FILE)
    stats = pstats.Stats(PROFILE_FILE)
    stats.sort_stats('cumulative').print_stats(30)


if __name__ == '__main__':
    profile()

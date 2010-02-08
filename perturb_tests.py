import unittest

import numpy
import numpy.testing.utils
import scipy.sparse
import scipy.linalg

from cmepy.solver import CmeSolver, CmeSolverMatrix
from cmepy.recorder import CmeRecorder

import perturb

def verify_does_not_contain_nan(x):
    x_flat = numpy.ravel(x)
    assert not numpy.logical_or.reduce(numpy.isnan(x_flat))

def verify_model_structure(model):
    assert model is not None
    assert type(model) is dict
    assert 'propensities' in model
    assert 'np' in model
    assert 'offset_vectors' in model
    assert 'doc' in model
    assert 'species counts' in model
    assert 'species' in model
    assert type(model['doc']) is str
    assert len(model['propensities']) == len(model['offset_vectors'])

class TestPerturb(unittest.TestCase):
    
    def testGetSimpleModel(self):
        model = perturb.get_simple_model()
        verify_model_structure(model)
        assert model['offset_vectors'][0] == (-1, 1)
        assert model['offset_vectors'][1] == (0, -1)
    
    def testSolveSimpleModelFull(self):
        model = perturb.get_simple_model()
        solver = CmeSolver(model)
        N_STEPS = 10
        time_steps = numpy.linspace(0.0, float(N_STEPS), N_STEPS+1)
        for t in time_steps:
            solver.step(t)
    
    def testComputeSystemMatrix(self):
        model = perturb.get_simple_model()
        matrix = perturb.compute_system_matrix(model)
        assert matrix is not None
        assert type(matrix) is scipy.sparse.coo_matrix
    
    def testDecomposeModel(self):
        full_model = perturb.get_simple_model()
        slow_reactions = (0, )
        fast_reactions = (1, )
        reaction_subsets = (slow_reactions, fast_reactions)
        (slow_model, fast_model) = perturb.decompose_model(full_model,
                                                           reaction_subsets)
        verify_model_structure(slow_model)
        verify_model_structure(fast_model)
        
        assert (slow_model['propensities'][0]
                == full_model['propensities'][0])
        assert (slow_model['offset_vectors'][0]
                == full_model['offset_vectors'][0])
        assert (fast_model['propensities'][0]
                == full_model['propensities'][1])
        assert (fast_model['offset_vectors'][0]
                == full_model['offset_vectors'][1])

    def testDecomposeModelThenComputeMatrices(self):
        full_model = perturb.get_simple_model()
        slow_reactions = (0, )
        fast_reactions = (1, )
        reaction_subsets = (slow_reactions, fast_reactions)
        (slow_model, fast_model) = perturb.decompose_model(full_model,
                                                           reaction_subsets)
        
        epsilon = 0.01
        slow_matrix = perturb.compute_system_matrix(slow_model)
        fast_matrix = perturb.compute_system_matrix(fast_model)*epsilon
    
    def testComputeAsymptoticMatrix(self):
        full_model = perturb.get_simple_model()
        slow_reactions = (0, )
        fast_reactions = (1, )
        reaction_subsets = (slow_reactions, fast_reactions)
        (slow_model, fast_model) = perturb.decompose_model(full_model,
                                                           reaction_subsets)
        
        epsilon = 0.01
        slow_matrix = perturb.compute_system_matrix(slow_model)
        fast_matrix = perturb.compute_system_matrix(fast_model)*epsilon
        
        # m := limit of exp(fast_matrix *t) as t --> +ive infty
        # approximate limit by using a large T
        
        T_INFTY = 100000.0
        dense_fast_matrix = fast_matrix.todense()
        m_matrix = scipy.linalg.expm(dense_fast_matrix*T_INFTY)
        verify_does_not_contain_nan(m_matrix)
    
    def testTruncatedSVDAggregation(self):
        m_matrix = numpy.array([[1, 1, 0, 0],
                                [0, 0, 0, 0],
                                [0, 0, 1, 1],
                                [0, 0, 0, 0]])
        
        # compute truncated svd of 2 largest singular values
        u_bar, s_bar, v_bar = perturb.truncated_svd(m_matrix,
                                                    k = 2)
        
        # use svd to compute aggregation and disaggregation matrices
        f = u_bar
        e = numpy.dot(numpy.diag(s_bar), v_bar)
        
        # recover m_matrix
        m_bar_matrix = numpy.dot(f, e)
        
        # verify that matrix is fully recovered
        numpy.testing.utils.assert_array_almost_equal(m_matrix, m_bar_matrix)
    
    def testComputeZerothOrderApproxSystem(self):
        full_model = perturb.get_simple_model()
        slow_reactions = (0, )
        fast_reactions = (1, )
        reaction_subsets = (slow_reactions, fast_reactions)
        (slow_model, fast_model) = perturb.decompose_model(full_model,
                                                           reaction_subsets)
        
        epsilon = 0.01
        slow_matrix = perturb.compute_system_matrix(slow_model)
        fast_matrix = perturb.compute_system_matrix(fast_model)*epsilon
        
        # m := limit of exp(fast_matrix *t) as t --> +ive infty
        # approximate limit by using a large T
        
        T_INFTY = 100000.0
        dense_fast_matrix = fast_matrix.todense()
        m_matrix = scipy.linalg.expm(dense_fast_matrix*T_INFTY)
        
        # compute truncated svd of 2 largest singular values
        u_bar, s_bar, v_bar = perturb.truncated_svd(m_matrix,
                                                    k = 2)
        
        # use svd to compute aggregation and disaggregation matrices
        f = u_bar
        e = numpy.dot(numpy.diag(s_bar), v_bar)
        
        # define initial distribution as system with
        # 1 copy of first species and 0 copies of second species
        p_0_fat = numpy.zeros(full_model['np'])
        p_0_fat[1, 0] = 1.0
        
        p_0 = numpy.ravel(p_0_fat)
        
        a_hat = numpy.dot(e, slow_matrix*f)
        a_hat_sparse = scipy.sparse.csr_matrix(a_hat)
        p_0_hat = numpy.dot(e, p_0)
        
        # create system of ODEs for aggregated system...
        np_hat = (numpy.size(p_0_hat),)
        
        # XXX TODO cme solver matrix shouldn't actually need ANY model since
        # the matrix and np are being explicitly supplied ...
        solver = CmeSolverMatrix(full_model, a_hat_sparse)
        solver.set_solver_params(np = np_hat)
        solver.set_initial_values(p0 = p_0_hat)
        
        
        # XXX TODO GENERALISE THIS AND MAKE IT NICER
        # boilerplate interface code:
        # trick recorder into thinking that aggregated solver
        # is actually the full solver by passing it this thing.
        # recorder only needs solver.model, solver.t, solver.get_p(),
        
        class WrappedSolver():
            def __init__(self, model, solver, unwrap):
                self.model = model
                self.solver = solver
                self.unwrap = unwrap
            def _gett(self):
                return solver.t
            def _getp(self):
                return self.unwrap(solver.get_p())
            t = property(_gett, None)
            p = property(_getp, None)
            def get_p(self):
                return self.p
        
        # define a map translating from the aggregated approx solution p_bar
        # back to the approx distribution p of the full model
        def unwrap(p_bar):
            return numpy.reshape(numpy.dot(f, p_bar), full_model['np'])
        
        wrapped_solver = WrappedSolver(full_model,
                                       solver,
                                       unwrap)
        
        recorder = CmeRecorder(wrapped_solver)
        recorder.add_target(output = ['expectation', 'std_dev'],
                        species = full_model['species'])
        
        n_time_steps = 100
        t_final = 10.0
        time_steps = numpy.linspace(0.0,
                                    t_final,
                                    n_time_steps + 1)
        for t in time_steps:
            solver.step(t)
            recorder.take_measurements()
        
        import pylab
        pylab.figure()
        for s_info in recorder.measurements('species'):
            pylab.plot(s_info.times, s_info.expectation, label = s_info.name)
        pylab.legend()
        pylab.title('species count expectation values')
        
        # plot std dev, over time, of all species counts
        pylab.figure()
        for s_info in recorder.measurements('species'):
            pylab.plot(s_info.times, s_info.std_dev, label = s_info.name)
        pylab.legend()
        pylab.title('species count std dev')
        pylab.show()
        
if __name__ == '__main__':
    unittest.main()

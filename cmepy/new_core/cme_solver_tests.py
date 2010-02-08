import unittest
from test import test_support

import itertools
import math

import numpy
from numpy.testing.utils import assert_almost_equal

import cmepy.solver
import cmepy.models

import cmepy.new_core.ode_solver as ode_solver
import cmepy.new_core.cme_solver as cme_solver


ALL_MODELS = [cmepy.models.A2B2A,
              cmepy.models.A2B2C,
              cmepy.models.BURR08,
              cmepy.models.DSMTS_001_01,
              cmepy.models.GOU07_A,
              cmepy.models.GOU07_B,
              cmepy.models.GOU07_C,
              cmepy.models.MM_SIMPLE,
              cmepy.models.MUNK08,
              cmepy.models.MUNK08_A]

FAST_MODELS = [cmepy.models.A2B2A,
               cmepy.models.A2B2C,
               cmepy.models.GOU07_A,
               cmepy.models.GOU07_B,
               cmepy.models.MM_SIMPLE,
               cmepy.models.MUNK08_A]

# some of the models are pretty slow and annoying to run, if they're
# solved during every unit test. So we only solve the faster ones

TEST_MODELS = FAST_MODELS

class CmeSolverTests(unittest.TestCase):
    
    def setUp(self):
        self.models = TEST_MODELS
        
    def testNewCmeSolverFunctions1(self):
        for model in self.models:
            inflated_shape = model['np']
            flattened_shape = (numpy.product(inflated_shape), )
            
            default_p_0 = cme_solver.create_default_p_0(model)
            flux_data = cme_solver.create_flux_data(model)
            dp_dt = cme_solver.create_diff_eqs(flux_data)
            
            assert default_p_0 is not None
            assert numpy.shape(default_p_0) == inflated_shape
            assert_almost_equal(numpy.ravel(default_p_0)[0], 1.0)
            assert_almost_equal(numpy.ravel(default_p_0)[1:],
                                numpy.zeros(flattened_shape)[1:])
            
            assert flux_data is not None
            assert len(flux_data) == len(model['propensities'])
            
            assert dp_dt is not None
            evaluated_dp_dt = dp_dt(0.0, default_p_0)
            assert evaluated_dp_dt is not None
            assert numpy.shape(evaluated_dp_dt) == numpy.shape(default_p_0)
    
    def testNewCmeSolverFunctions2(self):
        for model in self.models:
            inflated_shape = model['np']
            flattened_shape = (numpy.product(inflated_shape), )
            
            pack, unpack = cme_solver.create_packing_functions(model)
            default_p_0 = cme_solver.create_default_p_0(model)
            
            assert numpy.shape(pack(default_p_0)) == flattened_shape
            assert numpy.shape(unpack(pack(default_p_0))) == inflated_shape
            assert_almost_equal(unpack(pack(default_p_0)), default_p_0)
    
    def testNewCmeSolverFunctions3(self):
        for model in self.models:
            # compare manual construction of solver with
            # automatic construction ...
            
            # solver_a : `manually' constructed
            
            pack, unpack = cme_solver.create_packing_functions(model)
            default_p_0 = cme_solver.create_default_p_0(model)
            flux_data = cme_solver.create_flux_data(model)
            dp_dt = cme_solver.create_diff_eqs(flux_data)
            
            solver_a = ode_solver.Solver(dp_dt, default_p_0)
            solver_a.set_packing(pack, unpack)
            
            # solver_b : `automatically' constructed
            
            solver_b = cme_solver.create_cme_solver(model)
            
            assert_almost_equal(solver_a.t, solver_b.t)
            assert_almost_equal(solver_a.y, solver_b.y)
            assert_almost_equal(solver_a.dy_dt(0.0, default_p_0),
                                solver_b.dy_dt(0.0, default_p_0))
            
    def testCompareSolutionsWithExistingSolver(self):
        for model in self.models:
            # construct dependencies
            pack, unpack = cme_solver.create_packing_functions(model)
            default_p_0 = cme_solver.create_default_p_0(model)
            flux_data = cme_solver.create_flux_data(model)
            # XXX TODO INJECT TIME DEPENDENCE HERE [?]
            dp_dt = cme_solver.create_diff_eqs(flux_data)
            
            # initialise cme_solver as ode solver instance
            new_solver = ode_solver.Solver(dp_dt, default_p_0)
            new_solver.set_packing(pack, unpack)
            
            num_time_steps = 11
            time_steps = numpy.linspace(0.0, 0.5, num_time_steps)
            
            test_results = []
            
            for t in time_steps:
                new_solver.step(t)
                test_results.append(new_solver.y)
            
            # compute results for same model via old solver
            old_solver = cmepy.solver.CmeSolver(model)
            old_results = []
            for t in time_steps:
                old_solver.step(t)
                old_results.append(old_solver.get_p())
            
            # compare results and verify they agree
            for test_result, old_result in itertools.izip(test_results,
                                                          old_results):
                assert_almost_equal(test_result, old_result)
    
    def testTimeDependentFluxEvaluatorA(self):
        for model in self.models:
            pack, unpack = cme_solver.create_packing_functions(model)
            default_p_0 = cme_solver.create_default_p_0(model)
            flux_data = cme_solver.create_flux_data(model)
            
            # define time dependencies
            def create_phi(i):
                return lambda t : 0.5+0.5*math.sin(t*(i+1.0))
            
            num_reactions = len(model['propensities'])
            time_dependencies = [create_phi(i) for i in xrange(num_reactions)]
            flux_evaluator = cme_solver.create_time_dependent_flux_evaluator(time_dependencies)
            
            dp_dt = cme_solver.create_diff_eqs(flux_data, flux_evaluator)
            new_solver = ode_solver.Solver(dp_dt, default_p_0)
            new_solver.set_packing(pack, unpack)
            
            num_time_steps = 11
            time_steps = numpy.linspace(0.0, 0.5, num_time_steps)
            
            new_results = []
            for t in time_steps:
                new_solver.step(t)
                new_results.append(new_solver.y)
            
            model_old = dict(model)
            model_old['time_dependence'] = time_dependencies
            old_solver = cmepy.solver.CmeSolver(model_old)
            
            old_results = []
            for t in time_steps:
                old_solver.step(t)
                old_results.append(old_solver.get_p())
            
            # compare results and verify they agree
            for new_result, old_result in itertools.izip(new_results,
                                                          old_results):
                assert_almost_equal(new_result, old_result)

def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(CmeSolverTests)
    return suite

def main():
    test_support.run_unittest(CmeSolverTests)

if __name__ == '__main__':
    main()
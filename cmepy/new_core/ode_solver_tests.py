import unittest
from test import test_support

import numpy
import numpy.testing.utils

import cmepy.new_core.ode_solver as ode_solver

class OdeSolverTests(unittest.TestCase):
    
    def testInitialise(self):
        def dy_dt(t, y):
            return 0.0
        
        solver = ode_solver.Solver(dy_dt, y_0 = 1.0)
        assert solver is not None
        assert solver.dy_dt == dy_dt

    def testInitialiseSetInitialValuesPreventExternalModification(self):
        def dy_dt(t, y):
            return 0.0
        
        t_0 = 1.0
        y_0 = -3.0
        solver = ode_solver.Solver(dy_dt, y_0, t_0)
        assert solver is not None
        assert solver.dy_dt == dy_dt
        assert solver.t == t_0
        assert solver.y == y_0
        
        def externally_modify_solver_dy_dt():
            solver.dy_dt = None
        
        def externally_modify_solver_t():
            solver.t = 12345
        
        def externally_modify_solver_y():
            solver.y = 'banana'
        
        self.assertRaises(AttributeError, externally_modify_solver_dy_dt)    
        self.assertRaises(AttributeError, externally_modify_solver_t)
        self.assertRaises(AttributeError, externally_modify_solver_y)
    
    def testSolveScalarProblem1(self):
        def dy_dt(t, y):
            return 0.0
        
        y_0 = 11.0
        
        solver = ode_solver.Solver(dy_dt, y_0)
        
        time_steps = numpy.linspace(0.0, 10.0, 11)
        for t in time_steps:
            solver.step(t)
            assert solver.y == y_0
            assert solver.t == t
    
    def testSolveScalarProblem2(self):
        def dy_dt(t, y):
            return 1.0
        
        t_0 = 0.0
        y_0 = 15.0
        
        solver = ode_solver.Solver(dy_dt, y_0)
        
        time_steps = numpy.linspace(0.0, 10.0, 11)
        for t in time_steps:
            solver.step(t)
            numpy.testing.utils.assert_almost_equal(solver.y, y_0+t)
            assert solver.t == t
    
    def testPreventSolverSolutionCorruption(self):
        def dy_dt(t, y):
            return numpy.ones(numpy.shape(y))
        
        y_0 = numpy.zeros((4,))
        
        solver = ode_solver.Solver(dy_dt, y_0)
        
        time_steps = numpy.linspace(0.0, 10.0, 11)
        for t in time_steps:
            solver.step(t)
            assert solver.t == t
            y_current_1 = solver.y
            numpy.testing.utils.assert_almost_equal(y_current_1, y_0+t)
            # this operation should not corrupt the solver state ...
            y_current_1[:] = (-1.0e7, 1.0e4, 1.0e-4, 1.0e5)
            y_current_2 = solver.y
            numpy.testing.utils.assert_almost_equal(y_current_2, y_0+t)
            assert y_current_1 is not y_current_2
    
    def testPreventSolverBackwardsSteps(self):
        def dy_dt(t, y):
            return numpy.ones(numpy.shape(y))
        
        t_0 = 3.11
        y_0 = numpy.zeros((4,))
        
        solver = ode_solver.Solver(dy_dt, y_0, t_0)
        
        time_steps = numpy.linspace(0.0, 10.0, 11)
        for t in time_steps:
            if t < t_0:
                self.assertRaises(ValueError, solver.step, t)
            else:
                solver.step(t)
                assert solver.t == t
                numpy.testing.utils.assert_almost_equal(solver.y, y_0+(t-t_0))
    
    def testSolveWrappedProblem(self):
        
        shape = (10, 10)
        
        def dz_dt(t, z):
            return numpy.ones(shape)
        
        z_0 = -3.0*numpy.ones(shape)
        
        solver = ode_solver.Solver(dz_dt, z_0)
        
        def pack(z):
            return numpy.ravel(z)
        def unpack(y):
            return numpy.reshape(y, shape)
        
        solver.set_packing(pack, unpack)
        
        time_steps = numpy.linspace(0.0, 10.0, 11)
        for t in time_steps:
            solver.step(t)
            assert solver.t == t
            numpy.testing.utils.assert_almost_equal(solver.y, z_0 + t)
        
def suite():
    suite = unittest.TestLoader().loadTestsFromTestCase(OdeSolverTests)
    return suite

def main():
    test_support.run_unittest(OdeSolverTests)

if __name__ == '__main__':
    main()
"""
a simple test of the CmeSolverMatrix functionality vs CmeSolver
"""

import numpy

from cmepy.solver import CmeSolver, CmeSolverMatrix
import cmepy.models

def compare_solvers():
    """
    compare outputs of different solvers when applied to the same model
    """
    model = cmepy.models.A2B2A
    solvers = {'CmeSolver' : CmeSolver(model),
               'CmeSolverMatrix' : CmeSolverMatrix(model)}
    
    n_time_steps = 100
    time_steps = numpy.linspace(0.0, float(n_time_steps), n_time_steps+1)
    results = []
    for solver_name in solvers:
        print 'using solver: %s' % solver_name
        solver = solvers[solver_name]
        for t in time_steps:
            print 'stepping to t = %f' % t
            solver.step(t)
        print 'finished'
        results.append(solver.get_p())
    print 'comparing results...'
    numpy.testing.utils.assert_array_almost_equal(results[0],
                                                  results[1])
    print 'ok'
    
if __name__ == '__main__':    
    compare_solvers()

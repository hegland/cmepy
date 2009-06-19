"""
Example of CmeSolver set_initial_values usage.
"""


def main():
    import numpy
    from numpy.testing.utils import assert_array_almost_equal
    from cmepy.solver import CmeSolver
    from cmepy.models import A2B2C as model
    
    num_steps = 10
    time_steps = numpy.linspace(0.0, 10.0, num_steps)
    
    solver_a = CmeSolver(model)
    for t in time_steps:
        solver_a.step(t)
    
    p_a_final = solver_a.get_p()
    
    solver_b = CmeSolver(model)
    for t in time_steps:
        prev_p = solver_b.get_p()
        prev_t = solver_b.get_t()
        # rebuild solver, using output of previous solver
        # to specify the initial conditions
        solver_b = CmeSolver(model)
        solver_b.set_initial_values(t0 = prev_t,
                                    p0 = prev_p)
        solver_b.step(t)
    
    p_b_final = solver_b.get_p()

    
    print('verifying solutions are (almost) equal ...')
    assert_array_almost_equal(p_a_final, p_b_final)
    print('OK.')

if __name__ == '__main__':
    main()

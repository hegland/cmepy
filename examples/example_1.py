"""
Simple example illustrating use of CmeSolver.
"""

def main():
    """
    In this example we construct a model of a system of two Poisson processes,
    advance it via the CmeSolver, then plot output.
    """
    
    import numpy
    import pylab
    import cmepy.solver
    
    
    # initialise a cme solver instance for this model
    
    # -- specify the model, with reaction propensities
    
    #    The format each propensity function is
    #
    #        x1, x2, ..., xn -> prop ;
    #
    #    where x1, ..., xn are the counts of reactions 1, ..., n respectively,
    #    and prop is the propensity of the reaction, given these counts.
    #
    #    The propensity functions are packed into a tuple, so the first
    #    element is the propensity function for the first reaction,
    #    the second element is the propensity function corresponding to the
    #    second reaction, etc.
    
    model = {'propensities' : (lambda x1, x2 : 1.0,
                               lambda x1, x2 : 1.0,
                              )
            }
    
    # -- instantiate the solver
    solver = cmepy.solver.CmeSolver(model)
    
    # -- specify size of the (reaction count) state space
    np = (128, 128)
    solver.set_solver_params(np=np)
    
    # -- pass the initial time
    # -- the initial distribution can also be set - here we use the default
    #    initial distribution, which is both reaction counts zero with Pr = 1.
    solver.set_initial_values(t0=0.0)
    
    # initialise the plot
    pylab.figure()
    pylab.xlabel('reaction 1 count')
    pylab.ylabel('reaction 2 count')
    pylab.title('solutions of the CME for 2d Poisson process')
    
    time_steps = [0, 20, 40, 60]
    for time in time_steps:
        # advance the solution of the cme
        solver.step(time)
        
        # get the current probability distribution p
        p_current = solver.get_p()
        
        # draw a contour plot of the distribution p after each time step.
        # -- note we explicitly clamp all probabilities to be positive, in
        #    order to avoid the artifacts introduced when contour plotting
        #    very small negatives.
        pylab.contour(numpy.maximum(p_current.T, 0.0))
    
    # display the contour plots
    pylab.show()

if __name__ == '__main__':
    main()

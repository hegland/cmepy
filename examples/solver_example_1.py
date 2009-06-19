"""
Simple example of CmeSolver usage.
"""

def main():
    import numpy
    import cmepy.solver
    
    # specify a simple model, consisting of two pure birth reactions with constant
    # propensities, i.e.
    # * -> x1
    # * -> x2
    model = {'propensities' : (lambda x1, x2 : 1.0,  # propensity of first reaction
                               lambda x1, x2 : 2.0,),# propensity of second reaction
             'np' : (128, 128) } # dimensions of reaction-count state space
    
    # create a cme solver instance for this model
    solver = cmepy.solver.CmeSolver(model)
    
    distributions = []
    
    time_steps = numpy.linspace(0, 60, 10)
    for t in time_steps:
        # advance the solution to time t
        solver.step(t)
        # obtain and store the current distribution p
        distributions.append(solver.get_p())
        
    
    # plot distributions
    import pylab
    
    pylab.figure()
    for p in distributions:
        # negative probabilities of very small magnitude sometimes occur,
        # these introduce ugly artifacts when contour plotted
        p_threshold = numpy.maximum(p, 0.0)
        # n.b. we transpose the 2d array before contour plotting it
        # so that reaction 1 appears along the x axis.
        pylab.contour(p_threshold.T)
    
    pylab.xlabel('reaction 1 count')
    pylab.ylabel('reaction 2 count')
    pylab.title('solutions of the CME')
    pylab.show()

if __name__ == '__main__':
    main()

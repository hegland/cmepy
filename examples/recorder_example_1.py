"""
Example demonstrating the use of a CmeRecorder to compute common measurements.

For additional information, see the documentation for
    cmepy.recorder.cme_recorder
"""

def main():
    import numpy
    from cmepy.solver import CmeSolver
    from cmepy.recorder import CmeRecorder
    import cmepy.models
    
    model = cmepy.models.A2B2C
    solver = CmeSolver(model)
    
    recorder = CmeRecorder(solver)
    
    # we want to record expectation & marginal for all reactions
    recorder.add_target(output = ['expectation', 'marginal'],
                        reactions = model['reactions'])
    
    # we also want to record expectation & std dev for all species
    recorder.add_target(output = ['expectation', 'std_dev'],
                        species = model['species'])
    
    time_steps = numpy.linspace(0, 2.5, 51)
    for t in time_steps:
        # step the solver solution forward
        solver.step(t)
        # compute and store measurements, from solver, using current solution
        recorder.take_measurements()
    
    import pylab
    
    # plot expectation, over time, of all reaction counts
    pylab.figure()
    for r_info in recorder.measurements('reactions'):
        pylab.plot(r_info.times, r_info.expectation, label = r_info.name)
    pylab.legend()
    pylab.title('reaction count expectation values')
    
    # plot expectation, over time, of all species counts
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
    
    # get the info recorded for the first reaction count
    r_info = recorder.reactions[model['reactions'][0]]
    # look at the marginal distribution for the most recent recording
    marginal = r_info.marginal[-1]
    pylab.figure()
    pylab.plot(marginal)
    pylab.title('marginal for reaction '+str(r_info.name)+' at final time step')
    
    pylab.show()

if __name__ == '__main__':
    main()
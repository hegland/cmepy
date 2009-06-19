def main():
    import numpy
    from cmepy.solver import CmeSolver
    from cmepy.recorder import CmeRecorder
    
    # load a model
    import cmepy.models
    model = cmepy.models.MM_SIMPLE
    
    # create cme solver
    solver = CmeSolver(model)
    
    # create cme recorder, specifying that we want to record
    # the marginal of each reaction named in the model
    recorder = CmeRecorder(solver)
    recorder.add_target(output = ['marginal'],
                        reactions = model['reactions'])
    recorder.add_target(output = ['expectation'],
                        species = model['species'])
    
    # advance the solution, recording at each time
    time_steps = numpy.linspace(0, 2.5, 51)
    for t in time_steps:
        solver.step(t)
        recorder.take_measurements()
    
    
    import pylab
    # plot the recorded expectations, for each species count
    pylab.figure()
    for s_info in recorder.measurements('species'):
        pylab.plot(s_info.times, s_info.expectation, label = s_info.name)
    pylab.legend()
    pylab.title('Expectation values of species counts')
    
    # plot the recorded marginal distributions, for the final time step,
    # for each reaction count
    pylab.figure()
    for r_info in recorder.measurements('reactions'):
        final_marginal = r_info.marginal[-1]
        pylab.plot(final_marginal, label = r_info.name)
    
    pylab.legend()
    pylab.title('Final marginal distributions of reaction counts')
    
    pylab.show()

if __name__ == '__main__':
    main()

"""
A simple catalytic reaction example.

Demonstrates how the CmeSolver can be used to solve the chemical master
equation and display expected species counts, for the following catalytic
reaction:

    S_1 --k_1--> S_2
    S_2 --k_2--> S_3
    S_2 + S_4 --k_3--> S_2 + S_5

Here, the species S_2 acts as a catalyst for the third
reaction. The rate constants are

    k_1 = 1, k_2 = 1000, k_3 = 100;

while the initial species counts are

    15, 20, 0, 0, 0

for the species S_1, ..., S_5 respectively.

This reaction is taken from the paper 

    Mastny, E.A., Haseltine, E.L. and Rawlings, J.B.
    Two classes of quasi-steady-state model reductions for stochastic kinetics,
    The Journal of Chemical Physics, 2007, Vol 127 .

"""

def main():
    """
    Solves the CME for a simple catalytic reaction and displays time-stamped
    expected species count data.
    """
    
    import numpy
    from cmepy.solver import CmeSolver
    from cmepy.recorder import CmeRecorder
        
    # define the model
    model = {'doc' : 'simple catalytic reaction',
             'propensities' : (lambda *x : 1.0*(15-x[0]),
                               lambda *x : 1000.0*(x[0]-x[1]),
                               lambda *x : 100.0*(20-x[2])*(x[0]-x[1])),
             'reactions' : ('S1->S2', 'S2->S3', 'S2+S4->S2+S5' ),
             'species counts' : (lambda *x : 15-x[0],
                                 lambda *x : x[0]-x[1],
                                 lambda *x : x[1],
                                 lambda *x : 20-x[2],
                                 lambda *x : x[2]),
             'species' : ('S1', 'S2', 'S3', 'S4', 'S5'),
             'np' : (16, 16, 21), }
    
    # create cme solver
    solver = CmeSolver(model)
    
    # create cme recorder, specifying that we want to record
    # the expected count of each species named in the model
    recorder = CmeRecorder(solver)
    recorder.add_target(output = ['expectation'],
                        species = model['species'])
    
    # order species measurements by name
    species_measurements = list(recorder.measurements('species'))
    species_measurements.sort(key = lambda x : x.name)
    
    # print column headers
    print 't',
    for species_info in species_measurements:
        print ',\t'+species_info.name,
    print ''
    
    time_steps = numpy.linspace(0, 0.5, 101)
    
    for t in time_steps:
        solver.step(t)
        recorder.take_measurements()
        
        # print row containing time, followed by expected species counts
        print '%.3f' % t,
        for species_info in species_measurements:
            # display the most recent expected species count
            print ', % .5f' % (species_info.expectation[-1]),
        print ''

if __name__ == '__main__':
    main()

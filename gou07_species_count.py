import numpy

import cmepy.new_core.cme_solver as cme_solver
import cmepy.new_core.recorder as cme_recorder

import pylab

def create_gou07_c_model(max_p_trunc, max_q_trunc, s_0=10, d_0=2, vol=1.0):
    """
    Creates a species-count based model for the Gou07 c system of reactions:
    
        S->P
        D+P->D+P+P
        P+P->P+Q
        P+Q->Q+Q
        P->*
        Q->*
    
    The copy counts of the species S and D are assumed to be constant,
    with s_0 (default 10) copies of S and d_0 (default 2) copies of 2.
    """
    
    def nonneg(x):
        return numpy.maximum(x, 0.0)
    
    species = ('P', 'Q', 'S', 'D', )
    p = lambda *x : x[0]
    q = lambda *x : x[1]
    s = lambda *x : s_0
    d = lambda *x : d_0
    species_counts = (p, q, s, d)

    reactions = ('S->P',
                 'D+P->D+2P',
                 'P+P->P+Q',
                 'P+Q->2Q',
                 'P->*',
                 'Q->*',)
    
    props = (lambda *x : 0.002*(vol**-1)*s(*x),
             lambda *x : 0.001*(vol**-2)*d(*x)*p(*x),
             lambda *x : 0.005*(vol**-2)*p(*x)*nonneg(p(*x)-1)/2.0,
             lambda *x : 0.004*(vol**-2)*p(*x)*q(*x),
             lambda *x : 0.002*(vol**-1)*p(*x),
             lambda *x : 0.050*(vol**-1)*q(*x),)
    
    offsets = ((1, 0),
               (1, 0),
               (-1, 1),
               (-1, 1),
               (-1, 0),
               (0, -1),)
    
    np = (max_p_trunc+1, max_q_trunc+1)
    model = {'doc' : 'quadratic autocatalator (species count)',
             'reactions' : reactions,
             'propensities' : props,
             'offset_vectors' : offsets,
             'species' : species,
             'species counts' : species_counts,
             'np' : np}
    return model

def create_gou07_d_model(max_p_trunc, max_q_trunc, max_s_trunc, d_0=2, vol=1.0):
    """
    Creates a species-count based model for the Gou07 c system of reactions:
    
        S->P
        D+P->D+P+P
        P+P->P+Q
        P+Q->Q+Q
        P->*
        Q->*
    
    The copy count of the species D is assumed to be constant,
    with d_0 (default 2) copies of 2.
    """
    
    def nonneg(x):
        return numpy.maximum(x, 0.0)
    
    species = ('P', 'Q', 'S', 'D', )
    p = lambda *x : x[0]
    q = lambda *x : x[1]
    s = lambda *x : x[2]
    d = lambda *x : d_0
    species_counts = (p, q, s, d)

    reactions = ('S->P',
                 'D+P->D+2P',
                 'P+P->P+Q',
                 'P+Q->2Q',
                 'P->*',
                 'Q->*',)
    
    props = (lambda *x : 0.002*(vol**-1)*s(*x),
             lambda *x : 0.001*(vol**-2)*d(*x)*p(*x),
             lambda *x : 0.005*(vol**-2)*p(*x)*nonneg(p(*x)-1)/2.0,
             lambda *x : 0.004*(vol**-2)*p(*x)*q(*x),
             lambda *x : 0.002*(vol**-1)*p(*x),
             lambda *x : 0.050*(vol**-1)*q(*x),)
    
    offsets = ((1, 0, -1),
               (1, 0, 0),
               (-1, 1, 0),
               (-1, 1, 0),
               (-1, 0, 0),
               (0, -1, 0),)
    
    np = (max_p_trunc+1, max_q_trunc+1, max_s_trunc+1)
    model = {'doc' : 'quadratic autocatalator (species count)',
             'reactions' : reactions,
             'propensities' : props,
             'offset_vectors' : offsets,
             'species' : species,
             'species counts' : species_counts,
             'np' : np}
    return model

def display_plots(recorder, title):
    pylab.figure()
    for measurement in recorder.measurements('species'):
        pylab.plot(measurement.times,
                   measurement.expected_value,
                   label = measurement.name)
    pylab.legend()
    pylab.title(title+': species count expected value')
    
    pylab.figure()
    for measurement in recorder.measurements('species'):
        pylab.plot(measurement.times,
                   measurement.standard_deviation,
                   label = measurement.name)
    pylab.legend()
    pylab.title(title+': species count standard deviation')
    

def main_gou07_c():

    max_p_trunc = 20
    max_q_trunc = 20
    
    model = create_gou07_c_model(max_p_trunc, max_q_trunc)
    
    solver = cme_solver.create_cme_solver(model)
    recorder = cme_recorder.CmeRecorder(model)
    recorder.add_target('species',
                        ['expected value', 'standard deviation'],
                        model['species'][0:2],
                        model['species counts'][0:2])
    time_steps = numpy.linspace(0.0, 1000.0, 101)
    for t in time_steps:
        print ('stepping to t = %f' % t)
        solver.step(t)
        recorder.write(t, solver.y)
    
    display_plots(recorder, model['doc'])
    pylab.show()

def main_gou07_d():

    max_p_trunc = 30
    max_q_trunc = 20
    max_s_trunc = 20
    
    model = create_gou07_d_model(max_p_trunc, max_q_trunc, max_s_trunc)
    
    # we need to explicitly specify a non standard p_0
    p_0 = numpy.zeros(model['np'])
    # initial distribution concentrated at 0 copies of P, Q, maximal copies of S
    p_0[0, 0, -1] = 1.0
    solver = cme_solver.create_cme_solver(model, p_0)
    
    recorder = cme_recorder.CmeRecorder(model)
    recorder.add_target('species',
                        ['expected value', 'standard deviation'],
                        model['species'][0:3],
                        model['species counts'][0:3])
    time_steps = numpy.linspace(0.0, 1000.0, 101)
    for t in time_steps:
        print ('stepping to t = %f' % t)
        solver.step(t)
        recorder.write(t, solver.y)
    
    display_plots(recorder, model['doc'])
    pylab.figure()
    pylab.contourf(numpy.add.reduce(solver.y, axis=2))
    pylab.title('P & Q')
    pylab.figure()
    pylab.contourf(numpy.add.reduce(solver.y, axis=1))
    pylab.title('P & S')
    pylab.figure()
    pylab.contourf(numpy.add.reduce(solver.y, axis=0))
    pylab.title('Q & S')
    
    pylab.show()

def main():
    main_gou07_d()

if __name__ == '__main__':
    main()
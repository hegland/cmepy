import numpy
import pylab
from cmepy import model, solver, recorder

def main():
    initial_copies = 20
    
    m = model.create(
        propensities = [lambda x: initial_copies - x],
        transitions = [(1, )],
        shape = (initial_copies + 1, ),
        initial_state = (0, )
    )
    
    s = solver.create(
        model = m,
        sink = False
    )
    
    r = recorder.create(
        (('A->B', ), )
    )
    
    time_steps = numpy.linspace(0.0, 3.0, 6)
    for t in time_steps:
        s.step(t)
        r.write(t, s.y)
    
    pylab.figure()
    for t, d in zip(r['A->B'].times, r['A->B'].distributions):
        marginal = d.to_dense(m.shape)
        pylab.plot(marginal, label = 't = %.1f' % t)
    pylab.xlabel('Reaction count')
    pylab.ylabel('Probability')
    pylab.legend()
    pylab.savefig('simple_plot.png')

if __name__ == '__main__':
    main()
